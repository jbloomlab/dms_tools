"""Tests inference of preferences.

Written by Jesse Bloom."""


import sys
import unittest
import random
import Bio.Alphabet.IUPAC
aminoacids = Bio.Alphabet.IUPAC.IUPACProtein.letters
nts = Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters
codons = []
for nt1 in nts:
    for nt2 in nts:
        for nt3 in nts:
            codons.append('%s%s%s' % (nt1, nt2, nt3))
import numpy
import dms_tools.mcmc


class TestInferSitePreferences(unittest.TestCase):
    """Tests preference inference when priors are exactly correct."""

    def test_InferenceSitePreferences(self):
        """Inference with with correct priors."""
        seed = 1 # random number seed
        n_jobs = 2 # number of CPUs to use
        random.seed(seed)
        numpy.random.seed(seed)
        mutrate = (0.005, 0.02) # mutation rate drawn uniformly from this range after dividing by nchars
        errorratepre = (0.0001, 0.0002) # pre-selection error rate drawn uniformly from this range after dividing by nchars
        errorratepost = (0.0002, 0.0004) # pre-selection error rate drawn uniformly from this range after dividing by nchars
        pidist = (1e-5, 0.6) # pi drawn uniformly from this range (after correcting for proportionality)
        for (charactertype, characterlist) in [('DNA', nts), ('amino acid', aminoacids)]:#, ('codon', codons)]:
            sys.stderr.write('\n\nTesting preference inference at the %s level...\n' % charactertype) 
            sys.stderr.flush()
            wtchar = random.choice(characterlist)
            iwtchar = characterlist.index(wtchar)
            nchars = len(characterlist)
            mur = []
            pir = []
            epsilonr = []
            rhor = []
            for i in range(nchars):
                mur.append(random.uniform(mutrate[0] / float(nchars), mutrate[1] / float(nchars)))
                epsilonr.append(random.uniform(errorratepre[0] / float(nchars), errorratepre[1] / float(nchars)))
                rhor.append(random.uniform(errorratepost[0] / float(nchars), errorratepost[1] / float(nchars)))
                pir.append(random.uniform(pidist[0], pidist[1]))
            mur[iwtchar] = 1.0 - sum([mur[i] for i in range(nchars) if i != iwtchar])
            epsilonr[iwtchar] = 1.0 - sum([epsilonr[i] for i in range(nchars) if i != iwtchar])
            rhor[iwtchar] = 1.0 - sum([rhor[i] for i in range(nchars) if i != iwtchar])
            pir[iwtchar] = 1.0
            pir = [x / float(sum(pir)) for x in pir]
            mur = numpy.array(mur)
            pir = numpy.array(pir)
            epsilonr = numpy.array(epsilonr)
            rhor = numpy.array(rhor)
            deltar = numpy.zeros(nchars)
            deltar[iwtchar] = 1.0
            priors = {\
                    'pir_prior_params':dict([(char, pir[i] * nchars) for (i, char) in enumerate(characterlist)]),\
                    'mur_prior_params':dict([(char, mur[i] * nchars) for (i, char) in enumerate(characterlist)]),\
                    'epsilonr_prior_params':dict([(char, epsilonr[i] * nchars) for (i, char) in enumerate(characterlist)]),\
                    'rhor_prior_params':dict([(char, rhor[i] * nchars) for (i, char) in enumerate(characterlist)]),\
                    }
            difflist = []
            for (depth, maxdiffsum) in [(1e9, 0.01)]:
                sys.stderr.write("Testing with sequencing depth of %d...\n" % depth)
                sys.stderr.flush()
                # first with no errors
                nrpre = numpy.random.multinomial(depth, mur)
                nrpost = numpy.random.multinomial(depth, mur * pir / numpy.dot(mur, pir))
                counts = {\
                    'nrpre':dict([(char, nrpre[i]) for (i, char) in enumerate(characterlist)]),\
                    'nrpost':dict([(char, nrpost[i]) for (i, char) in enumerate(characterlist)]),\
                    }
                returntup = dms_tools.mcmc.InferSitePreferences(characterlist, wtchar, 'none', counts, priors, n_jobs=n_jobs)
                self.assertTrue(returntup, 'MCMC inference of preferences failed to converge for no-errors data for %s at depth of %d' % (charactertype, depth))
                if returntup:
                    (pi_means, pi_95credint) = returntup
                    inferred_pi = [pi_means[char] for char in characterlist]
                    diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
                    self.assertTrue(diffsum < maxdiffsum, 'MCMC inference of preferences had excessive large summed absolute difference of %g (max allowed is %g) from true values for no-errors data for %s at depth of %d' % (diffsum, maxdiffsum, charactertype, depth))
                sys.stderr.write('Converged to correct values for no-errors data (diffsum = %g) at depth %d for character type %s.\n' % (diffsum, depth, charactertype))
                sys.stderr.flush()
                # now with same errors pre and post
                nrpre = numpy.random.multinomial(depth, mur + epsilonr - deltar)
                nrpost = numpy.random.multinomial(depth, mur * pir / numpy.dot(mur, pir) + epsilonr - deltar)
                nrerr = numpy.random.multinomial(depth, epsilonr)
                counts = {\
                    'nrpre':dict([(char, nrpre[i]) for (i, char) in enumerate(characterlist)]),\
                    'nrpost':dict([(char, nrpost[i]) for (i, char) in enumerate(characterlist)]),\
                    'nrerr':dict([(char, nrerr[i]) for (i, char) in enumerate(characterlist)]),\
                    }
                returntup = dms_tools.mcmc.InferSitePreferences(characterlist, wtchar, 'same', counts, priors, n_jobs=n_jobs)
                self.assertTrue(returntup, 'MCMC inference of preferences failed to converge for same-errors data for %s at depth of %d' % (charactertype, depth))
                if returntup:
                    (pi_means, pi_95credint) = returntup
                    inferred_pi = [pi_means[char] for char in characterlist]
                    diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
                    self.assertTrue(diffsum < maxdiffsum, 'MCMC inference of preferences had excessive large summed absolute difference of %g (max allowed is %g) from true values for same-errors data for %s at depth of %d' % (diffsum, maxdiffsum, charactertype, depth))
                sys.stderr.write('Converged to correct values for same-errors data (diffsum = %g) at depth %d for character type %s.\n' % (diffsum, depth, charactertype))
                sys.stderr.flush()
                # now with different errors pre and post
                nrpre = numpy.random.multinomial(depth, mur + epsilonr - deltar)
                nrpost = numpy.random.multinomial(depth, mur * pir / numpy.dot(mur, pir) + rhor - deltar)
                nrerrpre = numpy.random.multinomial(depth, epsilonr)
                nrerrpost = numpy.random.multinomial(depth, rhor)
                counts = {\
                    'nrpre':dict([(char, nrpre[i]) for (i, char) in enumerate(characterlist)]),\
                    'nrpost':dict([(char, nrpost[i]) for (i, char) in enumerate(characterlist)]),\
                    'nrerrpre':dict([(char, nrerrpre[i]) for (i, char) in enumerate(characterlist)]),\
                    'nrerrpost':dict([(char, nrerrpost[i]) for (i, char) in enumerate(characterlist)]),\
                    }
                returntup = dms_tools.mcmc.InferSitePreferences(characterlist, wtchar, 'different', counts, priors, n_jobs=n_jobs)
                self.assertTrue(returntup, 'MCMC inference of preferences failed to converge for different-errors data for %s at depth of %d' % (charactertype, depth))
                if returntup:
                    (pi_means, pi_95credint) = returntup
                    inferred_pi = [pi_means[char] for char in characterlist]
                    diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
                    self.assertTrue(diffsum < maxdiffsum, 'MCMC inference of preferences had excessive large summed absolute difference of %g (max allowed is %g) from true values for different-errors data for %s at depth of %d' % (diffsum, maxdiffsum, charactertype, depth))
                sys.stderr.write('Converged to correct values for different-errors data (diffsum = %g) at depth %d for character type %s.\n' % (diffsum, depth, charactertype))
                sys.stderr.flush()
                # make sure sample with errors does NOT converge without error estimates
                returntup = dms_tools.mcmc.InferSitePreferences(characterlist, wtchar, 'none', counts, priors, n_jobs=n_jobs)
                self.assertFalse(returntup, 'MCMC inference of preferences with different-errors data converged even when using a no-errors model. Absolute summed difference of %g from true values for %s at depth of %d' % (charactertype, depth))
                sys.stderr.write('As expected, analyzing the different-errors data with the no-errors model did not converge at depth %d for character type %s. This is good.\n' % (depth, charactertype))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
