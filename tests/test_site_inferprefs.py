"""Tests inference of preferences.

Tests the Python functions on specific sites, **not** the
full ``dms_inferprefs`` program.

Written by Jesse Bloom."""


import sys
import os
import unittest
import random
import numpy
import dms_tools.mcmc
import dms_tools
aminoacids = dms_tools.aminoacids_nostop
nts = dms_tools.nts
codons = dms_tools.codons



class TestInferSitePreferences(unittest.TestCase):
    """Tests preference inference when priors are exactly correct.

    Makes sure inference converges to values that are very close
    to the correct one for all three error models:

        1) No errors

        2) Same errors pre- and post-selection

        3) Different errors pre- and post-selection

    Then also makes sure that use a no-error model for data with
    errors does **not** converge.

    Reads information about the random number seeds, character types,
    and number of CPUs from ``test_configuration.txt``.

    """

    def setUp(self):
        """Reads information from ``test_configuration.txt``."""
        fname = 'test_configuration.txt'
        if not os.path.isfile(fname):
            raise ValueError('Cannot find test configuration file %s in current directory.' % fname)
        lines = [line for line in open(fname) if not line.isspace()]
        self.seeds = None
        self.n_jobs = None
        self.charactertypes = None
        for line in lines:
            (key, value) = line.split(None, 1)
            if key == 'seeds':
                self.seeds = [int(x) for x in value.split(',')]
            elif key == 'charactertypes':
                self.charactertypes = set([x.strip() for x in value.split(',')])
                for x in self.charactertypes:
                    validvalues = ['DNA', 'amino acid', 'codon']
                    if x not in validvalues:
                        raise ValueError("Invalid charactertype of %s in %s. Valid values:\n%s" % (x, fname, ', '.join(validvalues)))
            elif key == 'n_jobs':
                self.n_jobs = int(value)
                if not (self.n_jobs == -1 or self.n_jobs >= 1):
                    raise ValueError('%s specifies invalid value of %d for n_jobs. Must be -1 or >= 1' % (fname, self.n_jobs))
            else:
                raise ValueError('Invalid line in %s:\n%s' % (fname, line))
        if self.seeds == None:
            raise ValueError('%s failed to specify seeds' % fname)
        if self.charactertypes == None:
            raise ValueError('%s failed to specify charactertypes' % fname)
        if self.n_jobs == None:
            raise ValueError('%s failed to specify n_jobs' % fname)


    def test_InferenceSitePreferences(self):
        """Inference with with correct priors."""
        mutrate = (0.005, 0.025) # mutation rate drawn uniformly from this range after dividing by nchars
        errorratepre = (0.0001, 0.0003) # pre-selection error rate drawn uniformly from this range after dividing by nchars
        errorratepost = (0.0002, 0.0005) # pre-selection error rate drawn uniformly from this range after dividing by nchars
        pidist = (1e-5, 0.6) # pi drawn uniformly from this range (after correcting for proportionality)
        for seed in self.seeds:
            random.seed(seed)
            numpy.random.seed(seed)
            for charactertype in self.charactertypes:
                characterlist = {'DNA':nts, 'amino acid':aminoacids, 'codon':codons}[charactertype]
                sys.stderr.write('\n\nTesting preference inference at the %s level for seed %d...\n' % (charactertype, seed)) 
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
                    (converged, pi_means, pi_95credint, logstring) = dms_tools.mcmc.InferSitePreferences(characterlist, wtchar, 'none', counts, priors, n_jobs=self.n_jobs)
                    self.assertTrue(converged, 'MCMC inference of preferences failed to converge for no-errors data for %s at depth of %d for seed %d' % (charactertype, depth, seed))
                    inferred_pi = [pi_means[char] for char in characterlist]
                    diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
                    self.assertTrue(diffsum < maxdiffsum, 'MCMC inference of preferences had excessive large summed absolute difference of %g (max allowed is %g) from true values for no-errors data for %s at depth of %d for seed %d' % (diffsum, maxdiffsum, charactertype, depth, seed))
                    sys.stderr.write('Converged to correct values for no-errors data (diffsum = %g) at depth %d for character type %s for seed %d.\n' % (diffsum, depth, charactertype, seed))
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
                    (converged, pi_means, pi_95credint, logstring) = dms_tools.mcmc.InferSitePreferences(characterlist, wtchar, 'same', counts, priors, n_jobs=self.n_jobs)
                    self.assertTrue(converged, 'MCMC inference of preferences failed to converge for same-errors data for %s at depth of %d for seed %d' % (charactertype, depth, seed))
                    inferred_pi = [pi_means[char] for char in characterlist]
                    diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
                    self.assertTrue(diffsum < maxdiffsum, 'MCMC inference of preferences had excessive large summed absolute difference of %g (max allowed is %g) from true values for same-errors data for %s at depth of %d for seed %d' % (diffsum, maxdiffsum, charactertype, depth, seed))
                    sys.stderr.write('Converged to correct values for same-errors data (diffsum = %g) at depth %d for character type %s for seed %d.\n' % (diffsum, depth, charactertype, seed))
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
                    (converged, pi_means, pi_95credint, logstring) = dms_tools.mcmc.InferSitePreferences(characterlist, wtchar, 'different', counts, priors, n_jobs=self.n_jobs)
                    self.assertTrue(converged, 'MCMC inference of preferences failed to converge for different-errors data for %s at depth of %d for seed %d' % (charactertype, depth, seed))
                    inferred_pi = [pi_means[char] for char in characterlist]
                    diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
                    self.assertTrue(diffsum < maxdiffsum, 'MCMC inference of preferences had excessive large summed absolute difference of %g (max allowed is %g) from true values for different-errors data for %s at depth of %d for seed %d' % (diffsum, maxdiffsum, charactertype, depth, seed))
                    sys.stderr.write('Converged to correct values for different-errors data (diffsum = %g) at depth %d for character type %s for seed %d.\n' % (diffsum, depth, charactertype, seed))
                    sys.stderr.flush()
                    # make sure sample with errors does NOT converge without error estimates
                    (converged, pi_means, pi_95credint, logstring) = dms_tools.mcmc.InferSitePreferences(characterlist, wtchar, 'none', counts, priors, n_jobs=self.n_jobs)
                    self.assertTrue(converged, 'MCMC inference of preferences failed to converge for different-errors data analyzed with no-errors model for %s at depth of %d for seed %d' % (charactertype, depth, seed))
                    inferred_pi = [pi_means[char] for char in characterlist]
                    diffsum = sum([abs(x - y) for (x, y) in zip(pir, inferred_pi)])
                    self.assertFalse(diffsum < maxdiffsum, 'MCMC inference of preferences with different-errors data converged (summed absolute difference of %g) even when using a no-errors model for %s at depth of %d for seed %d' % (diffsum, charactertype, depth, seed))
                    sys.stderr.write('As expected, analyzing the different-errors data with the no-errors model did not converge to the correct values (summed absolute difference of %g) at depth %d for character type %s for seed %d. This is good.\n' % (diffsum, depth, charactertype, seed))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
