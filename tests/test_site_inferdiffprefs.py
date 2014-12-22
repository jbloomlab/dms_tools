"""Tests inference of differential preferences.

Tests the Python functions on specific sites, **not** the
full ``dms_inferdiffprefs`` program.

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





class TestInferSiteDiffPreferences(unittest.TestCase):
    """Tests differential preference inference when priors are exactly correct.

    Makes sure inference converges to values that are very close
    to the correct ones both with and without seuqencing errors
    included.

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


    def test_InferenceSiteDiffPreferences(self):
        """Inference with with correct priors."""
        errorrate = (0.0001, 0.0003) # error rate drawn uniformly from this range after dividing by nchars
        pidist = (0.2, 0.8) # pi drawn uniformly from this range (after correcting for proportionality)
        frange = (0.05, 1.0)
        depth = 1e9
        maxdiffsum = 0.02
        for seed in self.seeds:
            random.seed(seed)
            numpy.random.seed(seed)
            for charactertype in self.charactertypes:
                characterlist = {'DNA':nts, 'amino acid':aminoacids, 'codon':codons}[charactertype]
                sys.stderr.write('\n\nTesting differential preference inference at the %s level for seed %d...\n' % (charactertype, seed)) 
                sys.stderr.flush()
                wtchar = random.choice(characterlist)
                iwtchar = characterlist.index(wtchar)
                nchars = len(characterlist)
                fr = []
                pir = []
                xir = []
                for i in range(nchars):
                    fr.append(random.uniform(frange[0], frange[1]) / float(nchars))
                    xir.append(random.uniform(errorrate[0] / float(nchars), errorrate[1] / float(nchars)))
                    pir.append(random.uniform(pidist[0], pidist[1]))
                fr[iwtchar] = 1.0 - sum([fr[i] for i in range(nchars) if i != iwtchar])
                xir[iwtchar] = 1.0 - sum([xir[i] for i in range(nchars) if i != iwtchar])
                pir[iwtchar] = 1.0
                pir = [x / float(sum(pir)) for x in pir]
                fr = numpy.array(fr)
                pir = numpy.array(pir)
                xir = numpy.array(xir)
                deltar = numpy.zeros(nchars)
                deltar[iwtchar] = 1.0
                priors = {\
                    'pirs1_prior_params':dict([(char, pir[i] * nchars) for (i, char) in enumerate(characterlist)]),\
                    'frstart_prior_params':dict([(char, fr[i] * nchars) for (i, char) in enumerate(characterlist)]),\
                    'xir_prior_params':dict([(char, xir[i] * nchars) for (i, char) in enumerate(characterlist)]),\
                    'deltapi_concentration':5.0,\
                    }
                for deltapir in ['zero', 'non-zero']:
                    sys.stderr.write('Testing with deltapi that is %s...\n' % deltapir)
                    if deltapir == 'zero':
                        deltapir = numpy.zeros(nchars)
                    else:
                        deltapir = numpy.random.dirichlet(pir * nchars * priors['deltapi_concentration']) - pir
                    # first with no error
                    nrstart = numpy.random.multinomial(depth, fr)
                    nrs1 = numpy.random.multinomial(depth, fr * pir / numpy.dot(fr, pir))
                    nrs2 = numpy.random.multinomial(depth, fr * (deltapir + pir) / numpy.dot(fr, deltapir + pir))
                    counts = {\
                        'nrstart':dict([(char, nrstart[i]) for (i, char) in enumerate(characterlist)]),\
                        'nrs1':dict([(char, nrs1[i]) for (i, char) in enumerate(characterlist)]),\
                        'nrs2':dict([(char, nrs2[i]) for (i, char) in enumerate(characterlist)]),\
                        }
                    (converged, deltapi_means, pr_deltapi_gt0, pr_deltapi_lt0, logstring) = dms_tools.mcmc.InferSiteDiffPrefs(characterlist, wtchar, 'none', counts, priors, n_jobs=self.n_jobs)
                    self.assertTrue(converged, 'MCMC inference of differential preferences failed to converge for no-errors data for %s at depth of %d for seed %d' % (charactertype, depth, seed))
                    inferred_deltapi = [deltapi_means[char] for char in characterlist]
                    diffsum = sum([abs(x - y) for (x, y) in zip(deltapir, inferred_deltapi)])
                    self.assertTrue(diffsum < maxdiffsum, 'MCMC inference of differential preferences had excessive large summed absolute difference of %g (max allowed is %g) from true values for no-errors data for %s at depth of %d for seed %d' % (diffsum, maxdiffsum, charactertype, depth, seed))
                    sys.stderr.write('Converged to correct values for no-errors data (diffsum = %g) at depth %d for character type %s for seed %d.\n' % (diffsum, depth, charactertype, seed))
                    sys.stderr.flush()
                    # now with errors
                    nrstart = numpy.random.multinomial(depth, fr + xir - deltar)
                    nrs1 = numpy.random.multinomial(depth, fr * pir / numpy.dot(fr, pir) + xir - deltar)
                    nrs2 = numpy.random.multinomial(depth, fr * (deltapir + pir) / numpy.dot(fr, deltapir + pir) + xir - deltar)
                    nrerr = numpy.random.multinomial(depth, xir)
                    counts = {\
                        'nrstart':dict([(char, nrstart[i]) for (i, char) in enumerate(characterlist)]),\
                        'nrs1':dict([(char, nrs1[i]) for (i, char) in enumerate(characterlist)]),\
                        'nrs2':dict([(char, nrs2[i]) for (i, char) in enumerate(characterlist)]),\
                        'nrerr':dict([(char, nrerr[i]) for (i, char) in enumerate(characterlist)]),\
                        }
                    (converged, deltapi_means, pr_deltapi_gt0, pr_deltapi_lt0, logstring) = dms_tools.mcmc.InferSiteDiffPrefs(characterlist, wtchar, 'same', counts, priors, n_jobs=self.n_jobs)
                    self.assertTrue(converged, 'MCMC inference of differential preferences failed to converge for same-errors data for %s at depth of %d for seed %d' % (charactertype, depth, seed))
                    inferred_deltapi = [deltapi_means[char] for char in characterlist]
                    diffsum = sum([abs(x - y) for (x, y) in zip(deltapir, inferred_deltapi)])
                    self.assertTrue(diffsum < maxdiffsum, 'MCMC inference of differential preferences had excessive large summed absolute difference of %g (max allowed is %g) from true values for same-errors data for %s at depth of %d for seed %d' % (diffsum, maxdiffsum, charactertype, depth, seed))
                    sys.stderr.write('Converged to correct values for same-errors data (diffsum = %g) at depth %d for character type %s for seed %d.\n' % (diffsum, depth, charactertype, seed))
                    sys.stderr.write('Now checking posterior probabilities of differential preferences greater or less than zero\n')
                    threshold = 0.05 # only worry about deltas at least this big
                    pr_cutoff = 0.95 # threshold for pr
                    for (ichar, char) in enumerate(characterlist):
                        prdeltapistring = 'For %s, deltapi = %g, pr_deltapi_gt0 = %g, and pr_deltapi_lt0 = %g' % (char, deltapir[ichar], pr_deltapi_gt0[char], pr_deltapi_lt0[char])
                        sys.stderr.write('%s\n'% prdeltapistring)
                        if deltapir[ichar] > threshold:
                            self.assertTrue(pr_deltapi_gt0[char] > pr_cutoff and pr_deltapi_lt0[char] < 1.0 - pr_cutoff, "Incorrect probabilities: %s" % prdeltapistring)
                        elif deltapir[ichar] < -threshold:
                            self.assertTrue(pr_deltapi_lt0[char] > pr_cutoff and pr_deltapi_gt0[char] < 1.0 - pr_cutoff, "Incorrect probabilities: %s" % prdeltapistring)
                    sys.stderr.flush()



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
