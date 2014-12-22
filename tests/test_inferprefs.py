"""Tests inference of preferences by ``dms_inferprefs`` program.

Uses data from ``../examples/Melnikov_et_al/infer_preferences_on_simulated_data``.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess



class TestInferPrefs(unittest.TestCase):
    """Runs ``dms_inferprefs`` on test data.
    
    Makes sure the program runs correctly, and gives values close
    to those used to simulate the counts.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.subdir = './test_inferprefs_files/'
        self.npre = '%s/simulated_errs_depth_2e+06_pre.txt' % self.subdir
        self.npost = '%s/simulated_errs_depth_2e+06_post.txt' % self.subdir
        self.nerrpre = '%s/simulated_errs_depth_2e+06_errpre.txt' % self.subdir
        self.nerrpost = '%s/simulated_errs_depth_2e+06_errpost.txt' % self.subdir
        self.actualprefs = '%s/actual_prefs.txt' % self.subdir
        for f in [self.npre, self.npost, self.nerrpre, self.nerrpost, self.actualprefs]:
            self.assertTrue(os.path.isfile(f), 'Cannot find required simulation file %s' % f)
        self.inferred = '%s/inferred_prefs.txt' % self.subdir
        self.corr = '%s/inferred_to_actual_corr.txt' % self.subdir
        for f in [self.inferred, self.corr]:
            if os.path.isfile(f):
                os.remove(f)


    def test_InferPrefs(self):
        """Runs ``dms_inferprefs``."""
        cmds = ['dms_inferprefs', self.npre, self.npost, self.inferred, '--errpre', self.nerrpre, '--errpost', self.nerrpost, '--ncpus', '-1', '--excludestop']
        self.assertFalse(os.path.isfile(self.inferred), '%s already exists' % self.inferred)
        sys.stderr.write('\nRunning the following command:\n%s\n' % ' '.join(cmds))
        subprocess.call(cmds)
        self.assertTrue(os.path.isfile(self.inferred), 'Failed to create file %s with command:\n%s\n' % (self.inferred, ' '.join(cmds)))
        self.assertFalse(os.path.isfile(self.corr), '%s already exists' % self.corr)
        cmds = ['dms_correlate', self.actualprefs, self.inferred, os.path.splitext(self.corr)[0], '--name1', 'actual', '--name2', 'inferred', '--corr_on_plot']
        subprocess.call(cmds)
        self.assertTrue(os.path.isfile(self.corr), 'Failed to create file %s with command:\n%s\n' % (self.corr, ' '.join(cmds)))
        with open(self.corr) as f:
            corr = float(f.readlines()[0].split('=')[1])
        self.assertTrue(corr > 0.95, 'The correlation of %g in %s is too low.' % (corr, self.corr))
        sys.stderr.write('\nTest passed. There is a good correlation of %g between the actual and inferred values.\n' % corr)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
