"""Tests inference of differential preferences by ``dms_inferdiffprefs`` program.

Uses data from ``../examples/Melnikov_et_al/infer_diffprefs_on_simulated_data``.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess



class TestInferDiffPrefs(unittest.TestCase):
    """Runs ``dms_inferdiffprefs`` on test data.
    
    Makes sure the program runs correctly, and gives values close
    to those used to simulate the counts.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.subdir = './test_inferdiffprefs_files/'
        self.nstart = '%s/simulated_depth_1e+07_start.txt' % self.subdir
        self.ns1 = '%s/simulated_depth_1e+07_s1.txt' % self.subdir
        self.ns2 = '%s/simulated_depth_1e+07_s2.txt' % self.subdir
        self.actualdiffprefs = '%s/actual_diff_prefs.txt' % self.subdir
        for f in [self.nstart, self.ns1, self.ns2, self.actualdiffprefs]:
            self.assertTrue(os.path.isfile(f), 'Cannot find required simulation file %s' % f)
        self.inferred = '%s/inferred_diffprefs.txt' % self.subdir
        self.ratioinferred = '%s/ratio_inferred_diffprefs.txt' % self.subdir
        self.corr = '%s/inferred_to_actual_corr.txt' % self.subdir
        self.ratiocorr = '%s/ratio_inferred_to_actual_corr.txt' % self.subdir
        for f in [self.inferred, self.ratioinferred, self.ratiocorr, self.corr]:
            if os.path.isfile(f):
                os.remove(f)


    def test_InferDiffPrefs(self):
        """Runs ``dms_inferdiffprefs``."""
        cmds = ['dms_inferdiffprefs', self.nstart, self.ns1, self.ns2, self.inferred, '--ncpus', '-1']
        self.assertFalse(os.path.isfile(self.inferred), '%s already exists' % self.inferred)
        sys.stderr.write('\nRunning the following command:\n%s\n' % ' '.join(cmds))
        subprocess.call(cmds)
        self.assertTrue(os.path.isfile(self.inferred), 'Failed to create file %s with command:\n%s\n' % (self.inferred, ' '.join(cmds)))

        self.assertFalse(os.path.isfile(self.corr), '%s already exists' % self.corr)
        cmds = ['dms_correlate', self.actualdiffprefs, self.inferred, os.path.splitext(self.corr)[0], '--name1', 'actual', '--name2', 'inferred', '--corr_on_plot']
        subprocess.call(cmds)
        self.assertTrue(os.path.isfile(self.corr), 'Failed to create file %s with command:\n%s\n' % (self.corr, ' '.join(cmds)))
        with open(self.corr) as f:
            corr = float(f.readlines()[0].split('=')[1])
        self.assertTrue(corr > 0.95, 'The correlation of %g in %s is too low.' % (corr, self.corr))
        sys.stderr.write('\nTest passed. There is a good correlation of %g between the actual and inferred values.\n' % corr) 

        self.assertFalse(os.path.isfile(self.ratioinferred), '%s already exists' % self.ratioinferred)
        cmds = ['dms_inferdiffprefs', self.nstart, self.ns1, self.ns2, self.ratioinferred, '--ncpus', '-1', '--ratio_estimation', '1']
        sys.stderr.write('\nRunning the following command:\n%s\n' % ' '.join(cmds))
        subprocess.call(cmds)
        self.assertTrue(os.path.isfile(self.ratioinferred), 'Failed to create file %s with command:\n%s\n' % (self.ratioinferred, ' '.join(cmds)))

        self.assertFalse(os.path.isfile(self.ratiocorr), '%s already exists' % self.ratiocorr)
        cmds = ['dms_correlate', self.actualdiffprefs, self.ratioinferred, os.path.splitext(self.corr)[0], '--name1', 'actual', '--name2', 'ratio-inferred', '--corr_on_plot']
        subprocess.call(cmds)
        self.assertTrue(os.path.isfile(self.ratiocorr), 'Failed to create file %s with command:\n%s\n' % (self.corr, ' '.join(cmds)))
        with open(self.ratiocorr) as f:
            corr = float(f.readlines()[0].split('=')[1])
        self.assertTrue(corr > 0.92, 'The correlation of %g in %s for ratio estimation is too low.' % (corr, self.ratiocorr))
        sys.stderr.write('\nTest passed. There is a good correlation of %g between the actual and ratio inferred values.\n' % corr)


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
