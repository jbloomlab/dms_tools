"""Tests merging of prefs by ``dms_merge`` program.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import dms_tools.file_io


class TestMergeCounts(unittest.TestCase):
    """Runs ``dms_merge`` on test data to average or rescale prefs.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.subdir = './test_merge_prefs_files/'
        self.prefs1 = '{0}/prefs1.txt'.format(self.subdir)
        self.prefs2 = '{0}/prefs2.txt'.format(self.subdir)
        self.expectedavg = '{0}/avg.txt'.format(self.subdir)
        self.expectedrescale = '{0}/prefs1_rescaled.txt'.format(self.subdir)
        for f in [self.prefs1, self.prefs2, self.expectedavg, self.expectedrescale]:
            self.assertTrue(os.path.isfile(f), 'Cannot find required file {0}'.format(f))


    def test_Avg(self):
        """Runs ``dms_merge`` to average prefs."""
        out = self.subdir + 'test_avg.txt'
        cmds = ['dms_merge', out, 'average', self.prefs1, self.prefs2]
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(' '.join(cmds)))
        subprocess.check_call(cmds)
        self.assertTrue(os.path.isfile(out), 'Failed to create file {0} with command:\n{1}\n'.format(out, ' '.join(cmds)))
        with open(out) as f_out, open(self.expectedavg) as f_expected:
            self.assertTrue(f_out.read() == f_expected.read(), "Failed to get expected result from averaging: {0} doesn't match {1}".format(out, self.expectedavg))
        sys.stderr.write('Successfully averaged prefs.')


    def test_Rescale(self):
        """Runs ``dms_merge`` to rescale prefs."""
        out = self.subdir + 'test_prefs1_rescaled.txt'
        cmds = ['dms_merge', out, 'rescale', self.prefs1, '--stringencyparameter', '2.0']
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(' '.join(cmds)))
        subprocess.check_call(cmds)
        self.assertTrue(os.path.isfile(out), 'Failed to create file {0} with command:\n{1}\n'.format(out, ' '.join(cmds)))
        with open(out) as f_out, open(self.expectedrescale) as f_expected:
            self.assertTrue(f_out.read() == f_expected.read(), "Failed to get expected result from averaging: {0} doesn't match {1}".format(out, self.expectedrescale))
        sys.stderr.write('Successfully rescaled prefs.')

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
