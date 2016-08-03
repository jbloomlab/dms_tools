"""Tests merging of differential selection by ``dms_merge`` program.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import numpy
import pandas


class TestMergeCounts(unittest.TestCase):
    """Runs ``dms_merge`` on test data to average differential selection values.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.subdir = './test_merge_diffsel_files/'
        self.diffsel1 = '{0}/diffsel1.txt'.format(self.subdir)
        self.diffsel2 = '{0}/diffsel2.txt'.format(self.subdir)
        self.diffsel3 = '{0}/diffsel3.txt'.format(self.subdir)
        self.expectedavg = '{0}/avg.txt'.format(self.subdir)
        self.expectedsiteavg = '{0}/siteavg.txt'.format(self.subdir)
        for f in [self.diffsel1, self.diffsel2, self.diffsel3, self.expectedavg, self.expectedsiteavg]:
            self.assertTrue(os.path.isfile(f), 'Cannot find required file {0}'.format(f))


    def test_Avg(self):
        """Runs ``dms_merge`` to average diffsel."""
        out = self.subdir + 'test_avg.txt'
        siteout = self.subdir + 'test_siteavg.txt'
        cmds = ['dms_merge', out, 'average', self.diffsel1, self.diffsel2, self.diffsel3, '--sitediffselfile', siteout]
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(' '.join(cmds)))
        subprocess.check_call(cmds)
        self.assertTrue(os.path.isfile(out), 'Failed to create file {0} with command:\n{1}\n'.format(out, ' '.join(cmds)))

        expect = pandas.read_csv(self.expectedavg)['diffsel'].values
        actual = pandas.read_csv(out)['diffsel'].values
        expect = expect[numpy.logical_not(numpy.isnan(expect))]
        actual = actual[numpy.logical_not(numpy.isnan(actual))]
        self.assertTrue(numpy.allclose(expect, actual), "Did not get expected average mutdiffsel values")
        sys.stderr.write('Averaging of mutdiffsel values successful\n')

        expect = pandas.read_csv(self.expectedsiteavg)['abs_diffsel'].values
        actual = pandas.read_csv(siteout)['abs_diffsel'].values
        expect = expect[numpy.logical_not(numpy.isnan(expect))]
        actual = actual[numpy.logical_not(numpy.isnan(actual))]
        self.assertTrue(numpy.allclose(expect, actual), "Did not get expected average absolute sitediffsel values")
        sys.stderr.write('Averaging of sitediffsel values successful\n')


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
