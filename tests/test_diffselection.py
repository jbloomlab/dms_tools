"""Tests ``dms_diffselectin`` program.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import numpy
import pandas


class TestDiffSelection(unittest.TestCase):
    """Runs ``dms_diffselection`` on test data.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.subdir = './test_diffselection_files/'
        self.mock = self.subdir + 'L1_mock_r1counts.txt'
        self.assertTrue(os.path.isfile(self.mock), 'Cannot find required simulation file {0}'.format(self.mock))
        self.concentrations = ['c1', 'c3']
        self.mincounts = [0, 10]
        self.selected = self.subdir + 'L1_H17L19_{concentration}_r1counts.txt'
        self.expecteddir = self.subdir + '/expected_output/'
        self.no_err_counts = self.subdir + 'no_err_counts.txt'
        self.err_counts = self.subdir + 'err_counts.txt'
        for c in self.concentrations:
            fname = self.selected.format(concentration=c)
            self.assertTrue(os.path.isfile(fname), "Cannot find required simulation file {0}".format(fname))
            for mincounts in self.mincounts:
                fprefix = self.expecteddir + '/{0}_mincounts{1}'.format(c, mincounts)
                for fsuffix in ['mutdiffsel.txt', 'sitediffsel.txt']:
                    fname = fprefix + fsuffix
                    self.assertTrue(os.path.isfile(fname), "Cannot find required simulation file {0}".format(fname))


    def test_DiffSelection(self):
        """Runs ``dms_diffselection``."""
        for c in self.concentrations:
            selected = self.selected.format(concentration=c)
            for mincounts in self.mincounts:
                outprefix = '{0}/{1}_mincounts{2}'.format(self.subdir, c, mincounts)
                cmds = ['dms_diffselection', self.mock, selected, outprefix, '--mincounts', str(mincounts)]
                sys.stderr.write('\nRunning the following command:\n%s\n' % ' '.join(cmds))
                subprocess.check_call(cmds)
                for (suffix, columntotest) in [
                        ('mutdiffsel.txt', 'diffsel'),
                        ('sitediffsel.txt', 'abs_diffsel')]:
                    f = outprefix + suffix
                    self.assertTrue(os.path.isfile(f), 'Failed to create {0} with command:\n{1}\n'.format(f, ' '.join(cmds)))
                    df = pandas.read_csv(f)
                    expected_f = '{0}/{1}_mincounts{2}{3}'.format(self.expecteddir, c, mincounts, suffix)
                    expected_df = pandas.read_csv(expected_f)
                    self.assertTrue(numpy.allclose(numpy.nan_to_num(df[columntotest].values), numpy.nan_to_num(expected_df[columntotest].values)), "Did not get expected results for comparing {0} and {1}".format(f, expected_f))


    def test_ErrorCorrection_NoErrors(self):
        """Tests ``--errorcontrolcounts`` option when there are no errors in control."""
        for c in self.concentrations:
            selected = self.selected.format(concentration=c)
            for mincounts in [0]:
                outprefix = '{0}/no_errors_{1}_mincounts{2}'.format(self.subdir, c, mincounts)
                cmds = ['dms_diffselection', self.mock, selected, outprefix, '--mincounts', str(mincounts), '--errorcontrolcounts', self.no_err_counts]
                sys.stderr.write('\nRunning the following command:\n%s\n' % ' '.join(cmds))
                subprocess.check_call(cmds)
                for (suffix, columntotest) in [
                        ('mutdiffsel.txt', 'diffsel'),
                        ('sitediffsel.txt', 'abs_diffsel')]:
                    f = outprefix + suffix
                    self.assertTrue(os.path.isfile(f), 'Failed to create {0} with command:\n{1}\n'.format(f, ' '.join(cmds)))
                    df = pandas.read_csv(f)
                    expected_f = '{0}/{1}_mincounts{2}{3}'.format(self.expecteddir, c, mincounts, suffix)
                    expected_df = pandas.read_csv(expected_f)
                    self.assertTrue(numpy.allclose(numpy.nan_to_num(df[columntotest].values), numpy.nan_to_num(expected_df[columntotest].values)), "Did not get expected results for comparing {0} and {1}".format(f, expected_f))

    def test_ErrorCorrection_WithErrors(self):
        """Tests ``--errorcontrolcounts`` option when there are errors in control."""
        for c in self.concentrations:
            selected = self.selected.format(concentration=c)
            for mincounts in [0]:
                outprefix = '{0}/errors_{1}_mincounts{2}'.format(self.subdir, c, mincounts)
                cmds = ['dms_diffselection', self.mock, selected, outprefix, '--mincounts', str(mincounts), '--errorcontrolcounts', self.err_counts]
                sys.stderr.write('\nRunning the following command:\n%s\n' % ' '.join(cmds))
                subprocess.check_call(cmds)
                for (suffix, columntotest) in [
                        ('mutdiffsel.txt', 'diffsel'),
                        ('sitediffsel.txt', 'abs_diffsel')]:
                    f = outprefix + suffix
                    self.assertTrue(os.path.isfile(f), 'Failed to create {0} with command:\n{1}\n'.format(f, ' '.join(cmds)))
                    df = pandas.read_csv(f)
                    expected_f = '{0}/errors_{1}_mincounts{2}{3}'.format(self.expecteddir, c, mincounts, suffix)
                    expected_df = pandas.read_csv(expected_f)
                    self.assertTrue(numpy.allclose(numpy.nan_to_num(df[columntotest].values), numpy.nan_to_num(expected_df[columntotest].values)), "Did not get expected results for comparing {0} and {1}".format(f, expected_f))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
