"""Tests merging of counts by ``dms_merge`` program.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import dms_tools.file_io


class TestMergeCounts(unittest.TestCase):
    """Runs ``dms_merge`` on test data to add counts.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.subdir = './test_merge_counts_files/'
        self.counts1 = '{0}/counts1.txt'.format(self.subdir)
        self.counts2 = '{0}/counts2.txt'.format(self.subdir)
        self.expected = '{0}/expected_merge.txt'.format(self.subdir)
        self.expectednorm = '{0}/expected_merge_normalized.txt'.format(self.subdir)
        for f in [self.counts1, self.counts2, self.expected, self.expectednorm]:
            self.assertTrue(os.path.isfile(f), 'Cannot find required file {0}'.format(f))
        self.testmerge = '{0}/test_merge.txt'.format(self.subdir)
        self.testmergenorm = '{0}/test_merge_normalized.txt'.format(self.subdir)
        for f in [self.testmerge, self.testmergenorm]:
            if os.path.isfile(f):
                os.remove(f)


    def test_Merge(self):
        """Runs ``dms_merge``."""
        for (out, expected, extracmd) in [
                (self.testmerge, self.expected, []),
                (self.testmergenorm, self.expectednorm, ['--normalize']),
                ]:
            cmds = ['dms_merge', out, 'sum', self.counts1, self.counts2] + extracmd
            sys.stderr.write('\nRunning the following command:\n{0}\n'.format(' '.join(cmds)))
            subprocess.call(cmds)
            self.assertTrue(os.path.isfile(out), 'Failed to create file {0} with command:\n{1}\n'.format(out, ' '.join(cmds)))
            actual = dms_tools.file_io.ReadDMSCounts(out, chartype='codon')
            expect = dms_tools.file_io.ReadDMSCounts(expected, chartype='codon')
            sites = set(actual.keys())
            self.assertTrue(sites == set(expect.keys()), 'Created output {0} does not match expected output {1} in terms of the sites keys.'.format(out, expected))
            for site in sites:
                self.assertTrue(actual[site] == expect[site], 'Created output {0} does not match expected output {1} at site {2}'.format(out, expected, site))
            sys.stderr.write('Successfully created expected output.')


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
