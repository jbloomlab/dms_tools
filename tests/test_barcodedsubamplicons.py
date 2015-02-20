"""Tests ``dms_barcodedsubamplicons``.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess



class TestBarcodedSubamplicons(unittest.TestCase):
    """Runs ``dms_barcodedsubamplicons`` on test data.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.subdir = './test_barcodedsubamplicons_files/'
        self.refseq = '%s/WSN-HA.fasta' % self.subdir
        self.r1 = '%s/R1.fastq.gz' % self.subdir
        self.r2 = '%s/R2.fastq.gz' % self.subdir
        self.alignspecs = '1,426,36,38 427,849,32,32 850,1275,31,37 1276,1698,46,45'
        self.correctcounts = '%s/correct_counts.txt' % self.subdir
        self.correctstats = '%s/correct_summarystats.txt' % self.subdir
        for f in [self.refseq, self.r1, self.r2, self.correctcounts, self.correctstats]:
            self.assertTrue(os.path.isfile(f), 'Cannot find required file %s' % f)
        self.prefix = '%s/alignment_' % self.subdir
        self.counts = '%scounts.txt' % self.prefix
        self.stats = '%ssummarystats.txt' % self.prefix
        for f in [self.counts, self.stats]:
            if os.path.isfile(f):
                os.remove(f)


    def test_BarcodedSubamplicons(self):
        """Runs ``dms_barcodedsubamplicons``."""
        cmds = ['dms_barcodedsubamplicons', self.prefix, self.refseq, self.r1, self.r2, self.alignspecs]
        self.assertFalse(os.path.isfile(self.counts), '%s already exists' % self.counts)
        self.assertFalse(os.path.isfile(self.stats), '%s already exists' % self.stats)
        sys.stderr.write('\nRunning the following command:\n%s\n' % ' '.join(cmds))
        os.system(' '.join(cmds))
        self.assertTrue(os.path.isfile(self.counts), 'Failed to create file %s with command:\n%s\n' % (self.counts, ' '.join(cmds)))
        self.assertTrue(os.path.isfile(self.stats), 'Failed to create file %s with command:\n%s\n' % (self.stats, ' '.join(cmds)))
        with open(self.counts) as f:
            counts = f.read()
        with open(self.correctcounts) as f:
            correctcounts = f.read()
        with open(self.stats) as f:
            stats = f.read()
        with open(self.correctstats) as f:
            correctstats = f.read()
        self.assertEqual(counts, correctcounts, "Counts did not match expected counts")
        self.assertEqual(stats, correctstats, "Counts did not match expected counts")


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
