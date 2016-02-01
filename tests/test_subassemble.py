"""Tests ``dms_subassemble`` and ``dms_matchbarcodedsubamplicons``.

Written by Jesse Bloom."""


import sys
import os
import unittest
import dms_tools.file_io


class TestSubasssemble(unittest.TestCase):
    """Runs ``dms_subassemble`` and ``dms_matchbarcodedsubamplicons`` on test data.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.subdir = './test_subassemble_files/'
        self.refseq = '%s/refseq.fasta' % self.subdir
        self.r1 = '%s/R1.fastq' % self.subdir
        self.r2 = '%s/R2.fastq' % self.subdir
        self.barcodes_r1 = '%s/barcodes_R1.fastq' % self.subdir
        self.barcodes_r2 = '%s/barcodes_R2.fastq' % self.subdir
        self.alignspecs = '1,3 61,1'
        self.barcodelength = 6
        self.correctsubassembly = '%s/actual_barcoded_variants.txt' % self.subdir
        self.correctstats = '%s/actual_stats.txt' % self.subdir
        self.correctallmutcounts = '%s/actual_allmutcounts.txt' % self.subdir
        self.correctsinglemutcounts = '%s/actual_singlemutcounts.txt' % self.subdir
        for f in [self.refseq, self.r1, self.r2, self.correctsubassembly, self.barcodes_r1, self.barcodes_r2, self.correctstats, self.correctallmutcounts, self.correctsinglemutcounts]:
            self.assertTrue(os.path.isfile(f), 'Cannot find required file %s' % f)
        self.prefix = '%s/test' % self.subdir
        self.subassembly = '%s/test_subassembled_variants.txt' % self.subdir
        self.stats = '%s/test_stats.txt' % self.subdir
        self.allmutcounts = '%s/test_allmutcounts.txt' % self.subdir
        self.singlemutcounts = '%s/test_singlemutcounts.txt' % self.subdir
        for f in [self.subassembly, self.stats, self.allmutcounts, self.singlemutcounts]:
            if os.path.isfile(f):
                os.remove(f)


    def test_Subassemble_and_MatchBarcodes(self):
        """Runs ``dms_subassemble`` and ``dms_matchbarcodedsubamplicons``."""
        # run dms_subassemble
        cmds = ['dms_subassemble', self.prefix, self.refseq, self.r1, self.r2, self.alignspecs, '--barcodelength', str(self.barcodelength)]
        self.assertFalse(os.path.isfile(self.subassembly), '%s already exists' % self.subassembly)
        sys.stderr.write('\nRunning the following command:\n%s\n' % ' '.join(cmds))
        os.system(' '.join(cmds))
        self.assertTrue(os.path.isfile(self.subassembly), 'Failed to create file %s with command:\n%s\n' % (self.subassembly, ' '.join(cmds)))
        with open(self.correctsubassembly) as f:
            correct = dict([line.split() for line in f if not line.isspace()])
        with open(self.subassembly) as f:
            results = dict([line.split()[ : 2] for line in f if not line.isspace()])
        barcodes_correct = set(correct.keys())
        barcodes_results = set(results.keys())
        self.assertFalse(barcodes_correct - barcodes_results, "Failed to subassemble all the expected barcodes. Did not get the following expected ones: %s" % ', '.join(barcodes_correct - barcodes_results))
        self.assertFalse(barcodes_results - barcodes_correct, "Subassembled unexpected barcodes. Got the following unexpected ones: %s" % ', '.join(barcodes_results - barcodes_correct))
        for barcode in barcodes_correct:
            self.assertEqual(len(correct[barcode]), len(results[barcode]), "Different length sequences for barcode %s" % barcode)
            self.assertEqual(correct[barcode], results[barcode], "Difference in subassembled sequences for barcode %s\nCorrect: %s\nActual: %s\nDifferences: %s" % (barcode, correct[barcode], results[barcode], ', '.join(['%s%d%s' % (correct[barcode][i], i + 1, results[barcode][i]) for i in range(len(correct[barcode])) if correct[barcode][i] != results[barcode][i]])))

        # now run dms_matchsubassembledbarcodes
        cmds = ['dms_matchsubassembledbarcodes', self.prefix, self.subassembly, self.barcodes_r1, self.barcodes_r2, '--r1start 3', '--r2end 2']
        for f in [self.stats, self.singlemutcounts, self.allmutcounts]:
            self.assertFalse(os.path.isfile(f), '%s already exists' % f)
        sys.stderr.write('\nRunning the following command:\n%s\n' % ' '.join(cmds))
        os.system(' '.join(cmds))
        for f in [self.stats, self.singlemutcounts, self.allmutcounts]:
            self.assertTrue(os.path.isfile(f), 'Failed to create file %s with command:\n%s\n' % (f, ' '.join(cmds)))
        stats = dms_tools.file_io.ReadSummaryStats(self.stats)[1]
        correctstats = dms_tools.file_io.ReadSummaryStats(self.correctstats)[1]
        self.assertTrue(all([stats[key] == correctstats[key] for key in correctstats.keys()]), 'Mismatched stats in %s and %s' % (self.stats, self.correctstats))
        for (f, actualf) in [(self.singlemutcounts, self.correctsinglemutcounts), (self.allmutcounts, self.correctallmutcounts)]:
            counts = dms_tools.file_io.ReadDMSCounts(f, 'codon')
            actualcounts = dms_tools.file_io.ReadDMSCounts(actualf, 'codon')
            self.assertTrue(counts == actualcounts, "Mismatched counts in %s and %s" % (f, actualf))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
