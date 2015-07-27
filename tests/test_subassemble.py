"""Tests ``dms_subassemble``.

Written by Jesse Bloom."""


import sys
import os
import unittest



class TestSubasssemble(unittest.TestCase):
    """Runs ``dms_subassemble`` on test data.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.subdir = './test_subassemble_files/'
        self.refseq = '%s/refseq.fasta' % self.subdir
        self.r1 = '%s/R1.fastq' % self.subdir
        self.r2 = '%s/R2.fastq' % self.subdir
        self.alignspecs = '1,3 61,1'
        self.barcodelength = 6
        self.correctsubassembly = '%s/actual_barcoded_variants.txt' % self.subdir
        for f in [self.refseq, self.r1, self.r2, self.correctsubassembly]:
            self.assertTrue(os.path.isfile(f), 'Cannot find required file %s' % f)
        self.prefix = '%s/test' % self.subdir
        self.results = '%s/test_subassembled_variants.txt' % self.subdir
        if os.path.isfile(self.results):
            os.remove(self.results)


    def test_Subassemble(self):
        """Runs ``dms_subassemble``."""
        cmds = ['dms_subassemble', self.prefix, self.refseq, self.r1, self.r2, self.alignspecs, '--barcodelength', str(self.barcodelength)]
        self.assertFalse(os.path.isfile(self.results), '%s already exists' % self.results)
        sys.stderr.write('\nRunning the following command:\n%s\n' % ' '.join(cmds))
        os.system(' '.join(cmds))
        self.assertTrue(os.path.isfile(self.results), 'Failed to create file %s with command:\n%s\n' % (self.results, ' '.join(cmds)))
        with open(self.correctsubassembly) as f:
            correct = dict([line.split() for line in f if not line.isspace()])
        with open(self.results) as f:
            results = dict([line.split()[ : 2] for line in f if not line.isspace()])
        barcodes_correct = set(correct.keys())
        barcodes_results = set(results.keys())
        self.assertFalse(barcodes_correct - barcodes_results, "Failed to subassemble all the expected barcodes. Did not get the following expected ones: %s" % ', '.join(barcodes_correct - barcodes_results))
        self.assertFalse(barcodes_results - barcodes_correct, "Subassembled unexpected barcodes. Got the following unexpected ones: %s" % ', '.join(barcodes_results - barcodes_correct))
        for barcode in barcodes_correct:
            self.assertEqual(len(correct[barcode]), len(results[barcode]), "Different length sequences for barcode %s" % barcode)
            self.assertEqual(correct[barcode], results[barcode], "Difference in subassembled sequences for barcode %s\nCorrect: %s\nActual: %s\nDifferences: %s" % (barcode, correct[barcode], results[barcode], ', '.join(['%s%d%s' % (correct[barcode][i], i + 1, results[barcode][i]) for i in range(len(correct[barcode])) if correct[barcode][i] != results[barcode][i]])))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
