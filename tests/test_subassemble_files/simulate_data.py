"""Running this script creates simulated test data for ``dms_subassemble``."""


import random
import os


def NTMutate(seq, mutrate):
    """Returns copy of *seq* with each nucleotide mutated at prob *mutrate*."""
    nts = 'ATCG'
    newseq = []
    for nt in seq:
        if random.random() < mutrate:
            choices = [x for x in nts if x != nt]
            mut = random.choice(choices)
            newseq.append(mut)
        else:
            newseq.append(nt)
    return ''.join(newseq)


def main():
    """Main body of script."""
    random.seed(1)
    nts = 'ATCG'
    ncodons = 40
    alignspecs = [(1, 3), (61, 1)]
    readlength = 60
    refseq = ''.join([random.choice(nts) for i in range(ncodons * 3)])
    with open('refseq.fasta', 'w') as f:
        f.write(">refseq\n%s" % refseq)
    barcodelength = 6
    nvariants = 8
    mutrate = 1.0 / ncodons
    errorrate = 1.0 / ncodons / 3.0
    nreadspersubamplicon = 40
    qscore = '9'
    r1header = '@M03100:54:000000000-AH22W:1:1101:15279:2347 1:N:0:1'
    r2header = '@M03100:54:000000000-AH22W:1:1101:15279:2347 2:N:0:1'
    r1file = open('R1.fastq', 'w')
    r2file = open('R2.fastq', 'w')
    barcodefile = open('actual_barcoded_variants.txt', 'w')
    for ivariant in range(nvariants):
        barcode = ''.join([random.choice(nts) for i in range(barcodelength)])
        variant = []
        nmuts = 0
        for icodon in range(ncodons):
            if random.random() < mutrate:
                variant.append(''.join([random.choice(nts) for i in range(3)]))
                nmuts += 1
            else:
                variant.append(refseq[icodon * 3 : icodon * 3 + 3])
        variant = ''.join(variant)
        barcodefile.write('%s %s\n' % (barcode, variant))
#        print "barcode %s is a variant with %d mutated codons out of %d codons total" % (barcode, nmuts, ncodons)
        assert len(variant) == len(refseq)
        for (refseqstart, r2start) in alignspecs:
            for iread in range(nreadspersubamplicon):
                r1 = barcode + refseq[-2 : ]
                r1 = NTMutate(r1, errorrate)
                r2 = ''.join([random.choice(nts) for i in range(r2start - 1)]) + NTMutate(variant[refseqstart - 1 : readlength + refseqstart - 1], errorrate) + ''.join([random.choice(nts) for i in range(7)])
                r1file.write('%s\n%s\n+\n%s\n' % (r1header, r1, qscore * len(r1)))
                r2file.write('%s\n%s\n+\n%s\n' % (r2header, r2, qscore * len(r2)))
    barcodefile.close()
    r1file.close()
    r2file.close()
#    os.system('dms_subassemble temp refseq.fasta R1.fastq R2.fastq 1,3 61,1 --barcodelength 6 --maxmuts 3')



if __name__ == '__main__':
    main()
