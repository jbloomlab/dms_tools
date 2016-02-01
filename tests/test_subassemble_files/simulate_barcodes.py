"""Creates simulated test data for ``dms_matchsubassemblebarcodes``."""


import random
import os
import Bio.SeqIO
import dms_tools
import dms_tools.cutils
import dms_tools.file_io


def GetCodonMuts(seq, refseq):
    """Returns comma-separated codon mutations in *seq* relative to *refseq*.
    
    Returns *no_mutations* if there aren't any."""
    assert len(seq) == len(refseq)
    assert len(seq) % 3 == 0
    mutlist = []
    for i in range(len(seq) // 3):
        wt = refseq[3 * i : 3 * i + 3]
        mut = seq[3 * i : 3 * i + 3]
        if wt != mut:
            mutlist.append('%s%d%s' % (wt, i + 1, mut))
    if mutlist:
        return ','.join(mutlist)
    else:
        return 'no_mutations'


def MakeReadPair(barcode, readlength, r1start, r2end, makemismatch, makelowq):
    """Makes read pair of *readlength* containing *barcode*.

    *r1start* and *r2end* indicate start / end of read
    in R1 / R2.

    *makemismatch* is Boolean switch indicating if there
    should be a mismatch.

    *makelowq* is Boolean switch indicating if there should
    be a low quality site.

    The return tuple is *(r1, r2, q1, q2)*.
    """
    highq = chr(35 + 33)
    lowq = chr(10 + 33)
    r1 = ''.join([random.choice(dms_tools.nts) for i in range(r1start - 1)]) + barcode
    while len(r1) < readlength:
        r1 += random.choice(dms_tools.nts)
    q1 = highq * readlength
    if makemismatch:
        i = random.choice(range(len(barcode)))
        mut = random.choice([nt for nt in dms_tools.nts if nt != barcode[i]])
        barcode = list(barcode)
        barcode[i] = mut
        barcode = ''.join(barcode)
    r2 = ''.join([random.choice(dms_tools.nts) for i in range(r2end)]) + barcode
    while len(r2) < readlength:
        r2 += random.choice(dms_tools.nts)
    r2 = dms_tools.cutils.ReverseComplement(r2)
    q2 = [random.choice([lowq, highq]) for i in range(readlength)]
    if makelowq:
        i = random.choice(range(len(barcode)))
        q1 = list(q1)
        q1[r1start - 1 + i] = lowq
        q1 = ''.join(q1)
        q2[readlength - len(barcode) - r2end + i] = lowq
    q2 = ''.join(q2)
    assert len(r1) == len(r2) == len(q1) == len(q2)
    return (r1, r2, q1, q2)


def main():
    """Main body of script."""
    random.seed(1)
    r1name = '@M03100:54:000000000-AH22W:1:1101:15279:2347 1:N:0:1'
    r2name = '@M03100:54:000000000-AH22W:1:1101:15279:2347 2:N:0:1'
    r1file = 'barcodes_R1.fastq'
    r2file = 'barcodes_R2.fastq'
    statsfile = 'actual_stats.txt'
    allmutcountsfile = 'actual_allmutcounts.txt'
    singlemutcountsfile = 'actual_singlemutcounts.txt'
    readlength = 20
    r1start = 3
    r2end = 2
    stats = {
        'nretained':25,
        'nunrecognized':4,
        'nmismatches':3,
        'nlowquality':2,
        }
    refseqfile = 'refseq.fasta'
    subassembledvariantsfile = 'actual_barcoded_variants.txt'
    assert all(map(os.path.isfile, [refseqfile, subassembledvariantsfile])), "Can't find expected files. Run simulate_data.py first"
    refseq = str(list(Bio.SeqIO.parse('refseq.fasta', 'fasta'))[0].seq)
    with open(subassembledvariantsfile) as f:
        bc_dict = dict([line.strip().split() for line in f])
    assert all([len(refseq) == len(seq) for seq in bc_dict.values()])
    bc_muts = dict([(bc, GetCodonMuts(seq, refseq)) for (bc, seq) in bc_dict.items()])
    bc_length = len(bc_dict.keys()[0])
    assert readlength - r1start - r2end + 1 >= bc_length
    bc_list = []
    read_list = []
    for i in range(stats['nretained']):
        bc = random.choice(bc_dict.keys())
        bc_list.append(bc)
        read_list.append(MakeReadPair(bc, readlength, r1start, r2end, False, False))
    for i in range(stats['nunrecognized']):
        bc = ''.join([random.choice(dms_tools.nts) for i in range(bc_length)])
        while bc in bc_dict:
            bc = ''.join([random.choice(dms_tools.nts) for i in range(bc_length)])
        read_list.append(MakeReadPair(bc, readlength, r1start, r2end, False, False))
    for i in range(stats['nmismatches']):
        bc = random.choice(bc_dict.keys())
        read_list.append(MakeReadPair(bc, readlength, r1start, r2end, True, False))
    for i in range(stats['nlowquality']):
        bc = random.choice(bc_dict.keys())
        read_list.append(MakeReadPair(bc, readlength, r1start, r2end, False, True))
    random.shuffle(read_list)
    for bc in bc_list:
        if bc_muts[bc] == 'no_mutations':
            key = 'n0mut'
        else:
            key = 'n%dmut' % (1 + bc_muts[bc].count(','))
        if key in stats:
            stats[key] += 1
        else:
            stats[key] = 1
    with open(r1file, 'w') as f1, open(r2file, 'w') as f2:
        for (r1, r2, q1, q2) in read_list:
            f1.write('%s\n%s\n+\n%s\n' % (r1name, r1, q1))
            f2.write('%s\n%s\n+\n%s\n' % (r2name, r2, q2))
    with open(statsfile, 'w') as f:
        f.write('\n'.join(['%s = %d' % tup for tup in stats.items()]))
    sites = [r for r in range(len(refseq) // 3)]
    allmutcounts = {}
    singlemutcounts = {}
    for r in sites:
        allmutcounts[str(r + 1)] = dict([(codon, 0) for codon in dms_tools.codons])
        singlemutcounts[str(r + 1)] = dict([(codon, 0) for codon in dms_tools.codons])
        allmutcounts[str(r + 1)]['WT'] = singlemutcounts[str(r + 1)]['WT'] = refseq[3 * r : 3 * r + 3]
    for bc in bc_list:
        muts = {}
        if bc_muts[bc] != 'no_mutations':
            for mutstring in bc_muts[bc].split(','):
                (wt, r, mut) = (mutstring[ : 3], int(mutstring[3 : -3]), mutstring[-3 : ])
                muts[r - 1] = mut
        for r in sites:
            if r in muts:
                allmutcounts[str(r + 1)][muts[r]] += 1
                if len(muts) == 1:
                    singlemutcounts[str(r + 1)][muts[r]] += 1
            else:
                allmutcounts[str(r + 1)][allmutcounts[str(r + 1)]['WT']] += 1
                if len(muts) == 0:
                    singlemutcounts[str(r + 1)][allmutcounts[str(r + 1)]['WT']] += 1
    dms_tools.file_io.WriteDMSCounts(allmutcountsfile, allmutcounts)
    dms_tools.file_io.WriteDMSCounts(singlemutcountsfile, singlemutcounts)


if __name__ == '__main__':
    main()
