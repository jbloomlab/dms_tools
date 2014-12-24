"""Converts input data files from Wu et al into dms_tools counts format."""


import dms_tools.file_io


def main():
    """Main body of script."""
    files_to_convert = [('IPT_NS1', 'input.txt'), ('CON_NS1', 'control.txt'), ('IFN_NS1', 'interferon.txt')]

    for (infile, outfile) in files_to_convert:
        print("Converting %s to %s" % (infile, outfile))
        counts = {}
        with open(infile) as f:
            for line in f:
                (wt, r, mut, n) = line.split()
                r = str(int(r) - 13) # renumber to start with one at start codon
                if int(r) > 693:
                    continue # coding sequence ends at nucleotide 693
                if r not in counts:
                    counts[r] = {'WT':wt, mut:int(n)}
                else:
                    assert counts[r]['WT'] == wt, "Mismatch at site %s: is wildtype %s or %s" (r, wt, counts[r]['WT'])
                    assert mut not in counts[r], "Duplicate count for %s at %s" % (mut, r)
                    counts[r][mut] = int(n)
        for r in list(counts.keys()):
            if all([n == 0 for n in counts[r].values() if isinstance(n, int)]):
                del counts[r] # delete sites where all counts are zero
        dms_tools.file_io.WriteDMSCounts(outfile, counts)
        

if __name__ == '__main__':
    main()
