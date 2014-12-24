"""Converts raw data from Melnikov et al into dms_tools counts files.

Written by Jesse Bloom."""


import os
import re
import glob
import subprocess
import dms_tools
import dms_tools.plot
import dms_tools.simulate


def Melnikov_to_DMSCounts(infile, outfile):
    """Converts deep mutational scanning counts files to ``dms_tools`` format.

    *infile* is existing file in ``*.aacounts.txt`` format of Melnikov et al.

    *outfile* is created file in ``dms_tools`` deep mutational scanning counts
    format.
    """
    with open(infile) as f:
        lines = f.readlines()
    assert len(lines) == 23, "Did not find 23 lines in %s" % infile
    counts = {}
    sites = lines[1].split()[1 : ]
    wts = lines[2].split()[1 : ]
    assert len(wts) == len(sites)
    counts_by_aa = []
    for (aa, line) in zip(dms_tools.aminoacids_nostop, lines[3 : ]):
        assert aa == line.split()[0]
        counts_by_aa.append([int(n) for n in line.split()[1 : ]])
        assert len(counts_by_aa[-1]) == len(sites)
    for (i, r) in enumerate(sites):
        rlist = [('WT', wts[i])]
        for (iaa, aa) in enumerate(dms_tools.aminoacids_nostop):
            rlist.append((aa, counts_by_aa[iaa][i]))
        counts[r] = dict(rlist)
        dms_tools.file_io.WriteDMSCounts(outfile, counts)



def main():
    """Main body of script."""

    # specify directory and file names
    raw_dir = './raw_data/' # raw data from Melnikov et al
    processed_dir = './processed_data/' # processed data from Melnikov et al
    if not os.path.isdir(processed_dir):
        os.mkdir(processed_dir)
    # list files with processed filename, raw filename
    files_to_convert = \
        [('%s/input_%d_dms_counts.txt' % (processed_dir, i), '%s/KKA2_Bkg%d.aacounts.txt' % (raw_dir, i)) for i in [1, 2]] + \
        [('%s/kanamycin_%d_dms_counts.txt' % (processed_dir, i), '%s/KKA2_S%d_Kan11_L%d.aacounts.txt' % (raw_dir, i, i)) for i in [1, 2]]

    # convert Melnikov et al raw data into dms_tools count files
    for (processed, raw) in files_to_convert:
        print('Processing %s to %s' % (raw, processed))
        Melnikov_to_DMSCounts(raw, processed)


if __name__ == '__main__':
    main()

