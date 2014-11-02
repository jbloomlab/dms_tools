"""
================================
``file_io`` module
================================

Module that reads and writes files.

Written by Jesse Bloom.

Functions in this module
---------------------------

* *ReadDMSCounts* : reads deep mutational scanning counts.

Function documentation
---------------------------

"""


import sys
import os
import cStringIO
import Bio.Alphabet.IUPAC

# lists of upper case DNA nucleotides and codons
aminoacids = Bio.Alphabet.IUPAC.IUPACProtein.letters
nts = Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters
nts = [nt.upper() for nt in nts]
codons = []
for nt1 in nts:
    for nt2 in nts:
        for nt3 in nts:
            codons.append('%s%s%s' % (nt1, nt2, nt3))


def ReadDMSCounts(f, chartype):
    """Reads deep mutational scanning counts.

    *f* is a readable file-like object that contains the counts, or 
    a string that gives the name of a file.

    *chartype* is the type of character. Valid values are the following strings:

        - *DNA* : DNA nucleotides, those listed in the *nts* variable
          defined in this module.

        - *codons* : DNA codons, those listed in the *codons* variable
          defined in this module.

    The counts are returned in the dictionary *counts*. The dictionary
    is keyed by **strings** giving the position (e.g. '1', '2', or '5A').
    For each position *r*, *strings[r]* is a dictionary keyed by the
    string *WT* and all characters (nucleotides or codons) as specified
    by *chartype*, in upper case. *counts[r]['WT']* is the wildtype identity
    at the site; *counts[r][x]* is the number of counts for character *x*.

    The specifications on the format of the contents of *f* are as follows:

        * Any lines beginning with *#* are treated as comments or headers, 
          and do not contain data.

        * Columns are space delimited. They can separated by one or more 
          spaces, or by one or more tabs.

        * Before the first data line (a line **not** beginning with *#*) 
          there must be a header line beginning with *# POSITION WT* and 
          then followed in arbitrary order by entries specifying valid 
          DNA nucleotides or valid DNA codons. The order of the entries 
          in this header line establish the order of the data in subsequent lines.

        * The case does matter, all nucleotides and codons **must be upper case**.

        * It is OK to intersperse columns that do not represent valid codons 
          or nucleotides as long as all of the codons or nucleotides are present. 
          This might be useful if you file is using *N* to denote ambiguous 
          nucleotides (which are ignored as not being definitive counts), or 
          is using lower-case letters (such as *a* to denote unsure calls).

        * For every data line, the first two entries will be *POSITION* and *WT*:

            - *POSITION* is the position number. These do not have to be 
              consecutive integers, so it is OK to skip sites (e.g. *1*, *2*, *5*) 
              or have negative site numbers (e.g. *-1*, *0*, *1*) or have sites 
              with letter suffixes (e.g. *1*, *2*, *2A*, *3*). All of these 
              options are useful when referring to numbering schemes that are 
              non-sequential, as is sometimes the case in gene numbering. Note 
              however that each position number must be unique.

            - *WT* must give a valid wildtype identity of the site as a 
               nucleotide or codon.

    Example of parsing nucleotide counts files with minimal required format and
    extra lines / columns.

    >>> f = cStringIO.StringIO()
    >>> f.write('# POSITION WT A T C G\\n')
    >>> f.write('1 A 310818 13 0 37\\n')
    >>> f.write('2 T 16 311782 37 4\\n')
    >>> f.write('2A G 27 27 11 312520\\n')
    >>> f.write('3 A 313647 13 4 31\\n')
    >>> f.seek(0)
    >>> f2 = cStringIO.StringIO()
    >>> f2.write('# extra comment line\\n')
    >>> f2.write('# POSITION WT N A T G C a\\n')
    >>> f2.write('1 A 1 310818\t13 37 0 15\\n')
    >>> f2.write('2 T 15 16 311782   4 37\\t 0\\n')
    >>> f2.write('# another comment line\\n')
    >>> f2.write('2A G 7 27 27 312520 11 1\\n')
    >>> f2.write('3 A 0 313647 13 31 4 2\\n')
    >>> f2.seek(0)
    >>> d = ReadDMSCounts(f, 'DNA')
    >>> d2 = ReadDMSCounts(f2, 'DNA')
    >>> d == d2
    True
    >>> d['1']['WT'] == 'A'
    True
    >>> d['2A']['WT'] == 'G'
    True
    >>> d['2']['T'] == 311782
    True
    >>> d['2']['C'] == 37
    True
    >>> 1 in d
    False
    >>> 'N' in d2['1']
    False
    >>> 'a' in d2['2']
    False
    >>> f.seek(0)
    >>> ReadDMSCounts(f, 'codon')
    Traceback (most recent call last):
        ...
    ValueError: file does not have a valid header for codon

    Example for codons

    >>> nts = ['A', 'C', 'G', 'T']
    >>> codons = []
    >>> for nt1 in nts:
    ...     for nt2 in nts:
    ...         for nt3 in nts:
    ...             codons.append('%s%s%s' % (nt1, nt2, nt3))
    >>> f = cStringIO.StringIO()
    >>> f.write('# POSITION WT %s\\n' % ' '.join(codons))
    >>> f.write('1 AAA %s\\n' % ' '.join([str(i) for i in range(len(codons))]))
    >>> f.seek(0)
    >>> d = ReadDMSCounts(f, 'codon')
    >>> d['1']['AAC'] == 1
    True
    >>> d['1']['WT'] == 'AAA'
    True
    >>> d['1']['WT'] == 'AAa'
    False
    >>> codons[5] = codons[5].lower()
    >>> f2 = cStringIO.StringIO()
    >>> f2.write('# POSITION WT %s\\n' % ' '.join(codons))
    >>> f2.write('1 AAA %s\\n' % ' '.join([str(i) for i in range(len(codons))]))
    >>> f2.seek(0)
    >>> d = ReadDMSCounts(f2, 'codon')
    Traceback (most recent call last):
        ...
    ValueError: file does not have a valid header for codon

    """
    if chartype.upper() == 'DNA':
        characters = nts
    elif chartype.upper() == 'CODON':
        characters = codons
    else:
        raise ValueError("Invalid chartype of %s" % chartype)
    characterindices = dict([(char, False) for char in characters])
    openedfile = False
    if isinstance(f, str):
        f = open(f)
        openedfile = True
    foundheader = False
    counts = {}
    for line in f:
        if line.isspace():
            continue
        elif line[0] == '#':
            if foundheader:
                continue
            entries = line[1 : ].split()
            if len(entries) >= 2 and entries[0].strip().upper() == 'POSITION' and entries[1].strip().upper() == 'WT':
                # appears to be header
                for i in range(2, len(entries)):
                    x = entries[i].strip()
                    if x in characterindices:
                        if characterindices[x] == False:
                            characterindices[x] = i
                        else:
                            raise ValueError("Duplicate entry for character %s in this line:\n%s" % (x, line))
                if not all([i != False for (char, i) in characterindices.iteritems()]):
                    raise ValueError('file does not have a valid header for %s' % chartype)
                foundheader = True
        else:
            if not foundheader:
                raise ValueError('Encountered data line before finding a valid header:\n%s' % line)
            # this is a data line
            entries = line.split()
            if len(entries) < 2 + len(characters):
                raise ValueError('Line has insufficient number of entries:\n%s' % line)
            r = entries[0].strip()
            wt = entries[1].strip()
            if r in counts:
                raise ValueError("Duplicate entry for position %s in line:\n%s" % (r, line))
            if wt not in characters:
                raise ValueError("Invalid wildtype character of %s is not a valid %s:\n%s" % (wt, chartype, line))
            counts[r] = {'WT':wt}
            for char in characters:
                i = characterindices[char]
                if len(entries) <= i:
                    raise ValueError("Not enough entries for counts for %s expected at index %d:\n%s" % (char, i, line))
                counts[r][char] = int(entries[i])
    if openedfile:
        f.close()
    return counts


if __name__ == '__main__':
    import doctest
    doctest.testmod()
