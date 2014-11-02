"""
================================
``io`` module
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
import StringIO
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

    *f* is a readable file-like object that contains the counts.

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

    >>>
    """

if __name__ == '__main__':
    import doctest
    doctest.testmod()
