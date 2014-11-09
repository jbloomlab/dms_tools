"""
================================
``file_io`` module
================================

Module that reads and writes files.

Written by Jesse Bloom.

Functions in this module
---------------------------

* *Versions* : Returns string giving versions of utilized packages.

* *ReadDMSCounts* : reads deep mutational scanning counts.

* *ReadMultipleDMSCountsFiles* : reads multiple counts files for the same sequence.

* *WritePreferences* : writes site-specific preferences.

* *ReadPreferences* : reads the preferences file written by *WritePreferences*.

Function documentation
---------------------------

"""


import sys
import os
import re
import time
import platform
import importlib
import cStringIO
import dms_tools
import dms_tools.utils



def Versions():
    """Returns a string with version information.

    You would call this function if you want a string giving detailed information
    on the version of ``dms_tools`` and the associated packages that it uses.
    This string might be useful if you want to print this information at the beginning
    of a script.
    """
    s = [\
            'Version information for dms_tools and associated packages.',
            '\tTime and date: %s' % time.asctime(),
            '\tPlatform: %s' % platform.platform(),
            '\tPython version: %s' % sys.version.replace('\n', ' '),
            '\tdms_tools version: %s' % dms_tools.__version__,
            ]
    for modname in ['Bio', 'numpy', 'matplotlib', 'cython', 'pystan']:
        try:
            v = importlib.import_module(modname).__version__
            s.append('\t%s version: %s' % (modname, v))
        except ImportError:
            s.append('\t%s cannot be imported into Python' % modname)
    return '\n'.join(s)


def ReadMultipleDMSCountsFiles(filelist, chartype):
    """Reads multiple deep mutational scanning counts for the same sequence.

    *filelist* is a list of names of files that can be read by *ReadDMSCounts*
    using *chartype*. 

    This function reads the counts for all files in *filelist*. It confirms
    that each file specifies counts for the same set of sites, and that
    the wildtype sequence specified in each file is the same. It then returns
    the following tuple: *(sites, wts, counts)*. In this tuple:

        - *sites* is a list of strings giving the sites, sorted so that the
          numeric portion of the string is interpreted as a number.

        - *wts* is a dictionary, with *wts[r]* giving the wildtype character
          for site *r* for each *r* in *sites*.

        - *counts* is a dictionary keyed by each file in *filelist*, with the value
          for each file name being the counts dictonary for that file in the format
          returned by *ReadDMSCounts*.
    """
    assert len(filelist) >= 1, "filelist must specify at least one file"
    counts = {}
    sites = None
    for fname in filelist:
        assert os.path.isfile(fname), "Cannot find DMS counts file %s" % fname
        counts[fname] = ReadDMSCounts(fname, chartype)
        if sites == None:
            sites = counts[fname].keys()
            assert sites, "%s does not specify information for any sites" % fname
            dms_tools.utils.NaturalSort(sites)
            wts = dict([(r, counts[fname][r]['WT']) for r in sites])
        else:
            if set(sites) != set(counts[fname].keys()):
                raise ValueError("Not all of the counts files have the same set of sites")
            for r in sites:
                if wts[r] != counts[fname][r]['WT']:
                    raise ValueError("The counts files are not consistent on the wildtype identity for site %s" % r)
    return (sites, wts, counts)



def ReadDMSCounts(f, chartype):
    """Reads deep mutational scanning counts.

    *f* is a readable file-like object that contains the counts, or 
    a string that gives the name of a file.

    *chartype* is the type of character. Valid values are the following strings:

        - *DNA* : DNA nucleotides, those listed in the *nts* variable
          defined in this module.

        - *codon* : DNA codons, those listed in the *codons* variable
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
        characters = dms_tools.nts
    elif chartype.upper() == 'CODON':
        characters = dms_tools.codons
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


def WritePreferences(f, sites, wts, pi_means, pi_95credint):
    """Writes site-specific preferences to a file.

    *f* is a writeable file-like object to which we write the counts,
     or a string giving the name of a file that we create.

    *sites* is a list of all site numbers (as strings) for which
     we have preferences.
     
    *wts* is a dictionary with *wts[r]* giving the wildtype character
     at site *r* for all *r* in *sites*.

    *pi_means* is a dictionary keyed by all sites (as strings). For each 
     site *r*, *pi_means[r]* is a dictionary where *pi_means[r][x]*
     is the preference of *r* for character *x*, where the
     characters are nucleotide, codons, or amino acids. 

    *pi_95credint* can either be *None* or a dictionary keyed by all
     sites that also key *pi_means, and has as values 2-tuples giving
     the upper and lower bounds of the 95% credible interval.

    The written file has the following format::

        # POSITION WT SITE_ENTROPY PI_A PI_C PI_G PI_T PI_A_95 PI_C_95 PI_G_95 PI_T_95
        1 G 1.61 0.21 0.09 0.58 0.12 0.17,0.28 0.04,0.12 0.51,0.69 0.09,0.15

    In this format, the first header line begins with ``#`` and indicates the
    characters as the suffix to ``PI_``. The 95% credible intervals (if given)
    are indicated by columns beginning with ``PI_`` and ending with ``_95``.
    The ``SITE_ENTROPY`` column gives the site entropy in bits (log base 2)
    as :math:`h_r = \sum_x \pi_{r,x} \\times \log_2 \pi_{r,x}`

    >>> sites = ['1']
    >>> pi_means = {'1':{'A':0.2, 'C':0.3, 'G':0.4, 'T':0.1}}
    >>> wts = {'1':'G'}
    >>> pi_95credint = {'1':{'A':(0.1, 0.3), 'C':(0.25, 0.35), 'G':(0.3, 0.5), 'T':(0.05, 0.15)}}
    >>> f = cStringIO.StringIO()
    >>> WritePreferences(f, sites, wts, pi_means, pi_95credint)
    >>> f.seek(0)
    >>> lines = f.readlines()
    >>> lines[0] == '# POSITION WT SITE_ENTROPY PI_A PI_C PI_G PI_T PI_A_95 PI_C_95 PI_G_95 PI_T_95\\n'
    True
    >>> f2 = cStringIO.StringIO()
    >>> WritePreferences(f2, sites, wts, pi_means, None)
    >>> f2.seek(0)
    >>> lines = f2.readlines()
    >>> lines[0] == '# POSITION WT SITE_ENTROPY PI_A PI_C PI_G PI_T\\n'
    True
    >>> f2.seek(0)
    >>> (sites2, wts2, pi_means2, pi_95credint2, h) = ReadPreferences(f2)
    >>> sites == sites2
    True
    >>> wts == wts2
    True
    >>> pi_95credint2 == None
    True
    >>> r = sites[0]
    >>> all([abs(pi_means[r][x] - pi_means2[r][x]) < 1.0e-6 for x in pi_means[r].keys()])
    True
    >>> print "%.3f" % h[r]
    1.846
    """
    openedfile = False
    if isinstance(f, str):
        f = open(f)
        openedfile = True
    assert sites and len(sites) >= 1 and not isinstance(sites, str), "sites must be a list with at least one entry"
    characters = pi_means[sites[0]].keys()
    characters.sort()
    f.write('# POSITION WT SITE_ENTROPY %s' % ' '.join(['PI_%s' % x for x in characters]))
    if pi_95credint:
        f.write(' %s\n' % ' '.join(['PI_%s_95' % x for x in characters]))
    else:
        f.write('\n')
    for r in sites:
        h = dms_tools.utils.SiteEntropy([pi_means[r][x] for x in characters])
        f.write('%s %s %g ' % (r, wts[r], h))
        f.write('%s' % ' '.join(['%g' % pi_means[r][x] for x in characters]))
        if pi_95credint:
            f.write(' %s\n' % ' '.join(['%g,%g' % pi_95credint[r][x] for x in characters]))
        else:
            f.write('\n')
    if openedfile:
        f.close()


def ReadPreferences(f):
    """Reads the site-specific preferences written by *WritePreferences*.

    *f* is the name of an existing file or a readable file-like object.

    The return value is the tuple: *(sites, wts, pi_means, pi_95credint, h)*
    where *sites*, *wts*, *pi_means*, and *pi_95credint* will all
    have the same values used to write the file with *WritePreferences*,
    and *h* is a dictionary with *h[r]* giving the site entropy (log base
    2) for each *r* in *sites*.

    See docstring of *WritePreferences* for example usage.
    """
    charmatch = re.compile('^PI_([A-z\*\-]+)$')
    if isinstance(f, str):
        open(f)
        lines = f.readlines()
        f.close()
    else:
        lines = f.readlines()
    characters = []
    sites = []
    wts = {}
    pi_means = {}
    pi_95credint = {}
    h = {}
    for line in lines:
        if line.isspace():
            continue
        elif line[0] == '#' and not characters:
            entries = line[1 : ].strip().split()
            if len(entries) < 4:
                raise ValueError("Insufficient entries in header:\n%s" % line)
            if not (entries[0] == 'POSITION' and entries[1] == 'WT' and entries[2] == 'SITE_ENTROPY'):
                raise ValueError("Not the correct first three header columns:\n%s" % line)
            i = 3
            while i < len(entries) and charmatch.search(entries[i]):
                characters.append(charmatch.search(entries[i]).group(1))
                i += 1
            if i  == len(entries):
                pi_95credint = None
                linelength = len(characters) + 3
            else:
                if not len(entries) - i == len(characters):
                    raise ValueError("Header line does not have valid credible interval format:\n%s" % line)
                if not all([entries[i + j] == 'PI_%s_95' % characters[j] for j in range(len(characters))]):
                    raise ValueError("mean and credible interval character mismatch in header:\n%s" % line)
                linelength = 2 * len(characters) + 3
        elif line[0] == '#':
            continue
        elif not characters:
            raise ValueError("Found data lines before encountering a valid header")
        else:
            entries = line.strip().split()
            if len(entries) != linelength:
                raise ValueError("Line does not have expected %d entries:\n%s" % (linelength, line))
            r = entries[0]
            assert r not in sites, "Duplicate site of %s" % r
            sites.append(r)
            wts[r] = entries[1]
            assert entries[1] in characters, "Character %s is not one of the valid ones in header. Valid possibilities: %s" % (entries[1], ', '.join(characters))
            h[r] = float(entries[2])
            pi_means[r] = dict([(x, float(entries[3 + i])) for (i, x) in enumerate(characters)])
            if pi_95credint != None:
                pi_95credint[r] = dict([(x, (float(entries[3 + len(characters) + i].split()[0], float(entries[3 + len(characters) + i].split()[1])))) for (i, x) in enumerate(characters)])
    return (sites, wts, pi_means, pi_95credint, h)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
