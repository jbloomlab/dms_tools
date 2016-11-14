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

* *WriteDMSCounts* : writes deep mutational scanning counts.

* *ReadMultipleDMSCountsFiles* : reads multiple counts files for the same sequence.

* *WritePreferences* : writes site-specific preferences.

* *ReadPreferences* : reads the preferences file written by *WritePreferences*.

* *WriteDiffPrefs* : writes site-specific differential preferences.

* *ReadDiffPrefs* : reads differential preferences file written by *WriteDiffPrefs*.

* *ReadMultiPrefOrDiffPref* : reads multiple preferences or differential preferences files.

* *IteratePairedFASTQ* : iterate over paired-read FASTQ files.

* *ReadSummaryStats* : reads alignment summary statistics.

* *ReadSubassembledVariants* : reads subassembled variants

Function documentation
---------------------------

"""


import sys
import os
import re
import time
import tempfile
import platform
import importlib
import math
import cStringIO
import subprocess
import gzip
import dms_tools
import dms_tools.utils
import pandas as pd


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
    for modname in ['Bio', 'numpy', 'scipy', 'matplotlib', 'cython', 'pystan', 'weblogolib', 'PyPDF2', 'pandas']:
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


def WriteDMSCounts(f, counts):
    """Writes deep mutational scanning counts file.

    This function is the inverse of *ReadDMSCounts*. 

    *f* is a writable file-like object, or a string that gives the name of
    a file that is created.

    *counts* is a dictionary of the same format returned by *ReadDMSCounts*.
    These counts are written to *f*.

    There is no *chartype* argument -- that is automatically determined from
    *counts*.

    >>> counts = {'1':{'A':2, 'C':3, 'G':4, 'T':1, 'WT':'A'}, '2':{'A':1, 'C':1, 'G':3, 'T':4, 'WT':'C'}}
    >>> counts['1']['COUNTS'] = 10
    >>> counts['2']['COUNTS'] = 9
    >>> counts['1']['F_A'] = counts['1']['A'] / float(counts['1']['COUNTS'])
    >>> counts['1']['F_C'] = counts['1']['C'] / float(counts['1']['COUNTS'])
    >>> counts['1']['F_G'] = counts['1']['G'] / float(counts['1']['COUNTS'])
    >>> counts['1']['F_T'] = counts['1']['T'] / float(counts['1']['COUNTS'])
    >>> counts['2']['F_A'] = counts['2']['A'] / float(counts['2']['COUNTS'])
    >>> counts['2']['F_C'] = counts['2']['C'] / float(counts['2']['COUNTS'])
    >>> counts['2']['F_G'] = counts['2']['G'] / float(counts['2']['COUNTS'])
    >>> counts['2']['F_T'] = counts['2']['T'] / float(counts['2']['COUNTS'])
    >>> f = cStringIO.StringIO()
    >>> WriteDMSCounts(f, counts)
    >>> f.seek(0)
    >>> readcounts = ReadDMSCounts(f, 'DNA')
    >>> counts == readcounts
    True

    """
    # if f is string, open this file
    if isinstance(f, str):
        fname = f
        if os.path.isfile(fname):
            os.remove(fname) # remove existing file
        fileopened = True
        f = open(f, 'w')
    else:
        fileopened = False

    sites = list(counts.keys())
    dms_tools.utils.NaturalSort(sites)
    assert sites, "No sites defined in counts"
    characters = set([x for x in counts[sites[0]].keys() if x != 'COUNTS' and x[ : 2] != 'F_'])
    assert 'WT' in characters, 'No wildtype specified by "WT" key'
    characters.remove('WT')
    for charset in [dms_tools.nts, dms_tools.codons, dms_tools.aminoacids_withstop, dms_tools.aminoacids_nostop]:
        if characters == set(charset):
            characters = charset
            break
    else:
        raise ValueError("Invalid character set:\n%s" % str(characters))
    
    try:
        f.write('# POSITION WT %s\n' % ' '.join(characters))
        for r in sites:
            assert set(characters + ['WT']).issubset(set(counts[r].keys())), "Not all the needed characters (or missing 'WT')"
            wt = counts[r]['WT']
            assert wt in characters, "Invalid wildtype of %s" % wt
            f.write('%s %s %s\n' % (r, wt, ' '.join(['%d' % counts[r][x] for x in characters])))
        if fileopened:
            f.close()
    except:
        if fileopened:
            try:
                f.close()
            except:
                pass
            if os.path.isfile(fname):
                os.remove(fname)
        raise


def ReadDMSCounts(f, chartype, translate_codon_to_aa=False, return_as_df=False):
    """Reads deep mutational scanning counts.

    *f* is a readable file-like object that contains the counts, or 
    a string that gives the name of a file.

    *chartype* is the type of character. Valid values are the following strings:

        - *DNA* : DNA nucleotides, those listed in the *nts* variable
          defined in this module.

        - *codon* : DNA codons, those listed in the *codons* variable
          defined in this module.

        - *aminoacids_nostop* : amino acids, not including stop codons.

        - *aminoacids_withstop* : amino acids, including stop codons (``*``).

    The counts are returned in the dictionary *counts*. The dictionary
    is keyed by **strings** giving the position (e.g. '1', '2', or '5A').
    For each position *r*, *counts[r]* is a dictionary keyed by the
    string *WT* and characters (nucleotides, codons, or amino-acids) specified
    by *chartype*, in upper case. *counts[r]['WT']* is the wildtype identity
    at the site; *counts[r][x]* is the number of counts for character *x*.
    In addition, *counts[r]['F_{0}'.format(x)]*
    is the frequency of counts of *x*, and *counts[r]['COUNTS']* is
    the total number of counts at site *r*. If there are no counts at a site,
    the frequency of counts for each character *x* at the site is set to 
    float(NaN).

    Alternatively, if *return_as_df* is set as *True*, the counts are returned
    as a pandas dataframe. 

    If *chartype* is *codon* and the option *translate_codon_to_aa* is set as
    True, then the returned counts will be for amino-acids as translated
    by the codon counts.

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
    elif chartype.upper() == 'AMINOACIDS_NOSTOP':
        characters = dms_tools.aminoacids_nostop
    elif chartype.upper() == 'AMINOACIDS_WITHSTOP':
        characters = dms_tools.aminoacids_withstop
    else:
        raise ValueError("Invalid chartype of %s" % chartype)
    if translate_codon_to_aa and chartype.upper() != 'CODON':
        raise ValueError("Can't use translate codon to aa option if chartype is not codon.")
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
            if not (wt in characters):
                raise ValueError("Invalid wildtype character of %s is not a valid %s:\n%s" % (wt, chartype, line))
            counts[r] = {'WT':wt}
            for char in characters:
                i = characterindices[char]
                if len(entries) <= i:
                    raise ValueError("Not enough entries for counts for %s expected at index %d:\n%s" % (char, i, line))
                try:
                    counts[r][char] = int(entries[i])
                except ValueError:
                    # number could be in scientific notation, handle this properly
                    if 'e' or 'E' in s:
                        n = float(entries[i])
                        if n % 1:
                            raise ValueError("Scientific notation number of %s is not an integer" % entries[1])
                        else:
                            counts[r][char] = int(n)
                    else:
                        raise
    if openedfile:
        f.close()

    # add a new key, 'COUNTS', to each site's counts dict. This is the total number of counts at the site.
    for site in counts.keys():
        total_counts = 0
        for c in characters:
            total_counts += counts[site][c]
        counts[site]['COUNTS'] = total_counts

    if translate_codon_to_aa:
        new_counts_dict = {}
        for site in counts.keys():
            siteaacountsdict = dms_tools.utils.SumCodonToAA(counts[site])
            siteaacountsdict['WT'] = dms_tools.codon_to_aa[counts[site]['WT']]
            siteaacountsdict['COUNTS'] = counts[site]['COUNTS']
            new_counts_dict[site] = siteaacountsdict
        counts =  new_counts_dict
        characters = dms_tools.aminoacids_withstop

    # add keys for the frequency of each character at each site, F_x:
    for site in counts.keys():
        for c in characters:
            try:
                counts[site]['F_%s' % c] = float(counts[site][c])/counts[site]['COUNTS']
            except ZeroDivisionError:
                counts[site]['F_%s' % c] = float('NaN')

    if return_as_df:
        custom_index = ['COUNTS', 'WT']
        [custom_index.append(c) for c in characters]
        [custom_index.append('F_%s' % c) for c in characters]
        series_dict = {}
        
        sites = counts.keys()
        dms_tools.utils.NaturalSort(sites)
        for site in sites:
            series_dict[site] = pd.Series(counts[site], index=custom_index)
        counts = pd.DataFrame(series_dict, columns=sites)

    return counts


def WritePreferences(f, sites, wts, pi_means, pi_95credint):
    """Writes site-specific preferences to a file.

    *f* is a writeable file-like object to which we write the preferences,
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
     sites that also key *pi_means*, and has as values 2-tuples giving
     the upper and lower bounds of the 95% credible interval. (If it is
     a dictionary with any value of *None*, that is the same as setting the
     overall variable to *None*).

    The written file has the following format::

        # POSITION WT SITE_ENTROPY PI_A PI_C PI_G PI_T PI_A_95 PI_C_95 PI_G_95 PI_T_95
        1 G 1.61 0.21 0.09 0.58 0.12 0.17,0.28 0.04,0.12 0.51,0.69 0.09,0.15

    In this format, the first header line begins with ``#`` and indicates the
    characters as the suffix to ``PI_``. The 95% credible intervals (if given)
    are indicated by columns beginning with ``PI_`` and ending with ``_95``.
    The ``SITE_ENTROPY`` column gives the site entropy in bits (log base 2)
    as :math:`h_r = \sum_x \pi_{r,x} \\times \log_2 \pi_{r,x}`.

    The ``POSITION`` header can also be ``SITE``.

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
    >>> f.seek(0)
    >>> tup = ReadPreferences(f)
    >>> tup[0] == sites
    True
    >>> tup[3] == pi_95credint
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
    if pi_95credint and None in pi_95credint.values():
        pi_95credint = None
    openedfile = False
    if isinstance(f, str):
        f = open(f, 'w')
        openedfile = True
    assert sites and len(sites) >= 1 and not isinstance(sites, str), "sites must be a list with at least one entry"
    characters = pi_means[sites[0]].keys()
    characters.sort()
    if '*' in characters:
        i = characters.index('*')
        characters = characters[ : i] + characters[i + 1 : ] + ['*']
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
        f = open(f)
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
            if not (entries[0] in ['POSITION', 'SITE'] and entries[1][ : 2] == 'WT' and entries[2] == 'SITE_ENTROPY'):
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
            assert entries[1] in characters or entries[1] == '?', "Character %s is not one of the valid ones in header. Valid possibilities: %s" % (entries[1], ', '.join(characters))
            h[r] = float(entries[2])
            pi_means[r] = dict([(x, float(entries[3 + i])) for (i, x) in enumerate(characters)])
            if pi_95credint != None:
                pi_95credint[r] = dict([(x, (float(entries[3 + len(characters) + i].split(',')[0]), float(entries[3 + len(characters) + i].split(',')[1]))) for (i, x) in enumerate(characters)])
    return (sites, wts, pi_means, pi_95credint, h)



def WriteDiffPrefs(f, sites, wts, deltapi_means, pr_deltapi_lt0, pr_deltapi_gt0):
    """Writes differential preferences to a file.

    *f* is a writeable file-like object to which we write the differential
     preferences, or a string giving the name of a file that we create.

    *sites* is a list of all site numbers (as strings) for which
     we have preferences.
     
    *wts* is a dictionary with *wts[r]* giving the wildtype character
     at site *r* for all *r* in *sites*.

    *deltapi_means* is a dictionary keyed by all sites (as strings). For each 
     site *r*, *deltapi_means[r]* is a dictionary where *deltapi_means[r][x]*
     is the preference of *r* for character *x*, where the
     characters are nucleotide, codons, or amino acids. 

    *pr_deltapi_gt0* and *pr_deltapi_lt0* can either both be *None*,
     or can both be dictionaries keyed by all sites in *deltapi_means*
     that have as their values the posterior probability that particular
     differential preference is greater than or less than zero, 
     respectively.

    The written file has the following format::

        # POSITION WT RMS_dPI dPI_A dPI_C dPI_G dPI_T Pr_dPI_A_lt0 Pr_dPI_A_gt0 Pr_dPI_C_lt0 Pr_dPI_C_gt0 Pr_dPI_G_lt0 Pr_dPI_G_gt0 Pr_dPI_T_lt0 Pr_dPI_T_gt0
        1 G 0.241 0.02 -0.4 0.2 0.18 0.41 0.59 0.90 0.10 0.2 0.8 0.3 0.7

    In this format, the first header line begins with ``#`` and indicates the
    characters as the suffix to ``dPI_``. The posterior probabilities
    are indicated by columns of the form ``Pr_dPI_A_gt0`` for greater than zero,
    and ``Pr_dPI_A_lt0`` for less than zero.
    The ``RMS_dPI`` column gives the root-mean-square differential preference
    as :math:`\\textrm{RMS}_{\Delta\pi_r} = \sqrt{\\frac{1}{\mathcal{N}_x}\sum_x \left(\Delta\pi_{r,x}\\right)^2}`

    >>> sites = ['1']
    >>> deltapi_means = {'1':{'A':-0.2, 'C':-0.3, 'G':0.4, 'T':0.1}}
    >>> wts = {'1':'G'}
    >>> f = cStringIO.StringIO()
    >>> WriteDiffPrefs(f, sites, wts, deltapi_means, None, None)
    >>> f.seek(0)
    >>> lines = f.readlines()
    >>> lines[0] == '# POSITION WT RMS_dPI dPI_A dPI_C dPI_G dPI_T\\n'
    True
    >>> entries = lines[1].split()
    >>> entries[0] == str(sites[0])
    True
    >>> entries[1] == wts[sites[0]]
    True
    >>> 1.0e-6 > abs(float(entries[2]) - math.sqrt(sum([x**2 for x in deltapi_means[sites[0]].values()]) / float(len(deltapi_means[sites[0]].values()))))
    True
    >>> 1.0e-6 > abs(float(entries[3]) - deltapi_means[sites[0]]['A'])
    True
    >>> f2 = cStringIO.StringIO()
    >>> pr_deltapi_gt0 = {'1':{'A':0.95, 'C':0.5, 'G':0.7, 'T':0.5}}
    >>> pr_deltapi_lt0 = {'1':{'A':0.05, 'C':0.5, 'G':0.3, 'T':0.5}}
    >>> WriteDiffPrefs(f2, sites, wts, deltapi_means, pr_deltapi_lt0, pr_deltapi_gt0)
    >>> f2.seek(0)
    >>> lines = f2.readlines()
    >>> lines[0] == '# POSITION WT RMS_dPI dPI_A dPI_C dPI_G dPI_T Pr_dPI_A_lt0 Pr_dPI_A_gt0 Pr_dPI_C_lt0 Pr_dPI_C_gt0 Pr_dPI_G_lt0 Pr_dPI_G_gt0 Pr_dPI_T_lt0 Pr_dPI_T_gt0\\n'
    True
    >>> f.seek(0)
    >>> tup = ReadDiffPrefs(f)
    >>> len(tup) == 6
    True
    >>> tup[0] == sites
    True
    >>> tup[1] == wts
    True
    >>> tup[2] == deltapi_means
    True
    >>> tup[3] == tup[4] == None
    True
    >>> 1.0e-6 > tup[5]['1'] - float(entries[2])
    True
    >>> f.close()
    >>> f2.seek(0)
    >>> tup = ReadDiffPrefs(f2)
    >>> len(tup) == 6
    True
    >>> tup[0] == sites
    True
    >>> tup[1] == wts
    True
    >>> tup[2] == deltapi_means
    True
    >>> tup[3] == pr_deltapi_lt0
    True
    >>> tup[4] == pr_deltapi_gt0
    True
    >>> 1.0e-6 > tup[5]['1'] - float(entries[2])
    True
    >>> f2.close()
    """
    openedfile = False
    if isinstance(f, str):
        f = open(f, 'w')
        openedfile = True
    assert sites and len(sites) >= 1 and not isinstance(sites, str), "sites must be a list with at least one entry"
    characters = deltapi_means[sites[0]].keys()
    characters.sort()
    if '*' in characters:
        i = characters.index('*')
        characters = characters[ : i] + characters[i + 1 : ] + ['*']
    f.write('# POSITION WT RMS_dPI %s' % ' '.join(['dPI_%s' % x for x in characters]))
    if (not pr_deltapi_gt0 and not pr_deltapi_lt0) or (all([not x for x in pr_deltapi_gt0.values()]) and (all([not x for x in pr_deltapi_lt0.values()]))):
        posterior_probs = False
        f.write('\n')
    else:
        f.write(' %s\n' % ' '.join(['Pr_dPI_%s_lt0 Pr_dPI_%s_gt0' % (x, x) for x in characters]))
        posterior_probs = True
    for r in sites:
        rms = dms_tools.utils.RMS([deltapi_means[r][x] for x in characters])
        f.write('%s %s %g ' % (r, wts[r], rms))
        f.write('%s' % ' '.join(['%g' % deltapi_means[r][x] for x in characters]))
        if posterior_probs:
            pr_lt0 = [pr_deltapi_lt0[r][x] for x in characters]
            pr_gt0 = [pr_deltapi_gt0[r][x] for x in characters]
            for (lt0, gt0) in zip(pr_lt0, pr_gt0):
                f.write(' %g %g' % (lt0, gt0))
        f.write('\n')
    if openedfile:
        f.close()


def ReadDiffPrefs(f):
    """Reads the differential preferences written by *WriteDiffPrefs*.

    *f* is the name of an existing file or a readable file-like object.

    The return value is the tuple: 
    *(sites, wts, deltapi_means, pr_deltapi_lt0, pr_deltapi_gt0, rms)*
    where *sites*, *wts*, *deltapi_means*, *pr_deltapi_lt0*, and 
    *pr_deltapi_gt0* all have the same values used to write the file
    with *WriteDiffPrefs* and *rms* is a dictionary with *rms[r]* giving
    the root-mean-square differential preference for each *r* in 
    *sites*.

    See docstring of *WriteDiffPrefs* for example usage.
    """
    charmatch = re.compile('^dPI_([A-z\*\-]+)$')
    if isinstance(f, str):
        f = open(f)
        lines = f.readlines()
        f.close()
    else:
        lines = f.readlines()
    characters = []
    sites = []
    wts = {}
    deltapi_means = {}
    pr_deltapi_lt0 = {}
    pr_deltapi_gt0 = {}
    rms = {}
    for line in lines:
        if line.isspace():
            continue
        elif line[0] == '#' and not characters:
            entries = line[1 : ].strip().split()
            if len(entries) < 4:
                raise ValueError("Insufficient entries in header:\n%s" % line)
            if not (entries[0] in ['POSITION', 'SITE'] and entries[1][ : 2] == 'WT' and entries[2] == 'RMS_dPI'):
                raise ValueError("Not the correct first three header columns:\n%s" % line)
            i = 3
            while i < len(entries) and charmatch.search(entries[i]):
                characters.append(charmatch.search(entries[i]).group(1))
                i += 1
            if i  == len(entries):
                pr_deltapi_lt0 = pr_deltapi_gt0 = None
                linelength = len(characters) + 3
            else:
                if not len(entries) - i == 2 * len(characters):
                    raise ValueError("Header line does not have posterior probabilities > and < format:\n%s" % line)
                if not all([entries[i + 2 * j] == 'Pr_dPI_%s_lt0' % characters[j] for j in range(len(characters))]): 
                    raise ValueError("Pr_dPI_?_lt0 mismatch in header:\n%s" % line)
                if not all([entries[i + 2 * j + 1] == 'Pr_dPI_%s_gt0' % characters[j] for j in range(len(characters))]):
                    raise ValueError("Pr_dPI_?_gt0 mismatch in header:\n%s" % line)
                linelength = 3 * len(characters) + 3
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
            assert entries[1] in characters or entries[1] == '?', "Character %s is not one of the valid ones in header. Valid possibilities: %s" % (entries[1], ', '.join(characters))
            rms[r] = float(entries[2])
            deltapi_means[r] = dict([(x, float(entries[3 + i])) for (i, x) in enumerate(characters)])
            if pr_deltapi_lt0 != None:
                pr_deltapi_lt0[r] = dict([(x, float(entries[3 + len(characters) + 2 * i])) for (i, x) in enumerate(characters)])
                pr_deltapi_gt0[r] = dict([(x, float(entries[3 + len(characters) + 2 * i + 1])) for (i, x) in enumerate(characters)])
    return (sites, wts, deltapi_means, pr_deltapi_lt0, pr_deltapi_gt0, rms)



def ReadMultiPrefOrDiffPref(infiles, removestop=False, sortsites=True):
    """Reads multiple preferences or differential preferences files.

    *infiles* is a list of preferences or differential preferences files in
    the formats read by *ReadPreferences* or *ReadDiffPrefs*. The files can
    be a mix of these two formats. Each of the files must have data for the
    same set of sites or a *ValueError* is raised. In addition, the
    files must use the same sets of characters or a *ValueError* is raised.
    Valid character sets are:

        * nucleotides

        * codons

        * amino acids without stop codons

        * amino acids with stop codons

        * a mix of amino acids with and without stop codons if *removestop*
          is *True*.

    *removestop* specifies whether we remove stop codons from amino acids. If
    this is done, all preferences are renormalized to sum to one by scaling
    them, and all differential preferences are renormalized to sum to zero
    by adding or subtracting a fixed amount to each.

    *sortsites* specifies that we do a natural sort on the site numbers
    in *infiles* rather than keeping them in their original order.

    The return value is as follows: 
    *(sites, characters, wts, databyfile)*. In this tuple:

        *sites* is a list of all site numbers (as strings) for which
        we have preferences or differential preferences.

        *characters* is a list of the characters that are being used
        (i.e. nucleotides, codons, or amino acids with or without
        stop codons).
     
        *wts* is a dictionary with *wts[r]* giving the wildtype character
        at site *r* for all *r* in *sites*. If the files do not all have
        the same wildtype at *r*, then *wts[r] = '?'*.

        *databyfile* is a dictionary keyed by each file name in *infiles*
        If that file is a preferences file, then the value is the dictionary
        *pi_means* that would be returned by *ReadPreferences*. If that file
        is a differential preferences file, then the value is the 
        dictionary *deltapi_means* that would be returned by *ReadDiffPrefs*.
    """
    assert len(infiles) >= 1, "infiles must specify at least one file"
    firstfileread = False
    for infile in infiles:
        assert os.path.isfile(infile), "Cannot find infile %s" % infile
        try:
            (isites, iwts, data, pi_95credint, h) = ReadPreferences(infile)
            if removestop:
                data = dms_tools.utils.RemoveStopFromPreferences(data)
        except:
            try:
                (isites, iwts, data, pr_deltapi_lt0, pr_deltapi_gt0, rms) = ReadDiffPrefs(infile)
                if removestop:
                    data = dms_tools.utils.RemoveStopFromDiffPrefs(data)
            except:
                raise IOError("infile %s is neither a valid preferences or differential preferences file" % infile)
        if sortsites:
            dms_tools.utils.NaturalSort(isites)
        if firstfileread:
            if isites != sites:
                raise ValueError("infile %s does not specify the same set of sites as at least one other file in: %s" % (infile, ', '.join(infiles)))
            for r in sites:
                if wts[r] != iwts[r]:
                    wts[r] = '?'
            databyfile[infile] = data
            if set(characters) != set(data[sites[0]].keys()):
                raise ValueError("Character sets do not match in all files. Expecting from the first file the following set:\n%s\nInstead got:\n%s\nin file %s" % (', '.join(characters), ', '.join(data[sites[0]].keys()), infile))
        else:
            firstfileread = True
            wts = iwts
            sites = isites
            databyfile = {infile:data}
            characters = data[sites[0]].keys()
            if set(characters) == set(dms_tools.nts):
                pass # characters are nucleotides
            elif set(characters) == set(dms_tools.codons):
                pass # characters are codons
            elif set(characters) == set(dms_tools.aminoacids_nostop):
                pass # characters are amino acids, no stop codon
            elif set(characters) == set(dms_tools.aminoacids_withstop):
                assert not removestop, "Failed to remove stop codon"
                pass # characters are amino acids, with stop codon
            else:
                raise ValueError("Invalid character set of: %s" % ', '.join(characters))
    return (sites, characters, wts, databyfile)



def IteratePairedFASTQ(r1files, r2files, gzipped, applyfilter, usegzip=False):
    """Iterates over FASTQ file pairs for paired-end sequencing reads.

    This function has been tested for its ability to read FASTQ files
    generated by the Illumina Casava 1.8 pipeline for paired-end reads.
    For FASTQ files generated by the Casava 1.8 pipeline it will ensure
    that the read directions are appropriate. If the sequence headers
    do NOT match the format of the Casava 1.8 pipeline, they will need
    to match the default format written by fastq-dump when extracting
    FASTQ files from SRA data, using the option '--split-files' to 
    generate r1 files and r2 files appropriately.

    CALLING VARIABLES:

    `r1files` : list of file name(s) containing the R1 reads.

    `r2files` : list of file name(s) containing the R2 reads. The files
    specified here must specify the exact same number of reads as
    those in `r1files`, and they must be in the same sequence (i.e.
    the first R1 read found going through `r1files` in order must
    match the first R2 read found going through `r2files` in order,
    etc).

    `gzipped` : Boolean switch specifying whether the FASTQ files
    specified by `r1files` and `r2files` are gzipped. If set to `True`,
    they are gzipped. This iterator will then read them without
    unzipping the files.

    `applyfilter` : Boolean switch specifying whether we remove read
    pairs in which one or more of the reads failed the Illumina
    chastity filter. If `True`, the iterator simply returns `False`
    each time it encounters a read pair where the chastity filter is
    failed. The chastity filter is indicated by Y or N in the
    sequence header -- values of Y indicate a failed filter. Note this
    option will only work if the sequence headers are the ones provided
    by the Illumina Casava 1.8 pipeline; for FASTQ files downloaded
    from the SRA using fastq-dump this filter information will not be
    present in the sequence header and the value of `applyfilter` will
    no longer be meaningful and all read pairs will be returned.

    `usegzip` : Boolean switch that is meaningful only if `gzipped`
    is `True`. This specifies that we use the Python ``gzip`` module to 
    unzip the reads. Otherwise we use the system ``gzip -cd command``,
    which may be faster. You should definitely specify
    `usegzip` to `True` if you are only reading a few lines, as with
    `usegzip` set to `False` the entire file is unzipped before anything happens.
    `False` by default.

    RESULT OF THIS FUNCTION:

    The function will iterate over all read pairs in the indicated
    files. It will first try to match sequence headers that follow the format
    of the Illumina Casava 1.8 pipeline, and if that fails it will match headers
    more permissively if they begin with the default structure written by fastq-dump:
    `@(samplename).(read number)`.
    
    For each iteration, it will return:

    * The return variable is `False` if ``applyfilter == True`` and one of the
      reads in the pair failed the Illumina chastity filter and the sequence 
      headers are of the type written by the Illumina Casava 1.8 pipeline.
    
    * Otherwise return variable is the following tuple:
      `(name, r1, r2, q1, q2)` where:

      - `name` : name of the read pair (a string)

      - `r1` : sequence of the first read R1 (a string)

      - `r2` : sequence of the second read R2 (a string)

      - `q1` : string of Q-scores for the R1 read 

      - `q2` : string of Q-scores for the R2 read

    """
    hm = re.compile('^\@(?P<name>[\-\w]+\:\d+\:[\w\-]+\:\d+\:\d+\:\d+\:\d+)'\
            + ' (?P<direction>[1,2])\:(?P<filtered>[Y,N])\:\d+\:'\
            + '\w+\\n$')
    hm_permissive = re.compile('^\@(?P<name>[\w\.]+).*\\n$')
    use_permissive = False # will switch to True if the Casava 1.8 header style is not matched
    if len(r1files) != len(r2files):
        raise IOError('r1files and r2files are of different length.')
    nfiles = len(r1files)
    ifile = 0
    while ifile < nfiles:
        if not os.path.isfile(r1files[ifile]):
            raise IOError("Can't find file %s" % r1files[ifile])
        if not os.path.isfile(r2files[ifile]):
            raise IOError("Can't find file %s" % r2files[ifile])
        if gzipped:
            # rather than using gzip module, we try using gzip -cd and then
            # read output. This is supposed to give better performance:
            # http://codebright.wordpress.com/2011/03/25/139/
            # unfortunately piping to subprocess.PIPE does not work
            # for the reason described in http://www.macaronikazoo.com/?p=607
            # so we have to use the following less efficient procedure of
            # first writing everything to a temporary file.
            try:
                if usegzip:
                    raise OSError # to make us go to except statement
                f1 = tempfile.TemporaryFile()
                p1 = subprocess.Popen(['gzip', '-cd', r1files[ifile]],\
                    stdout=f1.fileno())
                p1.wait()
                f1.seek(0)
                f2 = tempfile.TemporaryFile()
                p2 = subprocess.Popen(['gzip', '-cd', r2files[ifile]],\
                    stdout=f2.fileno())
                p2.wait()
                f2.seek(0)
            except OSError:
                # cannot use gzip -cd, try using gzip module
                if not usegzip:
                    warnings.warn("Failed to use gzip -cd, will have to use "+\
                            "the slower Python gzip module instead.",\
                            RuntimeWarning)
                f1 = gzip.open(r1files[ifile])
                f2 = gzip.open(r2files[ifile])
        else:
            f1 = open(r1files[ifile])
            f2 = open(r2files[ifile])
        try:
            while True:
                (h1, s1, x1, q1) = (f1.next(), f1.next(), f1.next(),\
                        f1.next())
                (h2, s2, x2, q2) = (f2.next(), f2.next(), f2.next(),\
                        f2.next())
                if use_permissive:
                    h1m = hm_permissive.search(h1)
                    h2m = hm_permissive.search(h2)
                else:
                    h1m = hm.search(h1)
                    h2m = hm.search(h2)
                    if (not h1m) and (not h2m): # neither header matched, so try the permissive mode from now on.
                        use_permissive=True
                        h1m = hm_permissive.search(h1)
                        h2m = hm_permissive.search(h2)
                if not h1m:
                    raise IOError("Failed to match head:\n%s" % h1)
                if not h2m:
                    raise IOError("Failed to match head:\n%s" % h2)
                if not use_permissive:
                    (n1, d1, cf1) = (h1m.group('name'),\
                        h1m.group('direction'), h1m.group('filtered'))
                    (n2, d2, cf2) = (h2m.group('name'),\
                        h2m.group('direction'), h2m.group('filtered'))
                else:
                    n1 = h1m.group('name') # the permissive matcher only finds 'name'; it doesn't assume that direction or filtered info is present.
                    n2 = h2m.group('name')
                if n1 != n2:
                    raise IOError("Name mismatch:\n%s\n%s" % (n1, n2))
                if not use_permissive:
                    if not (d1 == '1' and d2 == '2'):
                        raise IOError("Direction mismatch:\n%s\n%s" % (n1, n2))
                    if applyfilter and (cf1 == 'Y' or cf2 == 'Y'):
                        yield False # failed chastity filter
                s1 = s1.strip()
                s2 = s2.strip()
                q1 = q1.strip()
                q2 = q2.strip()
                yield (n1, s1, s2, q1, q2)
        except StopIteration:
            f1.close()
            f2.close()
            ifile += 1


def ReadSummaryStats(f):
    """Reads alignment summary statistics.

    *f* should be either a readable file-like object or a string giving
    the name of an existing file that contains alignment summary statistics.
    For instance, this might be the ``summarystats.txt`` file produced by
    ``dms_barcodedsubamplicons``. Here is example content from such
    a file::

        total read pairs = 4580258
        read pairs that fail Illumina filter = 0
        low quality read pairs = 915131
        barcodes randomly purged = 0
        discarded barcodes with 1 reads = 844768
        aligned barcodes with 2 reads = 466387
        discarded barcodes with 2 reads = 6457
        un-alignable barcodes with 2 reads = 40753
        aligned barcodes with 3 reads = 262149
        discarded barcodes with 3 reads = 5426
        un-alignable barcodes with 3 reads = 24891
        aligned barcodes with 4 reads = 120828
        discarded barcodes with 4 reads = 3281
        un-alignable barcodes with 4 reads = 7675
        aligned barcodes with 5 reads = 45342
        discarded barcodes with 5 reads = 1516
        un-alignable barcodes with 5 reads = 2863

    The return value is the 2-tuple *(keylist, stats)*.
    *keylist* is a list of all of the keys in *f* (these
    are the strings before the equals sign) in the same
    order as they appear in the file. *stats* is a dictionary
    keyed by each key in *keylist* and with the value being
    an integer giving the number of counts for that key.
    The reason for returning *keylist* as well as *stats*
    is that *keylist* preserves the order of the entries in *f*,
    whereas the dictionary *stats* will have an arbitrary order.

    Any lines in *f* that begin with ``#`` are treated
    as comments and are not processed into the return 
    variables.
    """
    if isinstance(f, str):
        if not os.path.isfile(f):
            raise IOError("Cannot find file %s" % f)
        openedfile = True
        f = open(f)
    keylist = []
    stats = {}
    for line in f:
        if line and not line.isspace() and line[0] != '#':
            tup = line.split('=')
            if len(tup) != 2:
                raise ValueError("Line does not have expected two entries separated by an equals sign: %s" % line)
            key = tup[0].strip()
            try:
                value = int(tup[1])
            except ValueError:
                raise ValueError("Second entry in line is not a valid integer.")
            keylist.append(key)
            stats[key] = value
    if openedfile:
        f.close()
    return (keylist, stats)


def ReadSubassembledVariants(infile, assumechartype='codon'):
    """Reads subassembled variants file produced by ``dms_subassemble``.

    *infile* : name of a ``*_subassembled_variants.txt`` file
    produced by ``dms_subassemble``.

    *assumechartype* is the type of nucleotide character assumed
    if it is **not** possible to auto-determine this from *infile*,
    which will be the case only if there are no mutations.

    The return value is *(bc_dict, bc_len, wtseq, chartype)*

    *bc_dict* is a dictionary keyed by each barcode, and with values
    being the list of mutations at that site.

    *bc_len* is the integer length of barcodes.

    *wtseq* is a string giving the wildtype sequence

    *chartype* indicates the character type of the sequence / 
    mutations as a string, which is auto-determined from
    *infile*. The *chartype* can be any of *codon*, *DNA*,
    or *aminoacid*. If all variants are unmutated and the sequence
    is a nucleotide sequence, then it is not possible to auto-
    determine whether the character type is *codon* or *DNA*. In that
    case, it is assumed to have the type of *assumechartype*.
    """
    bc_dict = {}
    mutmatch = re.compile('^(?P<wt>[ATCGatcg]+)(?P<pos>\d+)(?P<mut>[ATCGatcg]+)$')
    wtseq = len_char = bc_len = chartype = None
    with open(infile) as f:
        for line in f:
            cats = line.strip().split()
            assert len(cats) == 3, "Not 3 entries in subassembly line in %s:\n%s" % (infile, line)
            barcode = cats[0]
            if bc_len == None:
                bc_len = len(barcode)
            elif len(barcode) != bc_len:
                raise ValueError("Not all barcodes the same length in %s" % subassembled_variants_file)
            seq = cats[1].upper()
            mut_str = cats[2]
            if mut_str == "no_mutations":
                muts = []
                if wtseq == None:
                    wtseq = seq
                else:
                    assert wtseq == seq, "Mismatched wildtype sequence in %s:\n%s" % (infile, line)
            else:
                muts = mut_str.split(',')
                # the below checks wildtype sequence for match w/ muttions
                for mut in muts:
                    m = mutmatch.search(mut)
                    assert m, "Invalid mutation of %s in %s, line:\n%s" % (mut, infile, line)
                    (x, i, y) = (m.group('wt').upper(), int(m.group('pos')), m.group('mut').upper())
                    if len_char == None:
                        len_char = len(x)
                        if len_char == 3:
                            chartype = 'codon'
                            assert len(seq) % 3 == 0, "mutations indicate chartype of codon, but sequence length is not a multiple of 3 in %s" % infile
                            assert re.search('^[%s]+$' % ''.join(dms_tools.nts), seq), "mutations indicate chartype of codon, but sequence not all nucleotide characters in %s" % infile
                        elif len_char == 1:
                            if re.search('^[%s]+$' % ''.join(dms_tools.nts), seq):
                                chartype = 'DNA'
                            elif re.search('^[%s\*]+$' % ''.join(dms_tools.aminoacids_nostop), seq):
                                chartype = 'protein'
                            else:
                                raise ValueError("Unrecognized character type in sequence in %s:\n%s" % (nfile, seq))
                        else:
                            raise ValueError("Invalid character length in mutations in %s: %s" % (infile, mut))
                    assert len(x) == len(y) == len_char, "Invalid mutation of %s in %s. Expected characters of length %d. Line is:\n%s" % (mut, infile, len_char, line)
                    if chartype == 'codon':
                        i = (i - 1) * 3 + 1
                    if wtseq:
                        assert wtseq[i - 1 : i - 1 + len_char] == x, "Invalid wildtype in mutation %s in %s:\n%s" % (mut, infile, line)
                    assert seq[i - 1 : i - 1 + len_char] == y, "Invalid mutant in mutation %s in %s:\n%s" % (mut, infile, line)
                    seq = seq[ : i - 1] + x + seq[i - 1 + len_char : ]
                if wtseq == None:
                    wtseq = seq
                else:
                    assert wtseq == seq, "Mismatched wildtype sequence in %s:\n%s" % (infile, line)
            if barcode in bc_dict:
                raise ValueError("Duplicate barcode %s in %s" % (barcode, infile))
            bc_dict[barcode] = muts
    if chartype == None:
        if re.search('^[%s]+$' % ''.join(dms_tools.nts), wtseq):
            chartype = assumechartype
        elif re.search('^[%s\*]+$' % ''.join(dms_tools.aminoacids_nostop), wtseq):
            chartype = 'aminoacid'
        else:
            raise ValueError("Cannot auto-identify character type of wtseq:\n%s" % wtseq)
    return (bc_dict, bc_len, wtseq, chartype)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
