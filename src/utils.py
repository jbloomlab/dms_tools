"""
================================
``utils`` module
================================

Module with utilities for manipulating data related to ``dms_tools``.

Written by Jesse Bloom.

Functions in this module
---------------------------

* *NaturalSort* : Performs a natural sort (alphanumeric sort) on a list.

* *AvgMutRate* : Returns the average mutation rate across all sites.

* *SumCodonToAA* : Sums all codon entries for each amino acid.

* *SiteEntropy* : Computes site entropy.

* *RMS* : Computs root-mean-square value.

Function documentation
---------------------------

"""


import re
import math
import dms_tools


def RMS(xlist):
    """Computes root-mean-square value of entries in *xlist*.

    >>> xlist = [-0.2, 0.1, 0.03, 0.5]
    >>> print "%.5f" % RMS(xlist)
    0.27427
    """
    assert len(xlist) >= 1, "xlist must have lenght of at least one"
    return math.sqrt(sum([x**2 for x in xlist]) / float(len(xlist)))


def SiteEntropy(xlist, base=2.0):
    """Computes site entropy to log base *base*.

    *xlist* should be a list of numbers :math:`\ge 0` that
    sum to one.

    >>> print "%.3f" % SiteEntropy([1.0, 0.0])
    0.000
    >>> print "%.3f" % SiteEntropy([0.5, 0.5])
    1.000

    """
    assert abs(sum(xlist) - 1.0) < 1.0e5, "Sum of xlist is not close to one: %g" % sum(xlist)
    assert all([x >= 0 for x in xlist]), "xlist does not have all entries >= 0"
    h = 0.0
    for x in xlist:
        if x != 0:
            h -= x * math.log(x, base)
    return h


def NaturalSort(xlist):
    """Performs a natural sort (alphanumeric sort) in a list of strings.
   
    *xlist* is a list of strings. Sorts them, reading the integer portion
    as numbers. Currently works for integers, but not floats (anything 
    including or after a decimal is handled as a string, not number).

    This function simply returns *None*, but *xlist* is sorted in place.

    >>> xlist = ['1', '20', '100', '5A', '51A', '5']
    >>> NaturalSort(xlist)
    >>> print xlist
    ['1', '5', '5A', '20', '51A', '100']
    """
    nsre = re.compile('([0-9]+\.{0,1}[0-9]*)')
    def natural_sort_keys(s):
        return [int(text) if text.isdigit() else text.lower() for text in re.split(nsre, s)]
    xlist.sort(key=natural_sort_keys)


def AvgMutRate(counts, chartype):
    """Returns average mutation rate.
    
    *counts* is a dictionary of deep mutational scanning counts in the format returned
    by *file_io.ReadDMSCounts*.
    
    *chartype* is the character type of these counts, with the same meaning 
    as for *file_io.ReadDMSCounts*. So the allowable values are the strings
    *DNA* and *codon*.

    This function computes the average mutation rate over all counts for all sites for 
    mutations with each possible number of nucleotide changes. Specifically, we define the average
    rate of mutations to characters with :math:`m` nucleotide changes relative to wildtype as

    .. math:: 
       
       \overline{\epsilon_m} = \\frac{1}{L}\sum\limits_r \\frac{1}{N_r} \sum\limits_x n_{r,x} \\times \delta_{m,D_{x,\operatorname{wt}\left(r\\right)}}

    where :math:`r` ranges over all :math:`L` sites for which there are deep mutational scanning counts,
    :math:`N_r` is the total number of counts at site :math:`r`, :math:`n_{r,x}` is the number
    of counts that report character :math:`x` at site :math:`r`, :math:`D_{x,\operatorname{wt}\left(r\\right)}`
    is the number of nucleotide differences between character :math:`x` and the wildtype character
    :math:`\operatorname{wt}\left(r\\right)` at site :math:`r`, and :math:`\delta_{xy}` is the Kronecker
    delta. Note that this definition implies :math:`1 = \sum_m \overline{\epsilon_m}`. For nucleotides,
    there are just two values of :math:`m`: a value of :math:`m = 0` is no differences (the wildtype
    character) and a value of :math:`m = 1` is the a mutant nucleotide. For codons, :math:`m` can
    be 1, 2, or 3.

    The return value is the dictionary *mutrates* which is keyed by all values of :math:`m`
    for the *chartype* in question. The value of *mutrates[m]* is the mutation rate for
    characters with this number of differences.

    >>> counts = {'1':{'WT':'A', 'A':7, 'C':1, 'G':1, 'T':1},
    ...           '2':{'WT':'C', 'A':2, 'C':8, 'G':0, 'T':0}}
    >>> mutrates = AvgMutRate(counts, 'DNA')
    >>> print "%.3f" % mutrates[0]
    0.750
    >>> print "%.3f" % mutrates[1]
    0.250
    """
    if chartype == 'codon':
        characters = dms_tools.codons
        epsilon = dict([(m, 0.0) for m in range(4)])
    elif chartype == 'DNA':
        characters = dms_tools.nts
        epsilon = dict([(m, 0.0) for m in range(2)])
    else:
        raise ValueError("Invalid chartype of %s" % chartype)
    length = len(counts)
    for (site, sitecounts) in counts.iteritems():
        wt = sitecounts['WT']
        assert wt in characters, "wildtype of %s not in characters" % wt
        nr = sum([sitecounts[x] for x in characters])
        if not nr:
            continue # no counts for this site, doesn't contribute
        denom = float(length * nr)
        for x in characters:
            m = len([i for i in range(len(wt)) if wt[i] != x[i]])
            epsilon[m] += sitecounts[x] / denom
    assert abs(sum(epsilon.values()) - 1.0) < 1.0e-5, "Sum of mutrates is not close to one: %g" % (sum(epsilon.values()))
    return epsilon


def SumCodonToAA(codondict, includestop=True):
    """Sums all codon entries for each amino acid.

    *codondict* is a dictionary keyed by all codons (*dms_tools.codons*).
    The values should be numbers.

    *includestop* specifies whether we also sum values for stop codons to
    give an entry for ``*`` (the stop codon). Do this only if has default
    value of *True*.

    The return value is *aadict* which is keyed by all amino acids and
    has values equal to the sum of the entries for the encoding
    codons in *codondict*.

    >>> codondict = dict([(codon, 0) for codon in dms_tools.codons])
    >>> codondict['GGG'] = 0.1
    >>> codondict['GGA'] = 0.2
    >>> codondict['TGA'] = 0.3
    >>> codondict['CCC'] = 0.4
    >>> aadict = SumCodonToAA(codondict)
    >>> print "%.2f" % aadict['G']
    0.30
    >>> print "%.2f" % aadict['P']
    0.40
    >>> print "%.2f" % aadict['*']
    0.30
    >>> print "%.2f" % aadict['Y']
    0.00
    >>> aadict = SumCodonToAA(codondict, includestop=False)
    >>> print "%.2f" % aadict['G']
    0.30
    >>> print "%.2f" % aadict['P']
    0.40
    >>> print "%.2f" % aadict['*']
    Traceback (most recent call last):
        ...
    KeyError: '*'
    >>> print "%.2f" % aadict['Y']
    0.00
    """
    if includestop:
        aas = dms_tools.aminoacids_withstop
    else:
        aas = dms_tools.aminoacids_nostop
    aadict = dict([(aa, 0) for aa in aas])
    for codon in dms_tools.codons:
        aa = dms_tools.codon_to_aa[codon]
        if aa in aas:
            aadict[aa] += codondict[codon]
    return aadict


if __name__ == '__main__':
    import doctest
    doctest.testmod()
