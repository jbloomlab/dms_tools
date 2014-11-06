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

Function documentation
---------------------------

"""


import re
import dms_tools


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



if __name__ == '__main__':
    import doctest
    doctest.testmod()
