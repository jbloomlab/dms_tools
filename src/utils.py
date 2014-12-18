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

* *RMS* : Computes root-mean-square value.

* *RemoveStopFromPreferences* : removes stop codon

* *RemoveStopFromDiffPrefs* : removes stop codon

* *SumPrefsDiffPrefs* : sums preferences and differential preferences

* *Pref_or_DiffPref* : determines whether data represent preferences or differential preferences.

* *PrefsToEnrichments* : converts preferences to enrichment ratios.

* *AdjustErrorCounts* : adjust error counts so they don't exceed actual counts.

Function documentation
---------------------------

"""


import re
import math
import copy
import dms_tools


def PrefsToEnrichments(wts, pi_means, include_wt):
    """Converts site-specific preferences to enrichment ratios.

    *wts* and *pi_means* have the same meaning as the arguments by the
    same name returned by *file_io.ReadPreferences*.

    If *include_wt* is *True*, then the returned dictionary includes
    ratios for the wildtype characters (these are always one). If it is
    *False*, then ratios for the wildtype characters are not included.

    The return value is a dictionary *enrichments*, where *enrichments[r][x]*
    returns the enrichment ratio for character *x* at site *r* for all 
    sites and characters in *pi_means*.

    The enrichment ratio :math:`\phi_{r,x}` of site :math:`r` for 
    character :math:`x` is defined in terms of the preference
    :math:`\pi_{r,x}` as

    .. math::

        \phi_{r,x} = \\frac{\pi_{r,x}}{\pi_{r,\operatorname{wt}\left(r\\right)}}

    where :math:`\operatorname{wt}\left(r\\right)` is the wildtype residue
    at site :math:`r`.

    >>> wts = {'1':'A'}
    >>> pi_means = {'1':{'A':0.4, 'C':0.45, 'G':0.1, 'T':0.05}}
    >>> enrichments = PrefsToEnrichments(wts, pi_means, include_wt=True)
    >>> print "%.3f" % enrichments['1']['A']
    1.000
    >>> print "%.3f" % enrichments['1']['C']
    1.125
    >>> print "%.3f" % enrichments['1']['G']
    0.250
    >>> print "%.3f" % enrichments['1']['T']
    0.125
    >>> enrichments = PrefsToEnrichments(wts, pi_means, include_wt=False)
    >>> 'A' in enrichments
    False
    """
    enrichments = {}
    for (r, pir) in pi_means.items():
        wt = wts[r]
        wtpi = float(pir[wt])
        assert wtpi > 0, "wildtype preference is not > 0"
        enrichments[r] = dict([(x, pir[x] / wtpi) for x in pir.keys() if x != wt or include_wt])
    return enrichments


def Pref_or_DiffPref(data, tol=1.0e-3):
    """Determines if data represent preferences or differential preferences.

    *data* is a dictionary keyed by site, with each entry a dictionary
    keyed by character giving the preference / differential preference
    for that character.

    If values for each site are >= 0 and <= 1 and sum to one (within *tol*), then
    returns 'preferences'.

    If values are <= 1 and sum to zero (within *tol*), then returns
    'diffprefs'.

    Otherwise returns 'neither preferences nor diffprefs'

    >>> data = {'1':{'A':0.5, 'T':0.1, 'C':0.2, 'G':0.2}}
    >>> print Pref_or_DiffPref(data)
    preferences

    >>> data = {'1':{'A':-0.5, 'T':0.1, 'C':0.2, 'G':0.2}}
    >>> print Pref_or_DiffPref(data)
    diffprefs

    >>> data = {'1':{'A':-0.5, 'T':0.1, 'C':0.2, 'G':0.1}}
    >>> print Pref_or_DiffPref(data)
    neither preferences nor diffprefs

    """
    sums = []
    for (r, rdata) in data.iteritems():
        sums.append((sum(rdata.values())))
    if all([abs(rsum - 1.0) <= tol for rsum in sums]):
        return 'preferences'
    elif all([abs(rsum) <= tol for rsum in sums]):
        return 'diffprefs'
    else:
        return 'neither preferences nor diffprefs'


def AdjustErrorCounts(wt, counts, maxexcess):
    """Adjust error counts so that they don't exceed counts of interest.

    This deals with the case when the error control estimates a higher
    rate of a mutation than is found in the actual sample. This will confound
    the MCMC inference by ``dms_inferprefs`` and ``dms_inferdiffprefs``.
    So this script adjusts down the error counts until they are compatbile
    with the counts in the actual sample.

    Returns a dictionary like *counts* but with adjusted counts.

    This is a function utilized to clean counts data for ``dms_inferprefs``
    and ``dms_inferdiffprefs``.

    Essentially, it checks if the mutation rate in the error control
    is higher than that of the sample for which it is supposed to control.
    If so, it reduces the mutation counts for the error control so that
    they only exceed the controlled for sample by *maxexcess*.

    >>> wt = 'A'
    >>> counts = {'nrpre':{'A':500, 'C':10, 'G':40, 'T':20}, 'nrerrpre':{'A':250, 'C':1, 'G':30, 'T':10}, 'nrpost':{'A':500, 'C':20, 'G':20, 'T':20}, 'nrerrpost':{'A':1000, 'C':30, 'G':20, 'T':50}}
    >>> adjusted_counts = AdjustErrorCounts(wt, counts, maxexcess=1)
    >>> counts['nrpre'] == adjusted_counts['nrpre']
    True
    >>> counts['nrpost'] == adjusted_counts['nrpost']
    True
    >>> adjusted_counts['nrerrpre'] == {'A':250, 'C':1, 'G':21, 'T':10}
    True
    >>> adjusted_counts['nrerrpost'] == {'A':1000, 'C':30, 'G':20, 'T':41}
    True

    >>> wt = 'A'
    >>> counts = {'nrpre':{'A':500, 'C':10, 'G':40, 'T':20}, 'nrpost':{'A':500, 'C':20, 'G':20, 'T':20}, 'nrerr':{'A':1000, 'C':30, 'G':20, 'T':50}}
    >>> adjusted_counts = AdjustErrorCounts(wt, counts, maxexcess=1)
    >>> counts['nrpre'] == adjusted_counts['nrpre']
    True
    >>> counts['nrpost'] == adjusted_counts['nrpost']
    True
    >>> adjusted_counts['nrerr'] == {'A':1000, 'C':21, 'G':20, 'T':41}
    True

    >>> counts = {'nrstart':{'A':500, 'C':10, 'G':40, 'T':20}, 'nrs1':{'A':1000, 'C':20, 'G':80, 'T':40}, 'nrs2':{'A':500, 'C':20, 'G':20, 'T':20}, 'nrerr':{'A':1000, 'C':30, 'G':20, 'T':50}}
    >>> adjusted_counts = AdjustErrorCounts(wt, counts, maxexcess=1)
    >>> counts['nrstart'] == adjusted_counts['nrstart']
    True
    >>> counts['nrs1'] == adjusted_counts['nrs1']
    True
    >>> counts['nrs2'] == adjusted_counts['nrs2']
    True
    >>> adjusted_counts['nrerr'] == {'A':1000, 'C':21, 'G':20, 'T':41}
    True
    """
    returncounts = copy.deepcopy(counts)
    if 'nrpre' in counts and 'nrpost' in counts:
        # data for preferences
        mutcharacters = [x for x in counts['nrpre'].keys() if x != wt and x != 'WT']
        if len(counts) == 2:
            pass
        elif len(counts) == 3 and 'nrerr' in counts:
            # same errors
            for x in mutcharacters:
                returncounts['nrerr'][x] = min(counts['nrerr'][x], int(maxexcess + round(counts['nrpre'][x] / float(counts['nrpre'][wt]) * counts['nrerr'][wt])), int(maxexcess + round(counts['nrpost'][x] / float(counts['nrpost'][wt]) * counts['nrerr'][wt])))
        elif len(counts) == 4 and 'nrerrpre' in counts and 'nrerrpost' in counts:
            # different errors
            for x in mutcharacters:
                returncounts['nrerrpre'][x] = min(counts['nrerrpre'][x], int(maxexcess + round(counts['nrpre'][x] / float(counts['nrpre'][wt]) * counts['nrerrpre'][wt])))
                returncounts['nrerrpost'][x] = min(counts['nrerrpost'][x], int(maxexcess + round(counts['nrpost'][x] / float(counts['nrpost'][wt]) * counts['nrerrpost'][wt])))
        else:
            raise ValueError("counts does not have valid set of keys:\n%s" % str(counts.keys()))
    elif 'nrstart' in counts and 'nrs1' in counts and 'nrs2' in counts:
        # data for differential preferences
        mutcharacters = [x for x in counts['nrstart'].keys() if x != wt and x != 'WT']
        if len(counts) == 3:
            pass
        elif len(counts) == 4 and 'nrerr' in counts:
            for x in mutcharacters:
                returncounts['nrerr'][x] = min(counts['nrerr'][x], int(maxexcess + round(counts['nrstart'][x] / float(counts['nrstart'][wt]) * counts['nrerr'][wt])), int(maxexcess + round(counts['nrs1'][x] / float(counts['nrs1'][wt]) * counts['nrerr'][wt])), int(maxexcess + round(counts['nrs2'][x] / float(counts['nrs2'][wt]) * counts['nrerr'][wt])))
    else:
        raise ValueError("counts does not have valid set of keys:\n%s" % str(counts.keys()))
    return returncounts


def SumPrefsDiffPrefs(datalist, minus=[]):
    """Sums preferences and differential preferences.

    *datalist* is a list of preferences or differential preferences, which must be keyed
    by the same sets of sites and characters.

    *minus* is an optional list that is **subtracted** rather than added in the sum.

    Returns a dictionary giving the sums.

    >>> prefs = {'1':{'A':0.2, 'T':0.1, 'C':0.5, 'G':0.2}}
    >>> diffprefs = {'1':{'A':-0.1, 'T':0.3, 'C':-0.2, 'G':0.0}}
    >>> datasum = SumPrefsDiffPrefs([prefs, diffprefs])
    >>> print '%.2f' % datasum['1']['A']
    0.10
    >>> print '%.2f' % datasum['1']['T']
    0.40
    >>> print '%.2f' % datasum['1']['C']
    0.30
    >>> print '%.2f' % datasum['1']['G']
    0.20
    >>> minusprefs = {'1':{'A':0.1, 'T':0.1, 'C':-0.1, 'G':-0.1}}
    >>> datasum = SumPrefsDiffPrefs([prefs, diffprefs], minus=[minusprefs])
    >>> print '%.2f' % datasum['1']['A']
    0.00
    >>> print '%.2f' % datasum['1']['T']
    0.30
    >>> print '%.2f' % datasum['1']['C']
    0.40
    >>> print '%.2f' % datasum['1']['G']
    0.30

    """
    assert len(datalist) >= 1, "datalist must have at least one entry."
    sites = datalist[0].keys()
    assert len(sites) >= 1, "datalist must have at least one site"
    characters = datalist[0][sites[0]].keys()
    datasum = {}
    for r in sites:
        datasum[r] = dict([(x, 0.0) for x in characters])
        for idata in datalist:
            assert set(sites) == set(idata.keys()), "Not the same set of sites in all preferences / differential preferences"
            assert set(characters) == set(idata[r].keys()), "Not the same set of characters in all preferences / differential preferences"
            for x in characters:
                datasum[r][x] += idata[r][x]
        for idata in minus:
            assert set(sites) == set(idata.keys()), "Not the same set of sites in all preferences / differential preferences"
            assert set(characters) == set(idata[r].keys()), "Not the same set of characters in all preferences / differential preferences"
            for x in characters:
                datasum[r][x] -= idata[r][x]
    return datasum


def RemoveStopFromDiffPrefs(deltapi_means):
    """Removes stop codon and rescales differential preferences to sum to zero.

    Stop codons are denoted by ``*``.

    Essentially, the preference assigned to the stop codon is removed
    and all other preferences are scaled by an additive factor
    so they sum to zero.

    If there is not a stop codon present, then nothing is done.

    >>> aas = dms_tools.aminoacids_withstop
    >>> deltapi_means = {'1':dict([(aa, 1.0 / (len(aas) - 1)) for aa in aas if aa != '*'])}
    >>> deltapi_means['1']['*'] = -1.0
    >>> nostop_deltapi_means = RemoveStopFromDiffPrefs(deltapi_means)
    >>> deltapi_means == nostop_deltapi_means
    False
    >>> all(["%.3f" % nostop_deltapi_means['1'][aa] == '0.000' for aa in aas if aa != '*'])
    True
    >>> '*' in nostop_deltapi_means['1']
    False

    >>> aas = dms_tools.aminoacids_withstop
    >>> deltapi_means = {'1':dict([(aa, 1.0 / (len(aas) - 1)) for aa in aas[1 : ]])}
    >>> deltapi_means['1'][aas[0]] = -1.0
    >>> nostop_deltapi_means = RemoveStopFromDiffPrefs(deltapi_means)
    >>> deltapi_means == nostop_deltapi_means
    False
    >>> '*' in nostop_deltapi_means['1']
    False
    >>> print '%.5f' % nostop_deltapi_means['1'][aas[0]]
    -0.99750
    >>> print '%.5f' % nostop_deltapi_means['1'][aas[1]]
    0.05250

    >>> aas = dms_tools.aminoacids_nostop
    >>> deltapi_means = {'1':dict([(aa, 1.0 / (len(aas) - 1)) for aa in aas[1 : ]])}
    >>> deltapi_means['1'][aas[0]] = -1.0
    >>> nostop_deltapi_means = RemoveStopFromDiffPrefs(deltapi_means)
    >>> deltapi_means == nostop_deltapi_means
    True

    """
    nostop_deltapi_means = {}
    for (r, rdeltapi) in deltapi_means.iteritems():
        if '*' in rdeltapi:
            istopx = deltapi_means[r]['*']
            shift = istopx / float(len(rdeltapi) - 1)
            nostop_deltapi_means[r] = dict([(x, deltapix + shift) for (x, deltapix) in rdeltapi.iteritems() if x != '*'])
        else:
            nostop_deltapi_means[r] = rdeltapi
    return nostop_deltapi_means


def RemoveStopFromPreferences(pi_means):
    """Removes stop codon and renormalizes preference to sum to one.

    Stop codons are denoted by ``*``.

    Essentially, the preference assigned to the stop codon is removed
    and all other preferences are scaled by a multiplicative factor
    :math:`\ge 1` until they sum to one.

    If there is not a stop codon present, then nothing is done.

    >>> aas = dms_tools.aminoacids_withstop
    >>> pi_means = {'1':dict([(aa, 1.0 / len(aas)) for aa in aas]), '2':dict([(aas[0], 0.9)] + [(aa, 0.1 / (len(aas) - 1)) for aa in aas[1 : ]])}
    >>> nostop_pi_means = RemoveStopFromPreferences(pi_means)
    >>> pi_means == nostop_pi_means
    False
    >>> all(["%.3f" % nostop_pi_means['1'][aa] == '0.050' for aa in aas if aa != '*'])
    True
    >>> '*' in nostop_pi_means['1']
    False
    >>> '*' in nostop_pi_means['2']
    False
    >>> print '%.4f' % nostop_pi_means['2'][aas[0]]
    0.9045
    >>> print '%.6f' % nostop_pi_means['2'][aas[1]]
    0.005025

    >>> aas = dms_tools.aminoacids_nostop
    >>> pi_means = {'1':dict([(aa, 1.0 / len(aas)) for aa in aas]), '2':dict([(aas[0], 0.9)] + [(aa, 0.1 / (len(aas) - 1)) for aa in aas[1 : ]])}
    >>> nostop_pi_means = RemoveStopFromPreferences(pi_means)
    >>> pi_means == nostop_pi_means
    True
    """
    nostop_pi_means = {}
    for (r, rpi) in pi_means.iteritems():
        if '*' in rpi:
            istopx = pi_means[r]['*']
            assert 0 <= istopx < 1, "Preference for stop codon must be >= 0 and < 1, but got %g" % istopx
            scale = 1.0 / (1.0 - istopx)
            nostop_pi_means[r] = dict([(x, pix * scale) for (x, pix) in rpi.iteritems() if x != '*'])
        else:
            nostop_pi_means[r] = rpi
    return nostop_pi_means


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

    Any entries that are integers or floats are converted to strings.

    This function simply returns *None*, but *xlist* is sorted in place,
    with the integer/float to string conversion if relevant.

    >>> xlist = ['1', '20', '100', '5A', '51A', '5']
    >>> NaturalSort(xlist)
    >>> print xlist
    ['1', '5', '5A', '20', '51A', '100']

    >>> xlist = [1, '20', 100, '5A', '51A', '5']
    >>> NaturalSort(xlist)
    >>> print xlist
    ['1', '5', '5A', '20', '51A', '100']
    """
    for i in range(len(xlist)):
        if isinstance(xlist[i], (int, float)):
            xlist[i] = str(xlist[i])
    nsre = re.compile('([0-9]+\.{0,1}[0-9]*)')
    def natural_sort_keys(s):
        return [int(text) if text.isdigit() else text.lower() for text in re.split(nsre, s)]
    xlist.sort(key=natural_sort_keys)


def AvgMutRate(counts, chartype):
    """Returns average mutation rate.
    
    *counts* is a dictionary of deep mutational scanning counts in the format returned
    by *file_io.ReadDMSCounts*.
    
    *chartype* is the character type of these counts, with the same meaning 
    as for *file_io.ReadDMSCounts*. 

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
    elif chartype == 'aminoacids_withstop':
        characters = dms_tools.aminoacids_withstop
        epsilon = dict([(m, 0.0) for m in range(2)])
    elif chartype == 'aminoacids_nostop':
        characters = dms_tools.aminoacids_nostop
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
