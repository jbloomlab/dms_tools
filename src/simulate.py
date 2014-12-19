"""
================================
``simulate`` module
================================

Module for simulating test data for analysis with ``dms_tools``.

Written by Jesse Bloom.

Functions in this module
---------------------------

Function documentation
---------------------------

"""


import os
import re
import math
import random
import numpy.random
import dms_tools
import dms_tools.file_io


def SimulateLibraryCounts(depths_outfiles, prefs_file, avgmu, avgepsilon, avgrho, mu_concentration=300.0, epsilon_concentration=300.0, rho_concentration=300.0, seed=1):
    """Simulates counts data for deep mutational scanning experiment.

    This function can be used if you want to simulate the creation of a mutant library
    with a set of known site-specific preferences to then test how accurately you can
    infer these preferences.

    This function assumes mutations made at codon level, and that the site-specific preferences
    are at the amino-acid level.

    All random variables are drawn from Dirichlet distributions with the indicated
    concentration parameters.

    CALLING VARIABLES:

    * *depths_outfiles* : list of 2-tuples *(depth, outfile)*.

        - *depth* is an integer giving the per-site sequencing depth.

        - *outfile* is a string giving the prefix for the created
          deep mutational scanning counts file. These files are deleted
          if they already exist, and then re-created. The suffixes are:

            - ``_pre.txt`` for pre-selection counts

            - ``_post.txt`` for post-selection counts

            - ``_errpre.txt`` for pre-selection error-control counts,
              only created if *avgepsilon* is not zero.

            - ``_errpost.txt`` for post-selection error-control counts,
              only created if *avgrho* is not zero.

    * *prefs_file* : name of an existing file giving the actual site-specific amino-acid
      preferences, either with or without stop codons. Data are simulated using these preferences.
      If stop codons are not specified, the preferences for such codons are set to zero.

    * *avgmu* : average mutation rate at each site to **any** non-wildtype codon. The actual
      mutation rate to each possible codon is drawn from a Dirichlet.

    * *avgepsilon* : average error rate in pre-selection library. This is the average
      rate at which each nucleotide is mis-read as some other nucleotide. Set to zero
      if there are no errors.

    * *avgrho* : average error rate in post-selection library. This is the average
      rate at which each nucleotide is mis-read as some other nucleotide. Set to zero
      if there are no errors.

    * *mu_concentration* : mutation rates for any given codon drawn
      from Dirichlet with this concentration parameter.

    * *epsilon_concentration* : concentration parameter for pre-selection error rate.  

    * *rho_concentration* : concentratiion parameter for post-selection error rate.

    * *seed* is seed for random number generators

    RESULT OF THIS FUNCTION:

    The result of this function is creation of the files specified in *depths_outfiles*.
    """
    numpy.random.seed(seed)
    random.seed(seed)

    (sites, wts, pi_means, credints, h) = dms_tools.file_io.ReadPreferences(prefs_file)
    aas = list(pi_means[sites[0]].keys())
    if set(aas) == set(dms_tools.aminoacids_withstop):
        pass
    elif set(aas) == set(dms_tools.aminoacids_nostop):
        for site in sites:
            pi_means[site]['*'] = 0.0
            aas.append('*')
    else:
        raise ValueError("prefs_file must specify preferences for amino acids")
    codons = dms_tools.codons

    outfiles = {}
    assert len(depths_outfiles) == len(set([tup[1] for tup in depths_outfiles])), "Duplicate outfile names in depths_outfiles"
    for (depth, outfileprefix) in depths_outfiles:
        pre = '%s_pre.txt' % outfileprefix
        post = '%s_post.txt' % outfileprefix
        errpre = '%s_errpre.txt' % outfileprefix
        errpost = '%s_errpost.txt' % outfileprefix
        for f in (pre, post, errpre, errpost):
            if os.path.isfile(f):
                os.remove(f)
        outfiles[depth] = [open(pre, 'w'), open(post, 'w')]
        if avgepsilon:
            outfiles[depth].append(open(errpre, 'w'))
        else:
            outfiles[depth].append(None)
        if avgrho:
            outfiles[depth].append(open(errpost, 'w'))
        else:
            outfiles[depth].append(None)
        for f in outfiles[depth]:
            if f:
                f.write('# POSITION WT %s\n' % ' '.join(codons))

    for r in sites:
        wtaa = wts[r]
        wtcodon = random.choice([x for x in codons if dms_tools.codon_to_aa[x] == wtaa])
        iwtcodon = codons.index(wtcodon)

        mu = numpy.array([avgmu / (len(codons) - 1.0)] * len(codons))
        mu[iwtcodon] = 1.0 - avgmu
        mu = numpy.random.dirichlet(mu * len(codons) * mu_concentration)

        epsilon = []
        rho = []
        for (icodon, codon) in enumerate(codons):
            ndiffs = len([i for i in range(len(codon)) if codon[i] != wtcodon[i]])
            if ndiffs:
                epsilon.append((avgepsilon / 3.0)**ndiffs) # divide by 3 since there are 3 nucleotides
                rho.append((avgrho / 3.0)**ndiffs)
            else:
                epsilon.append(0.0)
                rho.append(0.0)
        epsilon[iwtcodon] = 1.0 - sum(epsilon)
        rho[iwtcodon] = 1.0 - sum(rho)
        epsilon = numpy.random.dirichlet(numpy.array(epsilon) * len(codons) * epsilon_concentration)
        rho = numpy.random.dirichlet(numpy.array(rho) * len(codons) * rho_concentration)

        delta = numpy.zeros(len(codons))
        delta[iwtcodon] = 1.0

        pi = numpy.array([pi_means[r][dms_tools.codon_to_aa[codon]] for codon in codons])

        for (depth, outfileprefix) in depths_outfiles:
            (f_pre, f_post, f_errpre, f_errpost) = outfiles[depth]
            f_pre.write('%s %s %s\n' % (r, wtcodon, ' '.join(["%g" % pix for pix in numpy.random.multinomial(depth, mu + epsilon - delta)])))
            if f_errpre:
                f_errpre.write('%s %s %s\n' % (r, wtcodon, ' '.join(["%g" % pix for pix in numpy.random.multinomial(depth, epsilon)])))
            if f_errpost:
                f_errpost.write('%s %s %s\n' % (r, wtcodon, ' '.join(["%g" % pix for pix in numpy.random.multinomial(depth, rho)])))
            f_post.write('%s %s %s\n' % (r, wtcodon, ' '.join(["%g" % pix for pix in numpy.random.multinomial(depth, pi * mu / numpy.dot(pi, mu) + rho - delta)])))

    for (depth, outfileprefix) in depths_outfiles:
        for f in outfiles[depth]:
            if f:
                f.close()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
