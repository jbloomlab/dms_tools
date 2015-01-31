"""
========================
``mcmc`` module
========================

Module that performs MCMC to infer preferences and differential preferences.

MCMC done using ``Stan`` (http://mc-stan.org/) via ``pystan``.

Written by Jesse Bloom.

Functions in this module
----------------------------------

* *InferSitePreferences* : use MCMC to infer site-specific preferences.

* *InferSitePreferencesFromEnrichmentRatios* :infers site-specific preferences directly from enrichment ratios.

* *InferSiteDiffPrefs* : use MCMC to infer site-specific differential preferences.

* *StanModelNoneErr* : ``pystan`` model used by *InferSitePreferences*.

* *StanModelSameErr* : ``pystan`` model used by *InferSitePreferences*.

* *StanModelDifferentErr* : ``pystan`` model used by *InferSitePreferences*.

* *StanModelDiffPrefNoErr* : ``pystan`` model used by *InferSiteDiffPrefs*.

* *StanModelDiffPrefSameErr* : ``pystan`` model used by *InferSiteDiffPrefs*.


Function documentation
----------------------------------

"""


import sys
import tempfile
import time
import math
import cPickle
import numpy
import numpy.random
import pystan


PRIOR_MIN_VALUE = 1.0e-7 # minimum value for Dirichlet prior elements


class StanModelDiffPrefNoErr(object):
    """PyStan model for *error_model* of 'none' for *InferSiteDiffPrefs*."""
    def __init__(self):
        self.pystancode =\
"""
data {
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=0> nrstart[Nchar]; // counts in starting library
    int<lower=0> nrs1[Nchar]; // counts in selection s1
    int<lower=0> nrs2[Nchar]; // counts in selection s2
    vector<lower=%g>[Nchar] pirs1_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] frstart_prior_params; // Dirichlet prior params
    real<lower=%g> deltapi_concentration; // concentration parameter for delta pi prior
}
parameters {
    simplex[Nchar] pirs1;
    simplex[Nchar] pirs2;
    simplex[Nchar] frstart;
}
transformed parameters {
    vector[Nchar] deltapir;
    simplex[Nchar] frs1;
    simplex[Nchar] frs2;
    frs1 <- frstart .* pirs1 / dot_product(frstart, pirs1);
    frs2 <- frstart .* pirs2 / dot_product(frstart, pirs2);
    deltapir <- pirs2 - pirs1;
}
model {
    pirs1 ~ dirichlet(pirs1_prior_params);
    frstart ~ dirichlet(frstart_prior_params);
    pirs2 ~ dirichlet(pirs1 * Nchar * deltapi_concentration);
    nrstart ~ multinomial(frstart);
    nrs1 ~ multinomial(frs1);
    nrs2 ~ multinomial(frs2);
}

""" % (PRIOR_MIN_VALUE, PRIOR_MIN_VALUE, PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode)


class StanModelDiffPrefSameErr(object):
    """PyStan model for *error_model* of 'same' for *InferSiteDiffPrefs*."""
    def __init__(self):
        self.pystancode =\
"""
data {
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=1, upper=Nchar> iwtchar; // index of wildtype character in 1, ... numbering
    int<lower=0> nrstart[Nchar]; // counts in starting library
    int<lower=0> nrs1[Nchar]; // counts in selection s1
    int<lower=0> nrs2[Nchar]; // counts in selection s2
    int<lower=0> nrerr[Nchar]; // counts in error control
    vector<lower=%g>[Nchar] pirs1_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] frstart_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] xir_prior_params; // Dirichlet prior params
    real<lower=%g> deltapi_concentration; // concentration parameter for delta pi prior
}
transformed data {
    simplex[Nchar] deltar;
    for (ichar in 1:Nchar) {
        deltar[ichar] <- 0.0;
    }
    deltar[iwtchar] <- 1.0;
}
parameters {
    simplex[Nchar] pirs1;
    simplex[Nchar] pirs2;
    simplex[Nchar] frstart;
    simplex[Nchar] xir;
}
transformed parameters {
    vector[Nchar] deltapir;
    simplex[Nchar] frs1;
    simplex[Nchar] frs2;
    simplex[Nchar] frstart_plus_err;
    frs1 <- frstart .* pirs1 / dot_product(frstart, pirs1) + xir - deltar;
    frs2 <- frstart .* pirs2 / dot_product(frstart, pirs2) + xir - deltar;
    frstart_plus_err <- frstart + xir - deltar;
    deltapir <- pirs2 - pirs1;
}
model {
    pirs1 ~ dirichlet(pirs1_prior_params);
    frstart ~ dirichlet(frstart_prior_params);
    xir ~ dirichlet(xir_prior_params);
    pirs2 ~ dirichlet(pirs1 * Nchar * deltapi_concentration);
    nrstart ~ multinomial(frstart_plus_err);
    nrs1 ~ multinomial(frs1);
    nrs2 ~ multinomial(frs2);
    nrerr ~ multinomial(xir);
}

""" % (PRIOR_MIN_VALUE, PRIOR_MIN_VALUE, PRIOR_MIN_VALUE, PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode)



class StanModelNoneErr(object):
    """PyStan model for *error_model* of 'none' for *InferSitePreferences*."""
    def __init__(self):
        self.pystancode =\
"""
data {
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=0> nrpre[Nchar]; // counts pre-selection
    int<lower=0> nrpost[Nchar]; // counts post-selection
    vector<lower=%g>[Nchar] pir_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] mur_prior_params; // Dirichlet prior params
}
parameters {
    simplex[Nchar] pir;
    simplex[Nchar] mur;
}
transformed parameters {
    simplex[Nchar] fr;
    fr <- pir .* mur / dot_product(pir, mur);
}
model {
    pir ~ dirichlet(pir_prior_params);
    mur ~ dirichlet(mur_prior_params);
    nrpre ~ multinomial(mur);
    nrpost ~ multinomial(fr);
}
""" % (PRIOR_MIN_VALUE, PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode)


class StanModelSameErr:
    """PyStan model for *error_model* of 'same' for *InferSitePreferences*."""
    def __init__(self):
        self.pystancode =\
"""
data {
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=1, upper=Nchar> iwtchar; // index of wildtype character in 1, ... numbering
    int<lower=0> nrpre[Nchar]; // counts pre-selection
    int<lower=0> nrpost[Nchar]; // counts post-selection
    int<lower=0> nrerr[Nchar]; // counts in error control
    vector<lower=%g>[Nchar] pir_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] mur_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] epsilonr_prior_params; // Dirichlet prior params
}
transformed data {
    simplex[Nchar] deltar;
    for (ichar in 1:Nchar) {
        deltar[ichar] <- 0.0;
    }
    deltar[iwtchar] <- 1.0;
}
parameters {
    simplex[Nchar] pir;
    simplex[Nchar] mur;
    simplex[Nchar] epsilonr;
}
transformed parameters {
    simplex[Nchar] fr_plus_err;
    simplex[Nchar] mur_plus_err;
    fr_plus_err <- pir .* mur / dot_product(pir, mur) + epsilonr - deltar;
    mur_plus_err <- mur + epsilonr - deltar;
}
model {
    pir ~ dirichlet(pir_prior_params);
    mur ~ dirichlet(mur_prior_params);
    epsilonr ~ dirichlet(epsilonr_prior_params);
    nrerr ~ multinomial(epsilonr);
    nrpre ~ multinomial(mur_plus_err);
    nrpost ~ multinomial(fr_plus_err);
}
""" % (PRIOR_MIN_VALUE, PRIOR_MIN_VALUE, PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode)


class StanModelDifferentErr:
    """PyStan model for *error_model* of 'different' for *InferSitePreferences*."""
    def __init__(self):
        self.pystancode =\
"""
data {
    int<lower=1> Nchar; // 64 for codons, 20 for amino acids, 4 for nucleotides
    int<lower=1, upper=Nchar> iwtchar; // index of wildtype character in 1, ... numbering
    int<lower=0> nrpre[Nchar]; // counts pre-selection
    int<lower=0> nrpost[Nchar]; // counts post-selection
    int<lower=0> nrerrpre[Nchar]; // counts in pre-selection error control
    int<lower=0> nrerrpost[Nchar]; // counts in post-selection error control
    vector<lower=%g>[Nchar] pir_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] mur_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] epsilonr_prior_params; // Dirichlet prior params
    vector<lower=%g>[Nchar] rhor_prior_params; // Dirichlet prior params
}
transformed data {
    simplex[Nchar] deltar;
    for (ichar in 1:Nchar) {
        deltar[ichar] <- 0.0;
    }
    deltar[iwtchar] <- 1.0;
}
parameters {
    simplex[Nchar] pir;
    simplex[Nchar] mur;
    simplex[Nchar] epsilonr;
    simplex[Nchar] rhor;
}
transformed parameters {
    simplex[Nchar] fr_plus_err;
    simplex[Nchar] mur_plus_err;
    fr_plus_err <- pir .* mur / dot_product(pir, mur) + rhor - deltar;
    mur_plus_err <- mur + epsilonr - deltar;
}
model {
    pir ~ dirichlet(pir_prior_params);
    mur ~ dirichlet(mur_prior_params);
    epsilonr ~ dirichlet(epsilonr_prior_params);
    rhor ~ dirichlet(rhor_prior_params);
    nrerrpre ~ multinomial(epsilonr);
    nrerrpost ~ multinomial(rhor);
    nrpre ~ multinomial(mur_plus_err);
    nrpost ~ multinomial(fr_plus_err);
}
""" % (PRIOR_MIN_VALUE, PRIOR_MIN_VALUE, PRIOR_MIN_VALUE, PRIOR_MIN_VALUE)
        self.model = pystan.StanModel(model_code=self.pystancode)


def _InitialValueDiffPrefs(error_model, nchains, iwtchar, nchars):
    """Gets valid initial values for PyStan differential preference inference.

    The values automatically initialized by PyStan frequently have invalid values.
    This function will generate random values that are valid
    for initialization and return a list that can be passed to the *StanModel*
    as the *init* argument.
    """
    initializationattempts = 10 # multiple attempts, as might fail for pathological random values
    nrescales = 10 # rescale non-wildtype values down this many times
    rescalefactor = 5.0 # rescale non-wildtype values down by this much
    deltar = numpy.zeros(nchars)
    deltar[iwtchar] = 1.0
    init = []
    concentrationparameters = [1.0 for ichar in range(nchars)]
    concentrationparameters[iwtchar] = 10.0 # make this element bigger, as it probably should be
    for chain in range(nchains):
        for iattempt in range(initializationattempts):
            chain_init = {\
                'pirs1':numpy.random.dirichlet(concentrationparameters),\
                'frstart':numpy.random.dirichlet(concentrationparameters),\
                'xir':numpy.random.dirichlet(concentrationparameters),\
                }
            chain_init['deltapir'] = numpy.random.dirichlet(concentrationparameters) - chain_init['frstart']
            irescale = 0
            while (irescale < nrescales) and (\
                    any(chain_init['deltapir'] + chain_init['pirs1'] < PRIOR_MIN_VALUE) or \
                    any(chain_init['deltapir'] + chain_init['pirs1'] > 1.0) or \
                    any(chain_init['pirs1'] < PRIOR_MIN_VALUE) or \
                    any(chain_init['frstart'] < PRIOR_MIN_VALUE) or \
                    any(chain_init['xir'] < PRIOR_MIN_VALUE) or \
                    any(chain_init['frstart'] + chain_init['xir'] - deltar < PRIOR_MIN_VALUE) or \
                    any(chain_init['frstart'] * chain_init['pirs1'] / numpy.dot(chain_init['frstart'], chain_init['pirs1']) + chain_init['xir'] - deltar < PRIOR_MIN_VALUE) or \
                    any(chain_init['frstart'] * (chain_init['pirs1'] + chain_init['deltapir']) / numpy.dot(chain_init['frstart'], chain_init['pirs1'] + chain_init['deltapir']) + chain_init['xir'] - deltar < PRIOR_MIN_VALUE)):
                irescale += 1
                chain_init['deltapir'] /= rescalefactor
                chain_init['xir'] /= rescalefactor
                chain_init['xir'][iwtchar] = 1.0 - sum([chain_init['xir'][ichar] for ichar in range(nchars) if ichar != iwtchar])
            if irescale < nrescales:
                break
        else:
            raise ValueError("Failed to initialize for same even after %d attempts" % initializationattempts)
        if error_model == 'none':
            del chain_init['xir']
        chain_init['pirs2'] = chain_init['deltapir'] + chain_init['pirs1']
        del chain_init['deltapir']
        init.append(chain_init)
    return init



def _InitialValuePreferences(error_model, nchains, iwtchar, nchars):
    """Gets valid initial values for PyStan preference inference.

    The values automatically initialized by PyStan frequently have invalid values.
    This function will generate random values that are valid
    for initialization and return a list that can be passed to the *StanModel*
    as the *init* argument.
    """
    initializationattempts = 10 # multiple attempts, as might fail for pathological random values
    nrescales = 10 # rescale non-wildtype values down this many times
    rescalefactor = 5.0 # rescale non-wildtype values down by this much
    deltar = numpy.zeros(nchars)
    deltar[iwtchar] = 1.0
    init = []
    concentrationparameters = [1.0 for ichar in range(nchars)]
    concentrationparameters[iwtchar] = 10.0 # make this element bigger, as it probably should be
    if error_model == 'none':
        for chain in range(nchains):
            while True:
                chain_init = {\
                    'pir':numpy.random.dirichlet(concentrationparameters),\
                    'mur':numpy.random.dirichlet(concentrationparameters),\
                    }
                if all(chain_init['pir'] > PRIOR_MIN_VALUE) and all(chain_init['mur'] > PRIOR_MIN_VALUE):
                    init.append(chain_init)
                    break
        return init
    elif error_model == 'same':
        posterrname = 'epsilonr'
    elif error_model == 'different':
        posterrname = 'rhor'
    else:
        raise ValueError("Invalid error_model of %s when getting initial values" % error_model)
    for chain in range(nchains):
        for iattempt in range(initializationattempts):
            chain_init = {\
                    'pir':numpy.random.dirichlet(concentrationparameters),\
                    'mur':numpy.random.dirichlet(concentrationparameters),\
                    'epsilonr':numpy.random.dirichlet(concentrationparameters),\
                    }
            if error_model == 'different':
                chain_init['rhor'] = numpy.random.dirichlet(concentrationparameters)
            irescale = 0
            while irescale < nrescales and (any(chain_init['pir'] * chain_init['mur'] / numpy.dot(chain_init['pir'], chain_init['mur']) + chain_init[posterrname] - deltar < PRIOR_MIN_VALUE) or any(chain_init['mur'] + chain_init['epsilonr'] - deltar < PRIOR_MIN_VALUE)): 
                irescale += 1
                chain_init['epsilonr'] /= rescalefactor
                chain_init['epsilonr'][iwtchar] = 1.0 - sum([chain_init['epsilonr'][ichar] for ichar in range(nchars) if ichar != iwtchar])
                if error_model == 'different':
                    chain_init['rhor'] /= rescalefactor
                    chain_init['rhor'][iwtchar] = 1.0 - sum([chain_init['rhor'][ichar] for ichar in range(nchars) if ichar != iwtchar])
            if all(chain_init['pir'] * chain_init['mur'] / numpy.dot(chain_init['pir'], chain_init['mur']) + chain_init[posterrname] - deltar > PRIOR_MIN_VALUE) and all(chain_init['mur'] + chain_init['epsilonr'] - deltar > PRIOR_MIN_VALUE):
                break
        else:
            raise ValueError("Failed to initialize for same even after %d attempts" % initializationattempts)
        init.append(chain_init)
    return init


def InferSitePreferencesFromEnrichmentRatios(characterlist, wtchar, error_model, counts, pseudocounts=1):
    """Infers site-specific preferences from enrichment ratios.

    This function mirrors the operations performed by *InferSitePreferences*, expect the preferences
    are calculated directly from enrichment ratios. 

    *characterlist*, *wtchar*, *error_model*, and *counts* have the same meaning as for *InferSitePreferences*.

    *pseudocounts* is a number > 0 giving the pseudocounts added to each count in *counts*.

    Briefly, for each character :math:`x`, we calculate the enrichment relative to the wildtype
    character :math:`\\rm{wt}` at this site as

    .. math::

        \phi_{x} = \\frac{\left(\\frac{\max\left(\mathcal{P}, n_x^{\\rm{post}} - n_x^{\\rm{err,post}} + \mathcal{P}\\right)}{\max\left(\mathcal{P}, n_{\\rm{wt}}^{\\rm{post}} - n_{\\rm{wt}}^{\\rm{err,post}} + \mathcal{P}\\right)}\\right)}{\left(\\frac{\max\left(\mathcal{P}, n_x^{\\rm{pre}} - n_x^{\\rm{err,pre}} + \mathcal{P}\\right)}{\max\left(\mathcal{P}, n_{\\rm{wt}}^{\\rm{pre}} - n_{\\rm{wt}}^{\\rm{err,pre}} + \mathcal{P}\\right)}\\right)}

    where :math:`\mathcal{P}` is the value of *pseudocounts*. When *error_model* is *none*, then 
    :math:`0 = n_x^{\\rm{post}} = n_x^{\\rm{pre}}` and the above equation reduces to a simple ratio
    of post- and pre-selection. With :math:`\mathcal{P} > 0`, this enrichment is always :math:`> 0` and 
    :math:`< \infty`.

    We then calculate the preference as

    .. math::

        \pi_x = \\frac{\phi_x}{\sum_y \phi_y}.

    The return value is: *(converged, pi_means, pi_95credint, logstring)*, where the tuple
    entries have the same meaning as for *InferSitePreferences* except that *pi_95credint* is
    *None* since no credible intervals can be estimated from direct enrichment ratio calculation
    as it is not a statistical model, and *converged* is *True* since this calculation
    always converges.
    """
    assert pseudocounts > 0, "pseudocounts must be greater than zero, invalid value of %g" % pseudocounts
    assert wtchar in characterlist, "wtchar %s not in characterlist %s" % (wtchar, str(characterlist))
    logstring = '\tComputed preferences directly from enrichment ratios.'
    psi = {}
    for x in characterlist:
        if error_model == 'none':
            nrerrpre = nrerrpost = 0.0
            nrerrprewt = nrerrpostwt = 0.0
        elif error_model == 'same':
            nrerrpre = nrerrpost = counts['nrerr'][x]
            nrerrprewt = nrerrpostwt = counts['nrerr'][wtchar]
        elif error_model == 'different':
            nrerrpre = counts['nrerrpre'][x]
            nrerrpost = counts['nrerrpost'][x]
            nrerrprewt = counts['nrerrpre'][wtchar]
            nrerrpostwt = counts['nrerrpost'][wtchar]
        else:
            raise ValueError("Invalid error_model of %s" % error_model)
        postratio = max(pseudocounts, counts['nrpost'][x] - nrerrpost + pseudocounts) / float(max(pseudocounts, counts['nrpost'][wtchar] - nrerrpostwt + pseudocounts))
        preratio = max(pseudocounts, counts['nrpre'][x] - nrerrpre + pseudocounts) / float(max(pseudocounts, counts['nrpre'][wtchar] - nrerrprewt + pseudocounts))
        psi[x] = postratio / preratio
    assert abs(psi[wtchar] - 1) < 1.0e-5, "wtchar does not have enrichment ratio of one: %g" % psi[wtchar]
    denom = sum(psi.values())
    pi = dict([(x, psi[x] / float(denom)) for x in characterlist])
    return (True, pi, None, logstring)



def InferSitePreferences(characterlist, wtchar, error_model, counts, priors, seed=1, niter=10000, increasetries=6, n_jobs=1, r_max=1.1, neff_min=100, nchains=4, increasefactor=2):
    """Infers site-specific preferences by MCMC for a specific site.

    Uses MCMC to infer the site-specific preferences :math:`\pi_{r,a}` for some site
    :math:`r` for each character :math:`a` by integrating over the posterior defined by the product 
    of the following priors and likelihoods, where for instance :math:`\\boldsymbol{\mathbf{\pi_r}}`
    is used to indicate the vector of :math:`\pi_{r,a}` values for all characters:

    .. math::

        \Pr\left(\\boldsymbol{\mathbf{\pi_r}}\\right) = \mathrm{Dirichlet}\left(\\boldsymbol{\mathbf{a_{\pi,r}}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{\mu_r}}\\right) = \mathrm{Dirichlet}\left(\\boldsymbol{\mathbf{a_{\mu,r}}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{\epsilon_r}}\\right) = \mathrm{Dirichlet}\left(\\boldsymbol{\mathbf{a_{\epsilon,r}}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{\\rho_r}}\\right) = \mathrm{Dirichlet}\left(\\boldsymbol{\mathbf{a_{\\rho,r}}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{n_r^{\\rm{pre}}}} \mid \\boldsymbol{\mathbf{\mu_r}}, \\boldsymbol{\mathbf{\epsilon_r}}\\right) = \mathrm{Multinomial}\left(\\boldsymbol{\mathbf{n_r^{\\rm{pre}}}}; \\boldsymbol{\mathbf{\mu_r}} + \\boldsymbol{\mathbf{\epsilon_r}} - \\boldsymbol{\mathbf{\delta_r}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{n_r^{\\rm{post}}}} \mid \\boldsymbol{\mathbf{\mu_r}}, \\boldsymbol{\mathbf{\epsilon_r}}\\right) = \mathrm{Multinomial}\left(\\boldsymbol{\mathbf{n_r^{\\rm{post}}}}; \\frac{\\boldsymbol{\mathbf{\mu_r}} \circ \\boldsymbol{\mathbf{\pi_r}}}{\\boldsymbol{\mathbf{\mu_r}} \cdot \\boldsymbol{\mathbf{\pi_r}}} + \\boldsymbol{\mathbf{\\rho_r}} - \\boldsymbol{\mathbf{\delta_r}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{n_r^{\\rm{err,pre}}}} \mid \\boldsymbol{\mathbf{\epsilon_r}}\\right) = \mathrm{Multinomial}\left(\\boldsymbol{\mathbf{n_r^{\\rm{err,pre}}}}; \\boldsymbol{\mathbf{\epsilon_r}}\\right)

        \Pr\left(\\boldsymbol{\mathbf{n_r^{\\rm{err,post}}}} \mid \\boldsymbol{\mathbf{\\rho_r}}\\right) = \mathrm{Multinomial}\left(\\boldsymbol{\mathbf{n_r^{\\rm{err,post}}}}; \\boldsymbol{\mathbf{\\rho_r}}\\right)

    where :math:`\\boldsymbol{\mathbf{\delta_r}}` is a vector that is zero except for a one at the element
    corresponding to the wildtype character.

    The MCMC automatically tries to guarantee convergence via the parameters specified
    by *r_max*, *neff_min*, *nchains*, *niter*, *increasefactor*, and *increasetries*.

    Here are the calling variables:    

        * *characterlist* is a list of all valid characters; entries must be unique. 
           Typically a list of all amino acids, codons, or nucleotides. Characters
           are case sensitive.

        * *wtchar* is a character that is in *characterlist* and defines the wildtype
          character at the site.

        * *error_model* specifies how the errors are estimated. It might
          be one of the following strings:

              - *none* : no errors (:math:`\\boldsymbol{\mathbf{\epsilon_r}} = \\boldsymbol{\mathbf{\\rho_r}} = \mathbf{0}`)

              - *same* : same error rates for pre and post-selection 
                (:math:`\\boldsymbol{\mathbf{\epsilon_r}} = \\boldsymbol{\mathbf{\\rho_r}}`).

              - *different* : different error rates for pre and post
                selection (:math:`\\boldsymbol{\mathbf{\epsilon_r}} \\ne \\boldsymbol{\mathbf{\\rho_r}}`).

          Alternatively, it could be an instance of *StanModelNoneErr*, 
          *StanModelSameErr*, or *StanModelDifferentErr*. These are the
          ``pystan`` models that implement the three strings described above.
          If you pass these model objects, then they will not need to be
          compiled by this function (strings require models to be compiled).
          So it is advantageous to pass the ``pystan`` models rather than
          the strings if calling this function repeatedly.

        * *counts* is a dictionary specifying the deep sequencing counts.
          Each string key should specify a dictionary keyed by all characters
          and with values giving the integer counts for that character. Keys:

            - *nrpre* : specifies :math:`\\boldsymbol{\mathbf{n_r^{\\rm{pre}}}}`

            - *nrpost* : specifies :math:`\\boldsymbol{\mathbf{n_r^{\\rm{post}}}}`

            - *nrerr* : only required if *error_model* is *same*, specifies
              :math:`\\boldsymbol{\mathbf{n_r^{\\rm{err,pre}}}} = \\boldsymbol{\mathbf{n_r^{\\rm{err,post}}}}`

            - *nrerrpre* and *nrerrpost* : only required if *error_model* is
              *different*, specify :math:`\\boldsymbol{\mathbf{n_r^{\\rm{err,pre}}}}` and
              :math:`\\boldsymbol{\mathbf{n_r^{\\rm{err,post}}}}`.

        * *priors* is a dictionary specifying the parameter vectors for the
          Dirichlet priors. Each string key should specify a dictionary keyed by
          all characters and with values giving the prior for that character. Keys:

            - *pir_prior_params* : specifies :math:`\\boldsymbol{\mathbf{a_{\pi,r}}}`

            - *mur_prior_params* : specifies :math:`\\boldsymbol{\mathbf{a_{\mu,r}}}`

            - *epsilonr_prior_params* : only required if *error_model* is *same*
              or *different*, specifies :math:`\\boldsymbol{\mathbf{a_{\epsilon,r}}}`

            - *rhor_prior_params* : only required if *error_model* is *different*,
              specifies :math:`\\boldsymbol{\mathbf{a_{\\rho,r}}}`

          Any values less than *PRIOR_MIN_VALUE* (a constant specified in this module)
          are automatically set to *PRIOR_MIN_VALUE*.

        * *seed* : integer seed for MCMC. Calling multiple times with same *seed*
          and data will give identical results provided architecture for floating
          point math on the computer is also the same.

        * *n_jobs* : The number of CPUs to use. If -1, use all available CPUs.

        * The following parameters specify how the MCMC tries to guarantee convergence.
          They all have reasonable default values, so you probably don't need to change them.
          The MCMC is considered to have converged if the mean Gelman-Rubin R
          statistic (http://www.jstor.org/stable/2246093) over all :math:`\pi_{r,x}`
          values is less than or equal to *r_max*, and the mean effective sample size
          is greater than or equal to *neff_min*. The MCMC first runs with *nchains*
          chains and *niter* iterations. If it fails to converge, it then increases
          the number of iterations by a factor of *increasefactor* and tries again,
          and repeats until it converges or until it has tried *increasetries*
          times to increase the number of iterations.

    The return value is: *(converged, pi_means, pi_95credint, logstring)*.
    The entries in this tuple have the following meanings:

        - *converged* is *True* if the MCMC converges and *False* otherwise.

        - *pi_means* is a dictionary keyed by all characters in *characterlist*
          with the value for each character giving the :math:`\pi_{r,a}` value.

        - *pi_95credint* is a dictionary keyed by all characters in *characterlist*
          with the value being a 2-tuple giving the lower and upper limits on the
          median-centered credible interval for :math:`\pi_{r,a}`.

        - *logstring* is a string that can be printed to give summary information
          about the MCMC run and its convergence.
    """
    logstring = ['\tBeginning MCMC at %s' % time.asctime()]
    numpy.random.seed(seed)
    assert nchains >= 2, "nchains must be at least two"
    assert niter >= 100, "niter must be at least 100"
    assert len(characterlist) == len(set(characterlist)), "characterlist contains a duplicate character:\n%s" % str(characterlist)
    assert wtchar in characterlist, "wtchar %s is not in characterlist %s" % (wtchar, str(characterlist))
    data = {'Nchar':len(characterlist), 'iwtchar':characterlist.index(wtchar) + 1}
    data['nrpre'] = [counts['nrpre'][char] for char in characterlist]
    data['nrpost'] = [counts['nrpost'][char] for char in characterlist]
    data['pir_prior_params'] = [max(PRIOR_MIN_VALUE, priors['pir_prior_params'][char]) for char in characterlist]
    data['mur_prior_params'] = [max(PRIOR_MIN_VALUE, priors['mur_prior_params'][char]) for char in characterlist]
    if error_model == 'none':
        sm = StanModelNoneErr().model
    elif isinstance(error_model, StanModelNoneErr):
        sm = error_model.model
        error_model = 'none'
    elif error_model == 'same' or isinstance(error_model, StanModelSameErr):
        if error_model == 'same':
            sm = StanModelSameErr().model
        else:
            sm = error_model.model
            error_model = 'same'
        data['nrerr'] = [counts['nrerr'][char] for char in characterlist]
        data['epsilonr_prior_params'] = [max(PRIOR_MIN_VALUE, priors['epsilonr_prior_params'][char]) for char in characterlist]
    elif error_model == 'different' or isinstance(error_model, StanModelDifferentErr):
        if error_model == 'different':
            sm = StanModelDifferentErr().model
        else:
            sm = error_model.model
            error_model = 'different'
        data['nrerrpre'] = [counts['nrerrpre'][char] for char in characterlist]
        data['nrerrpost'] = [counts['nrerrpost'][char] for char in characterlist]
        data['epsilonr_prior_params'] = [max(PRIOR_MIN_VALUE, priors['epsilonr_prior_params'][char]) for char in characterlist]
        data['rhor_prior_params'] = [max(PRIOR_MIN_VALUE, priors['rhor_prior_params'][char]) for char in characterlist]
    else:
        raise ValueError("Invalid error_model of %s" % error_model)
    ntry = 0
    while True: # run until converged or tries exhausted
        init = _InitialValuePreferences(error_model, nchains, characterlist.index(wtchar), len(characterlist))
        # this next command unfortunately creates a lot of output that I have not been able to re-direct...
        fit = sm.sampling(data=data, iter=niter, chains=nchains, seed=seed, n_jobs=n_jobs, refresh=-1, init=init)
        # extract output
        fitsummary = fit.summary()
        rownames = list(fitsummary['summary_rownames'])
        colnames = list(fitsummary['summary_colnames'])
        summary = fitsummary['summary']
        char_row_indices = dict([(char, rownames.index('pir[%d]' % characterlist.index(char))) for char in characterlist]) 
        rindex = colnames.index('Rhat')
        neffindex = colnames.index('n_eff')
        rlist = [summary[char_row_indices[char]][rindex] for char in characterlist]
        rhat_is_nan = [rhat for rhat in rlist if math.isnan(rhat)]
        rlist = [rhat for rhat in rlist if not math.isnan(rhat)]
        nefflist = [summary[char_row_indices[char]][neffindex] for char in characterlist]
        rmean = sum(rlist) / float(len(rlist))
        neffmean = sum(nefflist) / float(len(nefflist))
        logstring.append('\tAfter %d MCMC chains each of %d steps, mean R = %.2f and mean Neff = %d' % (nchains, niter, rmean, neffmean))
        if rhat_is_nan:
            logstring.append('\t\tThere are %d characters where R is nan' % len(rhat_is_nan))
        # allow convergence with extra stringent criteria when some sites have Rhat being nan
        # this is because pystan appears to sometimes give Rhat of nan for sites with low preference
        # turns out the pystan Rhat values of nan are a bug according to pystan developers 
        if ((not rhat_is_nan) and rmean <= r_max and neffmean >= neff_min) or (rhat_is_nan and len(rhat_is_nan) < 0.5 * len(characterlist) and rmean <= r_max and neffmean >= 3.0 * neff_min) or (rhat_is_nan and len(rhat_is_nan) <= 0.75 * len(characterlist) and rmean <= 1.0 + (r_max - 1.0) / 2.0 and neffmean >= 5.0 * neff_min):
            # converged
            logstring.append('\tMCMC is deemed to have converged at %s.' % time.asctime())
            meanindex = colnames.index('mean')
            lower95index = colnames.index('2.5%')
            upper95index = colnames.index('97.5%')
            pi_means = dict([(char, summary[char_row_indices[char]][meanindex]) for char in characterlist])
            pi_95credint = dict([(char, (summary[char_row_indices[char]][lower95index], summary[char_row_indices[char]][upper95index])) for char in characterlist])
            return (True, pi_means, pi_95credint, '\n'.join(logstring))
        else:
            # failed to converge
            if ntry < increasetries:
                ntry += 1
                niter = int(niter * increasefactor)
                logstring.append("\tMCMC failed to converge. Doing retry %d with %d iterations per chain." % (ntry, niter))
            else:
                with open('_no_converge_prefs_debugging.pickle', 'w') as f_debug:
                    cPickle.dump((counts, init, fitsummary), f_debug)
                logstring.append("\tMCMC FAILED to converge after all attempts at %s." % time.asctime())
                meanindex = colnames.index('mean')
                lower95index = colnames.index('2.5%')
                upper95index = colnames.index('97.5%')
                pi_means = dict([(char, summary[char_row_indices[char]][meanindex]) for char in characterlist])
                pi_95credint = dict([(char, (summary[char_row_indices[char]][lower95index], summary[char_row_indices[char]][upper95index])) for char in characterlist])
                return (False, pi_means, pi_95credint, '\n'.join(logstring))



def InferSiteDiffPrefs(characterlist, wtchar, error_model, counts, priors, seed=1, niter=10000, n_jobs=1, r_max=1.1, neff_min=100, nchains=4, increasefactor=2, increasetries=6):
    """Infers site-specific differntial preferences by MCMC for a specific site.

    The algorithm is explained in the documentation for the ``dms_tools`` package.

    The MCMC automatically tries to guarantee convergence via the parameters specified
    by *r_max*, *neff_min*, *nchains*, *niter*, *increasefactor*, and *increasetries*.

    Here are the calling variables:    

        * *characterlist* is a list of all valid characters; entries must be unique.
          Typically a list of all amino acids, codons, or nucleotides. Characters
          are case sensitive.

        * *wtchar* is a character that is in *characterlist* and defines the wildtype
          character at the site.

        * *error_model* specifies how the errors are estimated. It can
          be one of the following strings:

              - *none* : no errors 

              - *same* : same error rates for starting library and selections :math:`s1`
                and :math:`s2`.

          Alternatively, it could be an instance of *StanModelDiffPrefNoErr*
          or *StanModelDiffPrefSameErr*. These are the
          ``pystan`` models that implement the three strings described above.
          If you pass these model objects, then they will not need to be
          compiled by this function (strings require models to be compiled).
          So it is advantageous to pass the ``pystan`` models rather than
          the strings if calling this function repeatedly.

        * *counts* is a dictionary specifying the deep sequencing counts.
          Each string key should specify a dictionary keyed by all characters
          and with values giving the integer counts for that character. Keys:

            - *nrstart* : counts for starting library

            - *nrs1* : counts for selection :math:`s1`

            - *nrs2* : counts for selection :math:`s2`

            - *nrerr* : only required if *error_model* is
              *same*, specifis counts for error control.

        * *priors* is a dictionary specifying the parameter vectors for the
          Dirichlet priors. Each string key should specify a dictionary keyed by
          all characters and with values giving the prior for that character. Keys:

            - *pirs1_prior_params* : specifies Dirichlet prior vector for :math:`\pi_r^{s1}`

            - *frstart_prior_params* : specifies Dirichlet prior vector for 
              :math:`f_r^{\\textrm{start}}`

            - *deltapi_concentration* : specifies floating number giving concentration
              parameter for :math:`\Delta\pi_r`

            - *xir_prior_params* : only required if *error_model* is *same*,
              specified Dirichlet prior vector for error rates :math:`\\xi_r`

          Any values less than *PRIOR_MIN_VALUE* (a constant specified in this module)
          are automatically set to *PRIOR_MIN_VALUE*.

        * *seed* : integer seed for MCMC. Calling multiple times with same *seed*
          and data will give identical results provided architecture for floating
          point math on the computer is also the same.

        * *n_jobs* : The number of CPUs to use. If -1, use all available CPUs.

        * The following parameters specify how the MCMC tries to guarantee convergence.
          They all have reasonable default values, so you probably don't need to change them.
          The MCMC is considered to have converged if the mean Gelman-Rubin R
          statistic (http://www.jstor.org/stable/2246093) over all :math:`\pi_{r,x}`
          values is less than or equal to *r_max*, and the mean effective sample size
          is greater than or equal to *neff_min*. The MCMC first runs with *nchains*
          chains and *niter* iterations. If it fails to converge, it then increases
          the number of iterations by a factor of *increasefactor* and tries again,
          and repeats until it converges or until it has tried *increasetries*
          times to increase the number of iterations.

    The return value is: *(converged, deltapi_means, pr_deltapi_gt0, pr_deltapi_lt0, logstring)*.
    The entries in this tuple have the following meanings:

        - *converged* is *True* if the MCMC converges and *False* otherwise.

        - *deltapi_means* is a dictionary keyed by all characters in *characterlist*
          with the value for each character giving the :math:`\Delta\pi_{r,a}` value.

        - *pr_deltapi_gt0* is a dictionary keyed by all characters in *characterlist*
          with the value being the posterior probability that :math:`\Delta\pi_{r,a} > 0`.

        - *pr_deltapi_lt0* is a dictionary keyed by all characters in *characterlist*
          with the value being the posterior probability that :math:`\Delta\pi_{r,a} < 0`.

        - *logstring* is a string that can be printed to give summary information
          about the MCMC run and its convergence.
    """
    logstring = ['\tBeginning MCMC at %s' % time.asctime()]
    numpy.random.seed(seed)
    assert nchains >= 2, "nchains must be at least two"
    assert niter >= 100, "niter must be at least 100"
    assert len(characterlist) == len(set(characterlist)), "characterlist contains a duplicate character:\n%s" % str(characterlist)
    assert wtchar in characterlist, "wtchar %s is not in characterlist %s" % (wtchar, str(characterlist))
    data = {'Nchar':len(characterlist), 'iwtchar':characterlist.index(wtchar) + 1}
    data['nrstart'] = [counts['nrstart'][char] for char in characterlist]
    data['nrs1'] = [counts['nrs1'][char] for char in characterlist]
    data['nrs2'] = [counts['nrs2'][char] for char in characterlist]
    data['pirs1_prior_params'] = [max(PRIOR_MIN_VALUE, priors['pirs1_prior_params'][char]) for char in characterlist]
    data['frstart_prior_params'] = [max(PRIOR_MIN_VALUE, priors['frstart_prior_params'][char]) for char in characterlist]
    data['deltapi_concentration'] = max(PRIOR_MIN_VALUE, priors['deltapi_concentration'])
    if error_model == 'none':
        sm = StanModelDiffPrefNoErr().model
    elif isinstance(error_model, StanModelDiffPrefNoErr):
        sm = error_model.model
        error_model = 'none'
    elif error_model == 'same' or isinstance(error_model, StanModelDiffPrefSameErr):
        if error_model == 'same':
            sm = StanModelDiffPrefSameErr().model
        else:
            sm = error_model.model
            error_model = 'same'
        data['nrerr'] = [counts['nrerr'][char] for char in characterlist]
        data['xir_prior_params'] = [max(PRIOR_MIN_VALUE, priors['xir_prior_params'][char]) for char in characterlist]
    else:
        raise ValueError("Invalid error_model of %s" % error_model)
    ntry = 0
    while True: # run until converged or tries exhausted
        init = _InitialValueDiffPrefs(error_model, nchains, characterlist.index(wtchar), len(characterlist))
        # this next command unfortunately creates a lot of output that I have not been able to re-direct...
        fit = sm.sampling(data=data, iter=niter, chains=nchains, seed=seed, n_jobs=n_jobs, refresh=-1, init=init)
        # extract output
        fitsummary = fit.summary()
        rownames = list(fitsummary['summary_rownames'])
        colnames = list(fitsummary['summary_colnames'])
        summary = fitsummary['summary']
        char_row_indices = dict([(char, rownames.index('deltapir[%d]' % characterlist.index(char))) for char in characterlist]) 
        rindex = colnames.index('Rhat')
        neffindex = colnames.index('n_eff')
        rlist = [summary[char_row_indices[char]][rindex] for char in characterlist]
        rhat_is_nan = [rhat for rhat in rlist if math.isnan(rhat)]
        rlist = [rhat for rhat in rlist if not math.isnan(rhat)]
        nefflist = [summary[char_row_indices[char]][neffindex] for char in characterlist]
        rmean = sum(rlist) / float(len(rlist))
        neffmean = sum(nefflist) / float(len(nefflist))
        logstring.append('\tAfter %d MCMC chains each of %d steps, mean R = %.2f and mean Neff = %d' % (nchains, niter, rmean, neffmean))
        if rhat_is_nan:
            logstring.append('\t\tThere are %d characters where R is nan' % len(rhat_is_nan))
        # allow convergence with extra stringent criteria for a few sites with Rhat being nan
        # this is because pystan appears to sometimes give Rhat of nan for sites with low preference
        # turns out the pystan Rhat values of nan are a bug according to pystan developers 
        if ((not rhat_is_nan) and rmean <= r_max and neffmean >= neff_min) or (rhat_is_nan and len(rhat_is_nan) < 0.5 * len(characterlist) and rmean <= r_max and neffmean >= 3.0 * neff_min) or (rhat_is_nan and len(rhat_is_nan) <= 0.75 * len(characterlist) and rmean <= 1.0 + (r_max - 1.0) / 2.0 and neffmean >= 5.0 * neff_min):
            # converged
            logstring.append('\tMCMC is deemed to have converged at %s.' % time.asctime())
            meanindex = colnames.index('mean')
            deltapi_means = dict([(char, summary[char_row_indices[char]][meanindex]) for char in characterlist])
            deltapirsamples = fit.extract('deltapir')['deltapir'] 
            nsamples = niter * nchains // 2
            assert deltapirsamples.shape == (nsamples, len(characterlist)), "Unexpected shape of deltapirsamples. Got %s, expected (%d, %d)" % (str(deltapirsamples.shape), niter * nchains // 2, len(characterlist))
            pr_deltapi_gt0 = dict([(char, n / float(nsamples)) for (char, n) in zip(characterlist, sum(deltapirsamples > 0))])
            pr_deltapi_lt0 = dict([(char, n / float(nsamples)) for (char, n) in zip(characterlist, sum(deltapirsamples < 0))])
            return (True, deltapi_means, pr_deltapi_gt0, pr_deltapi_lt0, '\n'.join(logstring))
        else:
            # failed to converge
            if ntry < increasetries:
                ntry += 1
                niter = int(niter * increasefactor)
                logstring.append("\tMCMC failed to converge. Doing retry %d with %d iterations per chain." % (ntry, niter))
            else:
                with open('_no_converge_diffprefs_debugging.pickle', 'w') as f_debug:
                    cPickle.dump((counts, init, fitsummary), f_debug)
                logstring.append("\tMCMC FAILED to converge after all attempts at %s." % time.asctime())
                meanindex = colnames.index('mean')
                deltapi_means = dict([(char, summary[char_row_indices[char]][meanindex]) for char in characterlist])
                deltapirsamples = fit.extract('deltapir')['deltapir'] 
                nsamples = niter * nchains // 2
                assert deltapirsamples.shape == (nsamples, len(characterlist)), "Unexpected shape of deltapirsamples. Got %s, expected (%d, %d)" % (str(deltapirsamples.shape), niter * nchains // 2, len(characterlist))
                pr_deltapi_gt0 = dict([(char, n / float(nsamples)) for (char, n) in zip(characterlist, sum(deltapirsamples > 0))])
                pr_deltapi_lt0 = dict([(char, n / float(nsamples)) for (char, n) in zip(characterlist, sum(deltapirsamples < 0))])
                return (False, deltapi_means, pr_deltapi_gt0, pr_deltapi_lt0, '\n'.join(logstring))



if __name__ == '__main__':
    import doctest
    doctest.testmod()
