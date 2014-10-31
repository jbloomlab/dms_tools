"""
========================
``mcmc`` module
========================

Module that performs MCMC to infer preferences and differential preferences.

MCMC done using ``Stan`` (http://mc-stan.org/) via ``pystan``.

The first time some of the functions in this module are run, the ``pystan`` code
will be compiled, which will take some time.

Written by Jesse Bloom.

Functions in this module
----------------------------------

* *InferSitePreferences* : use MCMC to infer site-specific preferences.


Function documentation
----------------------------------

"""


import sys
import tempfile
import numpy
import numpy.random
import pystan


PRIOR_MIN_VALUE = 1.0e-7 # minimum value for Dirichlet prior elements


# StanModel when the error rates are zero
class _StanModelNoneErr:
    pystancode =\
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
    stanmodel = None
    @classmethod
    def Model(cls):
        if not cls.stanmodel:
            cls.stanmodel = pystan.StanModel(model_code=cls.pystancode)
        return cls.stanmodel

# StanModel when same error rate for pre and post-selection libraries
class _StanModelSameErr:
    pystancode =\
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
    stanmodel = None
    @classmethod
    def Model(cls):
        if not cls.stanmodel:
            cls.stanmodel = pystan.StanModel(model_code=cls.pystancode)
        return cls.stanmodel

# StanModel when different error rates for pre and post-selection libraries
class _StanModelDifferentErr:
    pystancode =\
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
    stanmodel = None
    @classmethod
    def Model(cls):
        if not cls.stanmodel:
            cls.stanmodel = pystan.StanModel(model_code=cls.pystancode)
        return cls.stanmodel


def _InitialValuePreferences(error_model, nchains, iwtchar, nchars):
    """Gets valid initial values for pir_inference_same_err_sm.

    The values automatically initialized by stypan frequently have invalid *mur_plus_err*
    or *fr_plus_err* due to negative values from the *mur*, *pir*, *epsilonr*, *rhor*
    choices. This function will generate random values that are valid
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
        raise ValueError("Invalid error_model of %s" % error_model)
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



def InferSitePreferences(characterlist, wtchar, error_model, counts, priors, n_jobs=1, seed=1, r_max=1.05, neff_min=200, nchains=3, niter=2000, increasefactor=2, increasetries=6):
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

        * *error_model* is a string specifying how the errors are estimated:

              - *none* : no errors (:math:`\\boldsymbol{\mathbf{\epsilon_r}} = \\boldsymbol{\mathbf{\\rho_r}} = \mathbf{0}`)

              - *same* : same error rates for pre and post-selection 
                (:math:`\\boldsymbol{\mathbf{\epsilon_r}} = \\boldsymbol{\mathbf{\\rho_r}}`).

              - *different* : different error rates for pre and post
                selection (:math:`\\boldsymbol{\mathbf{\epsilon_r}} \\ne \\boldsymbol{\mathbf{\\rho_r}}`).

        * *counts* is a dictionary specifying the deep sequencing counts.
          Each string key should specify a dictionary keyed by all characters
          and with values giving the integer counts for that character. Keys:

            - *npre* : specifies :math:`\\boldsymbol{\mathbf{n_r^{\\rm{pre}}}}`

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

        * *n_jobs* : The number of CPUs to use. If -1, use all available CPUs.

        * *seed* : integer seed for MCMC. Calling multiple times with same *seed*
          and data will give identical results provided architecture for floating
          point math on the computer is also the same.

        * The following parameters specify how the MCMC tries to guarantee convergence.
          They all have reasonable default values, so you probably don't need to change them.
          The MCMC is considered to have converged if for all :math:`\pi_{r,a}` values,
          both the Gelman-Rubin R statistic (http://www.jstor.org/stable/2246093)
          has a value less than or equal to *r_max* and the effective sample size
          is greater than or equal to *neff_min*. The MCMC first runs with *nchains*
          chains and *niter* iterations. If it fails to converge, it then increases
          the number of iterations by a factor of *increasefactor* and tries again,
          and repeats this process until it converges or until it has tried *increasetries*
          times to increase the number of iterations.

    The return values are as follows:

        * If the MCMC **fails** to converge, then the return value is *False*.

        * If the MCMC successfully converges, then the return value is the 2-tuple
          *(pi_means, pi_95credint)* where

            - *pi_means* is a dictionary keyed by all characters in *characterlist*
              with the value for each character giving the :math:`\pi_{r,a}` value.

            - *pi_95credint* is a dictionary keyed by all characters in *characterlist*
              with the value being a 2-tuple giving the lower and upper limits on the
              median-centered credible interval for :math:`\pi_{r,a}`.
    """
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
        sm = _StanModelNoneErr.Model()
    elif error_model == 'same':
        sm = _StanModelSameErr.Model()
        data['nrerr'] = [counts['nrerr'][char] for char in characterlist]
        data['epsilonr_prior_params'] = [max(PRIOR_MIN_VALUE, priors['epsilonr_prior_params'][char]) for char in characterlist]
    elif error_model == 'different':
        sm = _StanModelDifferentErr.Model()
        data['nrerrpre'] = [counts['nrerrpre'][char] for char in characterlist]
        data['nrerrpost'] = [counts['nrerrpost'][char] for char in characterlist]
        data['epsilonr_prior_params'] = [max(PRIOR_MIN_VALUE, priors['epsilonr_prior_params'][char]) for char in characterlist]
        data['rhor_prior_params'] = [max(PRIOR_MIN_VALUE, priors['rhor_prior_params'][char]) for char in characterlist]
    else:
        raise ValueError("Invalid error_model of %s" % error_model)
    ntry = 0
    while True: # run until converged or tries exhausted
        # run MCMC
        # this next command unfortunately creates a lot of output that I have not been able to re-direct...
        init = _InitialValuePreferences(error_model, nchains, characterlist.index(wtchar), len(characterlist))
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
        nefflist = [summary[char_row_indices[char]][neffindex] for char in characterlist]
        if (max(rlist) <= r_max) and (min(nefflist) >= neff_min):
            # converged
            #print "DEBUGGING: Converged with %d iterations, max R = %g and min Neff = %g" % (niter, max(rlist), min(nefflist))
            meanindex = colnames.index('mean')
            lower95index = colnames.index('2.5%')
            upper95index = colnames.index('97.5%')
            pi_means = dict([(char, summary[char_row_indices[char]][meanindex]) for char in characterlist])
            pi_95credint = dict([(char, (summary[char_row_indices[char]][lower95index], summary[char_row_indices[char]][upper95index])) for char in characterlist])
            return (pi_means, pi_95credint)
        else:
            # failed to converge
            #print "DEBUGGING: Failed to converge with %d iterations, max R = %g and min Neff = %g" % (niter, max(rlist), min(nefflist))
            if ntry < increasetries:
                ntry += 1
                niter = int(niter * increasefactor)
            else:
                return False # failed to converge after trying step increases

