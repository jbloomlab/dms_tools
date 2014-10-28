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


import pystan


PRIOR_MIN_VALUE = 1.0e-7 # minimum value for Dirichlet prior elements


# Code and StanModel when the error rates are zero
pir_inference_no_err_code =\
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
pir_inference_no_err_sm = None # global variable compiled to StanModel when first needed

# Code and StanModel when same error rate for pre and post-selection libraries
pir_inference_same_err_code =\
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
pir_inference_same_err_sm = None # global variable compiled to StanModel when first needed

# Code and StanModel when different error rates for pre and post-selection libraries
pir_inference_different_err_code =\
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
pir_inference_different_err_sm = None # global variable compiled to StanModel when first needed


def InferSitePreferences(characterlist, wtchar, error_model, counts, priors, use_all_cpus=True, seed=1):
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

        * *use_all_cpus* : If *True*, we use all available CPUs. If *False*,
            just use one CPU.

        * *seed* : integer seed for MCMC.
    """
    niter = 1000
    nchains = 4
    assert len(characterlist) == len(set(characterlist)), "characterlist contains a duplicate character:\n%s" % str(characterlist)
    assert wtchar in characterlist, "wtchar %s is not in characterlist %s" % (wtchar, str(characterlist))
    if use_all_cpus:
        n_jobs = -1 # pystan code for use all CPUs
    else:
        n_jobs = 1 # pystan code for one CPUs
    data = {'Nchar':len(characterlist), 'iwtchar':characterlist.index(wtchar) + 1}
    data['nrpre'] = [counts['nrpre'][char] for char in characterlist]
    data['nrpost'] = [counts['nrpost'][char] for char in characterlist]
    data['pir_prior_params'] = [max(PRIOR_MIN_VALUE, priors['pir_prior_params'][char]) for char in characterlist]
    data['mur_prior_params'] = [max(PRIOR_MIN_VALUE, priors['mur_prior_params'][char]) for char in characterlist]
    if error_model == 'none':
        if pir_inference_no_err_sm == None: 
            # compile global StanModel variable if this has not already been done
            pir_inference_no_err_sm = pystan.StanModel(model_code=pir_inference_no_err_code)
        sm = pir_inference_no_err_sm
    elif error_model == 'same':
        data['nrerr'] = [counts['nrerr'][char] for char in characterlist]
        data['epsilonr_prior_params'] = [max(PRIOR_MIN_VALUE, priors['epsilonr_prior_params'][char]) for char in characterlist]
        if pir_inference_same_err_sm == None: 
            # compile global StanModel variable if this has not already been done
            pir_inference_same_err_sm = pystan.StanModel(model_code=pir_inference_same_err_code)
        sm = pir_inference_same_err_sm
    elif error_model == 'different':
        data['nrerrpre'] = [counts['nrerrpre'][char] for char in characterlist]
        data['nrerrpost'] = [counts['nrerrpost'][char] for char in characterlist]
        data['epsilonr_prior_params'] = [max(PRIOR_MIN_VALUE, priors['epsilonr_prior_params'][char]) for char in characterlist]
        data['rhor_prior_params'] = [max(PRIOR_MIN_VALUE, priors['rhor_prior_params'][char]) for char in characterlist]
        if pir_inference_different_err_sm == None: 
            # compile global StanModel variable if this has not already been done
            pir_inference_different_err_sm = pystan.StanModel(model_code=pir_inference_different_err_code)
        sm = pir_inference_different_err_sm
    else:
        raise ValueError("Invalid error_model of %s" % error_model)
    fit = sm.sampling(data=data, iter=niter, chains=nchains, seed=seed, n_jobs=n_jobs)
    print fit.summary()['summary']
    print fit
        


