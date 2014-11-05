.. _dms_inferprefs:

==========================================
``dms_inferprefs``
==========================================

.. contents::

Overview
-------------
``dms_inferprefs`` is a program included with the `dms_tools`_ package. It infers the site-specific preferences :math:`\pi_{r,a}` for each character :math:`a` (an amino acid, codon, or nucleotide) at each site :math:`r`.

For a detailed description of the algorithm implemented by this program, see :ref:`inferprefs_algorithm`.

After you install `dms_tools`_, this program will be available to run at the command line.

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: InferPrefsParser
   :prog: dms_inferprefs

   n_pre
    This file gives the :math:`n_{r,x}^{\textrm{pre}}` values described in :ref:`inferprefs_algorithm`.
    
    See :ref:`dms_counts` for full specification of the format for this file.

   n_post
    This file gives the :math:`n_{r,x}^{\textrm{post}}` values described in :ref:`inferprefs_algorithm`.
    
    See :ref:`dms_counts` for full specification of the format for this file.

   outfile
    In addition to the inferred preferences :math:`\pi_{r,x}`, this file also gives the site entropies calculated from these preferences and the median-centered 95% credible intervals for each preference. 

    See :ref:`preferences_file` for full specification of the format for this file.

   \-\-errmodel
    Setting this option to "none" corresponds to setting to zero the error rates :math:`\epsilon_{r,x}` and :math:`\rho_{r,x}` described in :ref:`inferprefs_algorithm`.

    Setting this option to "same n_err" corresponds to setting :math:`\epsilon_{r,x} = \rho_{r,x}` for the error rates described in :ref:`inferprefs_algorithm`. In this case, the file specified by "nerr" should have the format specified by :ref:`dms_counts`.

    Setting this option to "same n_errpre n_errpost" correspdonds to allowing distinct values for the error rates :math:`\epsilon_{r,x}` and :math:`\rho_{r,x}` described in :ref:`inferprefs_algorithm`. In this case, the files specified by "n_errpre" and "n_errpost" (which given the :math:`n_{r,x}^{\textrm{err,pre}}` and :math:`n_{r,x}^{\textrm{err,post}}` values) should have the format specified by :ref:`dms_counts`.

   \-\-pi_alpha
    This is the parameter :math:`\alpha_{\pi}` described in :ref:`inferprefs_algorithm`. 
    
    Smaller values favor distributions of preferences that are more unequal among the amino acids at each site; larger values favor distributions of preferences more equal among the amino acids at each site.

   \-\-mu_alpha
    This is the parameter :math:`\alpha_{\mu}` described in :ref:`inferprefs_algorithm`.

   \-\-err_alpha
    This gives the parameters :math:`\alpha_{\epsilon}` and :math:`\alpha_{\rho}` described in :ref:`inferprefs_algorithm`; current these two parameters must be equal.

   \-\-chartype
    A value of "codon_to_aa" corresponds to using codon characters for the deep sequencing data, but inferring preferences for amino acids, as described in :ref:`chartype_codon_to_aa` This would be the sensible approach if you assume that selection is operating on the amino-acid level. Note that the priors over the mutagenesis rates assume that all codon mutations are made at equal frequencies (``NNN`` libraries).

    A value of "DNA" corresponds to using DNA characters for both the deep sequencing data and inferred preferences, as described in :ref:`chartype_DNA`

    A value of "codon" corresponds to using codon characters for both the deep sequencing data and inferred preferences, as described in :ref:`chartype_codon` You would prefer this option over "codon_to_aa" if you thought that there was different selection on different synonymous mutations. Note that the priors over the mutagenesis rates assume that all codon mutations are made at equal frequencies (``NNN`` libraries).

   \-\-includestop
    A value of "True" means that when using ``--chartype codon_to_aa``, we infer preferences for 21 amino acids, with stop codons (denoted by ``*``) one of the possibilities. A value of "False" means that we constrain the preference for stop codons to be zero regardless of whether or not there are counts for these codons in the data, and so only infer preferences for the 20 non-stop amino acids.

Output
----------
The inferred preferences are in the :ref:`preferences_file` specified by ``outfile``. A summary of the progress is in the file specified by ``logfile``. In addition, some information will be written to standard output. You should **not** be concerned if this information includes some warnings about convergence -- the program automatically tests for convergence as descrbed in :ref:`MCMC_inference`, and will stop with an error if the inference cannot converge. Unless that happens, the warnings don't indicate a problem, as the MCMC will continue to add steps until there is convergence.

Examples
-----------
Imagine we have data from codon-level mutagenesis of an influenza gene followed by functional sequencing and deep sequencing. We want to determine amino acid preferences. Let's say the pre-selection mutant library codon counts are in a file called ``mutDNA_codoncounts.txt`` and the post-selection codon counts are in a file called ``mutvirus_codoncounts.txt``. If we are **not** trying to account for error rates, we would use the command::

    dms_inferprefs mutDNA_codoncounts.txt mutvirus_codoncounts.txt amino_acid_preferences_no_errors.txt

The inferred preferences would then be in the created file ``amino_acid_preferences_no_errors.txt``. The log of the output would be in the file ``amino_acid_preferences_no_errors.log``.

The command above would include preferences for 21 different possible amino acids, as stop codons are considered a possible amino acid by default. If you do **not** want stop codons to be considered a possible amino acid (essentially constraining the preference for stop codons to zero), then use::

    dms_inferprefs --includestop False mutDNA_codoncounts.txt mutvirus_codoncounts.txt amino_acid_preferences_no_errors.txt

If we instead wanted to account for error rates, and we had sequenced unmutated plasmid to generate the codon counts file ``DNA_codoncounts.txt`` and had sequenced unmutated virus to generate the codon counts file ``virus_codoncounts.txt``, we would use the command::

    dms_inferprefs --errmodel different DNA_codoncounts.txt virus_codoncounts.txt mutDNA_codoncounts.txt mutvirus_codoncounts.txt amino_acid_preferences_corrected_for_errors.txt

If we ran that command and it was taking too long, we might want to use more CPUs if we have multi-core processors. If we want to use 12 CPUs, we would run::

    dms_inferprefs --ncpus 12 --errmodel different DNA_codoncounts.txt virus_codoncounts.txt mutDNA_codoncounts.txt mutvirus_codoncounts.txt amino_acid_preferences_no_errors.txt

`MCMC`_ is a stochastic algorithm. If we wanted to re-run it with a different random number seed to ensure we get the same results, we could use a different random number seed::

    dms_inferprefs --seed 2 --ncpus 12 --errmodel different DNA_codoncounts.txt virus_codoncounts.txt mutDNA_codoncounts.txt mutvirus_codoncounts.txt amino_acid_preferences_no_errors.txt

.. include:: weblinks.txt
