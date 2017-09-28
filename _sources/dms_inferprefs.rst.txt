.. _dms_inferprefs:

==========================================
``dms_inferprefs``
==========================================

.. contents::

Overview
-------------
``dms_inferprefs`` is a program included with the `dms_tools`_ package. It infers the site-specific preferences :math:`\pi_{r,x}` for each character :math:`x` (an amino acid, codon, or nucleotide) at each site :math:`r`. The inference is done by `MCMC`_ using `PyStan`_. 

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

   \-\-errpre
    Setting this option to "none" corresponds to setting to zero the error rates :math:`\epsilon_{r,x}` described in :ref:`inferprefs_algorithm`.

    If you set this option to the name of a file, that file is assumed to give the counts for sequencing the pre-selection error control (:math:`n_{r,x}^{\textrm{err,pre}}`) and is used to estimate the error rate :math:`\epsilon_{r,x}` described in :ref:`inferprefs_algorithm`. This file should have the format specified by :ref:`dms_counts`.

   \-\-errpost
    Setting this option to "none" corresponds to setting to zero the error rates :math:`\rho_{r,x}` described in :ref:`inferprefs_algorithm`.

    If you set this option to the name of a file, that file is assumed to give the counts for sequencing the post-selection error control (:math:`n_{r,x}^{\textrm{err,post}}`) and is used to estimate the error rate :math:`\rho_{r,x}` described in :ref:`inferprefs_algorithm`. This file should have the format specified by :ref:`dms_counts`.

   \-\-pi_alpha
    This is the parameter :math:`\alpha_{\pi}` described in :ref:`inferprefs_algorithm`. 
    
    Smaller values favor distributions of preferences that are more unequal among the possible characters at each site; larger values favor distributions of preferences more equal among the characters at each site.

   \-\-mu_alpha
    This is the parameter :math:`\alpha_{\mu}` described in :ref:`inferprefs_algorithm`.

   \-\-err_alpha
    This gives the parameters :math:`\alpha_{\epsilon}` and :math:`\alpha_{\rho}` described in :ref:`inferprefs_algorithm`; currently these two parameters must be equal.

   \-\-chartype
    Note that the **default** value is "codon_to_aa".

    A value of "codon_to_aa" corresponds to using codon characters for the deep sequencing data, but inferring preferences for amino acids, as described in :ref:`chartype_codon_to_aa` This would be the sensible approach if you assume that selection is operating on the amino-acid level. Note that the priors over the mutagenesis rates assume that all codon mutations are made at equal frequencies (``NNN`` libraries).

    A value of "DNA" corresponds to using DNA characters for both the deep sequencing data and inferred preferences, as described in :ref:`chartype_DNA`

    A value of "codon" corresponds to using codon characters for both the deep sequencing data and inferred preferences, as described in :ref:`chartype_codon` You would prefer this option over "codon_to_aa" if you thought that there was different selection on different synonymous mutations. Note that the priors over the mutagenesis rates assume that all codon mutations are made at equal frequencies (``NNN`` libraries).

    A value of "aa" corresponds to using amino-acid characters for both the deep sequence data and inferred preferences, as described in :ref:`chartype_aa`. You should **only** use this option if your data make it impossible to use "codon_to_aa".  Note that the priors for this option assume that all amino-acid mutations are made at equal frequencies.

   \-\-excludestop
    By default, when using ``--chartype codon_to_aa`` or ``--chartype aa``, we infer preferences for 21 amino acids, with stop codons (denoted by ``*``) one of the possibilities. If you specify the ``--excludestop`` option, then we constrain the preference for stop codons to be zero regardless of whether or not there are counts for these codons in the data, and so only infer preferences for the 20 non-stop amino acids.

   \-\-ratio_estimation
    A suggested value for *mincounts* is one. In general, larger values of *mincounts* conceptually corresponds to a stronger prior favoring all preferences being equal.    
    
    Preferences are computed using ratio estimation as follows. For each site, we set the enrichment ratio of the wildtype character :math:`\rm{wt}` equal to one
    
    .. math::
       
        \phi_{wt} = 1
    
    Next, for each non-wildtype character :math:`x`, we calculate the enrichment ratio relative to :math:`\rm{wt}` as

    .. math::
       
        \phi_{x} = 
        \frac{
          \left(\frac{
            \max\left(\frac{\mathcal{P}}{N_r^{\rm{post}}}, \frac{n_{r,x}^{\rm{post}}}{N_r^{\rm{post}}} - \frac{n_{r,x}^{\rm{errpost}}}{N_r^{\rm{errpost}}}\right)
          }{
            \frac{n_{r,\rm{wt}}^{\rm{post}}}{N_r^{\rm{post}}} + \delta - \frac{n_{r,\rm{wt}}^{\rm{errpost}}}{N_r^{\rm{errpost}}}
          }\right)
        }
        {
          \left(\frac{
            \max\left(\frac{\mathcal{P}}{N_r^{\rm{pre}}}, \frac{n_{r,x}^{\rm{pre}}}{N_r^{\rm{pre}}} - \frac{n_{r,x}^{\rm{errpre}}}{N_r^{\rm{errpre}}}\right)
          }{
            \frac{n_{r,\rm{wt}}^{\rm{pre}}}{N_r^{\rm{pre}}} + \delta - \frac{n_{r,\rm{wt}}^{\rm{errpre}}}{N_r^{\rm{errpre}}}
          }\right)
        }

    where :math:`\mathcal{P}` is the value of *mincounts*. When *error_model* is *none*, then all terms involving the error corrections (with superscript *err*) are ignored and :math:`\delta` is set to zero; otherwise :math:`\delta` is one.

    Next, for each character :math:`x`, including :math:`\rm{wt}`, we calculate the preference for :math:`x` as

    .. math::

        \pi_x = \frac{\phi_x}{\sum_y \phi_y}.


Output
----------
The inferred preferences are in the :ref:`preferences_file` specified by ``outfile``. A summary of the progress is in the file specified by ``logfile``. In addition, some information will be written to standard output. `Do not worry about Informational Messages`_. 

Because standard output will include lots of informational messages that obscure what is happening, it is most informative to track the progress by watching ``logfile``.

Do not worry about Informational Messages
---------------------------------------------------
When you run ``dms_inferprefs``, you may get a lot of ``Informational Message`` output on the screen. These messages look like this::

    Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
    stan::prob::dirichlet_log(N4stan5agrad3varE): prior sample sizes[1]  is 0:0, but must be > 0!
    If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
    but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

You do **not** need to worry about these messages. They are automatically produced by `PyStan`_. Although it seems like there are a lot of these messages, running ``dms_inferprefs`` on a 500-residue gene will invoke at least :math:`4 \times 500 = 2000` different `MCMC`_ chains, all of which may sporadically produce this warning.

``dms_inferprefs`` automatically tests for convergence and described in `Runtime and convergence`_ and :ref:`MCMC_inference` If the program actually fails to converge for any site, you will not get an ``outfile`` and the program will terminate with an error message.

Runtime and convergence
---------------------------
Because ``dms_inferprefs`` uses `MCMC`_ to infer the preferences as described in :ref:`inferprefs_algorithm`, it may take the program a while to run. The program will run fastest if you use multiple CPUs (such as by using the option ``--ncpus -1`` to use all available CPUs). It should take anywhere from a few minutes to a few hours to complete on a multi-processor machine, depending on the size of the protein, which characters are being used (nucleotides are fastest, codons are slowest), and whether error controls are being included in the analysis (via ``--errpre`` and ``--errpost``).

The program automatically checks for `MCMC`_ convergence using the criteria described in :ref:`MCMC_inference` The program will terminate with an error if it fails to converge; otherwise it converged properly.

If the program is taking too long, you can also just use the ``--ratio_estimation`` option to speed things up dramatically. This method is somewhat less accurate when the sequencing counts are fairly low (see `Software for the analysis and visualization of deep mutational scanning data (Bloom, 2015)`_), however it still should give fairly good results.

Examples
-----------
Imagine we have data from codon-level mutagenesis of an influenza gene followed by functional sequencing and deep sequencing. We want to determine amino acid preferences. Let's say the pre-selection mutant library codon counts are in a file called ``mutDNA_codoncounts.txt`` and the post-selection codon counts are in a file called ``mutvirus_codoncounts.txt``. If we are **not** trying to account for error rates, we would use the command::

    dms_inferprefs mutDNA_codoncounts.txt mutvirus_codoncounts.txt amino_acid_preferences_no_errors.txt

The inferred preferences would then be in the created file ``amino_acid_preferences_no_errors.txt``. The log of the output would be in the file ``amino_acid_preferences_no_errors.log``.

The command above would include preferences for 21 different possible amino acids, as stop codons are considered a possible amino acid by default. If you do **not** want stop codons to be considered a possible amino acid (essentially constraining the preference for stop codons to zero), then use::

    dms_inferprefs --excludestop mutDNA_codoncounts.txt mutvirus_codoncounts.txt amino_acid_preferences_no_errors.txt

If we instead wanted to account for error rates, and we had sequenced unmutated plasmid to generate the codon counts file ``DNA_codoncounts.txt`` and had sequenced unmutated virus to generate the codon counts file ``virus_codoncounts.txt``, we would use the command::

    dms_inferprefs --errpre DNA_codoncounts.txt --errpost virus_codoncounts.txt mutDNA_codoncounts.txt mutvirus_codoncounts.txt amino_acid_preferences_corrected_for_errors.txt

If we ran that command and it was taking too long, we might want to use more CPUs if we have multi-core processors. If we want to use 12 CPUs, we would run::

    dms_inferprefs --ncpus 12 --errpre DNA_codoncounts.txt --errpost virus_codoncounts.txt mutDNA_codoncounts.txt mutvirus_codoncounts.txt amino_acid_preferences_no_errors.txt

`MCMC`_ is a stochastic algorithm. If we wanted to re-run it with a different random number seed to ensure we get the same results, we could use a different random number seed::

    dms_inferprefs --seed 2 --ncpus 12 --errpre DNA_codoncounts.txt --errpost virus_codoncounts.txt mutDNA_codoncounts.txt mutvirus_codoncounts.txt amino_acid_preferences_no_errors.txt

.. include:: weblinks.txt
