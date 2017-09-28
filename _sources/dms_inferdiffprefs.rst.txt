.. _dms_inferdiffprefs:

==========================================
``dms_inferdiffprefs``
==========================================

.. contents::

Overview
-------------
``dms_inferdiffprefs`` is a program included with the `dms_tools`_ package. It infers the differential preference :math:`\Delta\pi_{r,x}` for each character :math:`x` (an amino acid, codon, or nucleotide) at each site :math:`r` in an alternative selection versus a control selection. Values of :math:`\Delta\pi_{r,x} > 0` indicate that :math:`x` is favored in the alternative selection versus the control selection at site :math:`r`.

The inference is done by `MCMC`_ using `PyStan`_. 

For a detailed description of the algorithm implemented by this program, see :ref:`inferdiffprefs_algorithm`.

After you install `dms_tools`_, this program will be available to run at the command line.

Should you use this program?
-----------------------------
Note that we are no longer certain if the method implemented in this program is superior to the much simpler :ref:`dms_diffselection`. So consider trying that program as well.

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: InferDiffPrefsParser
   :prog: dms_inferdiffprefs

   n_start
    This file gives the :math:`n_{r,x}^{\textrm{start}}` values described in :ref:`inferdiffprefs_algorithm`.
    
    See :ref:`dms_counts` for full specification of the format for this file.

   n_s1
    This file gives the :math:`n_{r,x}^{s1}` values described in :ref:`inferdiffprefs_algorithm`. These are for selection :math:`s1`, which is taken as the "control" selection.
    
    See :ref:`dms_counts` for full specification of the format for this file.

   n_s2
    This file gives the :math:`n_{r,x}^{s2}` values described in :ref:`inferdiffprefs_algorithm`. These are for selection :math:`s2`, which is taken as the "alternative" selection to the "control" selection :math:`s1`.
    
    See :ref:`dms_counts` for full specification of the format for this file.

   outfile
    In addition to the inferred differential preferences :math:`\Delta\pi_{r,x}` for each character :math:`x` in selection :math:`s2` versus selection :math:`s1`, this file also gives the root-mean-square (RMS) differential preference for each site, and the posterior probabilities that these preferences are greater than or less than zero.
    
    Note that the preferences are for :math:`s2` relative to :math:`s1`, so :math:`\Delta\pi_{r,x} > 0` means that :math:`x` is preferred more at site :math:`r` in selection :math:`s2` than selection :math:`s1`.

    See :ref:`diffpreferences_file` for full specification of the format for this file.

   \-\-err
    Setting this option to "none" corresponds to setting to zero the error rates :math:`\xi_{r,x}` described in :ref:`inferdiffprefs_algorithm`.

    If you set this option to the name of a file, that file is assumed to give the counts for sequencing the error control (:math:`n_{r,x}^{\textrm{err}}`) and is used to estimate the error rate :math:`\xi_{r,x}` described in :ref:`inferdiffprefs_algorithm`. This file should have the format specified by :ref:`dms_counts`.

   \-\-alpha_start
    This is the parameter :math:`\alpha_{\textrm{start}}` described in :ref:`inferdiffprefs_algorithm`. 
    
    Smaller values favor distributions of starting frequencies that are more unequal among the possible characters at each site; larger values favor distributions of frequencies more equal among the characters at each site.

   \-\-alpha_pis1
    This is the parameter :math:`\alpha_{\pi_{s1}}` described in :ref:`inferdiffprefs_algorithm`.

    Smaller values favor distributions of preferences that are more unequal among the possible characters at each site; larger values favor distributions of preferences more equal among the characters at each site.


   \-\-alpha_err
    This gives the parameter :math:`\alpha_{\xi}` described in :ref:`inferdiffprefs_algorithm`.

    Smaller values favor distributions of error rates that are more unequal among the possible characters at each site; larger values favor distributions of errors more equal among the characters at each site.

   \-\-alpha_deltapi  
    This gives the parameter :math:`\alpha_{\Delta\pi}` described in :ref:`\inferdiffprefs_algorithm`.

    Larger values correspond to a stronger prior expectation that the differential preference is zero. You might therefore consider using a value :math:`> 1` in order to enforce a prior expectation that there is no difference between the selections.

   \-\-chartype
    Note that the **default** value is "codon_to_aa".

    A value of "codon_to_aa" corresponds to using codon characters for the deep sequencing data, but inferring differential preferences for amino acids, as described in :ref:`chartype_codon_to_aa` This would be the sensible approach if you assume that selection is operating on the amino-acid level. Note that the priors over the mutagenesis rates assume that all codon mutations are made at equal frequencies (``NNN`` libraries).

    A value of "DNA" corresponds to using DNA characters for both the deep sequencing data and inferred preferences, as described in :ref:`chartype_DNA`

    A value of "codon" corresponds to using codon characters for both the deep sequencing data and inferred preferences, as described in :ref:`chartype_codon` You would prefer this option over "codon_to_aa" if you thought that there was different selection on different synonymous mutations. Note that the priors over the mutagenesis rates assume that all codon mutations are made at equal frequencies (``NNN`` libraries).

    A value of "aa" corresponds to using amino-acid characters for both the deep sequence data and inferred preferences,
    as described in :ref:`chartype_aa`. You should **only** use this option if your data make it impossible to use "codon_to_aa". Note that the priors for this option assume that all amino-acid mutations are made at equal frequencies.

   \-\-excludestop
    By default, when using ``--chartype codon_to_aa`` or ``--chartype aa``, we infer differential preferences for 21 amino acids, with stop codons (denoted by ``*``) one of the possibilities. If you specify the ``--excludestop`` option, then we constrain the differential preference for stop codons to be zero regardless of whether or not there are counts for these codons in the data, and so only infer differential preferences for the 20 non-stop amino acids.

   \-\-ratio_estimation
    A suggested value for *mincounts* is one. In general, larger values of *mincounts* conceptually corresponds to a stronger prior favoring all differential preferences being zero.
    
    Differential preferences are computed using ratio estimation as follows. For each site, we set the enrichment ratio of the wildtype character :math:`\rm{wt}` equal to one. We do this for each of the two selections :math:`s1` and :math:`s2`:
    
    .. math::
       
        \phi_{wt}^{\rm{s1}} = \phi_{wt}^{\rm{s2}} = 1
    
    Next, for each non-wildtype character :math:`x`, we calculate the enrichment ratio relative to :math:`\rm{wt}` for :math:`s1` and :math:`s2` as:

    .. math::
       
        \phi_{x}^{s1} = 
        \frac{
          \left(\frac{
            \max\left(\frac{\mathcal{P}}{N_r^{\rm{s1}}}, \frac{n_{r,x}^{\rm{s1}}}{N_r^{\rm{s1}}} - \frac{n_{r,x}^{\rm{err}}}{N_r^{\rm{err}}}\right)
          }{
            \frac{n_{r,\rm{wt}}^{\rm{s1}}}{N_r^{\rm{s1}}} + \delta - \frac{n_{r,\rm{wt}}^{\rm{err}}}{N_r^{\rm{err}}}
          }\right)
        }
        {
          \left(\frac{
            \max\left(\frac{\mathcal{P}}{N_r^{\rm{pre}}}, \frac{n_{r,x}^{\rm{pre}}}{N_r^{\rm{pre}}} - \frac{n_{r,x}^{\rm{err}}}{N_r^{\rm{err}}}\right)
          }{
            \frac{n_{r,\rm{wt}}^{\rm{pre}}}{N_r^{\rm{pre}}} + \delta - \frac{n_{r,\rm{wt}}^{\rm{err}}}{N_r^{\rm{err}}}
          }\right)
        }

    and
    
    .. math::
       
        \phi_{x}^{s2} = 
        \frac{
          \left(\frac{
            \max\left(\frac{\mathcal{P}}{N_r^{\rm{s2}}}, \frac{n_{r,x}^{\rm{s2}}}{N_r^{\rm{s2}}} - \frac{n_{r,x}^{\rm{err}}}{N_r^{\rm{err}}}\right)
          }{
            \frac{n_{r,\rm{wt}}^{\rm{s2}}}{N_r^{\rm{s2}}} + \delta - \frac{n_{r,\rm{wt}}^{\rm{err}}}{N_r^{\rm{err}}}
          }\right)
        }
        {
          \left(\frac{
            \max\left(\frac{\mathcal{P}}{N_r^{\rm{pre}}}, \frac{n_{r,x}^{\rm{pre}}}{N_r^{\rm{pre}}} - \frac{n_{r,x}^{\rm{err}}}{N_r^{\rm{err}}}\right)
          }{
            \frac{n_{r,\rm{wt}}^{\rm{pre}}}{N_r^{\rm{pre}}} + \delta - \frac{n_{r,\rm{wt}}^{\rm{err}}}{N_r^{\rm{err}}}
          }\right)
        }

    where :math:`\mathcal{P}` is the value of *mincounts*. When *error_model* is *none*, then all terms involving the error corrections (with superscript *err*) are ignored and :math:`\delta` is set to zero; otherwise :math:`\delta` is one.

    Next, for each character :math:`x`, including :math:`\rm{wt}`, we compute the preferences
    as normalized enrichment ratios for both :math:`s1` and :math:`s2` as:

        .. math::

            \pi_{x}^{s1} = \frac{\phi_x^{s1}}{\sum_y \phi_y^{s1}}

    and

        .. math::

            \pi_{x}^{s2} = \frac{\phi_x^{s2}}{\sum_y \phi_y^{s2}}

    and finally calculate the differential preference as 

        .. math::

            \Delta\pi_{x} = \pi_{x}^{s2} - \pi_{x}^{s1}
    

Output
----------
The inferred differential preferences are in the :ref:`diffpreferences_file` specified by ``outfile``. A summary of the progress is in the file specified by ``logfile``. In addition, some information will be written to standard output. `Do not worry about Informational Messages`_.

Becase standard output will include lots of informational messages that obscure what is happening, it will be more informative to track the progress by watching ``logfile``.

Do not worry about Informational Messages
---------------------------------------------------
When you run ``dms_inferdiffprefs``, you may get a lot of ``Informational Message`` output on the screen. These messages look like this::

    Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
    stan::prob::dirichlet_log(N4stan5agrad3varE): prior sample sizes[1]  is 0:0, but must be > 0!
    If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
    but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

You do **not** need to worry about these messages. They are automatically produced by `PyStan`_. Although it seems like there are a lot of these messages, running ``dms_inferdiffprefs`` on a 500-residue gene will invoke at least :math:`4 \times 500 = 2000` different `MCMC`_ chains, all of which may sporadically produce this warning.

``dms_inferdiffprefs`` automatically tests for convergence and described in `Runtime and convergence`_ and :ref:`MCMC_inference` If the program actually fails to converge for any site, you will not get an ``outfile`` and the program will terminate with an error message.

Runtime and convergence
---------------------------
Because ``dms_inferdiffprefs`` uses `MCMC`_ to infer the preferences as described in :ref:`inferdiffprefs_algorithm`, it may take the program a while to run. The program will run fastest if you use multiple CPUs (such as by using the option ``--ncpus -1`` to use all available CPUs). It should take anywhere from a few minutes to a few hours to complete on a multi-processor machine, depending on the size of the protein, which characters are being used (nucleotides are fastest, codons are slowest), and whether error controls are being included in the analysis (via ``--err``).

The program automatically checks for `MCMC`_ convergence using the criteria described in :ref:`MCMC_inference` The program will terminate with an error if it fails to converge; otherwise it converged properly.

You can make the program run dramatically faster by using the ``--ratio_estimation`` option. This option forgoes the Bayesian estimation and regularizing prior over the differential preferences, but will still be reasonably accurate if you have a fairly large number of counts.

Examples
-----------
Imagine we have created a library of codon mutants of an influenza gene and selected for functional variants. We deep sequence this starting pool of functional variants and put the codon counts in a file called ``mutvirus_codoncounts.txt``. We subject the starting pool of functional mutant variants to immune selection by passaging in the presence of antibody, and also do a control passage in the absence of immune selection. The codon counts from deep sequencing these pools are in the files ``immuneselected_codoncounts.txt`` and ``controlselected_codoncounts.txt``, respectively. We also do control sequencing of unmutated virus to quantify the error rates, and put the codon counts in the file ``unmutatedcontrol_codoncounts.txt``. We could then compute the differential preference :math:`\Delta\pi_{r,a}` for each amino acid :math:`a` at each site :math:`r` in the immune selected viruses versus the controls using::

    dms_inferdiffprefs --err unmutatedcontrol_codoncounts.txt mutvirus_codoncounts.txt controlselected_codoncounts.txt immuneselected_codoncounts.txt differentialpreferences.txt

The differential preferences are then in the created file called ``differentialpreferences.txt``, with :math:`\Delta\pi_{r,a} > 0` indicating that amino acid :math:`a` is preferentially favored at :math:`r` in the immune selected versus the control selected sample. A log file summarizing output will be written to ``differentialpreferences.log``.

.. include:: weblinks.txt
