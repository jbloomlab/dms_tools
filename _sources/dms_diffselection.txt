.. _dms_diffselection:

==========================================
``dms_diffselection``
==========================================

.. contents::

Overview
-------------
``dms_diffselection`` is a program included with the `dms_tools`_ package. 

It is designed to assess the differential selection on each mutation between a library that is mock-treated and subjected to some selective condition.

After you install `dms_tools`_, this program will be available to run at the command line.

Differential selection
++++++++++++++++++++++++++++++++++++++++

Say you have a deep mutational scanning library that is mock-treated  and selected. For instance, this might be a virus library that is passaged in the absence and presence of immune selection. For each site :math:`r` in the gene, let :math:`n_{r,x}^{\rm{mock}}` be the number of observed counts of character :math:`x` (might be amino acid, nucleotide, or codon) at site :math:`r` in the mock-treated condition, and let :math:`n_{r,x}^{\rm{selected}}` be the number of observed counts of :math:`x` at :math:`r` in the selected condition. 
Let :math:`\operatorname{wt}\left(r\right)` be denote the wildtype character at :math:`r`.
Then the relative enrichment of the mutant relative to the wildtype after selection is

.. math::
   :label: E_rx
  
   E_{r,x} = \frac{\left(n_{r,x}^{\rm{selected}} + f_{r, \rm{selected}} \times P\right) / \left(n_{r,\operatorname{wt}\left(r\right)}^{\rm{selected}} + f_{r, \rm{selected}} \times P\right)}{\left(n_{r,x}^{\rm{mock}} + f_{r, \rm{mock}} \times P\right) / \left(n_{r,\operatorname{wt}\left(r\right)}^{\rm{mock}} + f_{r, \rm{mock}} \times P\right)}

where :math:`P > 0` is a pseudocount that is added to each observed count (specified by ``--pseudocount`` option), and :math:`f_{r, \rm{selected}}` and :math:`f_{r, \rm{mock}}` are defined so that the pseudocount is scaled up for the library (*mock* or *selected*) with higher depth at site :math:`r`:

.. math::
   :label: f_rselected

   f_{r, \rm{selected}} = \max\left[1, \left(\sum_x n_{r,x}^{\rm{selected}}\right) / \left(\sum_x n_{r,x}^{\rm{mock}}\right)\right]

.. math::
   :label: f_rmock

   f_{r, \rm{mock}} = \max\left[1, \left(\sum_x n_{r,x}^{\rm{mock}}\right) / \left(\sum_x n_{r,x}^{\rm{selected}}\right)\right].

The reason for scaling the pseudocount by library depth is that the *mock* and *selected* libraries may often be sequenced at different depths. In that case, if the same pseudocount is added to both, then estimates of :math:`E_{r,x}` will tend to be systematically different than one even if the relative counts for the wildtype and mutant amino acid are the same in both two conditions. Scaling the pseudocounts by the ratio of depths fixes this problem. If you do not want to do this scaling, see the ``--no-scale-pseudocounts`` option.

Note that by definition, :math:`E_{r,\operatorname{wt}\left(r\right)}` is always one.

We then quantify the differential selection :math:`s_{r,x}` for :math:`x` at :math:`r` in the selected versus control condition as:

.. math::
   :label: s_rx

   s_{r,x} = \log_2 E_{r,x}

It is these :math:`s_{r,x}` values that are reported by ``dms_diffselection``. Note that the choice of logarithm base is arbitrary; we have chosen 2 here. Note that by definition, :math:`s_{r,\operatorname{wt}\left(r\right)}` is always zero.

When :math:`s_{r,x} > 0`, then :math:`x` is more favored in the selection versus mock condition.
When :math:`s_{r,x} < 0`, then :math:`x` is less favored in the selection versus mock condition.

Larger values of the pseudocount :math:`P` will avoid spuriously estimating strong selection when in fact you just have a lot of statistical noise due to small counts.

Error correction
++++++++++++++++++++++++
You can optionally correct for the potential inflation of some counts by sequencing errors by using the ``--errorcontrolcounts`` option to ``dms_diffselection``.

In this case, the counts defined in Equation :eq:`E_rx` are adjusted as follows. 
Let :math:`n_{r,x}` be the counts for character :math:`x` at site :math:`r` in the *mock* or *selected* sample (these are the numbers in the ``mockcounts`` or ``selectecounts`` files). 
Let :math:`n^{\rm{err}}_{r,x}` be the number of counts of :math:`x` at site :math:`r` in the error control specified by ``--errorcontrolcounts``.
Define 

.. math::
   :label: epsilon

   \epsilon_{r,x} = \left(n^{\rm{err}}_{r,x}\right) / \left(\sum_y n^{\rm{err}}_{r,y}\right).

When :math:`x \ne \operatorname{wt}\left(r\right)` then :math:`\epsilon_{r,x}` is the rate of errors to :math:`x` at site :math:`r`, and when :math:`x = \operatorname{wt}\left(r\right)` then :math:`\epsilon_{r,x}` is one minus the rate of errors away from the wildtype at site :math:`r`.

We then adjust observed counts :math:`n_{r,x}` to the error-corrected counts :math:`\hat{n}_{r,x}` by

.. math::
   :label: n_rx_adjusted

   \hat{n}_{r,x} = \begin{cases}
   \max\left[\left(\sum_y n_{r,y}\right) \left(\frac{n_{r,x}}{\sum_y n_{r,y}} - \epsilon_{r,x}\right), 0\right] & \mbox{if } x \ne \operatorname{wt}\left(r\right) \\
   n_{r,x} / \epsilon_{r,x} & \mbox{if } x = \operatorname{wt}\left(r\right) \\
   \end{cases}

These adjusted counts are then used in place of the un-adjusted counts in Equation :eq:`E_rx` above.

Note that if you are using ``--chartype`` of ``codon_to_aa``, the corrections are done on the codon counts, and these corrected which are then aggregated into the amino acid counts before computing the differential selection.

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: DiffSelectionParser
   :prog: dms_diffselection

   mockcounts
    This file gives the :math:`n_{r,x}^{\rm{mock}}` values described in `Differential selection`_.
       
    Should be in the format of a :ref:`dms_counts` file.

   selectedcounts
    This file gives the :math:`n_{r,x}^{\rm{mock}}` values described in `Differential selection`_.
       
    Should be in the format of a :ref:`dms_counts` file.

   outprefix
    For instance, if ``outprefix`` is ``antibodyselection_`` then the output files are ``antibodyselection_mutdiffsel.txt`` and ``antibodyselection_sitediffsel.txt``. See `Output`_ for an explanation of the contents of these files.

   \-\-errorcontrolcounts
    You can use this option if you have estimated your sequencing error rate by sequencing a "wildtype" sample where you expect all errors come from sequencing (or other related processes such as PCR or reverse transcription). The counts will then be adjusted to correct for the error rates estimated from this control.

    This file gives the :math:`n_{r,x}^{\rm{err}}` values described in `Error correction`_.

    Should be in the format of a :ref:`dms_counts` file.

   \-\-pseudocount
    This is the :math:`P` parameter described in the Equation :eq:`E_rx`. Note that by default this pseudocount is for the library with greater depth at that site, and the pseudocount for the other library is scaled up by the ratio of the relative depths as described by Equations :eq:`f_rselected` and :eq:`f_rmock`.

   \-\-mincounts
    If you set this option to anything other than its default of 0, then you will get differential selection (:math:`s_{r,x}` values) of ``NaN`` (not a number) for mutations that have very low counts both before and after selection. In some cases, this may be preferable than estimating a value for a mutation with very little data.

    Note that if you keep ``--mincounts`` at its default of 0 but have a reasonable value of ``--pseudocount``, then mutations with few or no counts will have differential selection values of 0 estimated for them.

   \-\-chartype
    A value of ``codon_to_aa`` means we have codon characters for the deep sequencing data, but infer selection for amino acids. Essentially, we sum the counts of all codons for each amino acid, and then do the calculation described in `Overview`_ on these amino-acid counts.

   \-\-no-scale-pseudocounts
    If you use this option, the scaling :math:`f_{r, \rm{selected}}` and :math:`f_{r, \rm{mock}}` in Equations :eq:`f_rselected` and :eq:`f_rmock` are always set to one. Note that this can give biased estimates if the depth is unequal in the two libraries.

   \-\-includestop
    If this is true, treat stop codons (denoted by ``*``) as a possible amino acid. Otherwise simply ignore any counts for stop codons.

Example usage
--------------
Let's say that ``mock_counts.txt`` has the codon counts in a mock-selected library, and ``antibody_counts`` has these same counts in a library selected with an antibody. To determine selection on the amino acids, you could then run::

    dms_diffselection mock_counts.txt antibody_counts.txt antibodyselection_

The created files would be ``antibodyselection_mutdiffsel.txt`` and ``antibodyselection_sitediffsel.txt``.

Output
----------
Running ``dms_diffselection`` will cause some descriptive output to be printed to standard output.

In addition, two files with the prefix ``outprefix`` are created (these files are overwritten if they already exist). These files are:


 * The ``*mutdiffsel.txt`` file gives the :math:`s_{r,x}` values in a CSV file. These are sorted by :math:`s_{r,x}` value. The site number and wildtype character are also indicated. Here are the first few lines of an example file::

    site,wt,mut,diffsel
    156,G,S,7.15912988065
    146,N,D,6.85510232461
    157,K,S,6.78492241017
    159,S,D,6.61445818152
    153,S,K,6.55644791926
    153,S,I,6.45786904581
    158,S,A,6.31170638356

   In these files, the entry for :math:`x` equal to the wildtype character is always ``NaN`` (not a number). So for instance, given the data above, we would also have entries like::

    156,G,G,NaN
    146,N,N,NaN

   In addition, if ``--mincounts`` is set to something other than zero, then there may be mutations for which an :math:`s_{r,x}` value cannot be computed. In that case, that value is also reported as ``NaN`` (not a number) in this output file.

 * ``*sitediffsel.txt`` file summarizes total differential selection on each site. It provides three summary statistics:

        - *abs_diffsel* is the sum of the absolute values of the selection on all mutations at the site, and so is a measure of the total selection at this site. It is computed as :math:`\sum_x \left|s_{r,x}\right|`.

        - *positive_diffsel* is a sum of the values of positively selected mutations at a site, and so is a measure of the total selection for favorable mutations at this site. It is computed as :math:`\sum_x \max\left(0, s_{r,x}\right)`.

        - *negative_diffsel* is a sum of the values of negatively selected mutations at a site, and so is a measure of the total selection against mutations at this site. It is computed as :math:`\sum_x \min\left(0, s_{r,x}\right)`.

   The file gives all three of these quantities for each site, sorted by *abs_diffsel*. Here are the first few lines of an example file::

    site,abs_diffsel,positive_diffsel,negative_diffsel
    157,92.5132004464,92.5132004464,0.0
    153,85.7023127899,84.7159706476,-0.986342142299
    136,56.4005369298,56.3664732104,-0.0340637194023
    158,56.3773123754,56.3773123754,0.0
    156,56.0690809464,56.0690809464,0.0
    175,31.899455678,3.79784216691,-28.1016135111

   Note that if ``--mincounts`` is not 0, then some of the :math:`s_{r,x}` values may be ``NaN`` (not a number). In that case, any characters :math:`x` that do not have a value for the :math:`s_{r,x}` value do **not** contribute to the sums used to define *abs_diffsel*, *positive_diffsel*, and *negative_diffsel* above (in other words, they are counted as if the differential selection for character :math:`x` is 0). If **all** characters have a differential selection of ``NaN``, then the three summed values shown in this output file are ``NaN`` as well.

These output files are in a form where they can directly read into `pandas`_ data frames via that package's ``read_csv`` function.

You can visualize the ``*mutdiffsel.txt`` files with :ref:`dms_logoplot`.

.. include:: weblinks.txt
