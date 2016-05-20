.. _dms_diffselection:

==========================================
``dms_diffselection``
==========================================

.. contents::

Overview
-------------
``dms_diffselection`` is a program included with the `dms_tools`_ package. 

It is designed to assess the differential selection on each mutation between a library that is mock treated and subjected to some selective condition.

Specifically, imagine that you have a deep mutational scanning library that is mock selected and selected. For instance, this might be a virus library that is passaged in the absence and presence of immune selection. For each site :math:`r` in the gene, let :math:`n_{r,x}^{\rm{mock}}` be the number of observed counts of character :math:`x` (might be amino acid, nucleotide, or codon) at site :math:`r` in the mock-selected condition, and let :math:`n_{r,x}^{\rm{selected}}` be the number of observed counts of :math:`x` at :math:`r` in the selected condition. 
Let :math:`\operatorname{wt}\left(r\right)` be denote the wildtype character at :math:`r`.
Then the relative enrichment of the mutant relative to the wildtype after selection is

.. math::
  
   E_{r,x} = \frac{\left(n_{r,x}^{\rm{selected}} + P\right) / \left(n_{r,\operatorname{wt}\left(r\right)}^{\rm{selected}} + P\right)}{\left(n_{r,x}^{\rm{mock}} + P\right) / \left(n_{r,\operatorname{wt}\left(r\right)}^{\rm{mock}} + P\right)}

where :math:`P > 0` is a pseudocount that is added to each observed count. 
Note that by definition, :math:`E_{r,\operatorname{wt}\left(r\right)}` is always one.

We then quantify the differential selection :math:`s_{r,x}` for :math:`x` at :math:`r` in the selected versus control condition as:

.. math::

   s_{r,x} = \log_2 E_{r,x}

It is these :math:`s_{r,x}` values that are reported by ``dms_diffselection``. Note that the choice of logarithm base is arbitrary; we have chosen 2 here. Note that by definition, :math:`s_{r,\operatorname{wt}\left(r\right)}` is always zero.

When :math:`s_{r,x} > 0`, then :math:`x` is more favored in the selection versus mock condition.
When :math:`s_{r,x} < 0`, then :math:`x` is less favored in the selection versus mock condition.

Larger values of the pseudocount :math:`P` will avoid spuriously estimating strong selection when in fact you just have a lot of statistical noise due to small counts.

After you install `dms_tools`_, this program will be available to run at the command line.

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: DiffSelectionParser
   :prog: dms_diffselection

   mockcounts
    This file gives the :math:`n_{r,x}^{\rm{mock}}` values described in the `Overview`_.
       
    Should be in the format of a :ref:`dms_counts` file.

   selectedcounts
    This file gives the :math:`n_{r,x}^{\rm{mock}}` values described in the `Overview`_.
       
    Should be in the format of a :ref:`dms_counts` file.

   outprefix
    For instance, if ``outprefix`` is ``antibodyselection_`` then the output files are ``antibodyselection_mutdiffsel.txt`` and ``antibodyselection_sitediffsel.txt``. See `Output`_ for an explanation of the contents of these files.

   pseudocount
    This is the :math:`P` parameter described in the `Overview`_.

   \-\-chartype
    A value of ``codon_to_aa`` means we have codon characters for the deep sequencing data, but infer selection for amino acids. Essentially, we sum the counts of all codons for each amino acid, and then do the calculation described in `Overview`_ on these amino-acid counts.

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


These output files are in a form where they can directly read into `pandas`_ data frames via that package's ``read_csv`` function.

.. include:: weblinks.txt