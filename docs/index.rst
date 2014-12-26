.. dms_tools documentation master file, created by
   sphinx-quickstart on Tue Oct 28 08:12:23 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for ``dms_tools``
=====================================

Overview
----------

`dms_tools`_ is a software package for analyzing and visualizing **d**\eep **m**\utational **s**\canning data.

Deep mutational scanning is a technique to profile the effects of all single mutations to a gene; see `Fowler and Fields (2014)`_ for a review. 

The :ref:`programs` that are installed as part of `dms_tools`_ allow you to easily analyze and visualize deep mutational scanning data; `dms_tools`_ also comes with a :ref:`pythonapi` that you can leverage to implement additional analyses.

One strength of `dms_tools`_ is that it analyzes deep mutational scanning data using a principled Bayesian framework that avoids the `statistical bias inherent in ratio estimation`_ (the most common strategy for analyzing deep mutational scanning data). A second strength of `dms_tools`_ is that it enables easy construction of logo plots (built via extension of `weblogo`_) that provide a compact and intuitive visualization of the results.

`dms_tools`_ is particularly suited for the following analyses:

    1) You have measured the effects of all single codon (or amino-acid or nucleotide) mutations to a gene, and wish to quantify the *preference* of each site for each possible identity. You can do this using :ref:`dms_inferprefs`, and visualize the results using :ref:`dms_logoplot`.

    2) You have subjected a library of mutant genes to selection under two different conditions, and wish to identify mutations that are favored under one condition versus the other. You can do this using :ref:`dms_inferdiffprefs`, and visualize the results using :ref:`dms_logoplot`.

    3) You have performed several `biological replicates`_ of a deep mutational scanning experiment (such replicates are always a good idea as they are the only sure way to quantify error), and wish to compare and combine the results. You can do this using :ref:`dms_correlate` and :ref:`dms_merge`.

`dms_tools`_ does **not** currently perform the actual alignment of sequencing reads and the counting of mutations; you must do that with some other program to create :ref:`dms_counts`\s that can be used as input to `dms_tools`_. 

Contents
------------------------

.. toctree:: 
   :maxdepth: 1

   installation
   programs
   fileformats
   examples
   algorithms
   pythonapi
   citations



Indices and tables
--------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. include:: weblinks.txt
