.. dms_tools documentation master file, created by
   sphinx-quickstart on Tue Oct 28 08:12:23 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for ``dms_tools``
=====================================

Quick overview
-----------------

`dms_tools`_ is a software package for analyzing and visualizing **d**\eep **m**\utational **s**\canning data.

The :ref:`programs` installed as part of `dms_tools`_ allow you to easily analyze and visualize deep mutational scanning data; `dms_tools`_ also comes with a :ref:`pythonapi`.

One strength of `dms_tools`_ is that it analyzes deep mutational scanning data using a principled Bayesian framework that avoids the `statistical bias inherent in ratio estimation`_ is that it enables easy construction of logo plots (built via extension of `weblogo`_) that provide a compact and intuitive visualization of the results.

`dms_tools`_ is suited for the following analyses:

    1) You have measured the effects of all single codon (or amino-acid or nucleotide) mutations to a gene, and wish to quantify the *preference* of each site for each possible identity. You can do this using :ref:`dms_inferprefs`, and visualize the results using :ref:`dms_logoplot`.

    2) You have subjected a library of mutant genes to selection under two different conditions, and wish to identify mutations that are favored under one condition versus the other. You can do this using :ref:`dms_inferdiffprefs`, and visualize the results using :ref:`dms_logoplot`.

    3) You have performed several `biological replicates`_ of a deep mutational scanning experiment, and wish to compare and combine the results. You can do this using :ref:`dms_correlate` and :ref:`dms_merge`.

    4) You have performed deep mutational scanning by sequencing barcoded subamplicons of your gene of interest, and wish to process the FASTQ files to count mutations. You can do this using :ref:`dms_barcodedsubamplicons`.

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
