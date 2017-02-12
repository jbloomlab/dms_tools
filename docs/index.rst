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

`dms_tools`_ is suited for the following analyses:

    1) If you have sequenced a gene using barcoded subamplicon sequencing, you can process the FASTQ files to count the mutations using :ref:`dms_barcodedsubamplicons`. If you are using subassembly, you can process the FASTQ files using :ref:`dms_subassemble` and then count the frequencies of variants / mutations using :ref:`dms_matchsubassembledbarcodes`. In either case, you can summarize the results for multiple samples with :ref:`dms_summarizealignments`.

    2) You have subjected a library of mutant genes to selection under two conditions, and wish to identify mutations that are favored under one versus the other. You can do this using :ref:`dms_diffselection` or :ref:`dms_inferdiffprefs`, and visualize the results using :ref:`dms_logoplot`.

    3) You have measured the effects of all single codon (or amino-acid or nucleotide) mutations to a gene, and wish to quantify the *preference* of each site for each identity. You can do this using :ref:`dms_inferprefs`, and visualize the results using :ref:`dms_logoplot`.

    4) You have subjected a library of mutant genes to selection under two conditions, and wish to identify mutations that are favored under one versus the other. You can do this using :ref:`dms_inferdiffprefs`, and visualize the results using :ref:`dms_logoplot`.

    5) You have performed several `biological replicates`_ of deep mutational scanning, and wish to compare and combine the results. You can do this using :ref:`dms_correlate` and :ref:`dms_merge`.

The `dms_tools source code`_ is freely available on GitHub; however, you will probably have an easier time you just follow the :ref:`installation` instructions rather than building from this source.

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
