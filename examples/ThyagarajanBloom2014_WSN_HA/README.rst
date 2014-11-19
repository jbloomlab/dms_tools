==========================================================================
`Thyagarajan and Bloom (2014)`_'s scanning of WSN influenza HA
==========================================================================

.. contents::

Overview
----------
This directory contains a re-analysis of the deep mutational scanning data of the hemagglutinin (HA) from the WSN strain of influenza as described by `Thyagarajan and Bloom (2014)`_.

That data was originally analyzed with `mapmuts`_ by `Thyagarajan and Bloom (2014)`_. That analysis is available in the `mapmuts example`_ and is described in the `mapmuts example documentation`_.

In the original analysis, `mapmuts`_ was used both to align the deep sequencing data to count codon mutations at each site, and then to infer the site-specific amino-acid preferences. Here we use the same codon counts from the `mapmuts`_ alignment, but re-infer the site-specific amino-acid preferences using `dms_tools`_. This allows us to compare the preferences from the two programs.

Input files
------------
Here are the files used for this analysis:

* The subdirectory ``./codoncounts/`` has the codon counts files from the previous `mapmuts`_ analysis of the deep mutational scanning data. These are for the three replicates of the experiment (*replicate_1*, *replicate_2*, *replicate_3*) for the *mutDNA* (mutant library), *mutvirus* (selected mutant viruses), *DNA* (error control for *mutDNA*) and *virus* (error control for *mutvirus*) samples. These are the files with the same names available in the various subdirectories of the `mapmuts example`_.


.. _`Thyagarajan and Bloom (2014)`: http://elifesciences.org/content/3/e03300
.. _`mapmuts`: http://jbloom.github.io/mapmuts/
.. _`mapmuts example`: https://github.com/jbloom/mapmuts/tree/master/examples/WSN_HA_2014Analysis
.. _`mapmuts example documentation`: http://jbloom.github.io/mapmuts/example_WSN_HA_2014Analysis.html
.. _`dms_tools`: https://github.com/jbloom/dms_tools
