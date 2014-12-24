==========================================================================
Influenza WSN HA deep mutational scanning by Thyagarajan and Bloom (2014)
==========================================================================

.. contents::

Overview
----------
This directory contains a re-analysis of the deep mutational scanning data of the hemagglutinin (HA) from the WSN strain of influenza as described by `Thyagarajan and Bloom (2014)`_.

That data was originally analyzed with `mapmuts`_ (see the `mapmuts example`_). Here we re-analyze the data using `dms_tools`_ to infer the site-specific preference for each amino acid at each site.

Input files
------------
The subdirectory ``./dms_counts/`` has the codon counts files from the previous `mapmuts`_ analysis of the deep mutational scanning data. These are for the three replicates of the experiment (*replicate_1*, *replicate_2*, *replicate_3*) for the *mutDNA* (mutant library), *mutvirus* (selected mutant viruses), *DNA* (error control for *mutDNA*) and *virus* (error control for *mutvirus*) samples. These are the files with the same names available in the various subdirectories of the `mapmuts example`_.

The subdirectory ``./PDB_structure/`` has the relative solvent accessibilities (``RSAs.txt``) and secondary structures (``SSs.txt``) as computed from the PDB structure `1RVX`_.

Analysis of site-specific preferences
---------------------------------------
The shell script ``infer_prefs`` performs the analysis. It can be run with::

    ./infer_prefs

This ``infer_prefs`` script performs the following steps:

1) It uses the `dms_tools`_ program ``dms_inferprefs`` to infer the site-specific preferences from each of the three experiment replicates into the files ``./inferred_preferences/prefs_replicate_1.txt``, etc.

2) It uses the `dms_tools`_ program ``dms_correlate`` to compute the pairwise correlations between the preferences for the three replicates, creating the plots ``./correlations/corr_1_vs_2.pdf``, etc. These plots are shown below.

3) It uses the `dms_tools`_ program ``dms_merge`` to average the preferences for the three replicates to create the file ``./inferred_preferences/prefs_avg.txt``.

4) It uses the `dms_tools`_ program ``dms_logoplot`` to create a visual display of the average site-specific preferences along with overlays showing the secondary structure and relative solvent accessibility. This plot is in the file ``prefs_logoplot.pdf`` shown below.

.. figure:: ./correlations/corr_1_vs_2.pdf
   :align: center
   :width: 20%
   :alt: ./correlations/corr_1_vs_2.pdf

   Correlation of preferences from replicates 1 and 2 (``./correlations/corr_1_vs_2.pdf``).

.. figure:: ./correlations/corr_1_vs_3.pdf
   :align: center
   :width: 20%
   :alt: ./correlations/corr_1_vs_3.pdf

   Correlation of preferences from replicates 1 and 3 (``./correlations/corr_1_vs_3.pdf``).

.. figure:: ./correlations/corr_2_vs_3.pdf
   :align: center
   :width: 20%
   :alt: ./correlations/corr_2_vs_3.pdf

   Correlation of preferences from replicates 2 and 3 (``./correlations/corr_2_vs_3.pdf``).

.. figure:: prefs_logoplot.pdf
   :align: center
   :width: 90%
   :alt: prefs_logoplot.pdf

   Logo plot showing the site-specific preferences (``prefs_logoplot.pdf``).

.. _`Thyagarajan and Bloom (2014)`: http://elifesciences.org/content/3/e03300
.. _`mapmuts`: http://jbloom.github.io/mapmuts/
.. _`mapmuts example`: https://github.com/jbloom/mapmuts/tree/master/examples/WSN_HA_2014Analysis
.. _`mapmuts example documentation`: http://jbloom.github.io/mapmuts/example_WSN_HA_2014Analysis.html
.. _`dms_tools`: https://github.com/jbloom/dms_tools
.. _`1RVX`: http://www.rcsb.org/pdb/explore.do?structureId=1rvx
