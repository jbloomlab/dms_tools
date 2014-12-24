=================================================================================
Tn5 deep mutational scanning by Melnikov et al (2014)
=================================================================================

.. contents::

Overview
------------
`Melnikov et al (2014)`_ perform deep mutational scanning on a Tn5 transposon. In this directory, `dms_tools`_ is used to analyze their data. 

This analysis demonstrates the use of `dms_tools`_ to infer and visualize site-specific preferences.

The analysis then tests the performance of `dms_tools`_ preference and differential inferences on data simulated based on data from `Melnikov et al (2014)`_. 

Input data
--------------------------
The following input data files are used by the scripts in this directory to generate all output:

    * ``./raw_data/`` : this directory contains the raw deep mutational scanning counts data from `Melnikov et al (2014)`_ as provided in Supplementary Data 2 of that paper (that supplementary data is a ZIP file, this directory simply contains the contents of the unzipped file).

    * ``./PDB_structure/`` contains information about the relative solvent accessibilities and secondary structures of the Tn5 transposon based on an analysis of the PDB structure with ``DSSP``. Specifically:

        - ``SSs.txt`` classifies the secondary structures.

        - ``RSAs.txt`` classifies the relative solvent accessibilities.

Inference of site-specific amino-acid preferences
-------------------------------------------------------------------------
Within this directory, the shell script ``infer_prefs`` infers the site-specific preferences using `dms_tools`_. The script can simply be run at the command line with::

    ./infer_prefs

The script performs the following steps:

    1) Use the custom ``Python`` script ``process_raw_to_dms_counts.py`` to convert the raw data files of `Melnikov et al (2014)`_ found in ``./raw_data/`` to deep mutational scanning counts files in `dms_tools`_ format. The counts files are in the directory ``./processed_data/``, and are for replicates 1 and 2 of selection on kanamycin.

    2) Use the `dms_tools`_ program ``dms_inferprefs`` to infer the site-specific preferences for each replicate of the experiment. The results are in ``./processed_data/prefs_1.txt`` and ``./processed_data/prefs_2.txt``.

    3) Use the `dms_tools`_ program ``dms_merge`` to average the site-specific preferences from the two replicates to create ``./processed_data/prefs.txt``.

    4) Use the `dms_tools`_ program ``dms_logoplot`` to visually display the site-specific preferences with overlays showing secondary structure and relative solvent accessibility. This final visualization of the preferences is in the file ``prefs_logoplot.pdf``.

.. figure:: prefs_logoplot.pdf
   :width: 75%
   :align: center
   :alt: prefs_logoplot.pdf

   The ``prefs_logoplot.pdf`` visualization of the site-specific preferences from `Melnikov et al (2014)`_ for selection on kanamycin.

Testing of inference of site-specific preferences on simulated data
----------------------------------------------------------------------
The subdirectory ``./infer_prefs_on_simulated_data/`` simulates a deep mutational scanning experiment using the actual preferences for Tn5 inferred above. The deep mutational scanning experiment is simulated at a variety of per-site sequencing depths, and both with and without sequencing errors. The simulation is designed to be "realistic" in terms of the parameters. The simulated deep mutational scanning counts files are then analyzed to see how well the inferred preferences match those actually used in the simulation. The inferences are done in two ways:

    1) Using the statistical model implemented in `dms_tools`_ via the ``dms_inferprefs`` program.

    2) Simply by calculating enrichment ratios for mutations, and then rescaling these to preferences.

The analysis is run by the script ``test_inferprefs.py`` via::

    cd infer_prefs_on_simulated_data
    python test_inferprefs.py

The results are shown in the plot below. Briefly, inference with ``dms_inferprefs`` is superior to simply calculating enrichment ratios. The accuracy increases with increasing sequencing depth. Slightly more data is needed for accurate inferences when there are sequencing errors.

.. figure:: infer_prefs_on_simulated_data/correlations.pdf
   :align: center
   :width: 80%
   :alt: infer_prefs_on_simulated_data/correlations.pdf

   Accuracy of the inference of site-specific preferences on simulated data.

Testing of inference of differential preferences on simulated data
----------------------------------------------------------------------
The subdirectory ``./infer_diffprefs_on_simulated_data/`` simulates passaging an already functionally selected mutant library in two conditions: a control condition, and a condition that induces change in amino-acid preferences (differential preferences) at 20% of the residues. The `dms_tools`_ program ``dms_inferdiffprefs`` is then used to analyze the simulated data to see how accurately it is possible to infer the differential preferences. This analysis can be run via::

    cd infer_diffprefs_on_simulated_data
    python test_inferdiffprefs.py

The results are shown below. With increasing read depth, the inferred differential preferences converge to the true values.

.. figure:: infer_diffprefs_on_simulated_data/correlations.pdf
   :align: center
   :width: 80% 
   :alt: infer_diffprefs_on_simulated_data/correlations.pdf

   Accuracy of inference of differential preferences on simulated data.

.. figure:: infer_diffprefs_on_simulated_data/inferred_diffprefs_logoplot.pdf
   :alt: infer_diffprefs_on_simulated_data/inferred_diffprefs_logoplot.pdf
   :width: 70%
   :align: center

   Inferred differential preferences visualized with ``dms_logoplot`` at a read depth of :math:`10^7`. Overlay bars show sites that actually had non-zero differential preferences in simulation.

.. _`Melnikov et al (2014)`: http://nar.oxfordjournals.org/content/42/14/e112
.. _`dms_tools`: https://github.com/jbloom/dms_tools
