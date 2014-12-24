=================================================================================
Influenza NS1 deep mutational scanning by Wu et al (2014)
=================================================================================

.. contents::

Overview
------------

This is an analysis of data from the study by `Wu et al (2014)`_. In this study, a library of influenza viruses carrying random nucleotide variants of NS1 are grown, and then are passaged in the presence and absence of interferon. The former conditions leads to differential preference for mutations that help reduce interferon sensitivity.

In this analysis, the ``dms_inferdiffprefs`` program from `dms_tools`_ is used to analyze and visualize this differential selection. Note that this analysis is performed at the nucleotide level since that is how the mutant library was made.

Input data
--------------------------
The input data is found in the subdirectory ``./counts_data/``. There are three files here that were obtained from the first author of `Wu et al (2014)`_ that give the counts of mutations in input library, the control selected library, and the IFN-selected library:

    - ``IPT_NS1``

    - ``CON_NS1``

    - ``IFN_NS1``

Running the ``Python`` script ``convert_to_dms_counts.py`` with::

    python convert_to_dms_counts.py

converts the counts files into formats compatible with `dms_tools`_:

    - ``input.txt``

    - ``control.txt``

    - ``interferon.txt``

Differential preference analysis
---------------------------------
The shell script ``infer_diffprefs`` runs `dms_tools`_ to compute the differential preferences in the interferon-selected treatment versus the control treatment. This script runs ``dms_inferdiffprefs`` to compute these differential preferences (in ``diffprefs.txt``), and then uses ``dms_logoplot`` to make the visualization in ``diffprefs_logoplot.pdf``.

Run the shell script with::

    ./infer_diffprefs

to create this image:

.. figure:: diffprefs_logoplot.pdf
   :width: 90%
   :align: center
   :alt: diffprefs_logoplot.pdf

   The file ``diffprefs_logoplot.pdf`` shows the differential preference for each nucleotide mutation in the interferon-selected versus the control passage of the NS1 mutant libraries.


.. _`Wu et al (2014)`: http://jvi.asm.org/content/88/17/10157
.. _`dms_tools`: https://github.com/jbloom/dms_tools
