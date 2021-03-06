Changelog
===========

1.2.2
---------
* Minor fixes to math formatting problems in docs

* Fixed bug caused by very high frequency mutations

1.2.1
---------
* Added ``merge_method`` of ``median`` to ``dms_merge``.

* Make compatible with ``weblogo`` version 3.5.

* Cleaned up some plot formatting

1.2.0
-------------
* Added ``dms_subassemble`` program

* Added ``dms_matchsubassembledbarcodes`` program

* Augmented ``dms_summarizealignments`` to handle output from ``dms_subassemble`` and ``dms_matchsubassembledbarcodes``

* Deleted loop counters that were previously still defined as package-level variables

* Added ``*singlemuttypes.pdf`` plot to ``dms_summarizealignments`` for *barcodedsubamplicons*.

* Added ``overlay_cmap`` option to ``dms_logoplot``

* Updated docs about installation

1.1.20
------------
* Added ``--mapmetric`` option to ``dms_logoplot`` to allow mapping of molecular weight, charge, or functional group to logoplot letter color.

* Added ``--markersize`` option to ``dms_correlate`` and turned clipping off to allow larger markers.

1.1.19
------------
* Some tweaks to ``dms_correlate`` to clarify parser and improve appearance of correlation plots.

* Added ``--restrictdiffsel`` option to ``dms_logoplot`` to only plot positive or negative site differential selection.

* Added ``--nosepline`` option to ``dms_logoplot`` to plot differential selection without black line separating positive and negative values.

1.1.18
----------
* Added ``--errorcontrolcounts`` to ``dms_diffselection``.

* Updated download URL in ``setup.py``.

* Some updates to documentation.

1.1.17
----------
* Print versions of all external packages (some were overlooked before)

* Make ``dms_merge`` handle the ``*mutdiffsel.txt`` files created by ``dms_diffselection``.

* Added ``translate_codon_to_aa`` and ``return_as_df`` options to dms_tools.file_io.ReadDMSCounts().

* Added ``dms_tools.utils.ParseNSMutFreqBySite()``.

* Make ``dms_correlate`` handle the ``mutdiffsel.txt`` and ``sitediffsel.txt``  files created by ``dms_diffselection``.

1.1.16
--------
* Added ``--mincounts`` option to ``dms_diffselection``, and make ``dms_logoplot`` properly plot differential selection with missing values.

* Added test (``./tests/test_diffselection.py``) for ``dms_diffselection``

* ``dms_barcodedsubamplicons`` now gives exactly reproducible output regardless of order of input reads. 

1.1.15
-------
* In ``dms_diffselection``, the pseudocounts are now scaled by the **smaller** depth library at each site.

1.1.14
------
* Added ``dms_diffselection`` program.

* Added ``--colormap`` and ``--overlay3`` options to ``dms_logoplot``

* Make ``dms_merge`` remove sites indicated as *None* when in *renumber* mode, added ``--skipfirstline`` option, and made commas an allowable separator.

* Increased number of colors allowed by *PlotMutCountFracs* (used by ``dms_summarizealignments``)

1.1.13
-----------
* Moved repository and docs from https://github.com/jbloom/dms_tools and http://jbloom.github.io/dms_tools to https://github.com/jbloomlab/dms_tools and http://jbloomlab.github.io/dms_tools

* Catch problems with ``dms_merge`` if multiple input files of the same name.

1.1.12
-------------
* Added option to sum count files in ``dms_merge``, with or without normalization of counts at each site. 

* Updated *dms_tools.file_io.IteratePairedFASTQ* to permit parsing of FASTQ files generated by ``fastq-dump``. 

1.1.11
-------------
* Fixed bug in ``dms_merge``, where normalization of prefs or diffprefs was not occurring when the user provided two or more identical pref or diffpref files

* Added y-axis label override option in ``plot.PlotDepth()`` and corrected ``dms_summarizealignments`` to label the y-axis as 'number of barcodes' when summarizing barcoded subamplicon data.

* Added ``merge_method`` choice ``rescale`` to ``dms_merge`` to rescale a preferences file by a ``--stringencyparameter``.

1.1.9
-------------
* Added option to trim reads R1 and R2 in ``dms_barcodedsubamplicons``

1.1.8
--------
* Fixed algorithm to estimate preferences and differential preferences using the ``--ratio_estimation`` option for both ``dms_inferprefs`` and ``dms_inferdiffprefs``.

* Fixed small bug in *dms_tools.file_io.WriteDiffPrefs* and in test suite.

1.1.7
--------
* ``dms_logoplot`` can now handle sites numbers that don't go sequentially starting with one or zero.

1.1.6
-------
* Added ``--stringencyparameter`` option to ``dms_logoplot``.

1.1.5
---------
* Added additional options for logo plots showing preferences and differential preferences 

1.1.4
----------
* Some minor documentation updates

* Added ``--ratio_estimation`` option to ``dms_inferdiffprefs``

1.1.3
--------
* Some minor documentation corrections; added BMC Bioinformatics (2015) citation

* Fixed minor bug in error handling by ``dms_editsites`` when there is an invalid file name.

* Added more detailed error message for invalid ``ydatamax`` in ``dms_logoplot``

1.1.2
--------
* Fixed bug in ``dms_summarizealignments`` processing of ``summarystats.txt`` files when ``--purgefracs`` is nonzero.

* Added ``--writemutfreqs`` option to ``dms_summarizealignments``

1.1.1
-------
* Fixed duplication of logging output by ``dms_inferprefs``

* Fixed bug that caused a crash during calculation of average mutation rate by ``dms_inferprefs`` when some sites have no counts data

* Fixed problems in MCMC convergence in ``dms_inferprefs`` and ``dms_inferdiffprefs`` when there are many Rhat values of ``nan`` by being more accommodating on the Rhat criterion if N_eff is sufficiently large

* Added test for ``dms_barcodedsubamplicons`` in ``./tests/``

* Some minor documentation updates

1.1.0
------
* Added ``dms_barcodedsubamplicons`` and ``dms_summarizealignments`` programs

1.0.1
--------
* Fixed bug in parsing *codon* option in ``dms_infeprefs`` and ``dms_inferdiffprefs``

* Relaxed convergence criterion for cases when *R* is ``nan`` for a few sites in MCMC

* Some minor documentation updates

1.0.0
--------
Initial release
