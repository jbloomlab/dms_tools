Changelog
===========

1.1.6
-------
Added ``--stringencyparameter`` option to ``dms_logoplot``.

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
* Added ``dms_barcodedsubamplicons`` and ``dms_summarizealignments` programs

1.0.1
--------
* Fixed bug in parsing *codon* option in ``dms_infeprefs`` and ``dms_inferdiffprefs``

* Relaxed convergence criterion for cases when *R* is ``nan`` for a few sites in MCMC

* Some minor documentation updates

1.0.0
--------
Initial release
