.. _dms_merge:

==========================================
``dms_merge``
==========================================

.. contents::

Overview
-------------
``dms_merge`` is a program included with the `dms_tools`_ package. It merges preferences, differential selection values, or differential preferences in different files by averaging them or adding / subtracting them. It also sums the counts from count files together. See `Examples`_ for illustrations of how you might do this.

After you install `dms_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: MergeParser
   :prog: dms_merge

   infile
    Files should be in the formats of a :ref:`preferences_file`, :ref:`diffpreferences_file`, or :ref:`dms_counts`.

   \-\-minus
    Files should be in the formats of a :ref:`preferences_file` or a :ref:`diffpreferences_file`. Currently, this option cannot be used with counts files.

   infiles 
    Note that for differential selection, the files can only be the ``*mutdiffsel.txt`` files created by :ref:`dms_diffselection`, **not** the ``*sitediffsel.txt`` files. However, you can create a ``*sitediffsel.txt`` file from the merged ``*mutdiffsel.txt`` files using the ``--sitediffselfile`` option.

   outfile
    For ``merge_method`` of "average", then this file will be of the same type as ``infiles`` (preferences, mutation-level differential selection values, or differential preferences). For merging differential selection, any values of *NaN* are ignored.

    For ``merge_method`` of "median", then infiles and outfile must be of the format of the ``*mutdiffsel.txt`` files created by :ref:`dms_diffselection`.

    For ``merge_method`` of "sum", this file may be either preferences, differential preferences, or counts:

        * If any of the files in ``infiles`` and ``--minus`` give preferences and the sum at each site is one, then ``outfile`` gives preferences.

        * If any of the files in ``infiles`` and ``--minus`` give preferences and the sum at each site is zero, then ``outfile`` gives differential preferences.

        * If all of the files in ``infiles`` and ``--minus`` are differential preferences, the ``outfile`` gives differential preferences.

        * If all of the files in ``infiles`` are counts, the ``outfile`` gives counts.

    For ``merge_method`` of "rescale", this created file will give preferences (after re-scaling).

    If ``outfile`` gives preferences, it will be in the :ref:`preferences_file` format. Note that this file will **not** have any columns giving 95% credible intervals, as these cannot be calculated when merging files.

    If ``outfile`` gives differential preferences, it will be in the :ref:`diffpreferences_file` format. Note that this file will **not** have any columns giving posterior probabilities, as these cannot be calculated when merging files.

    If ``outfile`` gives mutation-level differential selection values, it will be in the format of the ``*mutdiffsel.txt`` file created by :ref:`dms_diffselection`. To create a corresponding ``*sitediffsel.txt`` file, see the ``--sitediffselfile`` option.

    If ``outfile`` gives counts, it will be in the :ref:`dms_counts` format. 

   \-\-sitediffselfile 
    After merging the mutation-level differential selection values in ``*mutdiffsel.txt`` ``infiles``, you can create a new ``*sitediffsel.txt`` file (in the format of those created by :ref:`dms_diffselection` with this option.
    

Examples
-----------

Rescaling preferences by a stringency parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use::

    dms_merge rescaledprefs.txt rescale prefs.txt --stringencyparameter 2.5

to rescale the preferences in ``prefs.txt`` by a stringency parameter of :math:`beta = 2.5`. The created file ``rescaledprefs.txt`` will be in the :ref:`preferences_file` format.

Averaging differential selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use::

    dms_merge avgmutdiffsel.txt average mutdiffsel1.txt mutdiffsel2.txt

to average the two mutation differential selection files to create ``avgmutdiffsel.txt`` in the format of a ``*mutdiffsel.txt`` file from :ref:`dms_diffselection`. If you also want to create a ``*sitediffsel.txt`` file corresponding to ``avgmutdiffsel.txt``, use::

    dms_merge avgmutdiffsel.txt average mutdiffsel1.txt mutdiffsel2.txt --sitediffselfile avgsitediffsel.txt

to create ``avgsitediffsel.txt`` with the site differential selection values in the format created by :ref:`dms_diffselection`.

Averaging preferences
~~~~~~~~~~~~~~~~~~~~~~~
Use::

    dms_merge avgprefs.txt average prefs1.txt prefs2.txt

The created file ``avgprefs.txt`` will have the :ref:`preferences_file` format.

Adding and subtracting preferences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The command::

    dms_merge summedprefs.txt sum prefs1.txt prefs2.txt --minus prefs3.txt

would create an output file ``summedprefs.txt`` in the :ref:`preferences_file` format.

The command::

    dms_merge summedprefs.txt sum prefs1.txt prefs2.txt

is **invalid** since it does not create preferences that sum to one at each site (instead they sum to two).

The command::

    dms_merge diffs.txt sum prefs1.txt prefs2.txt --minus prefs3.txt prefs4.txt

would create a file ``diffs.txt`` in the :ref:`diffpreferences_file` format since the total for each site sums to zero.

Adding counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The command::

    dms_merge summedcounts.txt sum counts1.txt counts2.txt --chartype codon 

sums counts from the two provided counts files and writes them to the file ``summedcounts.txt`` in the :ref:`dms_counts` format.

Adding counts with normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The command::

    dms_merge summedcounts.txt sum counts1.txt counts2.txt --chartype codon --normalize

normalizes counts at each site to the minimum number of total counts observed at that site in ``counts1.txt`` and ``counts2.txt``. The normalized counts are then summed and written to the file ``summedcounts.txt`` in the :ref:`dms_counts` format.

.. include:: weblinks.txt
