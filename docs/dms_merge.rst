.. _dms_merge:

==========================================
``dms_merge``
==========================================

.. contents::

Overview
-------------
``dms_merge`` is a program included with the `dms_tools`_ package. It merges preferences or differential preferences in different files by averaging them or adding / subtracting them. It also sums the counts from count files together. See `Examples`_ for illustrations of how you might do this.

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

   outfile
    For ``merge_method`` of "average", then this file will be of the same type as ``infiles`` (either preferences or differential preferences).

    For ``merge_method`` of "sum", this file may be either preferences, differential preferences, or counts:

        * If any of the files in ``infiles`` and ``--minus`` give preferences and the sum at each site is one, then ``outfile`` gives preferences.

        * If any of the files in ``infiles`` and ``--minus`` give preferences and the sum at each site is zero, then ``outfile`` gives differential preferences.

        * If all of the files in ``infiles`` and ``--minus`` are differential preferences, the ``outfile`` gives differential preferences.

        * If all of the files in ``infiles`` are counts, the ``outfile`` gives counts.

    If ``outfile`` gives preferences, it will be in the :ref:`preferences_file` format. Note that this file will **not** have any columns giving 95% credible intervals, as these cannot be calculated when merging files.

    If ``outfile`` gives differential preferences, it will be in the :ref:`diffpreferences_file` format. Note that this file will **not** have any columns giving posterior probabilities, as these cannot be calculated when merging files.

    If ``outfile`` gives counts, it will be in the :ref:`dms_counts` format. 

Examples
-----------
These examples assume that you have used several experiments to measure counts and site-specific preferences in the files ``counts1.txt``, ``counts2.txt``, ``prefs1.txt``, ``prefs2.txt`` etc. You have also used several experiments to measure differential preferences in ``diffprefs1.txt``, ``diffprefs2.txt``, etc.

Averaging preferences
~~~~~~~~~~~~~~~~~~~~~~~
Use::

    dms_merge avgprefs.txt average prefs1.txt prefs2.txt

The created file ``avgprefs.txt`` will have the :ref:`diffpreferences_file` format.

Averaging differential preferences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use::

    dms_merge avgdiffprefs.txt average diffprefs1.txt diffprefs2.txt

The created file ``avgdiffprefs.txt`` will have the :ref:`preferences_file` format.

Invalid: averaging a mix of preferences and differential preferences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You cannot average a mix of preferences and differential preferences. So the following would be **invalid** and would cause an error::

    dms_merge mixedaverage.txt average prefs1.txt prefs2.txt diffprefs1.txt

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

Adding and subtracting differential preferences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The command::

    dms_merge summeddiffprefs.txt sum diffprefs1.txt diffprefs2.txt --minus diffprefs3.txt

creates ``summeddiffprefs.txt`` in the :ref:`diffpreferences_file` format.

Adding and subtracting a mix of preferences and differential preferences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The command::

    dms_merge sumprefs.txt prefs1.txt prefs2.txt diffprefs1.txt --minus prefs3.txt diffprefs2.txt

creates ``sumprefs.txt`` in :ref:`preferences_file` format, since the sum at each site is one.

In contrast, the command::

    dms_merge sumdiffprefs.txt prefs1.txt diffprefs1.txt --minus prefs2.txt

creates ``sumdiffprefs.txt`` in :ref:`diffpreferences_file` format, since the sum at each site is zero.

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
