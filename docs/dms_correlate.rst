.. _dms_correlate:

==========================================
``dms_correlate``
==========================================

.. contents::

Overview
-------------
``dms_correlate`` is a program included with the `dms_tools`_ package. It computes the correlations between pairs of preferences or differential preferences.

After you install `dms_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: CorrelateParser
   :prog: dms_correlate

   file1
    Should be in the formats of a :ref:`preferences_file` or a :ref:`diffpreferences_file`.

   file2
    Must be in the **same format** as "file1".

Examples
-----------
Imagine that you have used several experiments to measure site-specific preferences in the files ``prefs1.txt`` and ``prefs2.txt``. Then run the command::

    dms_correlate prefs1.txt prefs2.txt pref1_vs_pref2_corr --name1 "sample 1" --name2 "sample 2"

This will create the files ``pref1_vs_pref2_corr.txt`` and ``pref1_vs_pref2_corr.pdf``.

Here are the contents of ``pref1_vs_pref2_corr.txt``::

    R = 0.540113
    P = 0
    N = 11844

Here is an image of ``pref1_vs_pref2_corr.pdf``:

.. image:: pref1_vs_pref2_corr.pdf
   :width: 30%
   :alt: pref1_vs_pref2_corr.pdf
   :align: center



.. include:: weblinks.txt
