.. _dms_correlate:

==========================================
``dms_correlate``
==========================================

.. contents::

Overview
-------------
``dms_correlate`` is a program included with the `dms_tools`_ package. It computes the correlations between pairs of preferences, differential preferences, or differential selections.

After you install `dms_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: CorrelateParser
   :prog: dms_correlate

   file1
    Should be in the formats of a :ref:`preferences_file`, a :ref:`diffpreferences_file`, or a site or mutation differential selection file returned by :ref:`dms_diffselection`.

   file2
    Must be in the **same format** as "file1".

   \-\-enrichment
    The enrichment ratio :math:`\phi_{r,x}` for character :math:`x` at site :math:`r` is defined by

    .. math::

        \phi_{r,x} = \frac{\pi_{r,x}}{\pi_{r,\operatorname{wt}\left(r\right)}}

    where :math:`\pi_{r,x}` is the preference of :math:`r` for :math:`x`, and :math:`\operatorname{wt}\left(r\right)` is the wildtype identity at :math:`r`.

    Note that we do not include enrichment ratios for wildtype characters (which are one by definition) in the correlations or plots, and that the enrichment ratios are log transformed before plotting and computing correlations.

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
