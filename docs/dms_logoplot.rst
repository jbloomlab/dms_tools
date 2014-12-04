.. _dms_logoplot:

==========================================
``dms_logoplot``
==========================================

.. contents::

Overview
-------------
``dms_logoplot`` is a program included with the `dms_tools`_ package. It uses `weblogo`_ (via the ``weblogolib`` `Python`_ API) to make logo plots that visually display the preferences or differential preferences.

After you install `dms_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: LogoPlotParser
   :prog: dms_weblogo

   infile
    Should be in the formats of a :ref:`preferences_file` or a :ref:`diffpreferences_file`.

Examples
-----------

A preferences logo plot
---------------------------
Imagine that you have used inferred site-specific preferences into ``preferences.txt``, which has the format of a :ref:`preferences_file`. To display these preferences using a logo plot, run::

    dms_logoplot preferences.txt prefs_logoplot.pdf --nperline 81 --excludestop

This will create the file ``prefs_logoplot.pdf``, which will look something like the image below. Note that this logo plot does **not** include stop codons; if they were shown then they would be black ``*`` characters.

.. image:: prefs_logoplot.pdf
   :width: 90%
   :alt: prefs_logoplot.pdf
   :align: center



.. include:: weblinks.txt
