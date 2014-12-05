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

   logoplot
    See `Examples`_ for images of the types of plots that are created.

Examples
-----------

A preferences logo plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Imagine that you have used inferred site-specific preferences into ``preferences.txt``, which has the format of a :ref:`preferences_file`. To display these preferences using a logo plot, run::

    dms_logoplot preferences.txt prefs_logoplot.pdf --nperline 81 --excludestop

This will create the file ``prefs_logoplot.pdf``, which will look something like the image below. Note that this logo plot does **not** include stop codons; if they were shown then they would be black ``*`` characters. In the logo plot, the height of each letter is proportional to the preference for that amino acid. The letter heights sum to one, since :math:`\sum_a \pi_{r,a} = 1`.

.. image:: prefs_logoplot.pdf
   :width: 90%
   :alt: prefs_logoplot.pdf
   :align: center

A differential preferences logo plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Imagine that you have inferred differential preferences into ``diffprefs.txt``, which has the format of a :ref:`diffpreferences_file`. To display these differential preferences, run::

    dms_logoplot diffprefs.txt diffprefs_logoplot.pdf --nperline 81 --excludestop --diffprefheight 0.5

This will create the file ``diffprefs_logoplot.pdf``, which will look something like the image below. The height of letters above or below the center line are proportional to the differential preference for or against that amino acids. The overall negative and positive heights are equal since :math:`0 = \sum_a\Delta\pi_{r,a}`. For instance, in the plot below there is a strong differential preference for *M* at site 182, and a strong differential preference for *V* at position 90. Note that in the command above we used ``--diffprefheight 0.5``. This is because if we were not to rescale the y-axis (maximum height of the letter stacks) down, we would only use some of the dynamic range since the maximal differential preference is :math:`< 0`.

.. image:: diffprefs_logoplot.pdf
   :width: 90%
   :alt: diffprefs_logoplot.pdf
   :align: center

.. include:: weblinks.txt
