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
   :prog: dms_logoplot

   infile
    Should be in the format of a :ref:`preferences_file`, a :ref:`diffpreferences_file`, or one of the ``*mutdiffsel.txt`` files created by :ref:`dms_diffselection`.

   logoplot
    See `Examples`_ for images of the types of plots that are created.

   \-\-diffselheight
    So if you are using this option, you should specify a list of files in the format of the ``*mutdiffsel.txt`` files created by :ref:`dms_diffselection`.

   \-\-stringencyparameter
    Use this option if you have fit a stringency parameter :math:`\beta` and want to rescale the visualization of the preferences using this parameter. The preferences are rescaled so that each preference is proportional to :math:`\left(\pi_{r,a}\right)^{\beta}`, so values > 1 increase the weight of the preferences and values < 1 flatten them.

   \-\-overlay1
    For instance, contents of an example file specifying relative solvent accessibility for residues 2, 3, 4, and 6 is shown below. In the created ``logoplot``, residue 5 (assuming it exists in ``infile``) will not have any relative solvent accessibility shown::

        #SITE RSA
        2     0.3
        3     0.56
        4     0.02
        6     0.72

    A file specifying secondary structure for the same subset of residues is here::

        #SITE SS
        2     helix
        3     helix
        4     coil
        5     sheet

    For instance, if the above file was named ``secondary_structures.txt``, then you would use the option as ``--overlay1 secondary_structures.txt SS "secondary structure"``.

Examples
-----------

Preferences logo plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Imagine that you have used inferred site-specific preferences into ``preferences.txt``, which has the format of a :ref:`preferences_file`. To display these preferences using a logo plot, run::

    dms_logoplot preferences.txt prefs_logoplot.pdf --nperline 81 --excludestop

This will create the file ``prefs_logoplot.pdf``, which will look something like the image below. Note that this logo plot does **not** include stop codons; if they were shown then they would be black ``*`` characters. In the logo plot, the height of each letter is proportional to the preference for that amino acid. The letter heights sum to one, since :math:`\sum_a \pi_{r,a} = 1`.

.. image:: prefs_logoplot.pdf
   :width: 90%
   :alt: prefs_logoplot.pdf
   :align: center

Differential selection logo plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Imagine that you have inferred differential selection into ``mutdiffsel.txt`` using :ref:`dms_diffselection`.

You can visualize the differential selection (the :math:`s_{r,x}` values described in :ref:`dms_diffselection`) using::

    dms_logoplot mutdiffsel.txt diffsel_logoplot.pdf --nperline 115 --diffselheight mutdiffsel.txt mutdiffsel2.txt mutdiffsel3.txt

The ``--diffselheight`` is useful if you are making several such plots and want them to share a common y-axis.

Any mutations that have missing differential selection values (specified as ``NaN`` in the ``mutdiffsel.txt`` file) are treated as having a differential selection of zero in the plotting.

Here is an example of a created plot:

.. image:: diffsel_logoplot.pdf
   :width: 90%
   :alt: diffsel_logoplot.pdf
   :align: center

Differential preferences logo plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Imagine that you have inferred differential preferences into ``diffprefs.txt``, which has the format of a :ref:`diffpreferences_file`. To display these differential preferences, run::

    dms_logoplot diffprefs.txt diffprefs_logoplot.pdf --nperline 81 --excludestop --diffprefheight 0.5

This will create the file ``diffprefs_logoplot.pdf``, which will look something like the image below. The height of letters above or below the center line are proportional to the differential preference for or against that amino acids. The overall negative and positive heights are equal since :math:`0 = \sum_a\Delta\pi_{r,a}`. For instance, in the plot below there is a strong differential preference for *M* at site 182, and a strong differential preference for *V* at position 90. Note that in the command above we used ``--diffprefheight 0.5``. This is because if we were not to rescale the y-axis (maximum height of the letter stacks) down, we would only use some of the dynamic range since the maximal differential preference is :math:`< 0`.

.. image:: diffprefs_logoplot.pdf
   :width: 90%
   :alt: diffprefs_logoplot.pdf
   :align: center

Plot with an overlay
~~~~~~~~~~~~~~~~~~~~~~~~~
Now let's add an overlay to the `Preferences logo plot`_. We create files giving the secondary structure and the relative solvent accessibility. These files are ``SSs.txt``, which has the first few lines as follows::

    #SITE SS
    18 loop
    19 strand
    20 strand
    21 strand
    22 strand
    23 strand
    24 strand
    25 strand
    26 loop

and ``RSAs.txt``, which has the first few lines as follows::

    #SITE RSA
    18 0.170984455959
    19 0.168604651163
    20 0.0
    21 0.0778443113772
    22 0.00507614213198
    23 0.0
    24 0.0
    25 0.0803571428571
    26 0.0077519379845

We now use the command::

    dms_logoplot preferences.txt prefs_logoplot.pdf --nperline 81 --excludestop --overlay1 RSAs.txt RSA "relative solvent accessibility" --overlay2 SSs.txt SS "secondary structure"

to create the following image, which displays overlay bars with the secondary structure and solvent accessibility:

.. image:: prefs_logoplot_withoverlay.pdf
   :width: 90%
   :alt: prefs_logoplot_withoverlay.pdf
   :align: center


.. include:: weblinks.txt
