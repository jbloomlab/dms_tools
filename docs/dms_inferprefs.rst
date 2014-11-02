.. _dms_inferprefs:

==========================================
``dms_inferprefs``
==========================================

Overview
-------------
``dms_inferprefs`` is a program included with the `dms_tools`_ package. It infers the site-specific preferences :math:`\pi_{r,a}` for each character :math:`a` (an amino acid, codon, or nucleotide) at each site :math:`r`.

After you install `dms_tools`_, this program will be available to run at the command line.

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: InferPrefsParser
   :prog: dms_inferprefs

   n_pre
    See :ref:`dms_counts` for specification of the file format.

   n_post
    See :ref:`dms_counts` for specification of the file format.

   errmodel
    See :ref:`dms_counts` for specification of the file format for the error files "n_err* or "n_errpre" and "n_errpost".

Examples
-----------

.. include:: weblinks.txt
