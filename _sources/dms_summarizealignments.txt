.. _dms_summarizealignments:

==========================================
``dms_summarizealignments``
==========================================

.. contents::

Overview
-------------
``dms_summarizealignments`` is a program included with the `dms_tools`_ package. It can be used to make plots that summarize the results after examining several samples with :ref:`dms_barcodedsubamplicons`.

After you install `dms_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: SummarizeAlignmentsParser
   :prog: dms_summarizealignments

   outprefix
    See `Output files`_ for a detailed description of the created files and examples.

   alignment_type
    Currently, ``dms_summarizealignments`` assumes that the alignments were created using :ref:`dms_barcodedsubamplicons`.

   alignments
    Essentially, ``ALIGNPREFIX`` should be the value of ``--outprefix`` used when running :ref:`dms_barcodedsubamplicons`. This program expects to find the ``counts.txt`` and ``summarystats.txt`` files that are created by running :ref:`dms_barcodedsubamplicons` with ``--outprefix`` as ``ALIGNPREFIX``. ``NAME`` is simply the name given to each sample in the `Output files`_.

Example usage
--------------
Imagine that you have six samples, *DNA*, *mutDNA*, *virus-p1*, *mutvirus-p1*, *virus-p2*, and *mutvirus-p2*, and furthermore imagine that for each sample you ran :ref:`dms_barcodedsubamplicons` with its ``--outprefix`` option to values of ``DNA/DNA_``, ``mutDNA/mutDNA_``, etc. You could then summarize the results with the following command::

    dms_summarizealignments alignmentsummary_ barcodedsubamplicons DNA/DNA_,DNA mutDNA/mutDNA_,mutDNA virus-p1/virus-p1_,virus-p1 mutvirus-p1/mutvirus-p1_,mutvirus-p1 virus-p2/virus-p2_,virus-p2 mutvirus-p2/mutvirus-p2_,mutvirus-p2

Because this command uses ``outprefix`` of ``alignmentsummary_``, the `Output files`_ will have names like ``alignmentsummary_mutfreqs.pdf``, ``alignmentsummary_mutcounts_all.pdf``, etc. The ``alignments`` option specifies the prefix for each of the six samples, and then gives them the same names as above (*DNA*, *mutDNA*, etc).

Output files
---------------
Here are the output files, with the names that would be produced if running ``dms_summarizealignments`` as in the `Example usage`_.

``reads.pdf``
~~~~~~~~~~~~~~
This plot shows the number of read pairs that failed the Illumina chastity filter, were filtered by :ref:`dms_barcodedsubamplicons` as being of low quality, or were successfully processed and assigned to a barcode.

.. image:: alignmentsummary_reads.pdf
   :align: center
   :width: 60%
   :alt: alignmentsummary_reads.pdf

``barcodes.pdf``
~~~~~~~~~~~~~~~~~~~~~~
For all barcodes for which read were assigned, this plot shows the number barcodes with one, two, or at least three read pairs. For each category, it also indicates whether it was possible to build a consensus and align it to the reference sequence, whether the consensus was not alignable to the reference sequence, or whether the barcode simply had to be discarded due to insufficient matching reads.

.. image:: alignmentsummary_barcodes.pdf
   :align: center
   :width: 60%
   :alt: alignmentsummary_barcodes.pdf

``mutfreqs.pdf``
~~~~~~~~~~~~~~~~~
This plot shows the average per-codon mutation rate across the gene as determined from the aligned barcodes. Each bar is split: the left half categories mutations as to whether they are synonymous, nonsynonymous, or to stop codons; the right categorizes mutations by the number of nucleotide changes to that codon.

.. image:: alignmentsummary_mutfreqs.pdf
   :align: center
   :width: 60%
   :alt: alignmentsummary_mutfreqs.pdf

``depth.pdf``
~~~~~~~~~~~~~~
This plot shows the per-codon sequencing depth as a function of the codon position in the primary sequence.

.. image:: alignmentsummary_depth.pdf
   :align: center
   :width: 60%
   :alt: alignmentsummary_depth.pdf

``mutdepth.pdf``
~~~~~~~~~~~~~~~~~~~~~
This plot shows the per-codon mutation rate as a function of the codon position in the primary sequence.

.. image:: alignmentsummary_mutdepth.pdf
   :align: center
   :width: 60%
   :alt: alignmentsummary_mutdepth.pdf

``mutcounts_all.pdf``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This plot shows the number of times each of the possible codon mutations and each of the possible synonymous codon mutations was observed. For each number of counts, it indicates the fraction of all of these possible mutations that were observed at least that many times.

.. image:: alignmentsummary_mutcounts_all.pdf
   :align: center
   :width: 60%
   :alt: alignmentsummary_mutcounts_all.pdf

``mutcounts_multi_nt.pdf``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This plot shows the number of times each codon mutation was observed **only** for mutations that involve multi-nucleotide changes to the codon. This plot might be preferable to the one above in some cases because single-nucleotide codon changes can arise due to errors, but multi-nucleotide changes typically do not.

.. image:: alignmentsummary_mutcounts_multi_nt.pdf
   :align: center
   :width: 60%
   :alt: alignmentsummary_mutcounts_multi_nt.pdf


.. include:: weblinks.txt
