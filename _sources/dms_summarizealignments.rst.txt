.. _dms_summarizealignments:

==========================================
``dms_summarizealignments``
==========================================

.. contents::

Overview
-------------
``dms_summarizealignments`` is a program included with the `dms_tools`_ package. It can be used to make plots that summarize the results after examining several samples with any of the following programs:

    * :ref:`dms_barcodedsubamplicons` 
    
    * :ref:`dms_subassemble`

    * :ref:`dms_matchsubassembledbarcodes`

After you install `dms_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: SummarizeAlignmentsParser
   :prog: dms_summarizealignments

   outprefix
    See `Output files`_ for a detailed description of the created files and examples.

   alignments
    ``ALIGNPREFIX`` should be the value of ``--outprefix`` used when running :ref:`dms_barcodedsubamplicons`, :ref:`dms_subassemble`, or :ref:`dms_matchsubassembledbarcodes` (depending on the value of ``alignment_type``). This program expects to find the files that are created by running that program with ``--outprefix`` as ``ALIGNPREFIX``. ``NAME`` is simply the name given to each sample in the `Output files`_.

   \-\-groupbyfirst
    This option currently only applies when ``alignment_type`` is ``matchsubassembledbarcodes`. In this case, for the ``allmutdepth.pdf``, ``singlemutdepth.pdf``, ``allcumulfracs.pdf``, and ``singlecumulfracs.pdf`` `Output files`_, a separate output file is made by grouping samples with the same first word (prefix) in the name, where the first word is separated from the others by an underscore (``_``). This first name is then added directly before the plot extension in the file name.



Example usage
--------------
Imagine that you have six samples, *DNA*, *mutDNA*, *virus-p1*, *mutvirus-p1*, *virus-p2*, and *mutvirus-p2*, and furthermore imagine that for each sample you ran :ref:`dms_barcodedsubamplicons` with its ``--outprefix`` option to values of ``DNA/DNA_``, ``mutDNA/mutDNA_``, etc. You could then summarize the results with the following command::

    dms_summarizealignments alignmentsummary_ barcodedsubamplicons DNA/DNA_,DNA mutDNA/mutDNA_,mutDNA virus-p1/virus-p1_,virus-p1 mutvirus-p1/mutvirus-p1_,mutvirus-p1 virus-p2/virus-p2_,virus-p2 mutvirus-p2/mutvirus-p2_,mutvirus-p2

Because this command uses ``outprefix`` of ``alignmentsummary_``, the `Output files`_ will have names like ``alignmentsummary_mutfreqs.pdf``, ``alignmentsummary_mutcounts_all.pdf``, etc. The ``alignments`` option specifies the prefix for each of the six samples, and then gives them the same names as above (*DNA*, *mutDNA*, etc).

Alternatively, if you have subassembled samples *Lib1*, *Lib2*, *Lib3*, and *WT* with :ref:`dms_subassemble` with ``--outprefix`` as ``subassembly/Lib1``, ``subassembly/Lib2``, etc then you could summarize the results with::

    dms_summarizealignments subassemblysummary_ subassemble ./subassembly/Lib1,Lib1 ./subassembly/Lib2,Lib2 ./subassembly/Lib3,Lib3 ./subassembly/WT,WT

Because the command uses ``outprefix`` of ``subassemblysummary_``, the `Output files`_ will have names like ``subassemblysummary_reads.pdf``, ``subassemblysummary_depth.pdf``, and ``subassemblysummary_nsubassembled.pdf``.

Alternatively, if you have matched barcodes for samples after subassembly, you could summarize the results with::

    dms_summarizealignments matchbarcodesummary_ matchsubassembledbarcodes ./barcodes//Lib1_Input,Lib1_Input ./barcodes//Lib1_0-cip,Lib1_0-cip ./barcodes//Lib1_4-cip,Lib1_4-cip ./barcodes//Lib1_10-cip,Lib1_10-cip ./barcodes//Lib1_15-cip,Lib1_15-cip ./barcodes//Lib2_Input,Lib2_Input ./barcodes//Lib2_0-cip,Lib2_0-cip ./barcodes//Lib2_4-cip,Lib2_4-cip ./barcodes//Lib2_10-cip,Lib2_10-cip ./barcodes//Lib2_15-cip,Lib2_15-cip ./barcodes//Lib3_Input,Lib3_Input ./barcodes//Lib3_0-cip,Lib3_0-cip ./barcodes//Lib3_4-cip,Lib3_4-cip ./barcodes//Lib3_10-cip,Lib3_10-cip ./barcodes//Lib3_15-cip,Lib3_15-cip ./barcodes//WT_Input,WT_Input ./barcodes//WT_0-cip,WT_0-cip ./barcodes//WT_4-cip,WT_4-cip ./barcodes//WT_10-cip,WT_10-cip ./barcodes//WT_15-cip,WT_15-cip --groupbyfirst

Because the command uses ``outprefix`` of ``matchbarcodesummary``, the `Output files`_ will have names like ``matchbarcodesummary_reads.pdf``, etc.

Output files
---------------

For ``alignment_type`` of ``barcodedsubamplicons``
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

Here are the output files, with the names that would be produced if running ``dms_summarizealignments`` as in the `Example usage`_ for ``alignment_type`` of ``barcodedsubamplicons``:

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

If ``--writemutfreqs`` is used, a text file with the suffix ``mutfreqs.txt`` is created that has the raw data in this plot.

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

``singlemuttypes.pdf``
~~~~~~~~~~~~~~~~~~~~~~~~~~
This plot shows the frequency of different types of nucleotide mutations among **only** codons with single-nucleotide changes. This can be useful for assessing if oxidative damage is affecting your sample by inducing specific mutations. Note that this file shows a different set of samples than for the others above.

.. image:: alignmentsummary_singlemuttypes.pdf
   :align: center
   :width: 60%
   :alt: alignmentsummary_singlemutttypes.pdf



For ``alignment_type`` of ``subassemble``
++++++++++++++++++++++++++++++++++++++++++
Here are the output files, with the names that would be produced if running ``dms_summarizealignments`` as in the `Example usage`_ for ``alignment_type`` of ``subassemble``:

``reads.pdf``
~~~~~~~~~~~~~
The number of reads and whether or not they were part of a successful subassembly.

.. image:: subassemblysummary_reads.pdf
   :align: center
   :width: 60%
   :alt: subassemblysummary_reads.pdf


``depth.pdf``
~~~~~~~~~~~~~
The number of reads at each position, regardless of whether or not they were part of a successful subassembly.

.. image:: subassemblysummary_depth.pdf
   :align: center
   :width: 60%
   :alt: subassemblysummary_depth.pdf

``nsubassembled.pdf``
~~~~~~~~~~~~~~~~~~~~~~
The number of subassembled barcodes and how many mutations they have.

.. image:: subassemblysummary_nsubassembled.pdf
   :align: center
   :width: 60%
   :alt: subassemblysummary_nsubassembled.pdf

For ``alignment_type`` of ``matchsubassembledbarcodes``
++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Here are the output files, with the names that would be produced if running ``dms_summarizealignments`` as in the `Example usage`_ for ``alignment_type`` of ``matchsubassembledbarcodes``:

``reads.pdf``
~~~~~~~~~~~~~~
The number of read pairs for each sample, and their fate (only those indicated as *nretained* are actually successfully matched to a subassembled barcode).

.. image:: matchbarcodesummary_reads.pdf
   :align: center
   :width: 60%
   :alt: matchbarcodesummary_reads.pdf

``muts_per_matched_read.pdf``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For each read pair that matches a barcode (there are *nretained* such reads as indicated in the previous plot), this shows that number that match to a variant with the indicated number of mutations.

.. image:: matchbarcodesummary_muts_per_matched_read.pdf
   :align: center
   :width: 60%
   :alt: matchbarcodesummary_muts_per_matched_read.pdf

``allmutfreqs.pdf``
~~~~~~~~~~~~~~~~~~~~
This plot shows the average per-site mutation rate averaged over all sites and variants (unmutated, singly mutated, multiply mutated).

.. image:: matchbarcodesummary_allmutfreqs.pdf
   :align: center
   :width: 60%
   :alt: matchbarcodesummary_allmutfreqs.pdf

If ``--writemutfreqs`` is used, a text file with the suffix ``allmutfreqs.txt`` is created that has the raw data in this plot.

``singlemutfreqs.pdf``
~~~~~~~~~~~~~~~~~~~~~~~~
This plot shows the fraction of variants that have a mutation of the indicated type **only** among the unmutated and singly mutated variants.

.. image:: matchbarcodesummary_singlemutfreqs.pdf
   :align: center
   :width: 60%
   :alt: matchbarcodesummary_singlemutfreqs.pdf

If ``--writemutfreqs`` is used, a text file with the suffix ``singlemutfreqs.txt`` is created that has the raw data in this plot.

``allmutdepth.pdf``
~~~~~~~~~~~~~~~~~~~~
This plot shows the average per-site mutation rate as a function of primary sequence for all variants (unmutated, singly mutated, multiply mutated). Note that since we have used ``--groupbyfirst``, the `Example usage`_ above would create four plots. Here we show the first, ``allmutdepth_Lib1.pdf``.

.. image:: matchbarcodesummary_allmutdepth_Lib1.pdf
   :align: center
   :width: 60%
   :alt: matchbarcodesummary_allmutdepth_Lib1.pdf

``singlemutdepth.pdf``
~~~~~~~~~~~~~~~~~~~~~~
This plot shows the average per-site mutation rate at which each site has a singly mutated variant versus the number of unmutated variants. Note that since we have used ``--groupbyfirst``, the `Example usage`_ above would create four plots. Here we show the first, ``singlemutdepth_Lib1.pdf``.

.. image:: matchbarcodesummary_singlemutdepth_Lib1.pdf
   :align: center
   :width: 60%
   :alt: matchbarcodesummary_singlemutdepth_Lib1.pdf

``singlemutdepth`` and ``allmutdepth`` plots for **each sample**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A plot is also made for **each** sample that shows the mutation depth divided into nonsynonymous, synonymous, and stop-codon mutations. These plots have names like ``matchbarcodesummary_allmutdepth_Lib1_Input.pdf``. Note that the synonymous and stop-codon mutations frequencies are multiplied by a scale factor as indicated in the legend so that they are on a comparable scale to the nonsynonymous mutation frequencies.

Plots like this are made for *all* variants and just for *single* mutant variants.

.. image:: matchbarcodesummary_allmutdepth_Lib1_Input.pdf
   :align: center
   :width: 60%
   :alt: matchbarcodesummary_allmutdepth_Lib1_Input.pdf


``allcumulcounts.pdf``
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The number of times each mutation is observed among all variants (singly and multiply mutated). Since we used the ``--groupbyfirst``, the `Example usage`_ above would create four plots. here is the first, ``allcumulcounts_Lib1.pdf``.

.. image:: matchbarcodesummary_allcumulcounts_Lib1.pdf
   :align: center
   :width: 60%
   :alt: matchbarcodesummary_allcumulcounts_Lib1.pdf

``singlecumulcounts.pdf``
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The number of times each mutation is aboserved among singly mutated variants. Since we used the ``--groupbyfirst``, the `Example usage`_ above would create four plots. here is the first, ``singlecumulcounts_Lib1.pdf``.

.. image:: matchbarcodesummary_singlecumulcounts_Lib1.pdf
   :align: center
   :width: 60%
   :alt: matchbarcodesummary_singlecumulcounts_Lib1.pdf


.. include:: weblinks.txt
.. include:: weblinks.txt
