.. _dms_matchsubassembledbarcodes:

==========================================
``dms_matchsubassembledbarcodes``
==========================================

.. contents::

Overview
-------------
You would use this program if you have subassembled some sequences using :ref:`dms_subassemble`, and then have performed some selections on the subassembled sequences and sequenced the barcodes.

This program assumes that you have sequenced the barcodes with overlapping paired-end reads (the overlap helps with error correction). You give it the R1 and R2 FASTQ files as ``r1files`` and ``r2files`` , and it extracts the barcodes and uses them to determine the counts for individual variants / mutations.

This program assumes that the R1 files (``r1files``) are reading the barcodes in the same orientation as they are specified in the ``subassembled`` file created by :ref:`dms_subassemble`, and so that the R2 files (``r2files`` are reading the reverse complement of these barcodes). The options ``--r1start`` and ``--r2end`` determine where exactly in the reads the barcode is found.

The barcodes are only retained if all of the following conditions are met:

    * Both reads in the pair pass the Illumina filter.

    * At least one barcode has a quality score :math:`\ge` ``--minquality`` for each site.

    * Both read pairs agree on the identity of all positions in the barcode without any ``N`` nucleotides.



Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: MatchSubassembledBarcodesParser
   :prog: dms_matchsubassembledbarcodes

   outprefix
    See `Output files`_ for a list of the created files.

   subassembled
    This program matches the barcodes against barcodes that have been subassembled by :ref:`dms_subassemble`. You therefore need to provide it with a ``*_subassembled_variants.txt`` created by :ref:`dms_subassemble`. This should be in the :ref:`subassembly_file` format.

    The barcode length is determined from this file, since it will list the barcodes and so specify their length.

   \-\-r1start
    This argument specifies where in R1 the barcode starts. By default, we assume that it starts with the first nucleotide of the read.

    Must be an integer :math:`\ge 1`.

   \-\-r2end
    This argument specifies where in R2 the barcode ends. By default, we assume that it ends with the last nucleotide of the read, and so set a value of zero. If instead the barcode ends two nucleotides before the end of R2, you would set a value of 2.

    Must be an integer :math:`\ge 0`.

   \-\-purgefrac
    This option is useful if you are trying to figure out whether to sequence to greater depth. You can use it to randomly subsample your data, and so figure out if your results change if you have fewer reads (which is an indication if more reads will be helpful).

    The purging here is for the barcode read pairs in ``r1files`` and ``r2files``, **not** for the subassembled variants in ``subassembled``.

 

Output files
--------------
The following output files are created. Each has the prefix specified by ``outprefix`` with the following suffixes.

Single-mutant counts file
+++++++++++++++++++++++++++
This file has the prefix specified by ``outprefix`` followed by the suffix ``_singlemutcounts.txt``. This file is in the format of a :ref:`dms_counts` files. The mutation counts at each site are **only** for single-mutant variants; any subassembled variants with multiple mutations are ignored in this file. The wildtype counts for each site given by this file are **only** for fully wildtype sequence. 

In other words, if we have a sequence with a single mutation ``GGC31AGT``, it only contributes to the single count of the mutation to ``AGT`` at position 31. But a wildtype sequence contributes to the wildtype counts at all positions. Therefore, the statistic at each site give how that single mutant fared relative to the wildtype sequence.

This file is suitable for use by :ref:`dms_inferprefs` or :ref:`dms_inferdiffprefs` regardless of whether your library consisted of only singly mutated genes or singly and multiply mutated genes.

All-mutant counts file
+++++++++++++++++++++++++
This file has the prefix specified by ``outprefix`` followed by the suffix ``_allmutcounts.txt``. This file is in the format of a :ref:`dms_counts` files. The mutation counts here include all mutations at all sites. 

In other words, a sequence with a single mutation ``GGC31AGT`` contributes a count of a mutation to ``AGT`` at site 31 **and** a wildtype count at all other positions. A wildtype sequence contributes a wildtype count at all positions. A double mutant ``GAG40ATT,AAA46AAT`` contributes a count of a mutation to ``ATT`` at site 40, a mutation to ``AAT`` at site 46, and a wildtype count at all other positions.

This file is suitable for use by :ref:`dms_inferprefs` or :ref:`dms_inferdiffprefs` only if your library consists of singly and multiply mutated clones with a Poisson distribution of the number of mutations per clone. If your library is all single mutants, you should use the `Single-mutant counts file`_ instead.

Variant counts file
+++++++++++++++++++++
This file has the prefix specified by ``outprefix`` followed by the suffix ``_variantcounts.txt``. 

After a header line beginning with ``#``, each line has three space-delimited entries:

    * The number of counts for a barcode (including barcodes with no counts).

    * The barcode itself.

    * The mutations separated by a comma but no space, or *no_mutations* if there are no mutations.

Here are a few example lines::

    #COUNTS BARCODE MUTATIONS
    5 TTACCATG no_mutations
    3 GATACATG GGA31GCT
    1 ATACCGAT GGA10AAA,ATA22GTG
    0 CCATAGAT no_mutations

Read statistics file
+++++++++++++++++++++
This file has the prefix specified by ``outprefix`` followed by the suffix ``_stats.txt``.

It gives the fate of all the read pairs in ``r1files`` and ``r2files``.

It has the following format (although the lines can be in an arbitrary order)::

    nreads = 1000
    nfiltered = 10
    nlowquality = 47
    nmismatches = 50
    nunrecognized = 43
    nretained = 850
    n0mut = 250
    n1mut = 309
    n2mut = 191
    n3mut = 73
    n4mut = 20
    n5mut = 6
    n6mut = 1

The meanings are:

    * ``nreads`` is the total number of read pairs

    * ``nfiltered`` is the number removed by the Illumina quality filter

    * ``nlowquality`` is the number that fail to meet the quality filter specified by ``--minquality``

    * ``nmismatches`` is the number that fail because they are mismatched at one or more positions.

    * ``nunrecognized`` is the number of barcodes that don't match one of the subassembled ones present in ``subassembled``.

    * ``nretained`` is the number retained.

    * The retained reads are then further categorized by how many mutations they have, giving ``n0mut``, ``n1mut``, etc for reads that map to unmutated, singly mutated, etc variants.

    * If ``--purgefrac`` is nonzero, then there is also a key ``nrandomlypurged`` indicating the number of reads purged.

.. include:: weblinks.txt
