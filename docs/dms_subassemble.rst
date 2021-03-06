.. _dms_subassemble:

==========================================
``dms_subassemble``
==========================================

.. contents::

Overview
-------------
This program can process the FASTQ reads generated by `Subassembly sequencing`_ to subassemble gene variants identified by a unique random barcode sequence at their 3' end. It provides linkage information across sites.

After you install `dms_tools`_, this program will be available to run at the command line.

Subassembly sequencing
------------------------------------
Subassembly is a technique to use short read sequencing to sequence longer sequences (e.g. full-length genes) by attached a unique barcode (a string of ``N`` nucleotides) to each variant of the sequence, and then associating short reads that span the entire sequence with this barcode to build up the full sequence. The technique is originally described in `Hiatt et al (2010)`_. 

The ``dms_subassemble`` program currently works for the case in which the unique barcode is at the 3' end of the gene, and the fragments to be subassembled are generated by PCR at defined locations. The R1 read captures the barcode, and the R2 read begins at one of the defined PCR locations and captures part of the variant's sequence.

Here we give an example of the experimental workflow to subassemble *recA* from *E. coli*.

In order to subassemble the gene, it is necessary to generate appropriately sized subamplicons that contain the barcode sequence at their 3' ends (if the subamplicons are too large, then they won't form clusters efficiently). To do this, the plasmids were cut with XbaI & NdeI or XbaI & BstEII to remove 640 bp or 587 bp upstream of the barcode, respectively. The digested plasmids were then recircularized and used as template for the subamplicons spanning nucleotide positions -1 through 545. Uncut plasmid was used as template for the remaining subamplicons, which span nucleotide positions 481 through the terminator and barcode sequences downstream of the recA gene. Each subamplicon is amplified with a common reverse primer annealing to an Illumina adapter sequence that has been incorporated into the plasmid immediately downstream of the barcode sequence, and with one of 7 different forward primers that anneal to different regions of the recA gene. The forward primers were designed such that the distance from the 3'-end of one primer is :math:`\le 200` bp from the 3'-end of the next forward primer. The template and primer combinations used to generate the subamplicons are shown below:: 

                          Rnd1F-1: 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTccggcatgacaggagtaaa-3’                                                                                                        Rnd1F172: 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTcttggggcaggtggtctg                                                                                                                    Rnd1F354: 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTtcgacaacctgctgtgctc
   XbaI & NdeI digested recA with coding sequence in caps: 5’-ccggtattacccggcatgacaggagtaaaaATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGTTCGCTTTCACTGGATATCGCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAATCTTCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAGCGTGAAGGTAAAACCTGTGCGTTTATCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTCGATATCGACAACCTGCTGTGCTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCTGGCGCGTTCTGGCGCAGTAGACGTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTCATActagaNNNNNNNNNNNNNNNNNNagatcggaagagcgtcgtgtagggaaagagtgt-3’
   XbaI & NdeI digested recA with coding sequence in caps: 3’-ggccataatgggccgtactgtcctcattttTACCGATAGCTGCTTTTGTTTGTCTTTCGCAACCGCCGTCGTGACCCGGTCTAACTCTTTGTTAAACCATTTCCGAGGTAGTACGCGGACCCACTTCTGGCAAGGTACCTACACCTTTGGTAGAGATGGCCAAGCGAAAGTGACCTATAGCGCGAACCCCGTCCACCAGACGGCTACCCGGCATAGCAGCTTTAGATGCCTGGCCTTAGAAGGCCATTTTGGTGCGACTGCGACGTCCACTAGCGGCGTCGCGTCGCACTTCCATTTTGGACACGCAAATAGCTACGACTTGTGCGCGACCTGGGTTAGATGCGTGCATTTGACCCGCAGCTATAGCTGTTGGACGACACGAGGGTCGGCCTGTGGCCGCTCGTCCGTGACCTTTAGACACTGCGGGACCGCGCAAGACCGCGTCATCTGCAATAGCAGCAACTGAGGCACCGCCGTGACTGCGGCTTTCGCCTTTAGCTTCCGCTTTAGCCGCTGAGAGTATgatctNNNNNNNNNNNNNNNNNNtctagccttctcgcagcacatccctttctcaca-5’
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      Rnd1R: 3’-TCTAGCCTTCTCGCAGCACATCCCTTTCTCACA-5’

                          Rnd1F-1: 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTccggcatgacaggagtaaa-3’                                                                                                        Rnd1F172: 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTcttggggcaggtggtctg                                                                                                                    Rnd1F354: 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTtcgacaacctgctgtgctc
 XbaI & BstEII digested recA with coding sequence in caps: 5’-ccggtattacccggcatgacaggagtaaaaATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGTTCGCTTTCACTGGATATCGCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAATCTTCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAGCGTGAAGGTAAAACCTGTGCGTTTATCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTCGATATCGACAACCTGCTGTGCTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCTGGCGCGTTCTGGCGCAGTAGACGTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTCATATGGGCCTTGCGGCACGTATGATGAGCCAGGCGATGCGTAAGCTGGCGGGTAACctagaNNNNNNNNNNNNNNNNNNagatcggaagagcgtcgtgtagggaaagagtgt-3’
 XbaI & BstEII digested recA with coding sequence in caps: 3’-ggccataatgggccgtactgtcctcattttTACCGATAGCTGCTTTTGTTTGTCTTTCGCAACCGCCGTCGTGACCCGGTCTAACTCTTTGTTAAACCATTTCCGAGGTAGTACGCGGACCCACTTCTGGCAAGGTACCTACACCTTTGGTAGAGATGGCCAAGCGAAAGTGACCTATAGCGCGAACCCCGTCCACCAGACGGCTACCCGGCATAGCAGCTTTAGATGCCTGGCCTTAGAAGGCCATTTTGGTGCGACTGCGACGTCCACTAGCGGCGTCGCGTCGCACTTCCATTTTGGACACGCAAATAGCTACGACTTGTGCGCGACCTGGGTTAGATGCGTGCATTTGACCCGCAGCTATAGCTGTTGGACGACACGAGGGTCGGCCTGTGGCCGCTCGTCCGTGACCTTTAGACACTGCGGGACCGCGCAAGACCGCGTCATCTGCAATAGCAGCAACTGAGGCACCGCCGTGACTGCGGCTTTCGCCTTTAGCTTCCGCTTTAGCCGCTGAGAGTATACCCGGAACGCCGTGCATACTACTCGGTCCGCTACGCATTCGACCGCCCATTGgatctNNNNNNNNNNNNNNNNNNtctagccttctcgcagcacatccctttctcaca-5’
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           Rnd1R: 3’-TCTAGCCTTCTCGCAGCACATCCCTTTCTCACA-5’

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          Rnd1F481: 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTggaaatcgaaggcgaaatc-3’                                                                                   Rnd1F633: 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTgcaacccggaaaccactac-3’                                                                                                          Rnd1F809: 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTcctctacggcgaaggtatca-3’                                                              Rnd1F941: 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTtgcctggctgaaagataacc-3’
                E. coli recA with coding sequence in caps: 5’-ccggtattacccggcatgacaggagtaaaaATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGTTCGCTTTCACTGGATATCGCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAATCTTCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAGCGTGAAGGTAAAACCTGTGCGTTTATCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTCGATATCGACAACCTGCTGTGCTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCTGGCGCGTTCTGGCGCAGTAGACGTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTCATATGGGCCTTGCGGCACGTATGATGAGCCAGGCGATGCGTAAGCTGGCGGGTAACCTGAAGCAGTCCAACACGCTGCTGATCTTCATCAACCAGATCCGTATGAAAATTGGTGTGATGTTCGGCAACCCGGAAACCACTACCGGTGGTAACGCGCTGAAATTCTACGCCTCTGTTCGTCTCGACATCCGTCGTATCGGCGCGGTGAAAGAGGGCGAAAACGTGGTGGGTAGCGAAACCCGCGTGAAAGTGGTGAAGAACAAAATCGCTGCGCCGTTTAAACAGGCTGAATTCCAGATCCTCTACGGCGAAGGTATCAACTTCTACGGCGAACTGGTTGACCTGGGCGTAAAAGAGAAGCTGATCGAGAAAGCAGGCGCGTGGTACAGCTACAAAGGTGAGAAGATCGGTCAGGGTAAAGCGAATGCGACTGCCTGGCTGAAAGATAACCCGGAAACCGCGAAAGAGATCGAGAAGAAAGTACGTGAGTTGCTGCTGAGCAACCCGAACTCAACGCCGGATTTCTCTGTAGATGATAGCGAAGGCGTAGCAGAAACTAACGAAGATTTTTAAtcgtcttgtttgatacacaagggtcgcatctgcggcccttttgcttttttaagttgtaaggatatgccattctagaNNNNNNNNNNNNNNNNNNagatcggaagagcgtcgtgtagggaaagagtgt-3’
                E. coli recA with coding sequence in caps: 3’-ggccataatgggccgtactgtcctcattttTACCGATAGCTGCTTTTGTTTGTCTTTCGCAACCGCCGTCGTGACCCGGTCTAACTCTTTGTTAAACCATTTCCGAGGTAGTACGCGGACCCACTTCTGGCAAGGTACCTACACCTTTGGTAGAGATGGCCAAGCGAAAGTGACCTATAGCGCGAACCCCGTCCACCAGACGGCTACCCGGCATAGCAGCTTTAGATGCCTGGCCTTAGAAGGCCATTTTGGTGCGACTGCGACGTCCACTAGCGGCGTCGCGTCGCACTTCCATTTTGGACACGCAAATAGCTACGACTTGTGCGCGACCTGGGTTAGATGCGTGCATTTGACCCGCAGCTATAGCTGTTGGACGACACGAGGGTCGGCCTGTGGCCGCTCGTCCGTGACCTTTAGACACTGCGGGACCGCGCAAGACCGCGTCATCTGCAATAGCAGCAACTGAGGCACCGCCGTGACTGCGGCTTTCGCCTTTAGCTTCCGCTTTAGCCGCTGAGAGTATACCCGGAACGCCGTGCATACTACTCGGTCCGCTACGCATTCGACCGCCCATTGGACTTCGTCAGGTTGTGCGACGACTAGAAGTAGTTGGTCTAGGCATACTTTTAACCACACTACAAGCCGTTGGGCCTTTGGTGATGGCCACCATTGCGCGACTTTAAGATGCGGAGACAAGCAGAGCTGTAGGCAGCATAGCCGCGCCACTTTCTCCCGCTTTTGCACCACCCATCGCTTTGGGCGCACTTTCACCACTTCTTGTTTTAGCGACGCGGCAAATTTGTCCGACTTAAGGTCTAGGAGATGCCGCTTCCATAGTTGAAGATGCCGCTTGACCAACTGGACCCGCATTTTCTCTTCGACTAGCTCTTTCGTCCGCGCACCATGTCGATGTTTCCACTCTTCTAGCCAGTCCCATTTCGCTTACGCTGACGGACCGACTTTCTATTGGGCCTTTGGCGCTTTCTCTAGCTCTTCTTTCATGCACTCAACGACGACTCGTTGGGCTTGAGTTGCGGCCTAAAGAGACATCTACTATCGCTTCCGCATCGTCTTTGATTGCTTCTAAAAATTagcagaacaaactatgtgttcccagcgtagacgccgggaaaacgaaaaaattcaacattcctatacggtaagatctNNNNNNNNNNNNNNNNNNtctagccttctcgcagcacatccctttctcaca-5’
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      Rnd1R: 3’-TCTAGCCTTCTCGCAGCACATCCCTTTCTCACA-5’


After these subamplicons are generated in the Round 1 PCR, a second round of PCR is used to add the remaining Illumina adapter sequences and TruSeq index sequences for multiplexing. For example, here are two primers used in the Round 2 PCR reactions::

	>Rnd2IndexF: forward primer for the Round 2 PCR that adds the Illumina adaptor (uppercase) and a TruSeq index, shown here as 6 “n” nucleotides 
	5’-CAAGCAGAAGACGGCATACGAGATnnnnnngtgactggagttcagacgtgtgctcttcc-3’ 

	>Rnd2R: reverse primer for the Round 2 PCR that adds the Illumina adaptor (uppercase) 
	5’-AATGATACGGCGACCACCGAGATCTacactctttccctacacgacgctcttccgatct-3’ 


The products of the Round 2 PCR are sequenced with asymmetric paired-end reads with Read1 sequencing the barcode (:math:`\ge 18` nucleotides), and with Read2 of sufficient length to overlap with the next subamplicon (:math:`\ge 200` nucleotides). The sequencing data are analyzed with the standard Illumina pipeline to generate FASTQ files containing the R1 and R2 reads.

Subassembly algorithm
-----------------------
``dms_subassemble`` implements the following algorithm:

#. The paired *R1* and *R2* reads are read from ``r1files`` and ``r2files``, and any nucleotide that has a Q score < ``--minq`` is converted to an ``N`` (e.g. considered ambiguous).

#. The *R1* read is only processed to extract the barcode (the first ``--barcodelength`` nucleotides) and any remaining nucleotides in *R1* are discarded. The read pair is then discarded if:

    * There are any ``N`` nucleotides in the barcode.
    
    * The fraction of ``N`` nucleotides in the trimmed *R2* read is > ``--maxlowqfrac``.

#. An attempt is made to gaplessly align the trimmed *R2* read at each site specified by ``alignspecs``. Briefly, ``alignspecs`` specifies pairs of numbers *REFSEQSTART* and *R2START*. For *R2* read, we try to align the read starting at nucleotide *R2START* (1, 2, ... numbering) at nucleotide *REFSEQSTART* in ``refseq``. If ``--trimR2`` has its default value of *auto*, we only try to align up to where the next subamplicon would be (effectively trimming the unneeded part of *R2* from its 3' end). If the read gaplessly aligns with no more than ``--maxmuts`` mutations of character type ``--chartype``, the alignment is considered successful. Read that fail to align are written ``<outprefix>_unaligned.txt`` depending on the value of ``--no_write_unaligned`` and then discarded. Otherwise, reads are retained.

#. After collecting all alignable trimmed *R2* reads for a barcode, we then see if we have enough coverage to subassemble the gene. We call identities in the gene only if the character (which may be a codon character rather than a nucleotide character, see ``--chartype``) has at least ``--minreadspersite`` non-ambiguous (not ``N``) reads covering it and :math:`\ge` ``--minreadconcurrence`` of this reads concur on its identity.

Command-line usage
---------------------
.. argparse::
   :module: parsearguments
   :func: SubassembleParser
   :prog: dms_subassemble

   outprefix
    See `Output files`_ for a list of the created files.

   refseq
    This file should specify a valid in-frame coding sequence.

   alignspecs
    It is important to set ``alignspecs`` so that you don't count the part of the subamplicon that is in the primer binding site, since the nucleotide identities in this region come from the primers rather than the templates being sequences. Typically, *R2START* would be one greater than the length of the gene-binding region of the primer to avoid this.

    The alignments will fail if you don't set ``alignspecs`` exactly correctly, as the program only tries gapless alignment.

   \-\-no_write_barcode_reads
    If you set this option, then the program does **not** create the `All reads by barcode file`_. If you are not debugging, you may not want this file as it is very large.

Output files
--------------
The following output files are created. Each has the prefix specified by ``outprefix`` with the following suffixes.

Log file
++++++++++
This file has the suffix ``.log`` and tracks the progress of the program.

Subassembled variants file
+++++++++++++++++++++++++++
This file has the suffix ``_subassembled_variants.txt``. It is the primary output file, and lists all barcodes that can be subassembled according to the parameters passed to ``dms_subassemble``.

It is in the :ref:`subassembly_file` format. Note that for codon sequences (``--chartype`` of ``codon``), the mutations are numbered according to the **codon** position in 1, 2, ... numbering, not the nucleotide position. 

Here are a few example lines::

    TATTACATCTGCCCCCAA ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGTTCGCTTTCACTGGATATCGCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAATCTTCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAGCGTGAAGGTAAAACCTGTGCGTTTATCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGCAGTCGATATCGACAACCTGCTGTGCTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCTGGCGCGTTCTGGCGCAGTAGACGTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTCATATGGGCCTTGCGGCACGTATGATGAGCCAGGCGATGCGTAAGCTGGCGGGTAACCTGAAGCAGTCCAACACGCTGCTGATCTTCATCAACCAGATCCGTATGAAAATTGGTGTGATGTTCGGCAACCCGGAAACCACTACCGGTGGTAACGCGCTGAAATTCTACGCCTCTGTTCGTCTCGACATCCGTCGTATCGGCGCGGTGAAAGAGGGCGAAAACGTGGTGGGTAGCGAAACCCGCGTGAAAGTGGTGAAGAACAAAATCGCTGCGCCGTTTAAACAGGCTGAATTCCAGATCCTCTACGGCGAAGGTATCAACTTCTACGGCGAACTGGTTGACCTGGGCGTAAAAGAGAAGCTGATCGAGAAAGCAGGCGCGTGGTACAGCTACAAAGGTGAGAAGATCGGTCAGGGTAAAGCGAATGCGACTGCCTGGCTGAAAGATAACCCGGAAACCGCGAAAGAGATCGAGAAGAAAGTACGTGAGTTGCTGCTGAGCAACCCGAACTCAACGCCGGATTTCTCTGTAGATGATAGCGAAGGCGTAGCAGAAACTAACGAAGATTTTTAA GGC109GCA
    AACTCTGTGTTCCCATCA ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGTTCGCTTTCACTGGATATCGCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAATCTTCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAGCGTGAAGGTAAAACCTGTGCGTTTATCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTCGATATCGACAACCTGCTGTGCTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCTGGCGCGTTCTGGCGCAGTAGACGTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTCATATGGGCCTTGCGGCACGTATGATGAGCCAGGCGATGCGTAAGCTGGCGGGTAACCTGAAGCAGTCCAACACGCTGCTGATCTTCATCAACCAGATCCGTATGAAAATTGGTGTGATGTTCGGCAACCCGGAAACCACTACCGGTGGTAACGCGCTGAAATTCTACGCCTCTGTTCGTCTCGACATCCGTCGTATCGGCGCGGTGAAAGAGGGCGAAAACGTGGTGGGTAGCGAAACCCGCGTGAAAGTGGTGAAGAACAAAATCGCTGCGCCGTTTAAACAGGCTGAATTCCAGATCCTCTACGGCGAAGGTATCAACTTCTACGGCGAACTGGTTGACCTGGGCGTAAAAGAGAAGCTGATCGAGAAAGCAGGCGCGTGGTACAGCTACAAAGGTGAGAAGATCGGTCAGGGTAAAGCGAATGCGACTGCCTGGCTGAAAGATAACCCGGAAACCGCGAAAGAGATCGAGAAGAAAGTACGTGAGTTGCTGCTGAGCAACCCGAACTCAACGCCGGATTTCTCTGTAGATGATAGCGAAGGCGTAGCAGAAACTAACGAAGATTTTTAA no_mutations
    GTTAACCGATCAACGCAA ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGTTCGCTTGCACTGGATATCGCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAATCTTCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAGCGTGAAGGTAAAACCTGTGCGTTTATCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTCGATATCGACAACCTGCTGTGCTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCTGGCGCGTTCTGGCGCAGTAGACGTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTCATATGGGCCTTGCGGCACGTATGATGAGCCAGGCGATGCGTAAGCTGGCGGGTAACCTGAAGCAGTCCAACACGCTGCTGATCTTCATCAACCAGATCCGTATGAAAATTGGTGTGATGTTCGGCAACCCGGAAACCACTACCGGTGGTAACGCGCTGAAATTCTACGCCTCTGTTCGTCTCGACATCCGTCGTATCGGGGCGGTGAAAGAGGGCGAAAACGTGGTGGGTAGCGAAACCCGCGTGAAAGTGGTGAAGAACAAAATCGCTGCGCCGTTTAAACAGGCTGAATTCCAGATCCTCTACGGCGAAGGTATCAACTTCTACGGCGAACTGGTTGACCTGGGCGTAAAAGAGAAGCTGATCGAGAAAGCAGGCGCGTGGTACAGCTACAAAGGTGAGAAGATCGGTCAGGGTAAAGCGAATGCGACTGCCTGGCTGAAAGATAACCCGGAAACCGCGAAAGAGATCGAGAAGAAAGTACGTGAGTTGCTGCTGAGCAACCCGAACTCAACGCCGGATTTCTCTGTAGATGATAGCGAAGGCGTAGCAGAAACTAACGAAGATTTTTAA TCA47GCA,GGC230GGG

Each line lists a barcode, then the sequence subassembled for that barcode, and finally any mutations relative to ``refseq``.

All reads by barcode file
++++++++++++++++++++++++++
This file has the suffix ``_all_reads_by_barcode.txt``. It is a **very large** text file that lists each read that matches each barcode (both those successfully subassembled and those that aren't). It also explains why or why not a barcode was subassembled.

This file is **not** created if you set the ``--no_write_barcode_reads``. 

Summary statistics file
+++++++++++++++++++++++
This file has the suffix ``_summarystats.txt``. It lists summary statistics about the subassembly. Here is an example::

    barcodes (total) = 642641
    barcodes successfully subassembled = 98281
    barcodes with at least one alignable read = 615204
    read pairs (total) = 29611710
    read pairs aligned at site 1 = 2442839
    read pairs aligned at site 172 = 1944582
    read pairs aligned at site 354 = 2794256
    read pairs aligned at site 481 = 2352744
    read pairs aligned at site 633 = 2982037
    read pairs aligned at site 809 = 3024340
    read pairs aligned at site 941 = 2758838
    read pairs purged due to low quality = 8965508
    read pairs that are alignable = 18299636
    read pairs that are alignable and map to a subassembled barcode = 8068574
    read pairs that are unalignable = 2346566
    read pairs that fail Illumina filter = 0
    sites with insufficient concurrence due to mismatch between mutant and wildtype characters = 52936
    sites with insufficient concurrence due to mismatch between two mutant characters = 68862

Alignable reads per barcode file
++++++++++++++++++++++++++++++++++
This file has the suffix ``_alignablereadsperbarcode.txt``. It gives the distribution of the number of alignable reads per barcode. Here is an example of the first few lines::

    nreads  nbarcodes
    0   27437
    1   153500
    2   10580
    3   8714
    4   8670
    5   8825
    6   8991
    7   8958

Mutations among subassembled variants file
+++++++++++++++++++++++++++++++++++++++++++
This file has the suffix ``_nmuts_among_subassembled.txt``. It gives the distribution of the number of mutations per variant among subassembled variants. Here is an example::

    nmuts   nvariants
    0   35614
    1   32442
    2   17902
    3   7705
    4   3046
    5   1113
    6   316
    7   107
    8   33
    9   2
    10  1

Read start sites file
++++++++++++++++++++++
This file has the suffix ``_refseqstarts.txt``. It gives the number of reads that start at each of the positions in ``refseq`` specified in ``alignspecs``. Here is an example::

    refseqstart nreads
    1   2442839
    172 1944582
    354 2794256
    481 2352744
    633 2982037
    809 3024340
    941 2758838


.. include:: weblinks.txt
