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

   refseq
    This file should specify a valid in-frame coding sequence.

   alignspecs
    It is important to set ``alignspecs`` so that you don't count the part of the subamplicon that is in the primer binding site, since the nucleotide identities in this region come from the primers rather than the templates being sequences. Typically, *R2START* would be one greater than the length of the gene-binding region of the primer to avoid this.

    The alignments will fail if you don't set ``alignspecs`` exactly correctly, as the program only tries gapless alignment.

.. include:: weblinks.txt