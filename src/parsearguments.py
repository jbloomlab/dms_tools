"""
===========================
``parsearguments`` module
===========================

Module that parses command-line arguments for the scripts in ``dms_tools``.

Written by Jesse Bloom.

Defined in this module
----------------------------------

* *ArgumentParserNoArgHelp* : argument parser that prints help when no args

* *FloatGreaterThanZero* : parses whether string is float greater than zero.

* *FloatBetweenZeroAndOne* : parses whether is float between 0 and 1.

* *FloatBetweenHalfAndOne* : parses whether is float > 0.5 and <= 1.

* *NonNegativeInt* : parses whether string is non-negative integer.

* *IntGreaterEqual1* : parses whether string is integer >= 1.

* *IntGreaterEqual2* : parses whether string is integer >= 2.

* *ExistingFile* : parses whether string gives an existing file.

* *PDF_File* : parses whether string gives a PDF file name.

* *CommaSeparatedFASTQFiles* : parses list of FASTQ files.

* *AlignSpecs* : parses alignment specs for *BarcodedSubampliconsParser*

* *R2AlignSpecs* : parses alignment specs for *SubassembleParser*

* *AlignprefixName* : parses alignment prefixes/names for *SummarizeAlignmentsParser*

* *InferPrefsParser* : parser for ``dms_inferprefs``

* *InferDiffPrefsParser* : parser for ``dms_inferdiffprefs``

* *DiffSelectionParser* : parser for ``dms_diffselection``

* *MergeParser* : parser for ``dms_merge``

* *EditSitesParser* : parsers for ``dms_editsites``

* *CorrelateParser* : parser for ``dms_correlate``

* *LogoPlotParser* : parser for ``dms_logoplot``

* *BarcodedSubampliconsParser* : parser for ``dms_barcodedsubamplicons``

* *SubassembleParser* : parser for ``dms_subassemble``

* *MatchSubassembledBarcodesParser* : parser for ``dms_matchsubassembledbarcodes``

* *SummarizeAlignmentsParser* : parser for ``dms_summarizealignments``

* *SmartHelpFormatter* : new class for formatting argument helps.


Detailed documentation
----------------------------------

"""


import sys
import os
import re
import argparse
import dms_tools


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""
    def error(self, message):
        """Prints error message, then help."""
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


class SmartHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Raw text for help beginning with ``R|``, *ArgumentDefaultsHelpFormatter* otherwise.

    Based on *SmartFormatter* described at http://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text.

    This option allows you to specify exact help format by preceding ``R|``.
    """
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        else:
            return argparse.ArgumentDefaultsHelpFormatter._split_lines(self, text, width)


def IntGreaterEqual1(n):
    """If *n* is integer >= 1 returns it, otherwise an error."""
    if not isinstance(n, str):
        raise argparse.ArgumentTypeError('%r is not a string' % n)
    try:
        n = int(n)
    except:
        raise argparse.ArgumentTypeError('%s is not an integer' % n)
    if n < 1:
        raise argparse.ArgumentTypeError('%d is not greater or equal to one' % n)
    else:
        return n


def IntGreaterEqual2(n):
    """If *n* is integer >= 2 returns it, otherwise an error."""
    if not isinstance(n, str):
        raise argparse.ArgumentTypeError('%r is not a string' % n)
    try:
        n = int(n)
    except:
        raise argparse.ArgumentTypeError('%s is not an integer' % n)
    if n < 2:
        raise argparse.ArgumentTypeError('%d is not greater or equal to two' % n)
    else:
        return n


def IntGreaterEqual1(n):
    """If *n* is integer >= 1 returns it, otherwise an error."""
    if not isinstance(n, str):
        raise argparse.ArgumentTypeError('%r is not a string' % n)
    try:
        n = int(n)
    except:
        raise argparse.ArgumentTypeError('%s is not an integer' % n)
    if n < 1:
        raise argparse.ArgumentTypeError('%d is not greater or equal to one' % n)
    else:
        return n


def NonNegativeInt(n):
    """If *n* is non-negative integer returns it, otherwise an error.

    >>> print "%d" % NonNegativeInt('8')
    8

    >>> NonNegativeInt('8.1')
    Traceback (most recent call last):
       ...
    ArgumentTypeError: 8.1 is not an integer

    >>> print "%d" % NonNegativeInt('0')
    0

    >>> NonNegativeInt('-1')
    Traceback (most recent call last):
       ...
    ArgumentTypeError: -1 is not non-negative

    """
    if not isinstance(n, str):
        raise argparse.ArgumentTypeError('%r is not a string' % n)
    try:
       n = int(n)
    except:
        raise argparse.ArgumentTypeError('%s is not an integer' % n)
    if n < 0:
        raise argparse.ArgumentTypeError('%d is not non-negative' % n)
    else:
        return n


def FloatGreaterThanZero(x):
    """If *x* is string for float > 0, returns it, otherwise an error.

    Designed based on this: http://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin

    >>> print "%.3f" % FloatGreaterThanZero('0.1')
    0.100

    >>> FloatGreaterThanZero('0.0')
    Traceback (most recent call last):
        ...
    ArgumentTypeError: 0.0 not a float greater than zero

    >>> FloatGreaterThanZero('hi')
    Traceback (most recent call last):
        ...
    ValueError: could not convert string to float: hi
    """
    x = float(x)
    if x > 0:
        return x
    else:
        raise argparse.ArgumentTypeError("%r not a float greater than zero" % x)

def FloatBetweenZeroAndOne(x):
    """If *x* is float >= 0 or <= 1, returns it, otherwise error."""
    x = float(x)
    if 0 <= x <= 1:
        return x
    else:
        raise argparse.ArgumentTypeError("%r is not a float >= 0 and <= 1" % x)

def FloatBetweenHalfAndOne(x):
    """If *x* is float > 0.5 or <= 1, returns it, otherwise error."""
    x = float(x)
    if 0.5 < x <= 1:
        return x
    else:
        raise argparse.ArgumentTypeError("%r is not a float > 0.5 and <= 1" % x)

def ExistingFile(fname):
    """If *fname* is name of an existing file return it, otherwise an error.
    
    It is also acceptable for *fname* to be the string "none"."""
    if os.path.isfile(fname):
        return fname
    elif fname.lower() == 'none':
        return None
    else:
        raise argparse.ArgumentTypeError("%s does not specify a valid file name" % fname)


def PDF_File(fname):
    """If *fname* ends in PDF extension, returns it. Otherwise an error."""
    if os.path.splitext(fname)[1].lower() == '.pdf':
        return fname
    else:
        raise argparse.ArgumentTypeError('Not a valid PDF file name; should end in ".pdf":\n%s' % fname)


def CommaSeparatedFASTQFiles(files):
    """Returns name of comma-separated FASTQ files.

    Checks that all files exist. Returns error if files
    cannot be found or if they are not .fastq or .fastq.gz
    Otherwise returns list of files."""
    if re.search('\s', files):
        raise argparse.ArgumentTypeError("List of FASTQ files contains whitespace: %s" % files)
    if not files:
        raise argparse.ArgumentTypeError("Empty list of FASTQ files")
    files = files.split(',')
    for f in files:
        if not os.path.isfile(f):
            raise argparse.ArgumentTypeError("Cannot find FASTQ file %s" % f)
        if not (os.path.splitext(f)[1].lower() == '.fastq' or (os.path.splitext(os.path.splitext(f)[0])[1].lower() == '.fastq' and os.path.splitext(f)[1].lower() == '.gz')):
            raise argparse.ArgumentTypeError("FASTQ file %s does not end in extension .fastq or .fastq.gz")
    return files

def AlignprefixName(alignprefixname):
    """Parses prefix / name for *SummarizeAlignmentsParser*.

    >>> AlignprefixName('replicate_1/replicate_1_mutDNA_,mutDNA')
    ('replicate_1/replicate_1_mutDNA_', 'mutDNA')
    """
    tup = alignprefixname.split(',')
    if len(tup) != 2:
        raise argparse.ArgumentTypeError("Failed to find two comma-delimited strings in %s" % alignprefixname)
    return tuple(tup)


def AlignSpecs(alignspecs):
    """Parses alignment specs for *BarcodedSubampliconsParser*.

    >>> AlignSpecs('1,300,10,12')
    (1, 300, 10, 12)
    """
    aligntup = alignspecs.split(',')
    if len(aligntup) != 4:
        raise argparse.ArgumentTypeError("Failed to find four comma-delimited integers in alignspecs %s" % alignspecs)
    try:
        aligntup = [int(n) for n in aligntup]
    except ValueError:
        raise argparse.ArgumentTypeError("Non-integer in alignspecs %s" % alignspecs)
    if not all([n >= 1 for n in aligntup]):
        raise argparse.ArgumentTypeError('Not all integers >= 1 in alignspecs %s' % alignspecs)
    if not aligntup[0] < aligntup[1]:
        raise argparse.ArgumentTypeError("First entry must be less than second in alignspecs %s" % alignspecs)
    return tuple(aligntup)


def R2AlignSpecs(alignspecs):
    """Parses alignment specs for *SubasssembleParser*.

    >>> R2AlignSpecs('1,300')
    (1, 300)
    """
    aligntup = alignspecs.split(',')
    if len(aligntup) != 2:
        raise argparse.ArgumentTypeError("Failed to find two comma-delimited integers in alignspecs %s" % alignspecs)
    try:
        aligntup = [int(n) for n in aligntup]
    except ValueError:
        raise argparse.ArgumentTypeError("Non-integer in alignspecs %s" % alignspecs)
    if not all([n >= 1 for n in aligntup]):
        raise argparse.ArgumentTypeError('Not all integers >= 1 in alignspecs %s' % alignspecs)
    return tuple(aligntup)


def SummarizeAlignmentsParser():
    """Returns *argparse.ArgumentParser* for ``dms_summarizealignments``."""
    parser = ArgumentParserNoArgHelp(description='Makes plots that summarize alignments for one or more samples. Designed to be run after you have used another program to make your alignments, and now you want to visualize the results. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outprefix', help="Prefix for the output PDF plot files. Any existing files with these names are overwritten.")
    parser.add_argument('alignment_type', choices=['barcodedsubamplicons', 'subassemble', 'matchsubassembledbarcodes'], help="The type of data being summarized. Use 'barcodedsubamplicons' if output is from 'dms_barcodedsubamplicons', 'subassemble' if output is from 'dms_subassemble', 'matchsubassembledbarcodes' if output is from 'dms_matchsubassembledbarcodes'.")
    parser.add_argument('alignments', nargs='+', help="This argument is repeated to specify each alignment that is being summarized. Each repetition is two comma-delimited strings (no spaces): ALIGNPREFIX,NAME. ALIGNPREFIX is the value of '--outprefix' used when calling the alignment program and NAME is the name assigned to that sample in the summaries. We expect to find all of the alignment output files with ALIGNPREFIX; these are the input data.", type=AlignprefixName)
    parser.add_argument('--chartype', default='codon', choices=['codon'], help='Character type used for the alignments / mutation counting.')
    parser.add_argument('--maxmutcounts', type=NonNegativeInt, default=20, help="Maximum x-value for cumulative counts plots plots.")
    parser.add_argument('--maxperbarcode', type=NonNegativeInt, default=3, help="In 'barcodes.pdf' plot, group all barcodes with >= this many reads.")
    parser.set_defaults(writemutfreqs=False)
    parser.add_argument('--writemutfreqs', dest='writemutfreqs', action='store_true', help="Write a file 'mutfreqs.txt' that gives the numerical values plotted in 'mutfreqs.pdf'?")
    parser.set_defaults(groupbyfirst=False)
    parser.add_argument('--groupbyfirst', dest='groupbyfirst', action='store_true', help="Make separate 'mutdepth' and 'cumulfracs' plots for each name prefix (first word before underscore).")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def MatchSubassembledBarcodesParser():
    """Returns *argparse.ArgumentParser* for ``dms_matchsubassembledbarcodes``."""
    parser = ArgumentParserNoArgHelp(description='Extract barcodes after quality filter, match to subassembled variants, write variant counts files, make summary plots. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outprefix', help='Prefix for output files.')
    parser.add_argument('subassembled', type=ExistingFile, help="File output by 'dms_subassemble' that contains the barcoded subassembled variants.")
    parser.add_argument('r1files', type=CommaSeparatedFASTQFiles, help='Comma-separated list of R1 FASTQ files (no spaces) containing barcodes. Files can optionally be gzipped (extension .gz).')
    parser.add_argument('r2files', type=CommaSeparatedFASTQFiles, help="Like 'r1files' but for R2 read.")
    parser.add_argument('--minquality', type=int, default=25, help="Each site in barcode must have quality >= this for at least one read.")
    parser.add_argument('--r1start', type=IntGreaterEqual1, default=1, help="Barcode starts at this position in R1.")
    parser.add_argument('--r2end', type=NonNegativeInt, default=0, help="Barcode ends this many nucleotides before end of R2.")
    parser.add_argument('--purgefrac', type=FloatBetweenZeroAndOne, default=0, help="Randomly purge barcode read pairs in 'r1files' and 'r2files' with this probability (subsample the data).")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def SubassembleParser():
    """Returns *argparse.ArgumentParser* for ``dms_subassemble``."""
    parser = ArgumentParserNoArgHelp(description='Subassemble barcoded variants with linkage from short-read sequences. R1 has the barcode, R2 has the sequence fragment. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outprefix', help='Prefix for output files.')
    parser.add_argument('refseq', type=ExistingFile, help='Existing FASTA file containing wildtype gene we are subassembling.')
    parser.add_argument('r1files', type=CommaSeparatedFASTQFiles, help='Comma-separated list of R1 FASTQ files (no spaces). Files can optionally be gzipped (extension .gz).')
    parser.add_argument('r2files', type=CommaSeparatedFASTQFiles, help="Like 'r1files' but for R2. Must be same number of comma-separated entires as for 'r1files'.")
    parser.add_argument('alignspecs', nargs='+', help="This argument is repeated to specify each possible alignment location for R2. Each specification is two comma-delimited integers (no spaces): 'REFSEQSTART,R2START'. REFSEQSTART is nucleotide (1, 2, ... numbering) in 'refseq' where nucleotide R2START in R2 aligns.", type=R2AlignSpecs)
    parser.add_argument('--barcodelength', type=NonNegativeInt, default=18, help='Length of barcode (NNN...) which starts at beginning of R1.')
    parser.add_argument('--trimR2', choices=['auto', 'none'], default='auto', help="Trim R2 read? If 'auto', trim until start of next subamplicon specified by 'alignspecs'; if 'none' no trimming.")
    parser.add_argument('--minq', type=NonNegativeInt, default=20, help='Nucleotides with Q scores < this number are converted to N.')
    parser.add_argument('--maxlowqfrac', default=0.15, type=FloatBetweenZeroAndOne, help='Only retain if fraction of N nucleotides in R2 read <= this.')
    parser.add_argument('--maxmuts', type=NonNegativeInt, default=4, help='Only align read if <= this many mismatches with "refseq" counted in terms of "chartype".')
    parser.add_argument('--minreadspersite', default=2, type=IntGreaterEqual1, help='Call site only barcodes when >= this many reads give it a non-ambiguous identity.')
    parser.add_argument('--minreadconcurrence', default=0.75, type=FloatBetweenHalfAndOne, help="Only call sites when >= this fraction of reads concur.")
    parser.add_argument('--chartype', default='codon', choices=['codon'], help='Character for which we are counting mutations. Currently "codon" is only allowed value (in the future "nucleotide" might be added).')
    parser.set_defaults(no_write_barcode_reads=False)
    parser.add_argument('--no_write_barcode_reads', dest='no_write_barcode_reads', action='store_true', help="Don't write all barcodes and assigned reads to a file (saves space/time to not do this).")
    parser.add_argument('--purgefrac', type=FloatBetweenZeroAndOne, default=0, help='Randomly purge reads with this probability (subsample the data).')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def BarcodedSubampliconsParser():
    """Returns *argparse.ArgumentParser* for ``dms_barcodedsubamplicons`` script."""
    parser = ArgumentParserNoArgHelp(description='Gathers barcoded subamplicons, aligns to reference sequence, and counts mutations. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outprefix', help=\
        'Prefix for output files. ' +\
        'The suffixes of the created files are: ' +\
        '"counts.txt" - character counts at each site; ' +\
        '"summarystats.txt" - barcode summary stats; ' +\
        '".log" - file that logs progress; ' +\
        '"barcodeinfo.txt.gz" - lists all barcodes if using "--barcodeinfo". '+\
        'Existing files with these names are removed.')
    parser.add_argument('refseq', type=ExistingFile, help='Existing FASTA file containing gene to which we are aligning subamplicons and counting mutations.')
    parser.add_argument('r1files', type=CommaSeparatedFASTQFiles, help='Comma-separated list of R1 FASTQ files (no spaces). Files can optionally be gzipped (extension .gz).')
    parser.add_argument('r2files', type=CommaSeparatedFASTQFiles, help='Like "r1files" but R2 files. Must be same number of comma-separated entries as for "r1files".')
    parser.add_argument('alignspecs', nargs='+', help='This argument is repeated to specify each subamplicon. Each occurrence is four comma-delimited integers (no spaces): "REFSEQSTART,REFSEQEND,R1START,R2START". REFSEQSTART is nucleotide (1, 2, ... numbering) in "refseq" where nucleotide R1START in R1 aligns. REFSEQEND is nucleotide in "refseq" where nucleotide R2START in R2 aligns.', type=AlignSpecs)
    parser.add_argument('--barcodelength', type=NonNegativeInt, default=8, help='Length of barcodes (NNN... nucleotides) at the beginning of R1 and R2 reads.')
    parser.add_argument('--chartype', default='codon', choices=['codon'], help='Character for which we are counting mutations. Currently "codon" is only allowed value (in the future "nucleotide" might be added).')
    parser.add_argument('--maxmuts', type=NonNegativeInt, default=4, help='Only align subamplicons (consensus from a barcode) if <= this many mismatches with "refseq" counted in terms of "chartype".')
    parser.add_argument('--minq', type=NonNegativeInt, default=15, help='Only consider nucleotides with Q scores >= this number.')
    parser.add_argument('--maxlowqfrac', default=0.075, type=FloatBetweenZeroAndOne, help='Only retain read pairs if no "N" or Q < "minq" nucleotides in barcodes, and total fraction of such nucleotides is <= this number in each read individually and in the eventual subamplicon built from the barcodes.')
    parser.add_argument ('--minreadsperbarcode', type=IntGreaterEqual2, default=2, help='Retain only barcodes with >= this many reads that align gaplessly with >= "minreadidentity" identical high-quality nucleotides; should be >= 2.')
    parser.add_argument('--minreadidentity', default=0.9, type=FloatBetweenZeroAndOne, help='Retain only barcodes where all reads have >= this fraction of identical high-quality nucleotides that align gaplessly.')
    parser.add_argument('--minreadconcurrence', default=0.75, type=FloatBetweenHalfAndOne, help='For retained barcodes, only make mutation calls when >= this fraction of reads concur.')
    parser.add_argument('--maxreadtrim', type=NonNegativeInt, default=3, help="If R1 or R2 reads for same barcode are not all same length, trim up to this many nucleotides from 3' end; if still not same length then discard.")
    parser.add_argument('--R1_is_antisense', dest='R1_is_antisense', action='store_true', help='Is the R1 read in the antisense direction of "refseq"?')
    parser.set_defaults(R1_is_antisense=False)
    parser.add_argument('--R1trimlength', type=IntGreaterEqual1, help="Before performing any operations on the read pairs, trim nucleotides away from 3' end of R1 read to reach the read length specified here.")
    parser.add_argument('--R2trimlength', type=IntGreaterEqual1, help="Before performing any operations on the read pairs, trim nucleotides away from 3' end of R2 read to reach the read length specified here.")
    parser.add_argument('--barcodeinfo', dest='barcodeinfo', action='store_true', help='If you specify this option, create a file with suffix "barcodeinfo.txt.gz" containing information for each barcode. This file is quite large, and its creation will about double the program run time.');
    parser.add_argument('--purgereadfrac', type=FloatBetweenZeroAndOne, default=0, help='Randomly purge read pairs with this probability, thereby subsampling the data. You might want a value > 0 to estimate how the results depend on the sequencing depth. See also "--purgefrac".')
    parser.add_argument('--purgefrac', type=FloatBetweenZeroAndOne, default=0, help='Like "--purgereadfrac" but purges barcodes rather than read pairs. In deciding which option to use, remember that the number of unique barcodes sequenced does not necessarily increase linearly with read depth.')
    parser.add_argument('--seed', type=int, default=1, help='Random number seed used to select reads for purging when using "--purgereadfrac" or "--purgefrac".')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    parser.set_defaults(barcodeinfo=False)
    return parser


def LogoPlotParser():
    """Returns *argparse.ArgumentParser* for ``dms_logoplot`` script."""
    parser = ArgumentParserNoArgHelp(description='Make a logo plot visually displaying the preferences or differential preferences. Utilizes weblogo (https://code.google.com/p/weblogo/). This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', help='Existing file giving preferences, differential preferences, or differential selection. The program auto-detects which type of data is present. The values can be for DNA or amino acids (stop codons allowed, indicated by "*").', type=ExistingFile)
    parser.add_argument('logoplot', help='Name of created file containing the logo plot, must end in the extension ".pdf". Overwritten if it already exists.', type=PDF_File)
    parser.add_argument('--nperline', type=int, default=60, help='Put this many sites per line of the logo plot.')
    parser.add_argument('--numberevery', type=int, default=10, help='Number x-axis ticks for sites at this interval.')
    parser.add_argument('--restrictdiffsel', default='None', help="Specify 'positive' or 'negative' to restrict plotted differential selection to positive or negative selection. Only meaningful if 'infile' gives differential selection.", choices=['None', 'positive', 'negative'])
    parser.add_argument('--diffselheight', type=ExistingFile, nargs='*', help="Only meaningful if 'infile' gives differential selection. List other differential selection files, and y-axis limits will be set to the max and min of 'infile' and all files listed here. This is useful if you want to make multiple differential selection plots with the same y-axis limits.")
    parser.set_defaults(nosepline=False)
    parser.add_argument('--nosepline', dest='nosepline', action='store_true', help='Do not plot a black center line separating positive and negative values of differential selection or preferences.')
    parser.add_argument('--diffprefheight', type=float, default=1.0, help='This option is only meaningful if "infile" gives differential preferences. In that case, it gives the height of logo stacks (extends from minus this to plus this). Cannot be smaller than the maximum total differential preference range.')
    parser.add_argument('--excludestop', dest='excludestop', action='store_true', help='If we are using amino acids, do we remove stop codons (denoted by "*")? We only remove stop codons if this argument is specified. If this option is used, then data for stop codons is removed by re-normalizing preferences to sum to one, and differential preferences to sum to zero.')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--colormap', type=str, default='jet', help='Colormap for amino-acid hydrophobicity or molecular weight. Must specify a valid ``pylab`` color map. Only meaningful if "mapmetric" is "kd" or "mw".')
    parser.add_argument('--letterheight', default=1, type=FloatGreaterThanZero, help="Relative parameter indicating the height of the letter stacks in the logo plots.")
    parser.add_argument('--overlay1', default=None, nargs=3, metavar=('FILE', 'SHORTNAME', 'LONGNAME'), help='Specify an overlay bar above each line of the logo plot to illustrate a per-residue property such as relative solvent accessibility or secondary structure. Requires three arguments: FILE SHORTNAME LONGNAME. FILE is the name of an existing file. Except for comment lines beginning with "#", each line should have two whitespace delimited columns (additional columns are allowed but ignored). The first column gives the site number (matching that in "infile") and the second column giving the property for this site; properties must either all be non-whitespace strings giving a discrete category (such as secondary structure), or all be numbers (such as relative solvent accessibility). All listed sites must be in "infile", but not all sites in "infile" must be in FILE -- missing sites are assumed to lack a known value for the property and are shown in white. SHORTNAME is a short (3-5 character) name of the property, such as "RSA" for "relative solvent accessibility." LONGNAME is a longer name (such as "relative solvent accessibiity"), or the same as SHORTNAME if you do not have a separate long name.')
    parser.add_argument('--overlay2', default=None, nargs=3, metavar=('FILE', 'SHORTNAME', 'LONGNAME'), help='Specify a second overlay bar. Arguments have the same meaning as for "overlay1".')
    parser.add_argument('--overlay3', default=None, nargs=3, metavar=('FILE', 'SHORTNAME', 'LONGNAME'), help='Specify a third overlay bar. Arguments have the same meaning as for "overlay1".')
    parser.add_argument('--stringencyparameter', type=FloatGreaterThanZero, help="Scale preferences by this stringency parameter; only valid when 'infile' specifies preferences.")
    parser.add_argument('--mapmetric', default='kd', help="Specify the amino-acid metric used to map colors to amino-acids in the logoplot. 'kd' uses the Kyte-Doolittle hydrophobicity scale, 'mw' uses molecular weight, 'functionalgroup' divides the amino acids into seven functional groups, and 'charge' uses charge at neutral pH. When using 'kd' or 'mw', 'colormap' is used to map colors to the metric; when using 'charge', a black/red/blue colormapping is used for neutral/positive/negative; similarly 'functionalgroup' has its own colormap.", choices=['kd','mw','charge','functionalgroup'])
    parser.add_argument('--overlay_cmap', default=None, help="Specify color map for overlay bars. Should be name of valid matplotlib colormap such as 'jet' or 'OrRd' (http://matplotlib.org/users/colormaps.html)")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def CorrelateParser():
    """Returns *argparse.ArgumentParser* for ``dms_correlate`` script."""
    parser = ArgumentParserNoArgHelp(description='Determine and plot the Pearson correlation between pairs of preferences, differential selections, or differential preferences. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file1', help='Existing file giving first set of preferences, differential selections (either site or mutation), or differential preferences.', type=ExistingFile)
    parser.add_argument('file2', help='Existing file giving second set of data; must be the same type of data "file1".', type=ExistingFile)
    parser.add_argument('outfileprefix', help='Prefix for name of created output files (these files are overwritten if they already exist). The correlations are written to a file with this prefix and the suffix ".txt". Unless you use the option "--noplot", a scatter plot is written to a file with this prefix and the suffix ".pdf". In the correlations text file, the first line gives the Pearson correlation coefficient in the format: "R = 0.5312". The second line gives the P-value in the format: "P = 0.0000131". The third line gives the number of points in the format: "N = 4200".')
    parser.add_argument('--name1', default='data 1', help='Name used for the data in "file1" in the correlation plot. The string specified here uses LaTex formatting; names with spaces should be enclosed in quotes.')
    parser.add_argument('--name2', default='data 2', help='Name used for the data in "file2" in the correlation plot. The string specified here uses LaTex formatting; names with spaces should be enclosed in quotes.')
    parser.add_argument('--excludestop', dest='excludestop', action='store_true', help='If we are using amino acids, do we remove stop codons (denoted by "*")? We only remove stop codons if this argument is specified. If this option is used, then data for stop codons is removed by re-normalizing preferences to sum to one, and differential preferences to sum to zero.')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--noplot', dest='noplot', action='store_true', help='Normally this script creates a PDF scatter plot. If this option is specified, then no such plot will be created.')
    parser.set_defaults(noplot=False)
    parser.add_argument('--alpha', default=0.1, help='The transparency (alpha value) for the points on the scatter plot. A value of 1.0 correspond to no transparency; values close to zero give high transparency. Transparency (alpha < 1) might be helpful if you have many points on top of each other.', type=float)
    parser.add_argument('--markersize', default=4, help='The size of the marker for the points on the scatter plot.')
    parser.add_argument('--plot_title', default='None', help='Title put at the top of the correlation plot. The string specified here uses LaTex formatting; names with spaces should be enclosed in spaces. A value of "None" means no title.')
    parser.add_argument('--corr_on_plot', dest='corr_on_plot', action='store_true', help='If this option is used, then the correlation coefficient will be visually displayed on the plot.')
    parser.set_defaults(corr_on_plot=False)
    parser.add_argument('--r2', dest='r2', action='store_true', help='If this option is used, the correlation coefficient displayed on the plot when using "--corr_on_plot" will show R-squared rather than the R value.')
    parser.set_defaults(r2=False)
    parser.add_argument('--rms_dpi', dest='rms_dpi', action='store_true', help='If "file1" and "file2" specify differential preferences, this argument specifies that we compute the correlation between the root-mean-square (RMS) differential preference at each site rather than between the differential preferences themselves.')
    parser.set_defaults(rms_dpi=False)
    parser.add_argument('--pref_entropy', dest='pref_entropy', action='store_true', help='If "file1" and "file2" specify preferences, this argument specifies that we compute the correlation between the site entropy of the preferences (logarithm base 2, so bits) at each site rather than between the preferences themselves.')
    parser.set_defaults(pref_entropy=False)
    parser.set_defaults(enrichment=False)
    parser.add_argument('--enrichment', dest='enrichment', action='store_true', help='If this option is set, we plot the enrichment ratio for all mutations on a log scale rather than plotting the preferences. The computed correlations are also then for the log-transformed enrichment ratios. The enrichment ratio for the wildtype identity is always one, and so is not included.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    parser.add_argument('--restrictdiffsel', default='None', help='Specify "positive" or "negative" to restrict plotted correlation in site differential selection to positive or negative selection. Only meaningful if file1 and file2 are sitediffsel files.', choices=['None', 'positive', 'negative'])
    return parser


def EditSitesParser():
    parser = ArgumentParserNoArgHelp(description='Edits sites in a data file. Typically you would use this program if you wanted to renumber sites or remove certain sites. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', type=ExistingFile, help='Existing data file. This could be a deep mutational scanning counts file, a preferences file, a differential preferences file, or any other file with the following format: blank lines or lines beginning with "#" (comment lines) are ignored; every other line must begin with an entry giving a unique site number (such as "1" or "2A"). The line may then have an arbitrary number of other entries separated from the site number by whitespace. If the lines have no whitespace, then we look for comma separators.')
    parser.add_argument('outfile', help='The created output file in which the site editing has been performed on "infile". If this output file already exists, it is overwritten.')
    parser.add_argument('edit_method', choices=['renumber', 'remove', 'retain'], help='How to do the editing: renumber sites, remove specified sites, or retain only specified sites.')
    parser.add_argument('edit_file', type=ExistingFile, help='Existing file specifying how edits are made. If "edit_method" is "renumber", then all non-comment lines (those not beginning with "#") must have two space delimited entries specifying the existing site in "infile" and the new site number with which it is replaced; all sites must be specified, and if the new number is "None" then the site is removed in the created file. If "edit_method" is "remove", then each line should have as its first entry a site, and all of the listed sites are removed. If "edit_method" is "retain", then each line should have as its first entry a site, and only the listed sites are retained.')
    parser.set_defaults(skipfirstline=False)
    parser.add_argument('--skipfirstline', action='store_true', dest='skipfirstline', help='Skip the edit operation on the first site. This could be helpful if dealing with a CSV file in pandas format.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def MergeParser():
    """Returns *argparse.ArgumentParser* for ``dms_merge`` script."""
    parser = ArgumentParserNoArgHelp(description='Merge preferences, differential selection, or differential preferences by averaging or adding / subtracting the values in multiple files. Alternatively, sum counts by adding the counts for the characters at each site across multiple files. All files must specify the same character type: can be nucleotide, codon, or amino acid (see "--excludestop" if using amino acids). This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outfile', help='Created output file with merged values; removed if it already exists. If the "infiles" do not all have the same wildtype residue at a site, then the wildtype is indicated as "?" in "outfile".')
    parser.add_argument('merge_method', help='How to merge: If "average" then all "infiles" must specify preferences, differential selection on mutations (mutdiffsel.txt file), or differential preferences; these are then averaged (by taking the mean) in "outfile". If "median", then all "infiles" must specify differential selection on mutations (mutdiffsel.txt file); these are then averaged by taking the median. If "sum" then "infiles" can either all be count files, all be differential preferences, or can be a combination of preferences and differential preferences that (along with any additional files specified by "--minus") sum to a total preference or a total differential preference of zero at each site. If "rescale", then only one infile can be specified and it must be preferences; these preferences are then scaled by the value provided by "stringencyparameter"', choices=['average', 'median', 'sum', 'rescale'])
    parser.add_argument('infiles', nargs='+', help='Files to average or sum. Must all have the same sites and character type, but do not need to have the same wildtype residue at each site. Note that for differential selection, must be the *mutdiffsel.txt file.', type=ExistingFile)
    parser.add_argument('--excludestop', dest='excludestop', action='store_true', help='If we are using amino acids, do we remove stop codons (denoted by "*")? We only remove stop codons if this argument is specified. If this option is used, then any files with stop codons have these codons removed (re-normalizing preferences to sum to one, and differential preferences to sum to zero) before the merge.')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--minus', nargs='+', help='Files to subtract when summing. Can only be used if "merge_method" is "sum" and if files are either preferences or differential preferences.', type=ExistingFile)
    parser.add_argument('--sitediffselfile', help="If merging differential selection, 'infiles' and 'outfile' are the '*mutdiffsel.txt' file. If you also want to create a '*sitdiffsel.txt' file for the merged '*mutdiffsel.txt' file, specify the name of that file here.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    parser.add_argument('--stringencyparameter', help='Stringency parameter used to rescale preferences. Can only be used if "merge_method" is "rescale".', type=FloatGreaterThanZero)
    parser.add_argument('--normalize', dest='normalize', action='store_true', help='Whether to normalize the counts at each site to the minimum number of counts observed at that site across all provided count files, before summing the counts. Can only be used if "merge_method" is "sum".')
    parser.set_defaults(normalize=False)
    parser.add_argument('--chartype', default='codon', help='Characters for which counts are summed: "DNA" = counts for DNA; "codon" = counts for codons; "aa" = counts for amino acids (possibly including stop codons, see "--excludestop"). This option only needs to be specified when summing count files.', choices=['DNA', 'codon', 'aa'])
    return parser


def InferPrefsParser():
    """Returns *argparse.ArgumentParser* for ``dms_inferprefs`` script."""
    parser = ArgumentParserNoArgHelp(description='Infer site-specific preferences for amino acids, nucleotides, or codons. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('n_pre', help='Existing file with pre-selection counts from deep sequencing. For nucleotides, header line should be "# POSITION WT A C G T" and then subsequent lines give wildtype and counts for each site (i.e. "1 G 1013 23 19 47"); for codons just use the 64 codons instead (i.e. "AAA AAC AAG ...").', type=ExistingFile)
    parser.add_argument('n_post', help='Existing file with post-selection counts from deep sequencing. Same format as "n_pre".', type=ExistingFile)
    parser.add_argument('outfile', help='Created output file with site-specific preferences; overwritten if it already exists.')
    parser.add_argument('--chartype', default='codon_to_aa', help='Characters for which preferences are inferred: "codon_to_aa" = counts for codons and inferred preferences for amino acids; "DNA" = counts and inferred preference for DNA; "codon" = counts and inferred preferences for codons; "aa" = counts and inferred preferences for amino acids (possibly including stop codons, see "--excludestop").', choices=['codon_to_aa', 'DNA', 'codon', 'aa'],)
    parser.add_argument('--errpre', default='none', help='Control used to estimate errors in counts from "n_pre". If "none" then the counts in "n_pre" are assumed to have no errors; otherwise this should be a counts file with the same format as "n_pre" giving the counts for sequencing unmutated gene. Currently if this option is not "none" then --errpost also cannot be "none".', type=ExistingFile)
    parser.add_argument('--errpost', default='none', help='Control used to estimate errors in counts from "n_post". If "none" then the counts in "n_post" are assumed to have no errors; otherwise this should be a counts file with the same format as "n_pre" giving the counts for sequencing unmutated gene. If you have the same error control for "n_pre" and "n_post", then set that file name for both this option and --errpre. Currently if this option is not "none" then --errpre also cannot be "none".', type=ExistingFile)
    parser.add_argument('--excludestop', help='Exclude stop codons as a possible amino acid if using "--chartype" of "codon_to_aa" or "aa". Stop codons are only excluded if this option is specified; otherwise they are included, and are required to be specified if "--chartype" is "aa".', dest='excludestop', action='store_true')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--pi_alpha', help='Concentration parameter for Dirichlet prior over preferences (pi).', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--mu_alpha', help='Concentration parameter for Dirichlet prior over mutagenesis rate (mu).', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--err_alpha', help='Concentration parameter for Dirichlet priors over error rates (epsilon, rho).', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--logfile', help='Log progress to this file; overwritten if it already exists.', default='Base name of "outfile" with extension ".log"')
    parser.add_argument('--ncpus', default=1, help='Number of CPUs to use; set to -1 to use all available CPUs.', type=int)
    parser.add_argument('--seed', default=1, help='Random number seed.', type=int)
    parser.add_argument('--sites', default=None, nargs='+', help='Only perform the inference for the specified sites, which should be a space separated list such as "--sites 1 2 10". All of these sites must have data in the counts files.')
    parser.add_argument('--ratio_estimation', default=None, metavar='MINCOUNTS', type=FloatGreaterThanZero, help='Rather than use MCMC to estimate the preferences using a statistical model, we simply compute them by calculating enrichment ratios relative to wildtype post- and pre-selection, and then normalizing these ratios to sum to one at each site. The single argument for this option is MINCOUNTS, which is a number > 0. If the counts for a character are less than MINCOUNTS, either due to low counts or error correction, then the counts for that character are changed to MINCOUNTS to avoid estimating ratios of zero, less than zero, or infinity.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def DiffSelectionParser():
    """Returns *argparse.ArgumentParser* for ``dms_diffselection`` script."""
    parser = ArgumentParserNoArgHelp(description=\
        'Differential selection between mock-treated and selected samples. ' +
        'Quantified as log2[{(n_sel_mut + f_sel * P) / (n_sel_wt + f_sel * P)} / {(n_mock_mut + f_mock * P) / (n_mock_wt + f_mock * P)}] where P is pseudocount, f_sel = max(1, N_sel / N_mock), and f_mock = max(1, N_mock / N_sel) where N_mock and N_sel are the sequencing depth for mock and selected libraries at that site. ' +
        'This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mockcounts', type=ExistingFile, help='File with counts from mock-selected library.  For nucleotides, header line should be "# POSITION WT A C G T" and then subsequent lines give wildtype and counts for each site (i.e. "1 G 1013 23 19 47"); for codons use the 64 codons instead (i.e. "AAA AAC AAG ...").')
    parser.add_argument('selectedcounts', type=ExistingFile, help='File with counts from selected library')
    parser.add_argument('outprefix', help="Prefix for output files. Suffixes are 'mutdiffsel.txt' and 'sitediffsel.txt'")
    parser.add_argument('--errorcontrolcounts', type=ExistingFile, help='Counts file for error control (e.g., wildtype gene or virus) to correct mutation counts.')
    parser.add_argument('--pseudocount', type=FloatGreaterThanZero, default=10, help='Pseudocount added to each count for selected library. By default, pseudocounts is added to library with smaller depth, and pseudcount for other library is scaled by relative depth at that site (see --no-scale-pseudocounts).')
    parser.add_argument('--mincounts', type=NonNegativeInt, default=0, help='Only report differential selection for mutations for which at least one of mock or selected counts is >= this number. Mutations that do not meet this threshold have differential selection written as "NaN".')
    parser.add_argument('--no-scale-pseudocounts', dest='no_scale_pseudocounts', action='store_true', help="Use this option if you do NOT want to scale the pseudocounts for each library by relative depth, but instead want to use unscaled value specified by '--pseudocount'.")
    parser.set_defaults(no_scale_pseudocounts=False)
    parser.add_argument('--chartype', choices=['codon_to_aa', 'DNA', 'codon', 'aa'], default='codon_to_aa', help='Characters: "codon_to_aa" = counts for codons and selection for amino acids; "DNA" = counts and selection for DNA; "codon" = counts and selection for codons; "aa" = counts and selection for amino acids (possibly including stop codons, see "--includestop").')
    parser.set_defaults(includestop=False)
    parser.add_argument('--includestop', help='Include stop codons as a possible amino acid if using "--chartype" of "codon_to_aa" or "aa".', dest='includestop', action='store_true')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def InferDiffPrefsParser():
    """Returns *argparse.ArgumentParser* for ``dms_inferdiffprefs`` script."""
    parser = ArgumentParserNoArgHelp(description='Infer differential preferences for amino acids, nucleotides, or codons. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('n_start', help='Existing file with counts from deep sequencing of starting library. For nucleotides, header line should be "# POSITION WT A C G T" and then subsequent lines give wildtype and counts for each site (i.e. "1 G 1013 23 19 47"); for codons use the 64 codons instead (i.e. "AAA AAC AAG ...").', type=ExistingFile)
    parser.add_argument('n_s1', help='Existing file with counts from deep sequencing library subject to selection 1. Typically this is the control selection; differential preferences are for changes in selection 2 relative to selection 1. File has same format as "n_start".', type=ExistingFile)
    parser.add_argument('n_s2', help='Existing file with counts from deep sequencing library subject to selection 2. Typically this is the alternate selection; differential preferences are for changes in selection 2 relative to selection 1. File has same format as "n_start".', type=ExistingFile)
    parser.add_argument('outfile', help='Created output file with differential preferences; overwritten if it already exists.')
    parser.add_argument('--chartype', choices=['codon_to_aa', 'DNA', 'codon', 'aa'], default='codon_to_aa', help='Characters for which differential preferences are inferred: "codon_to_aa" = counts for codons and inferred differential preferences for amino acids; "DNA" = counts and inferred differential preference for DNA; "codon" = counts and inferred differential preferences for codons; "aa" = counts and inferred differential preferences for amino acids (possibly including stop codons, see "--excludestop").')
    parser.add_argument('--err', default='none', help='Control used to estimate errors in counts for all libraries ("n_start", "n_s1", "n_s2"). If "none" then the counts are assumed to have no errors; otherwise this should be a counts file with the same format as "n_start" giving the counts for sequencing unmutated gene.', type=ExistingFile)
    parser.add_argument('--excludestop', help='Exclude stop codons as a possible amino acid if using "--chartype" of "codon_to_aa" or "aa". Stop codons are only excluded if this option is specified; otherwise they are included, and are required to be specified if "--chartype" is "aa".', dest='excludestop', action='store_true')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--alpha_start', help='Concentration parameter for Dirichlet prior over starting frequencies.', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--alpha_pis1', help='Concentration parameter for Dirichlet prior over preferences for selection 1.', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--alpha_err', help='Concentration parameter for Dirichlet priors over error rate.', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--alpha_deltapi', help='Concentration parameter for Dirichlet priors over differential preferences. Larger values correspond to stronger expectation of differential preferences of zero.', default=2.0, type=FloatGreaterThanZero)
    parser.add_argument('--logfile', help='Log progress to this file; overwritten if it already exists.', default='Base name of "outfile" with extension ".log"')
    parser.add_argument('--ncpus', default=1, help='Number of CPUs to use; set to -1 to use all available CPUs.', type=int)
    parser.add_argument('--sites', default=None, nargs='+', help='Only perform the inference for the specified sites, which should be a space separated list such as "--sites 1 2 10". All of these sites must have data in the counts files.')
    parser.add_argument('--ratio_estimation', default=None, metavar='MINCOUNTS', type=FloatGreaterThanZero, help='Rather than use MCMC to estimate the differential preferences using a statistical model, we simply compute them by calculating enrichment ratios relative to wildtype for both selections, then normalizing these ratios to sum to one at each site, and finally taking the difference. The single argument for this option is MINCOUNTS, which is a number > 0. If the counts for a character are less than MINCOUNTS, either due to low counts or error correction, then the counts for that character are changed to MINCOUNTS to avoid estimating ratios of zero, less than zero, or infinity.')
    parser.add_argument('--seed', default=1, help='Random number seed.', type=int)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser



if __name__ == '__main__':
    import doctest
    doctest.testmod()
