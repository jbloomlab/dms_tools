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

* *ExistingFile* : parses whether string gives an existing file.

* *InferPrefsParser* : parser for ``dms_inferprefs``

* *InferDiffPrefsParser* : parser for ``dms_inferdiffprefs``

* *MergeParser* : parser for ``dms_merge``


Detailed documentation
----------------------------------

"""


import sys
import os
import argparse
import dms_tools


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""
    def error(self, message):
        """Prints error message, then help."""
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


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


def ExistingFile(fname):
    """If *fname* is name of an existing file return it, otherwise an error.
    
    It is also acceptable for *fname* to be the string "none"."""
    if os.path.isfile(fname) or fname.lower() == 'none':
        return fname
    else:
        raise argparse.ArgumentTypeError("%s does not specify a valid file name" % fname)


def MergeParser():
    """Returns *argparse.ArgumentParser* for ``dms_merge`` script."""
    parser = ArgumentParserNoArgHelp(description='Merge preferences or differential preferences by averaging or adding / subtracting the values in multiple files. All files must specify the same character type: can be nucleotide, codon, or amino acid (see "--excludestop" if using amino acids). This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outfile', help='Created output file with merged preferences / differential preferences; removed if it already exists. If the "infiles" do not all have the same wildtype residue at a site, then the wildtype is indicated as "?" in "outfile".')
    parser.add_argument('merge_method', help='How to merge: If "average" then all "infiles" must specify either preferences or differential preferences; these are then averaged in "outfile". If "sum" then "infiles" can either all be differential preferences, or can be a combination of preferences and differential preferences that (along with any additional files specified by "--minus") sum to a total preference or a total differential preference of zero at each site.', choices=['average', 'sum'])
    parser.add_argument('infiles', nargs='+', help='Files to average or sum. Must all have the same sites and character type, but do not need to have the same wildtype residue at each site.', type=ExistingFile)
    parser.add_argument('--excludestop', dest='excludestop', action='store_true', help='If we are using amino acids, do we remove stop codons (denoted by "*")? We only remove stop codons if this argument is specified. If this option is used, then any files with stop codons have these codons removed (re-normalizing preferences to sum to one, and differential preferences to sum to zero) before the merge.')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--minus', nargs='+', help='Files to subtract when summing. Can only be used if "merge_method" is "sum".', type=ExistingFile)
    return parser


def InferPrefsParser():
    """Returns *argparse.ArgumentParser* for ``dms_inferprefs`` script."""
    parser = ArgumentParserNoArgHelp(description='Infer site-specific preferences for amino acids, nucleotides, or codons. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('n_pre', help='Existing file with pre-selection counts from deep sequencing. For nucleotides, header line should be "# POSITION WT A C G T" and then subsequent lines give wildtype and counts for each site (i.e. "1 G 1013 23 19 47"); for codons just use the 64 codons instead (i.e. "AAA AAC AAG ...").', type=ExistingFile)
    parser.add_argument('n_post', help='Existing file with post-selection counts from deep sequencing. Same format as "n_pre".', type=ExistingFile)
    parser.add_argument('outfile', help='Created output file with site-specific preferences; overwritten if it already exists.')
    parser.add_argument('--chartype', default='codon_to_aa', help='Characters for which preferences are inferred: "codon_to_aa" = counts for codons and inferred preferences for amino acids; "DNA" = counts and inferred preference for DNA; "codon" = counts and inferred preferences for codons.')
    parser.add_argument('--errpre', default='none', help='Control used to estimate errors in counts from "n_pre". If "none" then the counts in "n_pre" are assumed to have no errors; otherwise this should be a counts file with the same format as "n_pre" giving the counts for sequencing unmutated gene. Currently if this option is not "none" then --errpost also cannot be "none".', type=ExistingFile)
    parser.add_argument('--errpost', default='none', help='Control used to estimate errors in counts from "n_post". If "none" then the counts in "n_post" are assumed to have no errors; otherwise this should be a counts file with the same format as "n_pre" giving the counts for sequencing unmutated gene. If you have the same error control for "n_pre" and "n_post", then set that file name for both this option and --errpre. Currently if this option is not "none" then --errpre also cannot be "none".', type=ExistingFile)
    parser.add_argument('--excludestop', help='Exclude stop codons as a possible amino acid if using "--chartype codon_to_aa". Stop codons are only excluded if this option is specified; otherwise they are included.', dest='excludestop', action='store_true')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--pi_alpha', help='Concentration parameter for Dirichlet prior over preferences (pi).', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--mu_alpha', help='Concentration parameter for Dirichlet prior over mutagenesis rate (mu).', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--err_alpha', help='Concentration parameter for Dirichlet priors over error rates (epsilon, rho).', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--logfile', help='Log progress to this file; overwritten if it already exists.', default='Base name of "outfile" with extension ".log"')
    parser.add_argument('--ncpus', default=1, help='Number of CPUs to use; set to -1 to use all available CPUs.', type=int)
    parser.add_argument('--seed', default=1, help='Random number seed.', type=int)
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def InferDiffPrefsParser():
    """Returns *argparse.ArgumentParser* for ``dms_inferdiffprefs`` script."""
    parser = ArgumentParserNoArgHelp(description='Infer differential preferences for amino acids, nucleotides, or codons. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('n_start', help='Existing file with counts from deep sequencing of starting library. For nucleotides, header line should be "# POSITION WT A C G T" and then subsequent lines give wildtype and counts for each site (i.e. "1 G 1013 23 19 47"); for codons use the 64 codons instead (i.e. "AAA AAC AAG ...").', type=ExistingFile)
    parser.add_argument('n_s1', help='Existing file with counts from deep sequencing library subject to selection 1. Typically this is the control selection; differential preferences are for changes in selection 2 relative to selection 1. File has same format as "n_start".', type=ExistingFile)
    parser.add_argument('n_s2', help='Existing file with counts from deep sequencing library subject to selection 2. Typically this is the alternate selection; differential preferences are for changes in selection 2 relative to selection 1. File has same format as "n_start".', type=ExistingFile)
    parser.add_argument('outfile', help='Created output file with differential preferences; overwritten if it already exists.')
    parser.add_argument('--chartype', default='codon_to_aa', help='Characters for which preferences are inferred: "codon_to_aa" = counts for codons and inferred preferences for amino acids; "DNA" = counts and inferred preference for DNA; "codon" = counts and inferred preferences for codons.')
    parser.add_argument('--err', default='none', help='Control used to estimate errors in counts for all libraries ("n_start", "n_s1", "n_s2"). If "none" then the counts are assumed to have no errors; otherwise this should be a counts file with the same format as "n_start" giving the counts for sequencing unmutated gene.', type=ExistingFile)
    parser.add_argument('--excludestop', help='Exclude stop codons as a possible amino acid if using "--chartype codon_to_aa". Stop codons are only excluded if this option is specified; otherwise they are included.', dest='excludestop', action='store_true')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--alpha_start', help='Concentration parameter for Dirichlet prior over starting frequencies.', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--alpha_pis1', help='Concentration parameter for Dirichlet prior over preferences for selection 1.', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--alpha_err', help='Concentration parameter for Dirichlet priors over error rate.', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--alpha_deltapi', help='Concentration parameter for Dirichlet priors over differential preferences. Larger values correspond to stronger expectation of differential preferences of zero.', default=10.0, type=FloatGreaterThanZero)
    parser.add_argument('--logfile', help='Log progress to this file; overwritten if it already exists.', default='Base name of "outfile" with extension ".log"')
    parser.add_argument('--ncpus', default=1, help='Number of CPUs to use; set to -1 to use all available CPUs.', type=int)
    parser.add_argument('--seed', default=1, help='Random number seed.', type=int)
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser



if __name__ == '__main__':
    import doctest
    doctest.testmod()
