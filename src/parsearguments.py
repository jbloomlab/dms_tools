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

* *InferPrefsParser* : parser for ``dms_inferprefs``

* *InferDiffPrefsParser* : parser for ``dms_inferdiffprefs``


Detailed documentation
----------------------------------

"""


import sys
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


def InferPrefsParser():
    """Returns *argparse.ArgumentParser* for ``dms_inferprefs`` script."""
    parser = ArgumentParserNoArgHelp(description='Infer site-specific preferences for amino acids, nucleotides, or codons. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('n_pre', help='Existing file with pre-selection counts from deep sequencing. For nucleotides, header line should be "# POSITION WT A C G T" and then subsequent lines give wildtype and counts for each site (i.e. "1 G 1013 23 19 47"); for codons just use the 64 codons instead (i.e. "AAA AAC AAG ...").')
    parser.add_argument('n_post', help='Existing file with post-selection counts from deep sequencing. Same format as "n_pre".')
    parser.add_argument('outfile', help='Created output file with site-specific preferences; overwritten if it already exists.')
    parser.add_argument('--chartype', default='codon_to_aa', help='Characters for which preferences are inferred: "codon_to_aa" = counts for codons and inferred preferences for amino acids; "DNA" = counts and inferred preference for DNA; "codon" = counts and inferred preferences for codons.')
    parser.add_argument('--errpre', default='none', help='Control used to estimate errors in counts from "n_pre". If "none" then the counts in "n_pre" are assumed to have no errors; otherwise this should be a counts file with the same format as "n_pre" giving the counts for sequencing unmutated gene. Currently if this option is not "none" then --errpost also cannot be "none".')
    parser.add_argument('--errpost', default='none', help='Control used to estimate errors in counts from "n_post". If "none" then the counts in "n_post" are assumed to have no errors; otherwise this should be a counts file with the same format as "n_pre" giving the counts for sequencing unmutated gene. If you have the same error control for "n_pre" and "n_post", then set that file name for both this option and --errpre. Currently if this option is not "none" then --errpre also cannot be "none".')
    parser.add_argument('--includestop', help='Are stop codons considered a possible amino acid if using "--chartype codon_to_aa"? Valid values are "True" or "False".', type=bool, default='True')
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
    parser.add_argument('n_start', help='Existing file with counts from deep sequencing of starting library. For nucleotides, header line should be "# POSITION WT A C G T" and then subsequent lines give wildtype and counts for each site (i.e. "1 G 1013 23 19 47"); for codons use the 64 codons instead (i.e. "AAA AAC AAG ...").')
    parser.add_argument('n_s1', help='Existing file with counts from deep sequencing library subject to selection 1. Typically this is the control selection; differential preferences are for changes in selection 2 relative to selection 1. File has same format as "n_start".')
    parser.add_argument('n_s2', help='Existing file with counts from deep sequencing library subject to selection 2. Typically this is the alternate selection; differential preferences are for changes in selection 2 relative to selection 1. File has same format as "n_start".')
    parser.add_argument('outfile', help='Created output file with differential preferences; overwritten if it already exists.')
    parser.add_argument('--chartype', default='codon_to_aa', help='Characters for which preferences are inferred: "codon_to_aa" = counts for codons and inferred preferences for amino acids; "DNA" = counts and inferred preference for DNA; "codon" = counts and inferred preferences for codons.')
    parser.add_argument('--err', default='none', help='Control used to estimate errors in counts for all libraries ("n_start", "n_s1", "n_s2"). If "none" then the counts are assumed to have no errors; otherwise this should be a counts file with the same format as "n_start" giving the counts for sequencing unmutated gene.')
    parser.add_argument('--includestop', help='Are stop codons considered a possible amino acid if using "--chartype codon_to_aa"? Valid values are "True" or "False".', type=bool, default='True')
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
