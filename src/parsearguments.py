"""
===========================
``parsearguments`` module
===========================

Module that parses command-line arguments for the scripts in ``dms_tools``.

Written by Jesse Bloom.

Defined in this module
----------------------------------

* *ArgumentParserNoArgHelp* : argument parser that prints help when no args

* *InferPrefsParser* : parser for ``dms_inferprefs``


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


def InferPrefsParser():
    """Returns *argparse.ArgumentParser* for ``dms_inferprefs`` script."""
    parser = ArgumentParserNoArgHelp(description='Infer site-specific preferences for amino acids, nucleotides, or codons. This script is part of %s (version %s) written by %s. The project webpage at %s provides detailed documentation.' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('n_pre', help='Existing file with pre-selection counts from deep sequencing.')
    parser.add_argument('n_post', help='Existing file with post-selection counts from deep sequencing.')
    parser.add_argument('outfile', help='Created output file with site-specific preferences; overwritten if it already exists.')
    parser.add_argument('--chartype', default='codon_to_aa', help='Characters for which preferences are inferred: "codon_to_aa" = counts for codons and inferred preferences for amino acids; "DNA" = counts and inferred preference for DNA; "codon" = counts and inferred preferences for codons.', choices=['codon_to_aa', 'DNA', 'codon'])
    parser.add_argument('--errmodel', default='none', help='Do we use controls to estimate errors in counts from "n_pre" and "n_post"? Options: "none" = counts assumed error-free; "same n_err" = same error rate for "n_pre" and "n_post" estimated from existing file "n_err" with error counts from deep sequencing; "different n_errpre n_errpost" = different error rates for "n_pre" and "n_post" estimated from existing files "n_errpre" and "n_errpost" with error counts from deep sequencing.')
    parser.add_argument('--pi_alpha', help='Concentration parameter for Dirichlet prior over preferences (pi).', default=1.0, type=float)
    parser.add_argument('--mu_alpha', help='Concentration parameter for Dirichlet prior over mutagenesis rate (mu).', default=1.0, type=float)
    parser.add_argument('--logfile', help='Log progress to this file; overwritten if it already exists.', default='Base name of "outfile" with extension ".log"')
    parser.add_argument('--err_alpha', help='Concentration parameter for Dirichlet priors over error rates (epsilon, rho).', default=1.0, type=float)
    parser.add_argument('--ncpus', default=1, help='Number of CPUs to use.', type=int)
    parser.add_argument('--seed', default=1, help='Random number seed.', type=int)
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser
