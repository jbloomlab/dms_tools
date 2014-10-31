"""
===========================
``parsearguments`` module
===========================

Module that parses command-line arguments for the scripts in ``dms_tools``.

Written by Jesse Bloom.

Functions in this module
----------------------------------

* *InferPreferencesParser* : parser for ``dms_inferpreferences``


Function documentation
----------------------------------

"""


import argparse


def InferPreferencesParser():
    """Returns *argparse.ArgumentParser* for ``dms_inferpreferences`` script."""
    parser = argparse.ArgumentParser(description='Inference of site-specific preferences for amino acids, nucleotides, or codons.')
    parser.add_argument('pre-selection_counts', help='Existing file with pre-selection counts from deep sequencing.')
    parser.add_argument('post-selection_counts', help='Existing file with post-selection counts from deep sequencing')
    parser.add_argument('outfile', help='Created output file with site-specific preferences')
    return parser
