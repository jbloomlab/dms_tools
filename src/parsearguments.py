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
    parser = argparse.ArgumentParser()
    parser.add_argument('foo', default=False, help='foo help')
    parser.add_argument('bar', default=True)
    return parser
