"""Setup script for ``dms_tools``.

Written by Jesse Bloom.
"""


import sys
import os
from distutils.core import setup
from distutils.core import Extension
from distutils.core import Command


# create a class to handle the 'test' command
class TestCommand(Command):
    """Run all of the tests for this package.

    This is an automatic test run class to make distutils perform the
    package testing. To run these tests, type:

    python setup.py build
    python setup.py test
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        """Run the test script tests/run_tests.py"""
        currdir = os.getcwd()
        os.chdir('tests')
        sys.path.insert(0, '')
        import run_tests
        run_tests.main()
        os.chdir(currdir)

# main setup command
setup(
    name = 'dms_tools', 
    version = '1.0', 
    author = 'Jesse D. Bloom', 
    author_email = 'jbloom@fredhutch.org', 
    url = 'https://github.com/jbloom/dms_tools', 
    description = 'Deep mutational scanning (DMS) analysis tools.',
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: Free for non-commercial use',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    platforms = 'Tested on Linux.',
    packages = ['dms_tools'],
    package_dir = {'dms_tools':'src'},
    cmdclass = {'test':TestCommand},
    scripts = [
            'scripts/dms_inferpreferences',
            ],
)
