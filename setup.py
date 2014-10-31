"""Setup script for ``dms_tools``.

Written by Jesse Bloom.
"""


import sys
import os
import re
from distutils.core import setup
from distutils.core import Extension
from distutils.core import Command

# get version
versionfile = 'src/_version.py'
versionline = open(versionfile).read()
versionstring = re.search("^__version__ = ['\"]([^'\"]+)['\"]", versionline)
if not versionstring:
    raise RuntimeError("Unable to parse version number from %s" % versionfile)
else:
    versionstring = versionstring.group(1)

# main setup command
setup(
    name = 'dms_tools', 
    version = versionstring, 
    author = 'Jesse D. Bloom', 
    author_email = 'jbloom@fredhutch.org', 
    url = 'https://github.com/jbloom/dms_tools', 
    description = 'Deep mutational scanning (DMS) analysis tools.',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: GPLv3',
        'Natural Language :: English',
        'Operating System :: Linux (and maybe also Mac OS X)',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    platforms = 'Tested on Linux.',
    packages = ['dms_tools'],
    package_dir = {'dms_tools':'src'},
    scripts = [
            'scripts/dms_inferpreferences',
            ],
)
