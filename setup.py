"""Setup script for ``dms_tools``.

Written by Jesse Bloom.
"""


import sys
import os
import re
from distutils.core import setup
from distutils.core import Extension
from distutils.core import Command

# get metadata
metadata = {}
lines = open('src/_metadata.py').readlines()
for dataname in ['version', 'author', 'author_email', 'url']:
    for line in lines:
        entries = line.split('=')
        assert len(entries) == 2, "Failed to parse metadata:\n%s" % line
        if entries[0].strip() == '__%s__' % dataname:
            if dataname in metadata:
                raise ValueError("Duplicate metadata for %s" % dataname)
            else:
                metadata[dataname] = entries[1].strip()[1 : -1]
    assert dataname in metadata, "Failed to find metadata for %s" % dataname

# main setup command
setup(
    name = 'dms_tools', 
    version = metadata['version'],
    author = metadata['author'],
    author_email = metadata['author_email'],
    url = metadata['url'],
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
            'scripts/dms_inferprefs',
            'scripts/dms_inferdiffprefs',
            'scripts/dms_merge',
            'scripts/dms_correlate',
            'scripts/dms_logoplot',
            ],
    package_data = {'dms_tools':['_weblogo_template.eps']}, # template from weblogo version 3.4
)
