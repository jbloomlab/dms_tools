"""Setup script for ``dms_tools``.

Written by Jesse Bloom.
"""

import sys
import os
import re
try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.core import Extension


if (sys.version_info[0], sys.version_info[1]) != (2, 7):
    raise RuntimeError('dms_tools is currently only compatible with Python 2.7.\nYou are using Python %d.%d' % (sys.version_info[0], sys.version_info[1]))

# get metadata, which is specified in another file
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


with open('README.rst') as f:
    readme = f.read()

# main setup command
setup(
    name = 'dms_tools', 
    version = metadata['version'],
    author = metadata['author'],
    author_email = metadata['author_email'],
    url = metadata['url'],
    download_url = 'https://github.com/jbloomlab/dms_tools/tarball/%s' % metadata['version'], # assumes appropriate tagged version is on GitHub
    description = 'Deep mutational scanning (DMS) analysis tools.',
    long_description = readme,
    license = 'GPLv3',
    install_requires = [
        'biopython>=1.6',
        'scipy>=0.13',
        'numpy>=1.8',
        'matplotlib>=1.5',
        'pystan>=2.5',
        'cython>=0.21',
        'weblogo>=3.4, <3.6',
        'PyPDF2>=1.23',
        'pandas>=0.18',
        ],
    platforms = 'Linux (and maybe also Mac OS X).',
    packages = ['dms_tools'],
    package_dir = {'dms_tools':'src'},
    ext_modules = [Extension('dms_tools.cutils', ['src/cutils.c'])],
    scripts = [
            'scripts/dms_inferprefs',
            'scripts/dms_inferdiffprefs',
            'scripts/dms_diffselection',
            'scripts/dms_merge',
            'scripts/dms_correlate',
            'scripts/dms_logoplot',
            'scripts/dms_editsites',
            'scripts/dms_barcodedsubamplicons',
            'scripts/dms_subassemble',
            'scripts/dms_matchsubassembledbarcodes',
            'scripts/dms_summarizealignments',
            ],
    package_data = {'dms_tools':['_weblogo_template.eps']}, # template from weblogo version 3.4
)
