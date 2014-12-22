===============================
Tests for `dms_tools`_
===============================

This directory contains tests for the `dms_tools`_ package. If you have modified the code and want to make sure the package is still working, you might run these tests.

Before running the tests, you need to install `dms_tools`_.

Each test is in a file of the format ``test_*.py``. The tests can be run individually with commands like::

    python test_inferprefs.py

More often, you might want to run all of the tests using the command::

    python run_tests.py

Some of the tests will take a long time (perhaps hours) to run.

You must ``cd`` into this ``./tests/`` subdirectory before running the tests.

The file ``test_configuration.txt`` in this directory specifies configuration options used by some of the tests. It should have three lines:

    1) The first line should begin with the word *seed* followed by a comma-delimited list of integer random number seeds used for the tests. Each test is run separately for each seed.

    2) The second line should begin with the word *charactertypes* followed by a comma-delimited list of one or more of the words *DNA*, *amino acid*, and *codon*. The tests are run for each of these character types. Tests for *codon* characters are particularly slow, so you might sometimes want to leave this out.

    3) The third line should begin with the word *n_jobs* and then specify the number of CPUs (at least one) to be used. Set to -1 to use all available CPUs.

For instance, here is an example of a valid ``test_configuration.txt`` file::

    seeds 0
    charactertypes DNA
    n_jobs 3

.. _`dms_tools`: https://github.com/jbloom/dms_tools
