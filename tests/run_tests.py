"""Runs tests on ``dms_tools``.

Tests every file in directory with name ``test_*.py``.

Also tests all modules in package with ``docutils``.

Written by Jesse Bloom."""


import os
import sys
import glob
import doctest
import unittest
import doctest
import pkgutil
import dms_tools


def main():
    """Run the tests."""

    sys.stderr.write('Running tests...\n')
    failurestrings = []

    # test all modules with doctest
    for (importer, modname, ispkg) in pkgutil.iter_modules(dms_tools.__path__):
        if (not ispkg) and modname[0] != '_':
            sys.stderr.write('\nTesting %s with doctest... ' % modname)
            module = __import__('dms_tools.%s' % modname, None, None, modname.split('.'))
            suite = doctest.DocTestSuite(module)
            del module
            result = unittest.TestResult()
            suite.run(result)
            if result.wasSuccessful():
                sys.stderr.write('all %d tests were successful.\n' % result.testsRun)
            else:
                sys.stderr.write('test FAILED!\n')
                for (testcase, failstring) in result.failures:
                    failurestrings.append(failstring)

    # all tests in files with names test_*.py
    for test in glob.glob('test_*.py'):
        test = os.path.splitext(test)[0]
        sys.stderr.write('\nRunning tests in %s...\n' % test)
        suite = unittest.TestLoader().loadTestsFromName(test)
        result = unittest.TestResult()
        suite.run(result)
        if result.wasSuccessful():
            sys.stderr.write('\nAll tests in %s were successful.\n' % test)
        else:
            sys.stderr.write('\nTest in %s FAILED!\n' % test)
            for (testcase, failstring) in result.failures + result.errors:
                failurestrings.append(failstring)

    # print summary of failures
    if not failurestrings:
        sys.stderr.write('\nTesting complete. All tests were passed successfully.\n')
    else:
        sys.stderr.write('\nTesting complete. Failed on the following tests:\n')
        for failstring in failurestrings:
            sys.stderr.write('\n*****************\n%s\n\n****************\n' % failstring)


if __name__ == '__main__':
    main() # run the program
