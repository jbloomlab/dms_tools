"""Runs tests on every file in directory with name ``test_*.py``.

Written by Jesse Bloom."""


import os
import sys
import glob
import doctest
import unittest


def main():
    """Run the tests."""

    sys.stderr.write('Running tests...\n')
    failurestrings = []

    # test all
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
        syst.stderr.write('\nTesting complete. Failed on the following tests:\n')
        for failstring in failurestrings:
            sys.stderr.write('\n*****************\n%s\n\n****************\n' % failstring)


if __name__ == '__main__':
    main() # run the program
