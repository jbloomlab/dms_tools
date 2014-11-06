"""
===========================
``dms_tools`` package
===========================
This package contains utilities for analyzing deep mutational scanning (DMS) data.

The package webpage is https://github.com/jbloom/dms_tools

The package is typically run via a series of command line scripts. However,
it also provides the following Python API.

Modules in the package
--------------------------
* *mcmc* : Perform MCMC inference

* *parsearguments* : Parses command line arguments to scripts.

Package-level variables
--------------------------
The following variables are defined for this package:

* *aminoacids* : list of the 20 amino acids in one-letter upper-case code.

* *aminoacids_nostop* : another name for *aminoacids*.

* *aminoacids_withstop* : like *aminoacids* but also includes the stop codon (``*``).

* *nts* : a list of the four upper-case DNA nucleotides.

* *codons* : a list of the 64 codons as upper-case three-letter strings.

"""
from _metadata import __version__
from _metadata import __author__
from _metadata import __author_email__
from _metadata import __url__
import Bio.Alphabet.IUPAC

# lists of upper case DNA nucleotides and codons
aminoacids = [aa.upper() for aa in Bio.Alphabet.IUPAC.IUPACProtein.letters]
aminoacids.sort()
aminoacids_nostop = aminoacids
aminoacids_withstop = aminoacids + ['*']
nts = Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters
nts = [nt.upper() for nt in nts]
nts.sort()
codons = []
for nt1 in nts:
    for nt2 in nts:
        for nt3 in nts:
            codons.append('%s%s%s' % (nt1, nt2, nt3))
