"""
===========================
``parsearguments`` module
===========================

Module that parses command-line arguments for the scripts in ``dms_tools``.

Written by Jesse Bloom.

Defined in this module
----------------------------------

* *ArgumentParserNoArgHelp* : argument parser that prints help when no args

* *FloatGreaterThanZero* : parses whether string is float greater than zero.

* *ExistingFile* : parses whether string gives an existing file.

* *PDF_File* : parses whether string gives a PDF file name.

* *InferPrefsParser* : parser for ``dms_inferprefs``

* *InferDiffPrefsParser* : parser for ``dms_inferdiffprefs``

* *MergeParser* : parser for ``dms_merge``

* *EditSitesParser* : parsers for ``dms_editsites``

* *CorrelateParser* : parser for ``dms_correlate``

* *LogoPlotParser* : parser for ``dms_logoplot``


Detailed documentation
----------------------------------

"""


import sys
import os
import argparse
import dms_tools


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""
    def error(self, message):
        """Prints error message, then help."""
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def FloatGreaterThanZero(x):
    """If *x* is string for float > 0, returns it, otherwise an error.

    Designed based on this: http://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin

    >>> print "%.3f" % FloatGreaterThanZero('0.1')
    0.100

    >>> FloatGreaterThanZero('0.0')
    Traceback (most recent call last):
        ...
    ArgumentTypeError: 0.0 not a float greater than zero

    >>> FloatGreaterThanZero('hi')
    Traceback (most recent call last):
        ...
    ValueError: could not convert string to float: hi
    """
    x = float(x)
    if x > 0:
        return x
    else:
        raise argparse.ArgumentTypeError("%r not a float greater than zero" % x)


def ExistingFile(fname):
    """If *fname* is name of an existing file return it, otherwise an error.
    
    It is also acceptable for *fname* to be the string "none"."""
    if os.path.isfile(fname) or fname.lower() == 'none':
        return fname
    else:
        raise argparse.ArgumentTypeError("%s does not specify a valid file name" % fname)


def PDF_File(fname):
    """If *fname* ends in PDF extension, returns it. Otherwise an error."""
    if os.path.splitext(fname)[1].lower() == '.pdf':
        return fname
    else:
        raise argparse.ArgumentTypeError('Not a valid PDF file name; should end in ".pdf":\n%s' % fname)


def LogoPlotParser():
    """Returns *argparse.ArgumentParser* for ``dms_logoplot`` script."""
    parser = ArgumentParserNoArgHelp(description='Make a logo plot visually displaying the preferences or differential preferences. Utilizes weblogo (https://code.google.com/p/weblogo/). This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', help='Existing file giving the preferences or differential preferences. The program auto-detects which type of data is present. The preferences can be for DNA or amino acids (stop codons allowed, indicated by "*").', type=ExistingFile)
    parser.add_argument('logoplot', help='Name of created file containing the logo plot, must end in the extension ".pdf". Overwritten if it already exists.', type=PDF_File)
    parser.add_argument('--nperline', type=int, default=60, help='Put this many sites per line of the logo plot.')
    parser.add_argument('--numberevery', type=int, default=10, help='Number x-axis ticks for sites at this interval.')
    parser.add_argument('--diffprefheight', type=float, default=1.0, help='This option is only meaningful if "infile" gives differential preferences. In that case, it gives the height of logo stacks (extends from minus this to plus this). Cannot be smaller than the maximum total differential preference range.')
    parser.add_argument('--excludestop', dest='excludestop', action='store_true', help='If we are using amino acids, do we remove stop codons (denoted by "*")? We only remove stop codons if this argument is specified. If this option is used, then data for stop codons is removed by re-normalizing preferences to sum to one, and differential preferences to sum to zero.')
    parser.add_argument('--overlay1', default=None, nargs=3, metavar=('FILE', 'SHORTNAME', 'LONGNAME'), help='Specify an overlay bar above each line of the logo plot to illustrate a per-residue property such as relative solvent accessibility or secondary structure. Requires three arguments: FILE SHORTNAME LONGNAME. FILE is the name of an existing file. Except for comment lines beginning with "#", each line should have two whitespace delimited columns (additional columns are allowed but ignored). The first column gives the site number (matching that in "infile") and the second column giving the property for this site; properties must either all be non-whitespace strings giving a discrete category (such as secondary structure), or all be numbers (such as relative solvent accessibility). All listed sites must be in "infile", but not all sites in "infile" must be in FILE -- missing sites are assumed to lack a known value for the property and are shown in white. SHORTNAME is a short (3-5 character) name of the property, such as "RSA" for "relative solvent accessibility." LONGNAME is a longer name (such as "relative solvent accessibiity"), or the same as SHORTNAME if you do not have a separate long name.')
    parser.add_argument('--overlay2', default=None, nargs=3, metavar=('FILE', 'SHORTNAME', 'LONGNAME'), help='Specify a second overlay bar. Arguments have the same meaning as for "overlay1".')
    parser.set_defaults(excludestop=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def CorrelateParser():
    """Returns *argparse.ArgumentParser* for ``dms_correlate`` script."""
    parser = ArgumentParserNoArgHelp(description='Determine and plot the Pearson correlation between pairs of preferences or pairs of differential preferences. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file1', help='Existing file giving first set of preferences or differential preferences.', type=ExistingFile)
    parser.add_argument('file2', help='Existing file giving second set of preferences or differential preferences. Must give the same type of data (preferences or differential preferences) for the same set of sites and same type of character as "file1".', type=ExistingFile)
    parser.add_argument('outfileprefix', help='Prefix for name of created output files (these files are overwritten if they already exist). The correlations are written to a file with this prefix and the suffix ".txt". Unless you use the option "--noplot", a scatter plot is written to a file with this prefix and the suffix ".pdf". In the correlations text file, the first line gives the Pearson correlation coefficient in the format: "R = 0.5312". The second line gives the P-value in the format: "P = 0.0000131". The third line gives the number of points in the format: "N = 4200".')
    parser.add_argument('--name1', default='data 1', help='Name used for the data in "file1" in the correlation plot. The string specified here uses LaTex formatting; names with spaces should be enclosed in quotes.')
    parser.add_argument('--name2', default='data 2', help='Name used for the data in "file2" in the correlation plot. The string specified here uses LaTex formatting; names with spaces should be enclosed in quotes.')
    parser.add_argument('--excludestop', dest='excludestop', action='store_true', help='If we are using amino acids, do we remove stop codons (denoted by "*")? We only remove stop codons if this argument is specified. If this option is used, then data for stop codons is removed by re-normalizing preferences to sum to one, and differential preferences to sum to zero.')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--noplot', dest='noplot', action='store_true', help='Normally this script creates a PDF scatter plot. If this option is specified, then no such plot will be created.')
    parser.set_defaults(noplot=False)
    parser.add_argument('--alpha', default=0.1, help='The transparency (alpha value) for the points on the scatter plot. A value of 1.0 correspond to no transparency; values close to zero give high transparency. Transparency (alpha < 1) might be helpful if you have many points on top of each other.', type=float)
    parser.add_argument('--plot_title', default='None', help='Title put at the top of the correlation plot. The string specified here uses LaTex formatting; names with spaces should be enclosed in spaces. A value of "None" means no title.')
    parser.add_argument('--corr_on_plot', dest='corr_on_plot', action='store_true', help='If this option is used, then the correlation coefficient will be visually displayed on the plot.')
    parser.set_defaults(corr_on_plot=False)
    parser.add_argument('--r2', dest='r2', action='store_true', help='If this option is used, the correlation coefficient displayed on the plot when using "--corr_on_plot" will show R-squared rather than the R value.')
    parser.set_defaults(r2=False)
    parser.add_argument('--rms_dpi', dest='rms_dpi', action='store_true', help='If "file1" and "file2" specify differential preferences, this argument specifies that we compute the correlation between the root-mean-square (RMS) differential preference at each site rather than between the differential preferences themselves.')
    parser.set_defaults(rms_dpi=False)
    parser.add_argument('--pref_entropy', dest='pref_entropy', action='store_true', help='If "file1" and "file2" specify preferences, this argument specifies that we compute the correlation between the site entropy of the preferences (logarithm base 2, so bits) at each site rather than between the preferences themselves.')
    parser.set_defaults(pref_entropy=False)
    parser.set_defaults(enrichment=False)
    parser.add_argument('--enrichment', dest='enrichment', action='store_true', help='If this option is set, we plot the enrichment ratio for all mutations on a log scale rather than plotting the preferences. The computed correlations are also then for the log-transformed enrichment ratios. The enrichment ratio for the wildtype identity is always one, and so is not included.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def EditSitesParser():
    parser = ArgumentParserNoArgHelp(description='Edits sites in a data file. Typically you would use this program if you wanted to renumber sites or remove certain sites. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', type=ExistingFile, help='Existing data file. This could be a deep mutational scanning counts file, a preferences file, a differential preferences file, or any other file with the following format: blank lines or lines beginning with "#" (comment lines) are ignored; every other line must begin with an entry giving a unique site number (such as "1" or "2A"). The line may then have an arbitrary number of other entries separated from the site number by whitespace.')
    parser.add_argument('outfile', help='The created output file in which the site editing has been performed on "infile". If this output file already exists, it is overwritten.')
    parser.add_argument('edit_method', choices=['renumber', 'remove', 'retain'], help='How to do the editing: renumber sites, remove specified sites, or retain only specified sites.')
    parser.add_argument('edit_file', type=ExistingFile, help='Existing file specifying how edits are made. If "edit_method" is "renumber", then all non-comment lines (those not beginning with "#") must have two space delimited entries specifying the existing site in "infile" and the new site number with which it is replaced; all sites must be specified. If "edit_method" is "remove", then each line should have as its first entry a site, and all of the listed sites are removed. If "edit_method" is "retain", then each line should have as its first entry a site, and only the listed sites are retained.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def MergeParser():
    """Returns *argparse.ArgumentParser* for ``dms_merge`` script."""
    parser = ArgumentParserNoArgHelp(description='Merge preferences or differential preferences by averaging or adding / subtracting the values in multiple files. All files must specify the same character type: can be nucleotide, codon, or amino acid (see "--excludestop" if using amino acids). This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outfile', help='Created output file with merged preferences / differential preferences; removed if it already exists. If the "infiles" do not all have the same wildtype residue at a site, then the wildtype is indicated as "?" in "outfile".')
    parser.add_argument('merge_method', help='How to merge: If "average" then all "infiles" must specify either preferences or differential preferences; these are then averaged in "outfile". If "sum" then "infiles" can either all be differential preferences, or can be a combination of preferences and differential preferences that (along with any additional files specified by "--minus") sum to a total preference or a total differential preference of zero at each site.', choices=['average', 'sum'])
    parser.add_argument('infiles', nargs='+', help='Files to average or sum. Must all have the same sites and character type, but do not need to have the same wildtype residue at each site.', type=ExistingFile)
    parser.add_argument('--excludestop', dest='excludestop', action='store_true', help='If we are using amino acids, do we remove stop codons (denoted by "*")? We only remove stop codons if this argument is specified. If this option is used, then any files with stop codons have these codons removed (re-normalizing preferences to sum to one, and differential preferences to sum to zero) before the merge.')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--minus', nargs='+', help='Files to subtract when summing. Can only be used if "merge_method" is "sum".', type=ExistingFile)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def InferPrefsParser():
    """Returns *argparse.ArgumentParser* for ``dms_inferprefs`` script."""
    parser = ArgumentParserNoArgHelp(description='Infer site-specific preferences for amino acids, nucleotides, or codons. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('n_pre', help='Existing file with pre-selection counts from deep sequencing. For nucleotides, header line should be "# POSITION WT A C G T" and then subsequent lines give wildtype and counts for each site (i.e. "1 G 1013 23 19 47"); for codons just use the 64 codons instead (i.e. "AAA AAC AAG ...").', type=ExistingFile)
    parser.add_argument('n_post', help='Existing file with post-selection counts from deep sequencing. Same format as "n_pre".', type=ExistingFile)
    parser.add_argument('outfile', help='Created output file with site-specific preferences; overwritten if it already exists.')
    parser.add_argument('--chartype', default='codon_to_aa', help='Characters for which preferences are inferred: "codon_to_aa" = counts for codons and inferred preferences for amino acids; "DNA" = counts and inferred preference for DNA; "codon" = counts and inferred preferences for codons; "aa" = counts and inferred preferences for amino acids (possibly including stop codons, see "--excludestop").', choices=['codon_to_aa', 'DNA', 'codon', 'aa'],)
    parser.add_argument('--errpre', default='none', help='Control used to estimate errors in counts from "n_pre". If "none" then the counts in "n_pre" are assumed to have no errors; otherwise this should be a counts file with the same format as "n_pre" giving the counts for sequencing unmutated gene. Currently if this option is not "none" then --errpost also cannot be "none".', type=ExistingFile)
    parser.add_argument('--errpost', default='none', help='Control used to estimate errors in counts from "n_post". If "none" then the counts in "n_post" are assumed to have no errors; otherwise this should be a counts file with the same format as "n_pre" giving the counts for sequencing unmutated gene. If you have the same error control for "n_pre" and "n_post", then set that file name for both this option and --errpre. Currently if this option is not "none" then --errpre also cannot be "none".', type=ExistingFile)
    parser.add_argument('--excludestop', help='Exclude stop codons as a possible amino acid if using "--chartype" of "codon_to_aa" or "aa". Stop codons are only excluded if this option is specified; otherwise they are included, and are required to be specified if "--chartype" is "aa".', dest='excludestop', action='store_true')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--pi_alpha', help='Concentration parameter for Dirichlet prior over preferences (pi).', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--mu_alpha', help='Concentration parameter for Dirichlet prior over mutagenesis rate (mu).', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--err_alpha', help='Concentration parameter for Dirichlet priors over error rates (epsilon, rho).', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--logfile', help='Log progress to this file; overwritten if it already exists.', default='Base name of "outfile" with extension ".log"')
    parser.add_argument('--ncpus', default=1, help='Number of CPUs to use; set to -1 to use all available CPUs.', type=int)
    parser.add_argument('--seed', default=1, help='Random number seed.', type=int)
    parser.add_argument('--sites', default=None, nargs='+', help='Only perform the inference for the specified sites, which should be a space separated list such as "--sites 1 2 10". All of these sites must have data in the counts files.')
    parser.add_argument('--ratio_estimation', default=None, metavar='PSEUDOCOUNTS', type=FloatGreaterThanZero, help='Rather than use MCMC to estimate the preferences using a statistical model, we simply compute by calculating enrichment ratios relative to wildtype post- and pre-selection, and then normalizing these ratios to sum to one at each site. The single argument for this option is PSEUDOCOUNTS, which is a number > 0 that specifies the pseudocount added to each count to avoiding estimating ratios of zero or infinity.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser


def InferDiffPrefsParser():
    """Returns *argparse.ArgumentParser* for ``dms_inferdiffprefs`` script."""
    parser = ArgumentParserNoArgHelp(description='Infer differential preferences for amino acids, nucleotides, or codons. This script is part of %s (version %s) written by %s. Detailed documentation is at %s' % (dms_tools.__name__, dms_tools.__version__, dms_tools.__author__, dms_tools.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('n_start', help='Existing file with counts from deep sequencing of starting library. For nucleotides, header line should be "# POSITION WT A C G T" and then subsequent lines give wildtype and counts for each site (i.e. "1 G 1013 23 19 47"); for codons use the 64 codons instead (i.e. "AAA AAC AAG ...").', type=ExistingFile)
    parser.add_argument('n_s1', help='Existing file with counts from deep sequencing library subject to selection 1. Typically this is the control selection; differential preferences are for changes in selection 2 relative to selection 1. File has same format as "n_start".', type=ExistingFile)
    parser.add_argument('n_s2', help='Existing file with counts from deep sequencing library subject to selection 2. Typically this is the alternate selection; differential preferences are for changes in selection 2 relative to selection 1. File has same format as "n_start".', type=ExistingFile)
    parser.add_argument('outfile', help='Created output file with differential preferences; overwritten if it already exists.')
    parser.add_argument('--chartype', choices=['codon_to_aa', 'DNA', 'codon', 'aa'], default='codon_to_aa', help='Characters for which differential preferences are inferred: "codon_to_aa" = counts for codons and inferred differential preferences for amino acids; "DNA" = counts and inferred differential preference for DNA; "codon" = counts and inferred differential preferences for codons; "aa" = counts and inferred differential preferences for amino acids (possibly including stop codons, see "--excludestop").')
    parser.add_argument('--err', default='none', help='Control used to estimate errors in counts for all libraries ("n_start", "n_s1", "n_s2"). If "none" then the counts are assumed to have no errors; otherwise this should be a counts file with the same format as "n_start" giving the counts for sequencing unmutated gene.', type=ExistingFile)
    parser.add_argument('--excludestop', help='Exclude stop codons as a possible amino acid if using "--chartype" of "codon_to_aa" or "aa". Stop codons are only excluded if this option is specified; otherwise they are included, and are required to be specified if "--chartype" is "aa".', dest='excludestop', action='store_true')
    parser.set_defaults(excludestop=False)
    parser.add_argument('--alpha_start', help='Concentration parameter for Dirichlet prior over starting frequencies.', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--alpha_pis1', help='Concentration parameter for Dirichlet prior over preferences for selection 1.', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--alpha_err', help='Concentration parameter for Dirichlet priors over error rate.', default=1.0, type=FloatGreaterThanZero)
    parser.add_argument('--alpha_deltapi', help='Concentration parameter for Dirichlet priors over differential preferences. Larger values correspond to stronger expectation of differential preferences of zero.', default=2.0, type=FloatGreaterThanZero)
    parser.add_argument('--logfile', help='Log progress to this file; overwritten if it already exists.', default='Base name of "outfile" with extension ".log"')
    parser.add_argument('--ncpus', default=1, help='Number of CPUs to use; set to -1 to use all available CPUs.', type=int)
    parser.add_argument('--sites', default=None, nargs='+', help='Only perform the inference for the specified sites, which should be a space separated list such as "--sites 1 2 10". All of these sites must have data in the counts files.')
    parser.add_argument('--seed', default=1, help='Random number seed.', type=int)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=dms_tools.__version__))
    return parser



if __name__ == '__main__':
    import doctest
    doctest.testmod()
