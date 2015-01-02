"""Master script to test differential preference inference on simulated deep mutational scanning data.

Written by Jesse Bloom."""


import os
import re
import glob
import subprocess
import dms_tools
import dms_tools.plot
import dms_tools.simulate



def SimulateCounts(simulated_dir, start_counts, fnameprefix, actual_diff_prefs):
    """Simulate deep mutational scanning counts. 

    Simulations done at several sequencing depths, and with and without
    sequencing errors.
    
    Returns list of all created files with the simulated data. For each simulation,
    there is a 4-tuple *(nstart, ns1, ns1, nerr)*. For simulations with
    no errors, the last entry is 'none'.
    
    The files are placed in *simulated_dir*, which should be an existing directory.

    The files have the prefix *fnameprefix* followed by a suffix giving information
    about the depth and whether the simulations include sequencing errors.
    
    The simulations are done using starting counts proportional to those in *start_counts*.

    *actual_diff_prefs* is created file with actual differential preferences.
    """    

    assert os.path.isdir(simulated_dir), "Can't find %s" % simulated_dir
    assert os.path.isfile(start_counts), "Can't find %s" % start_counts
        
    depths = [2e5, 5e5, 1e6, 2e6, 5e6, 1e7]

    filelist = []

    avgepsilon = 0.0

    fprefix = '%s/%s_depth' % (simulated_dir, fnameprefix)
    depths_outfileprefixes = [(depth, '%s_%.0e' % (fprefix, depth)) for depth in depths]
    dms_tools.simulate.SimulateDiffPrefCounts(depths_outfileprefixes, start_counts, actual_diff_prefs, avgepsilon)
    for (depth, outfileprefix) in depths_outfileprefixes:
        filetup = ['%s_start.txt' % outfileprefix, '%s_s1.txt' % outfileprefix, '%s_s2.txt' % outfileprefix]
        if avgepsilon:
            filetup.append('%s_err.txt' % outfileprefix)
        else:
            filetup.append('none')
        filelist.append(tuple(filetup))
    return filelist



def main():
    """Main body of script."""

    # specify directory and file names
    simulated_dir = './simulated_data/' # simulated data and its subsequent analysis
    if not os.path.isdir(simulated_dir):
        os.mkdir(simulated_dir)
    start_counts = 'start_counts.txt' # plausible starting counts based on Melnikov et al data
    actual_diff_prefs = './simulated_data/actual_diff_prefs.txt'
    one_over_mu = 264 # gene length, inverse of mutation rate

    # simulate data, then analyze
    simulatedfiles = SimulateCounts(simulated_dir, start_counts, 'simulated', actual_diff_prefs)
    corr_plots = [] 
    depths = set([])
    alphas = set([])
    for (nstart, ns1, ns2, nerr) in simulatedfiles:
        for alpha in ['2']:
            outfile = '%s_%s_inferred.txt' % ('_'.join(nstart.split('_')[ : -1]), alpha)
            depth = float(re.search('_depth_([\d\.eE\+]+)_', outfile).group(1))
            depths.add(depth)
            alphas.add(float(alpha))
            subprocess.call(['dms_inferdiffprefs', nstart, ns1, ns2, outfile, '--err', nerr, '--ncpus', '-1', '--alpha_deltapi', alpha])
            corrprefix = '%s_corr' % os.path.splitext(outfile)[0]
            corr_plots.append((depth, float(alpha), '%s.pdf' % corrprefix))
            subprocess.call(['dms_correlate', actual_diff_prefs, outfile, corrprefix, '--name1', 'actual', '--name2', 'inferred', '--corr_on_plot', '--r2'])

    # make merged summary PDFs
    depths = list(depths)
    alphas = list(alphas)
    depths.sort()
    alphas.sort()
    plotwidth = 2
    latexhead = '\n'.join(['\documentclass[12pt,letterpaper]{article}',\
                           '\usepackage{mathptmx}',\
                           '\usepackage{graphicx}',\
                           '\usepackage{setspace}',\
                           '\usepackage{rotating}',\
                           '\pagestyle{empty}',\
                           '\usepackage[paperheight=%.2fin,paperwidth=%.2fin,top=0.05in,bottom=0.05in,right=0.05in,left=0.05in]{geometry}' % (len(alphas) * plotwidth * 1.05 + 0.2, len(depths) * plotwidth * 1.03 + int(len(alphas) > 1) * 0.6),\
                           '\\begin{document}',\
                           '\\begin{figure}',\
                          ])
    corr_plots.sort()
    fname = 'correlations.tex'
    firstcolumn = True
    if len(alphas) == 1:
        firstcolumn = False
    with open(fname, 'w') as f:
        f.write(latexhead)
        for (idepth, depth) in enumerate(depths):
            iwidth = plotwidth
            if firstcolumn:
                iwidth += 0.6
            Nmu = str(int(depth / float(one_over_mu)))
            depth = dms_tools.plot.Base10Formatter(depth, 2, 0, 0)
            if firstcolumn:
                f.write('\n\\begin{minipage}{%.2fin}\n\centerline{\large \hspace{0.7in} $\mathbf{N\overline{\mu} = \\frac{%s}{%d} \\approx %s}$}\n\hspace{0.01in}\n' % (iwidth, depth, one_over_mu, Nmu))
            else:
                f.write('\n\\begin{minipage}{%.2fin}\n\centerline{\large \hspace{0.4in}$\mathbf{N\overline{\mu} = \\frac{%s}{%d} \\approx %s}$}\n\hspace{0.01in}\n' % (iwidth, depth, one_over_mu, Nmu))
            for (ialpha, alpha) in enumerate(alphas):
                plot = corr_plots[ialpha + len(alphas) * idepth][2]
                if firstcolumn:
                    f.write('\centerline{\parbox[b][%gin][c]{0.65in}{\large ${\\alpha_{\Delta\pi} = %d}$} \includegraphics[width=%.2fin]{%s}}\n' % (plotwidth, alpha, plotwidth, plot))
                else:
                    f.write('\centerline{\includegraphics[width=%.2fin]{%s}}\n' % (plotwidth, plot))
            f.write('\end{minipage}')
            firstcolumn = False
        f.write('\n\end{figure}\n\end{document}')
    subprocess.call(['pdflatex', fname])
    for f in glob.glob('%s.*' % os.path.splitext(fname)[0]):
        if os.path.splitext(f)[1] not in ['.pdf']:
            os.remove(f)

    # make logoplot
    actual_nonzero_diffpref = '%s/actual_nonzero_diffpref.txt' % simulated_dir
    f_out = open(actual_nonzero_diffpref, 'w')
    f_out.write('#POSITION\tNONZERO_DIFFPREF')
    with open(actual_diff_prefs) as f_in:
        for line in f_in:
            if line and line[0] != '#' and not line.isspace():
                if all([float(entry) == 0 for entry in line.split()[2 : ]]):
                    f_out.write('\n%s\tno' % line.split()[0])
                else:
                    f_out.write('\n%s\tyes' % line.split()[0])
    f_out.close()
    subprocess.call(['dms_logoplot', '%s/simulated_depth_1e+07_2_inferred.txt' % simulated_dir, 'inferred_diffprefs_logoplot.pdf', '--nperline', '53', '--overlay1', actual_nonzero_diffpref, '$\\ne 0$?', 'Is actual differential preference non-zero?', '--diffprefheight', '0.45'])


if __name__ == '__main__':
    main()

