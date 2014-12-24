"""Master script to test preference inference on simulated deep mutational scanning data.

Analyzes deep mutational scanning data of Melnikov et al to get preferences, then
simulates new data with same site-specific preferences to test
inference methods.

Written by Jesse Bloom."""


import os
import re
import glob
import subprocess
import dms_tools
import dms_tools.plot
import dms_tools.simulate



def SimulateCounts(simulated_dir, actualprefs, fnameprefix):
    """Simulate deep mutational scanning counts. 

    Simulations done at several sequencing depths, and with and without
    sequencing errors.
    
    Returns list of all created files with the simulated data. For each simulation,
    there is a 4-tuple *(nrpre, nrpost, nrerrpre, nrerrpost)*. For simulations with
    no errors, the last two entries are 'none'.
    
    The files are placed in *simulated_dir*, which should be an existing directory.

    The files have the prefix *fnameprefix* followed by a suffix giving information
    about the depth and whether the simulations include sequencing errors.
    
    The simulations are done using the actual preferences in *actualprefs*.
    """    

    assert os.path.isdir(simulated_dir), "Can't find %s" % simulated_dir
    assert os.path.isfile(actualprefs), "Can't find %s" % actualprefs
        
    depths = [1e5, 2e5, 5e5, 1e6, 2e6]
    avgmu = 1.0 / 264.0 # one mutation per gene

    filelist = []

    for (suffix, avgepsilon, avgrho) in [('', 0.0, 0.0), ('_errs', 5e-5, 8e-5)]:
        fprefix = '%s/%s%s_depth' % (simulated_dir, fnameprefix, suffix)
        depths_outfileprefixes = [(depth, '%s_%.0e' % (fprefix, depth)) for depth in depths]
        dms_tools.simulate.SimulateLibraryCounts(depths_outfileprefixes, actualprefs, avgmu, avgepsilon, avgrho)
        for (depth, outfileprefix) in depths_outfileprefixes:
            filetup = ['%s_pre.txt' % outfileprefix, '%s_post.txt' % outfileprefix]
            if avgepsilon:
                filetup.append('%s_errpre.txt' % outfileprefix)
            else:
                filetup.append('none')
            if avgrho:
                filetup.append('%s_errpost.txt' % outfileprefix)
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
    dms_avgprefs = '../processed_data/prefs.txt' # dms_tools inferred prefs
    assert os.path.isfile(dms_avgprefs), "You need to have already run the analysis in the parent directory to create %s" % dms_avgprefs

    # simulate data, then analyze
    simulatedfiles = SimulateCounts(simulated_dir, dms_avgprefs, 'simulated')
    corr_plots = [] 
    for (nrpre, nrpost, nrerrpre, nrerrpost) in simulatedfiles:
        for (infertype, inferargs) in [('preference', []), ('ratio', ['--ratio_estimation', '1.0'])]:
            outfile = '%s_%s-inferred.txt' % ('_'.join(nrpre.split('_')[ : -1]), infertype)
            depth = float(re.search('_depth_(\S+)_', outfile).group(1))
            subprocess.call(['dms_inferprefs', nrpre, nrpost, outfile, '--errpre', nrerrpre, '--errpost', nrerrpost, '--ncpus', '-1', '--excludestop'] + inferargs)
            corrprefix = '%s_corr' % os.path.splitext(outfile)[0]
            corr_plots.append((depth, '%s.pdf' % corrprefix))
            subprocess.call(['dms_correlate', dms_avgprefs, outfile, corrprefix, '--name1', 'actual', '--name2', 'inferred', '--corr_on_plot', '--r2'])

    # make merged summary PDFs
    plotwidth = 2
    latexhead = '\n'.join(['\documentclass[12pt,letterpaper]{article}',\
                           '\usepackage{mathptmx}',\
                           '\usepackage{graphicx}',\
                           '\usepackage{setspace}',\
                           '\usepackage{rotating}',\
                           '\pagestyle{empty}',\
                           '\usepackage[paperheight=%.2fin,paperwidth=%.2fin,top=0.05in,bottom=0.05in,right=0.05in,left=0.05in]{geometry}' % (4.5 * plotwidth, len(corr_plots) // 4 * plotwidth * 1.1 + 0.3),\
                           '\\begin{document}',\
                           '\\begin{figure}',\
                           '\\begin{minipage}{0.3in}',\
                           '{\Huge \\bf A}',\
                           '\\vspace{0.13in}\n',
                           '\\rotatebox{90}{\parbox{3.5in}{\\begin{spacing}{1.3}\centering \large \\bf simulated without sequencing errors \\\\ ratio inference \hspace{0.32in} preference inference \end{spacing}}}',\
                           '\\vspace{0.38in}\n',
                           '{\Huge \\bf B}\n',\
                           '\\vspace{0.13in}\n',
                           '\\rotatebox{90}{\parbox{3.5in}{\\begin{spacing}{1.3}\centering \large \\bf simulated with sequencing errors \\\\ ratio inference \hspace{0.32in} preference inference \end{spacing}}}',\
                           '\\vspace{0.32in}\n',
                           '\end{minipage}',\
                           '\hspace{0.1in}',\
                          ])
    corr_plots.sort()
    fname = 'correlations.tex'
    with open(fname, 'w') as f:
        f.write(latexhead)
        for i in range(len(corr_plots) // 4):
            depth = dms_tools.plot.Base10Formatter(corr_plots[4 * i][0], 2, 0, 0)
            f.write('\n\\begin{minipage}{%.2fin}\n\centerline{\large \hspace{0.4in}\\bf depth of $\mathbf{%s}$}\n\hspace{0.01in}\n\centerline{\includegraphics[width=%.2fin]{%s}}\n\centerline{\includegraphics[width=%.2fin]{%s}}\n\\vspace{0.2in}\n\n\centerline{\includegraphics[width=%.2fin]{%s}}\n\centerline{\includegraphics[width=%.2fin]{%s}}\n\end{minipage}' % (1.05 * plotwidth, depth, plotwidth, corr_plots[4 * i][1], plotwidth, corr_plots[4 * i + 1][1], plotwidth, corr_plots[4 * i + 2][1], plotwidth, corr_plots[4 * i + 3][1]))
        f.write('\n\end{figure}\n\end{document}')
    subprocess.call(['pdflatex', fname])
    for f in glob.glob('%s.*' % os.path.splitext(fname)[0]):
        if os.path.splitext(f)[1] not in ['.pdf']:
            os.remove(f)


if __name__ == '__main__':
    main()

