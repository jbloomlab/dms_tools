"""
========================
``plot`` module
========================
This module uses ``pylab`` and ``matplotlib`` to make plots. 
The ``pdf`` backend is used for ``matplotlib`` / ``pylab``. 


Written by Jesse Bloom.

Functions in this module
---------------------------

`PlotCorrelation` : plots correlation between two variables.

`Base10Formatter` : formats numbers in base 10

`PlotPairedMutFracs` : bar graph of mutation frequencies per sample

`PlotMutCountFracs` : plots fraction of mutations that occur >= a given number of times.

`PlotDepth` : plots per-site depth along primary sequence

`PlotReadStarts` : plots start of read distributions for subassembly.

`PlotSampleBargraphs` : plots bargraphs of categories for several samples.

Function documentation
-----------------------------
"""


import os
import math
import matplotlib
matplotlib.use('pdf')
import pylab
import dms_tools.file_io



def Base10Formatter(number, exp_cutoff, exp_decimal_digits, decimal_digits):
    """Converts a number into Latex formatting with scientific notation.

    Takes a number and converts it to a string that can be shown
    in LaTex using math mode. It is converted to scientific notation
    if the criteria specified by `exp_cutoff` is met.


    `number` : the number to be formatted, should be a float or integer.
    Currently only works for `number >= 0`.

    `exp_cutoff` : convert to scientific notation if ``abs(math.log10(number)) 
    >= exp_cutoff``.

    `exp_decimal_digits` : show this many digits after the decimal if number
    is converted to scientific notation.

    `decimal_digits` : show this many digits after the decimal if number
    is NOT converted to scientific notation.

    The returned value is a LaTex formatted string. If the number is zero, the
    returned string is simply '0'.

    EXAMPLES:

    >>> Base10Formatter(103, 3, 1, 1)
    '103.0'

    >>> Base10Formatter(103.0, 2, 1, 1)
    '1.0 \\\\times 10^{2}'

    >>> Base10Formatter(103.0, 2, 2, 1)
    '1.03 \\\\times 10^{2}'

    >>> Base10Formatter(2892.3, 3, 1, 1) 
    '2.9 \\\\times 10^{3}'

    >>> Base10Formatter(0.0, 3, 1, 1) 
    '0'

    >>> Base10Formatter(0.012, 2, 1, 1)
    '1.2 \\\\times 10^{-2}'
  
    >>> Base10Formatter(-0.1, 3, 1, 1)
    Traceback (most recent call last):
        ...
    ValueError: number must be >= 0
    """
    if number < 0:
        raise ValueError('number must be >= 0')
    if number == 0:
        return '0'
    exponent = int(math.log10(number))
    if math.log10(number) < exponent and number < 1:
        exponent -= 1
    if abs(exponent) >= exp_cutoff:
        x = number / (10.**exponent)
        formatstr = '%.' + '%d' % exp_decimal_digits + 'f \\times 10^{%d}'
        return formatstr % (x, exponent)
    else:
        formatstr = '%.' + '%d' % decimal_digits + 'f'
        return formatstr % number



def PlotCorrelation(xs, ys, plotfile, xlabel, ylabel, logx=False, logy=False,\
        corr=None, title=False, alpha=1.0, symmetrize=False, fixaxes=False, additionalxy=[], bigmargin=0.35, xsize=1.8, r2=False,
        marker_style='b.', additional_marker_style='r^', marker_size=4, additional_marker_size=3):
    """Plots the correlation between two variables as a scatter plot.
    
    The data is plotted as a scatter plot.
    This function uses ``pylab`` / ``matplotlib``. 
    The calling variables use LaTex format for strings. So for
    example, '$10^5$' will print the LaTex equivalent of this 
    string. Similarly, certain raw text strings (such as those
    including underscores) will cause problems if you do not
    escape the LaTex format meaning. For instance, 'x_label'
    will cause a problem since underscore is not valid
    outside of math mode in LaTex, so you would need to use
    'x\_label' to escape the underscore.

    CALLING VARIABLES:

    * *xs* and *ys* are lists of numbers, with the lists
      being of the same length. Entry *xs[i]* is plotted
      on the x-axis against entries *ys[i]* on the y-axis.

    * *plotfile* is a string giving the name of the plot PDF file
      that we create. It should end in the extension ``.pdf``.
      If this plot already exists, it is overwritten.

    * *xlabel* is a string giving the label placed on the x-axis.

    * *ylabel* is a string giving the label placed on the y-axis.

    * *logx* specifies that we log transform the data in *xs*.
      This is *False* by default; set to *True* if you want
      to log transform (base 10 logarithms) the data.

    * *logy* is like *logx*, but for the data in *ys*.

    * *corr* specifies that we include a correlation coefficient on the
      plot. If it has its default value of *None*, then no correlation
      coefficient is included. Otherwise, *corr = (r, p)* where *r* is
      the correlation coefficient (displayed with the label ``R``) and *p*
      is the P-value (displayed with the label ``P``).

    * *r2* : if *True* and using *corr*, then display :math:`R^2` rather
      than :math:`R`.

    * *title* is a string giving the title placed above the plot. 
      It can be *False* if no title is to be used. Otherwise, it should
      be the title string (using LaTex formatting, spaces are allowed).
      Is *False* by default.

    * *alpha* is the transparency of the plotted points. By default
      it is one, which means that there is no transparency. If you make
      the value closer to zero, then the points become partially transparent.
      The rationale is that if you have many overlapping points, making
      them partially transparent helps you better see the density.
      At any position on the plot, the intensity will be saturated
      when there are 1.0 / alpha points plotted. So a reasonable value
      of *alpha* might be something like 0.1.

    * *symmetrize* is an optional argument that is *False* by default.
      If *True*, we make the X and Y limits on the plot the same.

    * *fixaxes* is an optional argument that is *False* by default.
      If *True*, we fix both the X and Y axes to go from 0 to 1, 
      with ticks at 0, 0.5, and 1. If you set this option to *True*,
      then you must set *logx* and *logy* to *False*.
    
    * *additionalxy* is an optional list. By default, it is an empty list.
      However, you can use it to plot additional data points with a different
      color. The main data (specified by *xs* and *ys*) is plotted with
      blue circles. For each addition set of data that you want to plot,
      include a 2-tuple of lists of the form *(x2s, y2s)* in 
      *additionalxy*. Currently, only one such 2-tuple is allowed
      (so only one additional data set). These are plotted as red
      triangles, whereas the main data is plotted as blue
      circles. The same *alpha* value set by the *alpha* parameter applies
      to these points as well. 

    * *bigmargin* is an optional argument that is 0.26 by default. It is
      the fraction of the plot width taken up by the larger margin, which
      is the bottom and left. Make larger if you need more space for axis
      labels.

    * *xsize* is an optional argument that is the width of the plot in inches.
      It is 1.8 by default.

    * *marker_style* and *additional_marker_style* are optional arguments to 
      change the color/style of the marker for the main and additional data points,
      respectively. See the matplotlib documentation for a full explanation.

    * *marker_size* and *additional_marker_size* are optional arguments to
      set the size of the markers.
    """
    if not os.path.splitext(plotfile)[1].upper() == '.PDF':
        raise ValueError("plotfile does not end in PDF extension: %s " % plotfile)
    if not (len(xs) == len(ys) >= 2):
        raise ValueError("xs and ys do not specify lists of the same length with >= 2 entries")
    if fixaxes and (logy or logx):
        raise ValueError("Cannot use fixaxes with logx or logy")
    smallmargin = 0.04
    (lmargin, rmargin, bmargin, tmargin) = (bigmargin, smallmargin, bigmargin, smallmargin)
    plotmargin = 0.03 # add this much above and below the last data point
    logplotmargin = 2.0 # scale limits by this much if log scale
    if title:
        titlemargin = 0.09 * (0.8 + title.count('\n') + title.count('\\\\'))
        tmargin += titlemargin
    ysize = xsize * (1.0 - lmargin - rmargin) / (1.0 - tmargin - bmargin)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=9)
    matplotlib.rc('legend', fontsize=10)
    figure = pylab.figure(figsize=(xsize, ysize), facecolor='white')
    ax = pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
    pylab.plot(xs, ys, marker_style, markersize=marker_size, alpha=alpha, mew=0, clip_on=False)
    if additionalxy:
        if len(additionalxy) != 1:
            raise ValueError("Currently additionalxy only works for one additional set of data")
        (x2s, y2s) = additionalxy[0]
        assert len(x2s) == len(y2s) > 0, "additionalxy does not specify two non-empty data lists of equal length"
        xs += x2s # add for later correlation and limit calculations
        ys += y2s
        pylab.plot(x2s, y2s, additional_marker_style, markersize=additional_marker_size, alpha=alpha)
    (xmin, xmax, ymin, ymax) = (min(xs), max(xs), min(ys), max(ys))
    if fixaxes:
        xmin = ymin = 0.0
        xmax = ymax = 1.0
    elif symmetrize:
        xmin = ymin = min(xmin, ymin)
        xmax = ymax = max(xmax, ymax)
    if logy:
        pylab.gca().set_yscale('log')
        ax.set_ylim([ymin / logplotmargin, ymax * logplotmargin])
        ys = [math.log(y) for y in ys]
    else:
        ymargin = plotmargin * (ymax - ymin)
        ax.set_ylim([ymin - ymargin, ymax + ymargin])
    if logx:
        pylab.gca().set_xscale('log')
        ax.set_xlim([xmin / logplotmargin, xmax * logplotmargin])
        xs = [math.log(x) for x in xs]
    else:
        xmargin = plotmargin * (xmax - xmin)
        ax.set_xlim([xmin - xmargin, xmax + xmargin])
    pylab.xlabel(xlabel, size=10)
    pylab.ylabel(ylabel, size=10)
    if title:
        pylab.title(title, size=10)
    if corr:
        (r, p) = corr
        if r2:
            r = '$R^2 = %.2f$' % (r**2)
        else:
            r = '$R = %.2f$' % r
        if p < 1e-10:
            p = '$P < 10^{-10}$'
        else:
            p = '$P = %s$' % Base10Formatter(p, 2, 1, 2)
        text = '%s\n%s' % (r, p)
        pylab.text(0.05, 0.96, text, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, size=10)
    if logy:
        yticker = matplotlib.ticker.LogLocator(numticks=4)
    elif fixaxes:
        yticker = matplotlib.ticker.FixedLocator([0, 0.5, 1])
    else:
        yticker = matplotlib.ticker.MaxNLocator(4, integer=True)
    pylab.gca().yaxis.set_major_locator(yticker)
    if logx:
        xticker = matplotlib.ticker.LogLocator(numticks=4)
    elif fixaxes:
        xticker = matplotlib.ticker.FixedLocator([0, 0.5, 1])
    else:
        xticker = matplotlib.ticker.MaxNLocator(4, integer=True)
        
    spineOffset = {'left': 3, 'bottom': 3}    
    [spine.set_position(('outward',spineOffset[loc])) if loc in ['left','bottom'] else spine.set_color('none') for loc, spine in ax.spines.items() ] 
    ax.tick_params(axis='x', direction='out')
    ax.tick_params(axis='y', direction='out')
    
    pylab.gca().xaxis.set_major_locator(xticker)
    pylab.gca().get_xaxis().tick_bottom()
    pylab.gca().get_yaxis().tick_left()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotDepth(codon_counts, names, plotfile, mutdepth=False, y_axis_label=None, 
        separatemuttypes=False):
    """Plots per-site depth along primary sequence.

    `codon_counts` : a list of dictionaries giving the codon counts for
    each sample as read by `dms_tools.file_io.ReadDMSCounts`.

    `names` : a list of strings giving the names of the samples
    corresponding to each entry in `codon_counts`.

    `plotfile` : name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    `mutdepth` : Boolean switch, if *True* then rather than plotting
    sequencing depth, we plot per-site mutation rate.

    `y_axis_label` : a string, if specified, overrides the default y-axis
    label 'reads' or 'mutation frequency'.

    `separatemuttypes` : Boolean switch specifying that we plot a different
    line for each of nonsynonymous, synonymous, and stop codon mutations.
    Only can be used if *mutdepth* is *True*.
    """
    if os.path.splitext(plotfile)[1].lower() != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    if separatemuttypes and not mutdepth:
        raise ValueError("Can't use separatemuttypes without mutdepth")
    assert len(codon_counts) == len(names) > 0, "codon_counts and names are not non-empty lists of the same length"
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)
    matplotlib.rc('legend', fontsize=10)
    sites = list(codon_counts[0].iterkeys())
    codon_counts = [dms_tools.utils.ClassifyCodonCounts(counts) for counts in codon_counts]
    xs = list(range(len(sites)))
    dms_tools.utils.NaturalSort(sites)
    if separatemuttypes:
        nlegendcols = 2
        nlegendrows = int(math.ceil(3 * len(names) / float(nlegendcols)))
    else:
        nlegendcols = 4
        nlegendrows = int(math.ceil(len(names) / float(nlegendcols)))
    fig = pylab.figure(figsize=(5.5, 2.16 * (0.76 + 0.23 * nlegendrows)))
    (lmargin, rmargin, bmargin, tmargin) = (0.1, 0.02, 0.19, 0.01 + 0.11 * nlegendrows)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 - bmargin - tmargin])
    lines = []
    all_ys = []
    legendnames = []
    for (name, counts) in zip(names, codon_counts):
        if mutdepth:
            if separatemuttypes:
                muttypes = [(' nonsyn', 'N_NS', 1), (' syn X 20', 'N_SYN', 20), (' stop X 20', 'N_STOP', 20)]
            else:
                muttypes = [('', 'all', 1)]
            for (mtype, key, scalefac) in muttypes:
                ys = []
                for r in sites:
                    if counts[r]['COUNTS']:
                        if key == 'all':
                            ys.append((counts[r]['COUNTS'] - counts[r]['N_WT']) / float(counts[r]['COUNTS']))
                        else:
                            ys.append(scalefac * counts[r][key] / float(counts[r]['COUNTS']))
                    else:
                        ys.append(0)
                all_ys += ys
                line = pylab.plot(xs, ys, lw=1.2)
                lines.append(line[0])
                legendnames.append(name.replace('_', ' ') + mtype)
        else:
            ys = [counts[r]['COUNTS'] for r in sites]
            all_ys += ys
            line = pylab.plot(xs, ys, lw=1.2)
            lines.append(line[0])
            legendnames.append(name.replace('_', ' '))
    pylab.xlabel('codon position')
    if mutdepth:
        pylab.ylabel('mutation frequency')
    else:
        pylab.ylabel('number of reads')
    if y_axis_label:
        pylab.ylabel(y_axis_label)
    all_ys.sort()
    # if the top y value is excessively large, set a small ymax to avoid distorting the y-axis
    if all_ys[-1] >= 2.5 * all_ys[-len(legendnames) - 1]:
        pylab.gca().set_ylim([0, 2.5 * all_ys[-len(legendnames) - 1]])
    else:
        pylab.gca().set_ylim([0, pylab.gca().get_ylim()[1]])
    yticker = matplotlib.ticker.MaxNLocator(4)
    pylab.gca().yaxis.set_major_locator(yticker)
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-2, 3))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    pylab.gca().set_xlim([0, len(xs) - 1])
    if len(xs) <= 250:
        # xtick every 50
        xlocator = matplotlib.ticker.FixedLocator(list(range(0, len(xs), 50)))
        xformatter = matplotlib.ticker.FixedFormatter([sites[i] for i in range(0, len(xs), 50)])
    else:
        # xtick every 100
        xlocator = matplotlib.ticker.FixedLocator(list(range(0, len(xs), 100)))
        xformatter = matplotlib.ticker.FixedFormatter([sites[i] for i in range(0, len(xs), 100)])
    pylab.gca().xaxis.set_major_locator(xlocator)
    pylab.gca().xaxis.set_major_formatter(xformatter)
    pylab.legend(lines, legendnames, handlelength=1.2, handletextpad=0.3, columnspacing=0.9, bbox_to_anchor=(0.54, 1.03 + 0.14 * nlegendrows), loc='upper center', ncol=nlegendcols)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotPairedMutFracs(codon_counts, names, plotfile, ylabel='fraction'):
    """Makes paired bar graph of mutation fractions per codon.
    
    For each sample, calculates the fraction of codon counts
    that represent synonymous / nonsynoymous /
    stop codons, and one- / two- / three-nucleotide mutations. Makes
    paired stacked bar graphs showing these fractions. The fractions
    for all sites in the gene are aggregated together. 

    `codon_counts` : a list of dictionaries giving the codon counts for
    each sample as read by `dms_tools.file_io.ReadDMSCounts`.

    `names` : a list of strings giving the names of the samples
    corresponding to each entry in `codon_counts`.

    `plotfile` : name of the output plot file created by this method
    (such as 'plot.pdf'). The extension must be ``.pdf``.

    `ylabel` : the text placed on the ylabel.
    """
    if os.path.splitext(plotfile)[1].lower() != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    names = [name.replace('_', '\_') for name in names]
    if not (isinstance(codon_counts, list)) and codon_counts:
        raise ValueError("codon_counts does not specify a non-empty list")
    if len(codon_counts) != len(names):
        raise ValueError("codon_counts and names differ in length")
    bar1 = ['synonymous', 'nonsynonymous', 'stop codon']
    bar2 = ['1 nucleotide mutation', '2 nucleotide mutations', '3 nucleotide mutations']
    d = dict([(key, []) for key in bar1 + bar2])
    for counts in codon_counts:
        counts = dms_tools.utils.ClassifyCodonCounts(counts)
        denom = float(counts['TOTAL_COUNTS'])
        if not denom:
            raise ValueError("no counts for a sample")
        for key in bar1 + bar2:
            if key == 'nonsynonymous':
                d[key].append(counts['TOTAL_NS'] / denom)
            elif key == 'synonymous':
                d[key].append(counts['TOTAL_SYN'] / denom)
            elif key == 'stop codon':
                d[key].append(counts['TOTAL_STOP'] / denom)
            elif key == '1 nucleotide mutation':
                d[key].append(counts['TOTAL_N_1MUT'] / denom)
            elif key == '2 nucleotide mutations':
                d[key].append(counts['TOTAL_N_2MUT'] / denom)
            elif key == '3 nucleotide mutations':
                d[key].append(counts['TOTAL_N_3MUT'] / denom)
            else:
                raise ValueError("Invalid key of %s" % key)
    nsamples = len(names)
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=10)
    matplotlib.rc('legend', fontsize=10)
    matplotlib.rc('xtick', labelsize=10)
    matplotlib.rc('patch', linewidth=0.5)
    # make stacked bar graph
    fig = pylab.figure(figsize=(6.6, 3.75))
    (lmargin, rmargin, bmargin, tmargin) = (0.08, 0.01, 0.43, 0.13)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 -\
            bmargin - tmargin])
    bars = []
    for (ibar, keys, colors) in [(0, bar1, 'brg'), (1, bar2, 'myc')]:
        bottoms = [0] * nsamples
        barwidth = 0.35
        indices = [i - barwidth + ibar * barwidth for i in range(1, nsamples + 1)]
        totalheights = [0 for i in range(nsamples)]
        for (key, color) in zip(keys, colors):
            totalheights = [totalheights[i] + d[key][i] for i in range(nsamples)]
            b = pylab.bar(indices, d[key], width=barwidth, bottom=bottoms, color=color)
            bars.append(b)
            for i in range(nsamples):
                bottoms[i] += d[key][i]
    ymax = max(bottoms)
    pylab.gca().set_ylim([0, 1.08 * ymax])
    pylab.xticks([i for i in indices], names, rotation=90)
    pylab.gca().set_xlim([0.4, nsamples + 0.6])
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-2, 2))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    yticker = matplotlib.ticker.MaxNLocator(5)
    pylab.gca().yaxis.set_major_locator(yticker)
    barmarks = [b[0] for b in bars]
    barlabels = bar1 + bar2
    # reorder bar labels since pylab puts them down columns
    barmarks = [barmarks[0], barmarks[3], barmarks[1], barmarks[4], barmarks[2], barmarks[5]]
    barlabels = [barlabels[0], barlabels[3], barlabels[1], barlabels[4], barlabels[2], barlabels[5]]
    pylab.legend(barmarks, barlabels, handlelength=1.1,\
            bbox_to_anchor=(0.54, 1.31), loc='upper center', ncol=3)
    pylab.ylabel(ylabel, size=10)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotMutCountFracs(plotfile, title, names, all_cumulfracs, syn_cumulfracs, all_counts, syn_counts, legendloc, writecounts=True, nmax=None):
    """Plots fraction of mutations with >= a given number of counts.

    Does this for both all mutations and synonymous mutations. The plots
    are placed side by side.

    The plots show the fractions of mutations
    that are found >= *n* times for a range of values of *n*.

    CALLING VARIABLES:

    * *plotfile* : string giving name of the created plot file. Must end in
      the extension ``.pdf``.

    * *title* : string giving the title placed about the plot.

    * *names* : a list of strings giving the names of the samples to
      plot.

    * *all_cumulfracs* : a list of the same length as *names* giving
      the cumulative fractions for all mutations. Each entry should be a list, 
      and all of these lists should be the same length. *cumulfracs[i][n]*
      gives the fraction of mutations to sample *names[i]* that
      are found >= *n* times. The x-axis of the created plot
      will go from 0 to *len(all_cumulfracs[0]) - 1*.

    * *syn_cumulfracs* : a list like *all_cumulfracs* except for synonymous
      mutations.

    * *all_counts* : integer counts of all mutations (the total number of mutations
      used for *all_cumulfracs*), used to create plot title.

    * *syn_counts* : like *all_counts* but for synonymous mutations.

    * *legendloc* : specifies the location of the legend. Should be a string.
      Valid values are:

        - *bottom* : put legend at the bottom of the plot.

        - *right* : put legend at the right of the plot.

    * *writecounts* is a Boolean switch specifying whether we include the counts of all
      mutations (specified by *all_counts* and *syn_counts*) in the plot title. We do
      this if *writecounts* is *True*, and do not if it is *False*. 

    *nmax* if specified, should be an integer > 1 giving the x-axis maximum.
    """

    if os.path.splitext(plotfile)[1].lower() != '.pdf':
        raise ValueError('plotfile must end in .pdf: %s' % plotfile)

    # plot setup stuff
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=11)
    matplotlib.rc('legend', fontsize=11)
    if legendloc == 'bottom':
        ncol = 4 # number of columns
        legendrowheight = 0.2
        nrows = math.ceil(len(names) / float(ncol))
        xlegendmargin = 0.0
    elif legendloc == 'right':
        ncol = 1
        nrows = 0 # we specify by columns
        legendrowheight = 0
        xlegendmargin = 1.45
    else:
        raise ValueError("Invalid legendloc of %s" % legendloc)
    (xsize, ysize) = (4.75 + xlegendmargin, legendrowheight * nrows + 2.4)
    styles = ['k-.', 'y-.', 'b-', 'r:', 'g:', 'c--', 'm--', 'k-', 'y-', 'b:', 'r-', 'g-', 'c:', 'm:', 'k--', 'y:', 'b--', 'g--', 'c-', 'm-', 'k:']
    styles = ['k-.', 'y-.', 'b-', 'r:', 'g:', 'c--', 'm--', 'k:', 'y-', 'b-.', 'r-', 'g--', 'm:']
    lmargin = 0.11 * (xsize - xlegendmargin) / xsize # left margin for plot
    rmargin = 0.01 + xlegendmargin / xsize # right margin for plot
    centermargin = 0.02 # horizontal space between plots
    tmargin = 0.16 # top margin
    bmargin = 0.22 + (legendrowheight * nrows) / ysize # bottom margin
    plotwidth = (1.0 - lmargin - rmargin - centermargin) / 2.0

    assert 0 < len(syn_cumulfracs) == len(all_cumulfracs) == len(names) <= len(styles), "Must specify equal numbers of fractions and names, and no more than %d" % len(styles)

    # start plotting
    fig = pylab.figure(figsize=(xsize, ysize))
    ax_all = pylab.axes([lmargin, bmargin, plotwidth, 1 - bmargin - tmargin])
    ax_syn = pylab.axes([lmargin + plotwidth + centermargin, bmargin, plotwidth, 1 - bmargin - tmargin])
    for (cumulfracs, ax, write_ylabel, ax_title) in [(syn_cumulfracs, ax_syn, False, 'synonymous (%d total)' % syn_counts), (all_cumulfracs, ax_all, True, 'all (%d total)' % all_counts)]:
        if not writecounts:
            ax_title = ax_title.split()[0]
        pylab.axes(ax)
        if not nmax:
            nmax = max([len(x) for x in cumulfracs])
            assert nmax, "Length of entries in cumulfracs must be >= 1"
        else:
            assert nmax > 1, "nmax should be > 1"
        lines = []
        xs = [n for n in range(0, nmax)]
        for i in range(len(names)):
            i_cumulfracs = (cumulfracs[i] + [0] * (nmax - len(cumulfracs[i])))[ : nmax]
            plotline = pylab.plot(xs, i_cumulfracs, styles[i], lw=1.5)
            lines.append(plotline[0])
            pylab.xlabel("Mutation counts", size=11)
            if write_ylabel:
                pylab.ylabel("Frac. $\ge$ this many counts", size=11)
            else:
                yformatter = matplotlib.ticker.NullFormatter()
                pylab.gca().yaxis.set_major_formatter(yformatter)
        pylab.gca().set_ylim([-0.02, 1.02])
        yticker = matplotlib.ticker.FixedLocator([0.0, 0.5, 1.0])
        pylab.gca().yaxis.set_major_locator(yticker)
        pylab.gca().set_xlim([0, nmax - 1])
        xticker = matplotlib.ticker.MaxNLocator(4)
        pylab.gca().xaxis.set_major_locator(xticker)
        if title:
            pylab.title(ax_title, size=11)
    pylab.suptitle("{\\bf %s}" % title, size=11)
    if legendloc == 'bottom':
        fig.legend(lines, [name.replace('_', ' ') for name in names], handlelength=2.25, handletextpad=0.2, columnspacing=0.8, ncol=ncol, bbox_to_anchor=(0.5, -0.01), loc='lower center')
    elif legendloc == 'right':
        fig.legend(lines, [name.replace('_', ' ') for name in names], handlelength=2.25, handletextpad=0.2, columnspacing=0.8, ncol=ncol, bbox_to_anchor=(1.0, 0.52), loc='center right')
    else:
        raise ValueError("Invalid legendloc of %s" % legendloc)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotReadStarts(names, depthbystart, plotfile):
    """Plots distribution of read starts.

    *names* : list of sample names.

    *depthbystart* : dictionary keyed by each name in *names*. V
    Value is dictionary keyed by read start position, value is
    number of reads. All names must have same start positions.

    *plotfile* : name of created PDF plot.
    """
    if os.path.splitext(plotfile)[1].lower() != '.pdf':
        raise ValueError("plotfile must in in .pdf but got %s" % plotfile)
    assert len(names) == len(depthbystart) > 0, "names and depthbystart are not non-empty and of the same length"
    starts = list(depthbystart[0].keys())
    starts.sort()
    assert all([set(istarts.keys()) == set(starts) for istarts in depthbystart]), "Not same starts for all samples"
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('legend', fontsize=10)
    matplotlib.rc('font', size=10)
    fig = pylab.figure(figsize=(5, 3.25))
    (lmargin, rmargin, bmargin, tmargin) = (0.1, 0.02, 0.12, 0.21)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 - bmargin - tmargin])
    lines = []
    styles = ['bo-', 'rs-', 'gd-.', 'c^:', 'mx-', 'y*--', 'kv-.']
    assert len(names) <= len(styles), "too many names for specified styles"
    for (name, depth, style) in zip(names, depthbystart, styles):
        line = pylab.plot(starts, [depth[start] for start in starts], style)
        lines.append(line[0])
    pylab.legend(lines, [name.replace('_', ' ') for name in names], ncol=int(math.ceil(len(names) / 2.0)), loc='upper center', numpoints=1, bbox_to_anchor=(0.54, 1.25))
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-2, 4))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    pylab.gca().set_xlim([starts[0] - 0.03 * (starts[-1] - starts[0]), starts[-1] + 0.03 * (starts[-1] - starts[0])])
    pylab.xticks(starts)
    pylab.xlabel('read start position')
    pylab.ylabel('number of reads')
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


def PlotSampleBargraphs(names, categories, data, plotfile, ylabel, groupbyfirstword=True, ncolumns=3):
    """Plots bargraph of counts for different samples.

    *names* is a list of the names of the different samples.

    *categories* is a list of the different categories for each sample.

    *data* is a dictionary keyed by each name in *names*, and that value
    is in turn a dictionary keyed by each category in *categories*.

    *plotfile* is the name of the PDF file that we are creating.

    *ylabel* is the label placed on the y-axis.

    *groupbyfirstword* is a Boolean switch that specifies whether we group 
    consecutive categories with the same first word in the string to have
    the same color.

    *ncolumns* is the number of columns in each legend line.
    """
    if os.path.splitext(plotfile)[1].lower() != '.pdf':
        raise ValueError("plotfile must in in .pdf but got %s" % plotfile)
    assert len(names) == len(data) > 0, "names and data are not non-empty and of the same length"
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    hatches = ['', '//////', '*****', '.....', '*', 'o'] 
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('legend', fontsize=10)
    matplotlib.rc('font', size=10)
    matplotlib.rc('patch', linewidth=0.75)
    ncol = ncolumns # number of legend columns
    nlegendrows = int(math.ceil(len(categories) / float(ncol)))
    fig = pylab.figure(figsize=(6, 3.96 * (1.02 + 0.055 * nlegendrows)))
    (lmargin, rmargin, bmargin, tmargin) = (0.09, 0.01, 0.41, 0.0175 + 0.046 * nlegendrows)
    ax = pylab.axes([lmargin, bmargin, 1 - lmargin - rmargin, 1 - bmargin - tmargin])
    bottoms = [0] * len(names)
    bars = []
    barwidth = 0.5
    indices = [i - barwidth / 2. for i in range(1, len(names) + 1)]
    icolor = ihatch = 0
    previousfirstword = None
    for category in categories:
        if previousfirstword == None:
            pass
        elif groupbyfirstword:
            if previousfirstword == category.split()[0]:
                ihatch += 1
            else:
                ihatch = 0
                icolor += 1
        else:
            icolor += 1
        previousfirstword = category.split()[0]
        b = pylab.bar(indices, [data[name][category] for name in names], width=barwidth, bottom=bottoms, color=colors[icolor % len(colors)], hatch=hatches[ihatch % len(hatches)])
        bars.append(b)
        for (iname, name) in enumerate(names):
            bottoms[iname] += data[name][category]
    ymax = max(bottoms)
    pylab.gca().set_ylim([0, 1.04 * ymax])
    pylab.xticks([i + barwidth / 2. for i in indices], [name.replace('_', ' ') for name in names], rotation=90)
    pylab.gca().set_xlim([0.5, len(names) + 0.5])
    yformatter = pylab.ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-3, 3))
    pylab.gca().yaxis.set_major_formatter(yformatter)
    pylab.legend([b[0] for b in bars], categories, handletextpad=0.3, handlelength=1.3, columnspacing=1.1, bbox_to_anchor=(0.53, 1.037 + 0.11 * nlegendrows), loc='upper center', ncol=ncol)
    pylab.ylabel(ylabel)
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
