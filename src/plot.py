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

Function documentation
-----------------------------
"""


import os
import math
import matplotlib
matplotlib.use('pdf')
import pylab



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
        corr=None, title=False, alpha=1.0, symmetrize=False, fixaxes=False, additionalxy=[], bigmargin=0.29, xsize=1.8, r2=False):
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
    pylab.plot(xs, ys, 'b.', markersize=4, alpha=alpha)
    if additionalxy:
        if len(additionalxy) != 1:
            raise ValueError("Currently additionalxy only works for one additional set of data")
        (x2s, y2s) = additionalxy[0]
        assert len(x2s) == len(y2s) > 0, "additionalxy does not specify two non-empty data lists of equal length"
        xs += x2s # add for later correlation and limit calculations
        ys += y2s
        pylab.plot(x2s, y2s, 'r^', markersize=3, alpha=alpha)
    (xmin, xmax, ymin, ymax) = (min(xs), max(xs), min(ys), max(ys))
    if fixaxes:
        xmin = ymin = 0.0
        xmax = ymax = 1.0
    elif symmetrize:
        xmin = ymin = min(xmin, ymin)
        xmax = ymax = max(xmax, xmin)
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
        yticker = matplotlib.ticker.MaxNLocator(4)
    pylab.gca().yaxis.set_major_locator(yticker)
    if logx:
        xticker = matplotlib.ticker.LogLocator(numticks=4)
    elif fixaxes:
        xticker = matplotlib.ticker.FixedLocator([0, 0.5, 1])
    else:
        xticker = matplotlib.ticker.MaxNLocator(4)
    pylab.gca().xaxis.set_major_locator(xticker)
    pylab.gca().get_xaxis().tick_bottom()
    pylab.gca().get_yaxis().tick_left()
    pylab.savefig(plotfile)
    pylab.clf()
    pylab.close()



if __name__ == '__main__':
    import doctest
    doctest.testmod()
