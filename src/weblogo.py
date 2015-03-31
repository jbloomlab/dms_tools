"""
============================
``weblogo`` module
============================

Module for making sequence logos with the *weblogolib* package distributed with
``weblogo`` (http://weblogo.threeplusone.com/). This module interfaces with the 
*weblogolib* API, and so is only known to work with *weblogolib* version 3.4.

Written by Jesse Bloom.


Functions in this module
----------------------------------
* *LogoPlot* : plots a sequence logo of site-specific preferences or differential preferences.

* *LogoOverlay* : makes logo plot overlays, called by *LogoPlot*.

* *KyteDoolittleColorMapping* : maps amino-acid hydrophobicities to colors.

Function documentation
-----------------------------

"""


import os
import tempfile
import string
import math
import shutil
import numpy
import matplotlib
matplotlib.use('pdf')
import pylab
import PyPDF2
import dms_tools.utils
# the following are part of the weblogo 3.4 library
import weblogolib # weblogo library
if weblogolib.__version__ != '3.4':
    raise ImportError("This module only works with weblogolib version 3.4. You have version %s" % weblogolib.__version__)
import weblogolib.colorscheme # weblogo library
import corebio.matrix # weblogo library
import corebio.utils # weblogo library


def KyteDoolittleColorMapping(maptype='jet', reverse=True):
    """Maps amino-acid hydrophobicities to colors.

    Uses the Kyte-Doolittle hydrophobicity scale defined by::

        J. Kyte & R. F. Doolittle: 
        "A simple method for displaying the hydropathic character of a protein." 
        J Mol Biol, 157, 105-132

    More positive values indicate higher hydrophobicity, while more negative values
    indicate lower hydrophobicity.

    The returned variable is the 3-tuple *(cmap, mapping_d, mapper)*:

        * *cmap* is a ``pylab`` *LinearSegmentedColorMap* object.

        * *mapping_d* is a dictionary keyed by the one-letter amino-acid
          codes. The values are the colors in CSS2 format (e.g. #FF0000
          for red) for that amino acid. The value for a stop codon 
          (denoted by a ``*`` character) is black (#000000).

        * *mapper* is the actual *pylab.cm.ScalarMappable* object.

    The optional calling argument *maptype* should specify a valid ``pylab`` color map.

    The optional calling argument *reverse* specifies that we set up the color
    map so that the most hydrophobic residue comes first (in the Kyte-Doolittle
    scale the most hydrophobic comes last as it has the largest value). This option
    is *True* by default as it seems more intuitive to have charged residues red
    and hydrophobic ones blue.
    """
    d = {'A':1.8, 'C':2.5, 'D':-3.5, 'E':-3.5, 'F':2.8, 'G':-0.4,\
         'H':-3.2, 'I':4.5, 'K':-3.9, 'L':3.8, 'M':1.9, 'N':-3.5,\
         'P':-1.6, 'Q':-3.5, 'R':-4.5, 'S':-0.8, 'T':-0.7,\
         'V':4.2, 'W':-0.9, 'Y':-1.3}
    hydrophobicities = [d[aa] for aa in dms_tools.aminoacids_nostop]
    if reverse:
        hydrophobicities = [-1 * x for x in hydrophobicities]
    mapper = pylab.cm.ScalarMappable(cmap=maptype)
    mapper.set_clim(min(hydrophobicities), max(hydrophobicities))
    mapping_d = {'*':'#000000'}
    for (aa, h) in zip(dms_tools.aminoacids_nostop, hydrophobicities):
        tup = mapper.to_rgba(h, bytes=True)
        (red, green, blue, alpha) = tup
        mapping_d[aa] = '#%02x%02x%02x' % (red, green, blue)
        assert len(mapping_d[aa]) == 7
    cmap = mapper.get_cmap()
    return (cmap, mapping_d, mapper)



def LogoPlot(sites, datatype, data, plotfile, nperline, numberevery=10, allowunsorted=False, ydatamax=1.01, overlay=None):
    """Constructs a sequence logo showing amino-acid or nucleotide preferences.

    The heights of each letter is equal to the preference of
    that site for that amino acid or nucleotide.

    Note that stop codons may or may not be included in the logo
    depending on whether they are present in *pi_d*.  

    CALLING VARIABLES:

    * *sites* is a list of all of the sites that are being included
      in the logo, as strings or numbers. They must be in natural sort order (as
      is done by *dms_tools.utils.NaturalSort*) or an error will
      be raised **unless** *allowunsorted* is *True*. The sites
      in the plot are ordered in the same arrangement
      listed in *sites*. These should be **strings**, not integers.

    * *datatype* should be one of the two following strings depending on whether
      we are making a plot of preferences or differential preferences: 
      'prefs' or 'diffprefs'

    * *data* is a dictionary that has a key for every entry in
      *sites*. For every site *r* in *sites*, *sites[r][x]*
      is the preference or differential preference for character *x*. 
      Preferences must sum to one; differential preferences to zero.
      All sites must have the same set of characters. The characters
      must be the set of nucleotides (*dms_tools.nts*)
      or the set of amino acids with or without stop codons
      (*dms_tools.aminoacids_nostop* or *dms_tools.aminoacids_withstop*).

    * *plotfile* is a string giving the name of the created PDF file 
      of the logo plot.
      It must end in the extension ``.pdf``.

    * *nperline* is the number of sites per line. Often 40 to 80 are good values.

    * *numberevery* is specifies how frequently we put labels for the sites on
      x-axis.

    * *allowunsorted* : if *True* then we allow the entries in *sites* to 
      **not** be sorted. This means that the logo plot will **not** have
      sites in sorted order.

    * *ydatamax* : meaningful only if *datatype* is 'diffprefs'. In this case, it gives
      the maximum that the logo stacks extend in the positive and negative directions.
      Cannot be smaller than the maximum extent of the differential preferences.

    * *overlay* : this argument allows you to make overlay bars that indicated
      other properties for the sites. By default, this option is *None*, meaning that
      no overlay is created. If you set it to something else, it must be a list
      giving either one or two properties. Each property is a tuple:
      *(prop_d, shortname, longname)* where:

        - *prop_d* is a dictionary keyed by site numbers that are in *sites*.
          For each *r* in *sites*, *prop_d[r]* gives the value of the property,
          or if there is no entry in *prop_d* for *r*, then the property
          is undefined and is colored white. Properties can either be:

            * continuous: in this case, all of the values should be numbers.

            * discrete : in this case, all of the values should be strings.
              While in practice, if you have more than a few discrete
              categories (different strings), the plot will be a mess.

        - *shortname* : short name for the property; will not format well
          if more than 4 or 5 characters.

        - *longname* : longer name for property used on axes label. Can be the
          same as *shortname* if you don't need a different long name.

    """
    assert datatype in ['prefs', 'diffprefs']

    # check data, and get characters
    assert sites, "No sites specified"
    assert set(sites) == set(data.keys()), "There is not a complete match between sites and the keys of data"
    characters = data[sites[0]].keys()
    if set(characters) == set(dms_tools.nts):
        alphabet_type = 'nt'
    elif set(characters) == set(dms_tools.aminoacids_nostop) or set(characters) == set(dms_tools.aminoacids_withstop):
        alphabet_type = 'aa'
    else:
        raise ValueError("Invalid set of character keys in data. Do not specify either nucleotides or amino acids:\n%s" % str(characters))
    for r in sites:
        if set(data[r].keys()) != set(characters):
            raise ValueError("Not all sites in data have the same set of characters")

    firstblankchar = 'B' # character for first blank space for diffprefs
    assert firstblankchar not in characters, "firstblankchar in characters"
    lastblankchar = 'b' # character for last blank space for diffprefs
    assert lastblankchar not in characters, "lastblankchar in characters"
    separatorchar = '-' # separates positive and negative for diffprefs
    assert lastblankchar not in characters, "lastblankchar in characters"
    separatorheight = 0.02 # height of separator as fraction of total for diffprefs

    if os.path.splitext(plotfile)[1].lower() != '.pdf':
        raise ValueError("plotfile must end in .pdf: %s" % plotfile)
    if os.path.isfile(plotfile):
        os.remove(plotfile) # remove existing plot

    if not allowunsorted:
        sorted_sites = [r for r in sites]
        dms_tools.utils.NaturalSort(sorted_sites)
        if sorted_sites != sites:
            raise ValueError("sites is not properly sorted")

    # Following are specifications of weblogo sizing taken from its documentation
    stackwidth = 9.5 # stack width in points, not default size of 10.8, but set to this in weblogo call below
    barheight = 5.5 # height of bars in points if using overlay
    barspacing = 2.0 # spacing between bars in points if using overlay
    stackaspectratio = 4.4 # ratio of stack height:width, doesn't count part going over maximum value of 1
    if overlay:
        if not (1 <= len(overlay) <= 2):
            raise ValueError("overlay must be a list of one or two entries; instead it had %d entries" % len(overlay))
        ymax = (stackaspectratio * stackwidth + len(overlay) * (barspacing + barheight)) / float(stackaspectratio * stackwidth)
        aspectratio = ymax * stackaspectratio # effective aspect ratio for full range
    else:
        ymax = 1.0 
        aspectratio = stackaspectratio
    rmargin = 11.5 # right margin in points, fixed by weblogo
    stackheightmargin = 16 # margin between stacks in points, fixed by weblogo

    try:
        # write data into transfacfile (a temporary file)
        (fd, transfacfile) = tempfile.mkstemp()
        f = os.fdopen(fd, 'w')
        ordered_alphabets = {} # keyed by site index (0, 1, ...) with values ordered lists for characters from bottom to top
        if datatype == 'prefs':
            chars_for_string = characters
            f.write('ID ID\nBF BF\nP0 %s\n' % ' '.join(chars_for_string))
            for (isite, r) in enumerate(sites):
                f.write('%s %s\n' % (r, ' '.join([str(data[r][x]) for x in characters])))
                pi_r = [(data[r][x], x) for x in characters]
                pi_r.sort()
                ordered_alphabets[isite] = [tup[1] for tup in pi_r] # order from smallest to biggest
        elif datatype == 'diffprefs':
            chars_for_string = characters + [firstblankchar, lastblankchar, separatorchar]
            ydatamax *= 2.0 # maximum possible range of data, multiply by two for range
            f.write('ID ID\nBF BF\nP0 %s\n' % ' '.join(chars_for_string))
            for (isite, r) in enumerate(sites):
                positivesum = sum([data[r][x] for x in characters if data[r][x] > 0]) + separatorheight / 2.0
                negativesum = sum([data[r][x] for x in characters if data[r][x] < 0]) - separatorheight / 2.0
                if abs(positivesum + negativesum) > 1.0e-3:
                    raise ValueError("Differential preferences sum of %s is not close to zero for site %s" % (positivesum + negativesum, r))
                if 2.0 * positivesum > ydatamax:
                    raise ValueError("You need to increase ydatamax: the total differential preferences sum to more than the y-axis limits. Right now, ydatamax is %.3f while the total differential preferences are %.3f" % (ydatamax, 2.0 * positivesum))
                f.write('%s' % r)
                deltapi_r = []
                for x in characters:
                    deltapi_r.append((data[r][x], x))
                    f.write(' %s' % (abs(data[r][x]) / float(ydatamax)))
                deltapi_r.sort()
                firstpositiveindex = 0
                while deltapi_r[firstpositiveindex][0] < 0:
                    firstpositiveindex += 1
                ordered_alphabets[isite] = [firstblankchar] + [tup[1] for tup in deltapi_r[ : firstpositiveindex]] + [separatorchar] + [tup[1] for tup in deltapi_r[firstpositiveindex : ]] + [lastblankchar] # order from most negative to most positive with blank characters and separators
                f.write(' %g %g %g\n' % (0.5 * (ydatamax + 2.0 * negativesum) / ydatamax, 0.5 * (ydatamax + 2.0 * negativesum) / ydatamax, separatorheight)) # heights for blank charactors and separators
        else:
            raise ValueError("Invalid datatype of %s" % datatype)
        f.close()

        # create web logo
        charstring = ''.join(chars_for_string)
        assert len(charstring) == len(chars_for_string), "Length of charstring doesn't match length of chars_for_string. Do you have unallowable multi-letter characters?\n%s" % (str(chars_for_string))
        logoprior = weblogolib.parse_prior('equiprobable', charstring, 0)
        motif = _my_Motif.read_transfac(open(transfacfile), charstring)
        logodata = weblogolib.LogoData.from_counts(motif.alphabet, motif, logoprior)
        logo_options = weblogolib.LogoOptions()
        logo_options.fineprint = None
        logo_options.stacks_per_line = nperline
        logo_options.stack_aspect_ratio = aspectratio
        logo_options.stack_width = stackwidth
        logo_options.unit_name = 'probability'
        logo_options.show_yaxis = False
        logo_options.yaxis_scale = ymax 
        if alphabet_type == 'aa':
            (cmap, colormapping, mapper) = KyteDoolittleColorMapping()
        elif alphabet_type == 'nt':
            colormapping = {}
            colormapping['A'] = '#008000'
            colormapping['T'] = '#FF0000'
            colormapping['C'] = '#0000FF'
            colormapping['G'] = '#FFA500'
        else:
            raise ValueError("Invalid alphabet_type %s" % alphabet_type)
        colormapping[firstblankchar] = colormapping[lastblankchar] = '#000000' # black, but color doesn't matter as modified weblogo code replaces with empty space
        colormapping[separatorchar] = '#000000' # black
        color_scheme = weblogolib.colorscheme.ColorScheme()
        for x in chars_for_string:
            color_scheme.groups.append(weblogolib.colorscheme.ColorGroup(x, colormapping[x], "'%s'" % x))
        logo_options.color_scheme = color_scheme
        logo_options.annotate = [{True:r, False:''}[0 == isite % numberevery] for (isite, r) in enumerate(sites)]
        logoformat = weblogolib.LogoFormat(logodata, logo_options)
        # _my_pdf_formatter is modified from weblogo version 3.4 source code
        # to allow custom ordering of the symbols.
        pdf = _my_pdf_formatter(logodata, logoformat, ordered_alphabets) 
        open(plotfile, 'w').write(pdf)
        assert os.path.isfile(plotfile), "Failed to find expected plotfile %s" % plotfile
    finally:
        # close if still open
        try:
            f.close()
        except:
            pass
        # remove temporary file
        if os.path.isfile(transfacfile):
            os.remove(transfacfile)

    # now build the overlay
    if overlay:
        try:
            (fdoverlay, overlayfile) = tempfile.mkstemp(suffix='.pdf')
            (fdmerged, mergedfile) = tempfile.mkstemp(suffix='.pdf')
            foverlay = os.fdopen(fdoverlay, 'wb')
            foverlay.close() # close, but we still have the path overlayfile...
            fmerged = os.fdopen(fdmerged, 'wb')
            LogoOverlay(sites, overlayfile, overlay, nperline, sitewidth=stackwidth, rmargin=rmargin, logoheight=stackwidth * stackaspectratio + stackheightmargin, barheight=barheight, barspacing=barspacing)
            plotfile_f = open(plotfile, 'rb')
            plot = PyPDF2.PdfFileReader(plotfile_f).getPage(0)
            overlayfile_f = open(overlayfile, 'rb')
            overlay = PyPDF2.PdfFileReader(overlayfile_f).getPage(0)
            xshift = overlay.artBox[2] - plot.artBox[2]
            overlay.mergeTranslatedPage(plot, xshift, 0)
            overlay.compressContentStreams() 
            output = PyPDF2.PdfFileWriter()
            output.addPage(overlay)
            output.write(fmerged)
            fmerged.close()
            shutil.move(mergedfile, plotfile)
        finally:
            try:
                plotfile_f.close()
            except:
                pass
            try:
                overlayfile_f.close()
            except:
                pass
            try:
                foverlay.close()
            except:
                pass
            try:
                fmerged.close()
            except:
                pass
            for fname in [overlayfile, mergedfile]:
                if os.path.isfile(fname):
                    os.remove(fname)


#########################################################################
# The following code is modified from weblogo (version 3.4), which
# comes with the following license:
#
# -------------------------------- WebLogo --------------------------------

#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006-2011, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without 
#  modification, are permitted provided that the following conditions are met: 
#
#  (1) Redistributions of source code must retain the above copyright notice, 
#  this list of conditions and the following disclaimer. 
#
#  (2) Redistributions in binary form must reproduce the above copyright 
#  notice, this list of conditions and the following disclaimer in the 
#  documentation and or other materials provided with the distribution. 
#
#  (3) Neither the name of the University of California, Lawrence Berkeley 
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors 
#  may be used to endorse or promote products derived from this software 
#  without specific prior written permission. 
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#  POSSIBILITY OF SUCH DAMAGE. 

# Replicates README.txt

def _my_pdf_formatter(data, format, ordered_alphabets) :
    """ Generate a logo in PDF format.
    
    Modified from weblogo version 3.4 source code.
    """
    eps = _my_eps_formatter(data, format, ordered_alphabets).decode()
    gs = weblogolib.GhostscriptAPI()    
    return gs.convert('pdf', eps, format.logo_width, format.logo_height)


def _my_eps_formatter(logodata, format, ordered_alphabets) :
    """ Generate a logo in Encapsulated Postscript (EPS)
    
    Modified from weblogo version 3.4 source code. 

    *ordered_alphabets* is a dictionary keyed by zero-indexed
    consecutive sites, with values giving order of characters
    from bottom to top.
    """
    substitutions = {}
    from_format =[
        "creation_date",    "logo_width",           "logo_height",      
        "lines_per_logo",   "line_width",           "line_height",
        "line_margin_right","line_margin_left",     "line_margin_bottom",
        "line_margin_top",  "title_height",         "xaxis_label_height",
        "creator_text",     "logo_title",           "logo_margin",
        "stroke_width",     "tic_length",           
        "stacks_per_line",  "stack_margin",
        "yaxis_label",      "yaxis_tic_interval",   "yaxis_minor_tic_interval",
        "xaxis_label",      "xaxis_tic_interval",   "number_interval",
        "fineprint",        "shrink_fraction",      "errorbar_fraction",
        "errorbar_width_fraction",
        "errorbar_gray",    "small_fontsize",       "fontsize",
        "title_fontsize",   "number_fontsize",      "text_font",
        "logo_font",        "title_font",          
        "logo_label",       "yaxis_scale",          "end_type",
        "debug",            "show_title",           "show_xaxis",
        "show_xaxis_label", "show_yaxis",           "show_yaxis_label",
        "show_boxes",       "show_errorbars",       "show_fineprint",
        "rotate_numbers",   "show_ends",            "stack_height",
        "stack_width"
        ]
   
    for s in from_format :
        substitutions[s] = getattr(format,s)

    substitutions["shrink"] = str(format.show_boxes).lower()


    # --------- COLORS --------------
    def format_color(color):
        return  " ".join( ("[",str(color.red) , str(color.green), 
            str(color.blue), "]"))  

    substitutions["default_color"] = format_color(format.default_color)

    colors = []  
    for group in format.color_scheme.groups :
        cf = format_color(group.color)
        for s in group.symbols :
            colors.append( "  ("+s+") " + cf )
    substitutions["color_dict"] = "\n".join(colors)
        
    data = []
    
    # Unit conversion. 'None' for probability units
    conv_factor = None #JDB
    #JDB conv_factor = std_units[format.unit_name]
    
    data.append("StartLine")

    seq_from = format.logo_start- format.first_index
    seq_to = format.logo_end - format.first_index +1

    # seq_index : zero based index into sequence data
    # logo_index : User visible coordinate, first_index based
    # stack_index : zero based index of visible stacks
    for seq_index in range(seq_from, seq_to) :
        logo_index = seq_index + format.first_index 
        stack_index = seq_index - seq_from
        
        if stack_index!=0 and (stack_index % format.stacks_per_line) ==0 :
            data.append("")
            data.append("EndLine")
            data.append("StartLine")
            data.append("")
        
        data.append("(%s) StartStack" % format.annotate[seq_index] )

        if conv_factor: 
            stack_height = logodata.entropy[seq_index] * std_units[format.unit_name]
        else :
            stack_height = 1.0 # Probability

        # The following code modified by JDB to use ordered_alphabets
        # and also to replace the "blank" characters 'b' and 'B'
        # by spaces.
        s_d = dict(zip(logodata.alphabet, logodata.counts[seq_index]))
        s = []
        for aa in ordered_alphabets[seq_index]:
            if aa not in ['B', 'b']:
                s.append((s_d[aa], aa))
            else:
                s.append((s_d[aa], ' '))
#        s = [(s_d[aa], aa) for aa in ordered_alphabets[seq_index]]

        # Sort by frequency. If equal frequency then reverse alphabetic
        # (So sort reverse alphabetic first, then frequencty)
        # TODO: doublecheck this actual works
        #s = list(zip(logodata.counts[seq_index], logodata.alphabet))
        #s.sort(key= lambda x: x[1])
        #s.reverse()
        #s.sort(key= lambda x: x[0])
        #if not format.reverse_stacks: s.reverse()

        C = float(sum(logodata.counts[seq_index])) 
        if C > 0.0 :
            fraction_width = 1.0
            if format.scale_width :
                fraction_width = logodata.weight[seq_index] 
            # print(fraction_width, file=sys.stderr)
            for c in s:
                data.append(" %f %f (%s) ShowSymbol" % (fraction_width, c[0]*stack_height/C, c[1]) )

        # Draw error bar on top of logo. Replaced by DrawErrorbarFirst above.
        if logodata.entropy_interval is not None and conv_factor and C>0.0:

            low, high = logodata.entropy_interval[seq_index]
            center = logodata.entropy[seq_index]
            low *= conv_factor
            high *= conv_factor
            center *=conv_factor
            if high> format.yaxis_scale : high = format.yaxis_scale 

            down = (center - low) 
            up   = (high - center) 
            data.append(" %f %f DrawErrorbar" % (down, up) )
            
        data.append("EndStack")
        data.append("")
               
    data.append("EndLine")
    substitutions["logo_data"] = "\n".join(data)  


    # Create and output logo
    template = corebio.utils.resource_string( __name__, '_weblogo_template.eps', __file__).decode()
    logo = string.Template(template).substitute(substitutions)

    return logo.encode()
#
# End of code modified from weblogo    
#########################################################################

#########################################################################
# More code modified from weblogo version 3.4 by Jesse Bloom to allow non-
# alphabetic characters in motifs.

#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006 John Gilman

#  This software is distributed under the MIT Open Source License.
#  <http://www.opensource.org/licenses/mit-license.html>
#
#  Permission is hereby granted, free of charge, to any person obtaining a 
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
#  IN THE SOFTWARE.


class _my_Motif(corebio.matrix.AlphabeticArray) :
    """A two dimensional array where the second dimension is indexed by an 
    Alphabet. Used to represent sequence motifs and similar information.

    
    Attr:
    - alphabet     -- An Alphabet
    - array        -- A numpy array
    - name         -- The name of this motif (if any) as a string.
    - description  -- The description, if any.
    
    """
            
    def __init__(self, alphabet, array=None, dtype=None, name=None,
            description = None, scale=None) :
        corebio.matrix.AlphabeticArray.__init__(self, (None, alphabet), array, dtype)
        self.name = name
        self.description = description
        self.scale = scale    

    @property
    def alphabet(self):
        return self.alphabets[1]
        
    def reindex(self, alphabet) :
        return  _my_Motif(alphabet, corebio.matrix.AlphabeticArray.reindex(self, (None, alphabet)))
  
    # These methods alter self, and therefore do not return a value.
    # (Compare to Seq objects, where the data is immutable and therefore methods return a new Seq.)
    # TODO: Should reindex (above) also act on self?
    
    def reverse(self):
        """Reverse sequence data"""
#        self.array = na.array(self.array[::-1]) # This is a view into the origional numpy array.
        self.array = self.array[::-1] # This is a view into the origional numpy array.
    
    @staticmethod #TODO: should be classmethod?
    def read_transfac( fin, alphabet = None) :
        """ Parse a sequence matrix from a file. 
        Returns a tuple of (alphabet, matrix)
        """
   
        items = []

        start=True
        for line in fin :
            if line.isspace() or line[0] =='#' : continue
            stuff = line.split()
            if start and stuff[0] != 'PO' and stuff[0] != 'P0': continue
            if stuff[0]=='XX' or stuff[0]=='//': break
            start = False
            items.append(stuff)
        if len(items) < 2  :
            raise ValueError("Vacuous file.")

        # Is the first line a header line?
        header = items.pop(0)
        hcols = len(header)
        rows = len(items)
        cols = len(items[0])
        if not( header[0] == 'PO' or header[0] =='P0' or hcols == cols-1 or hcols == cols-2) :
            raise ValueError("Missing header line!")

        # Do all lines (except the first) contain the same number of items?
        cols = len(items[0])
        for i in range(1, len(items)) :
            if cols != len(items[i]) :
                raise ValueError("Inconsistant length, row %d: " % i)

        # Vertical or horizontal arrangement?
        if header[0] == 'PO' or header[0] == 'P0': header.pop(0)

        position_header = True    
        alphabet_header = True    
        for h in header :
            if not corebio.utils.isint(h) : position_header = False
#allow non-alphabetic            if not str.isalpha(h) : alphabet_header = False

        if not position_header and not alphabet_header :
            raise ValueError("Can't parse header: %s" % str(header))

        if position_header and alphabet_header :
            raise ValueError("Can't parse header")


        # Check row headers
        if alphabet_header :
            for i,r in enumerate(items) :
                if not corebio.utils.isint(r[0]) and r[0][0]!='P' : 
                    raise ValueError(
                        "Expected position as first item on line %d" % i)
                r.pop(0)
                defacto_alphabet = ''.join(header)
        else :
            a = []
            for i,r in enumerate(items) :
                if not ischar(r[0]) and r[0][0]!='P' : 
                    raise ValueError(
                        "Expected position as first item on line %d" % i)
                a.append(r.pop(0))
            defacto_alphabet = ''.join(a)                

        # Check defacto_alphabet
        defacto_alphabet = corebio.seq.Alphabet(defacto_alphabet)

        if alphabet :
            if not defacto_alphabet.alphabetic(alphabet) :
                raise ValueError("Incompatible alphabets: %s , %s (defacto)"
                                 % (alphabet, defacto_alphabet))
        else :            
            alphabets = (unambiguous_rna_alphabet,
                        unambiguous_dna_alphabet,                      
                        unambiguous_protein_alphabet,
                      )
            for a in alphabets :
                if defacto_alphabet.alphabetic(a) :
                    alphabet = a
                    break
            if not alphabet :
                alphabet = defacto_alphabet
   

        # The last item of each row may be extra cruft. Remove
        if len(items[0]) == len(header) +1 :
            for r in items :
                r.pop()

        # items should now be a list of lists of numbers (as strings) 
        rows = len(items)
        cols = len(items[0])
        matrix = numpy.zeros( (rows,cols) , dtype=numpy.float64) 
        for r in range( rows) :
            for c in range(cols):
                matrix[r,c] = float( items[r][c]) 

        if position_header :
            matrix.transpose() 

        return _my_Motif(defacto_alphabet, matrix).reindex(alphabet)

# End of code modified from weblogo version 3.4
#==============================================================


def LogoOverlay(sites, overlayfile, overlay, nperline, sitewidth, rmargin, logoheight, barheight, barspacing):
    """Makes overlay for *LogoPlot*.

    This function creates colored bars overlay bars showing up to two
    properties.
    The trick of this function is to create the bars the right
    size so they align when they overlay the logo plot. 

    CALLING VARIABLES:

    * *sites* : same as the variable of this name used by *LogoPlot*.

    * *overlayfile* is a string giving the name of created PDF file containing
      the overlay. It must end in the extension ``.pdf``.

    * *overlay* : same as the variable of this name used by *LogoPlot*.

    * *nperline* : same as the variable of this name used by *LogoPlot*.

    * *sitewidth* is the width of each site in points.

    * *rmargin* is the right margin in points.

    * *logoheight* is the total height of each logo row in points.

    * *barheight* is the total height of each bar in points.

    * *barspacing* is the vertical spacing between bars in points.

    * *cmap* is a ``pylab`` *LinearSegmentedColorMap* used for the bar coloring.
    """
    if not pylab.get_backend().lower() == 'pdf':
        raise ValueError("You cannot use this function without first setting the matplotlib / pylab backend to 'pdf'. Do this with: matplotlib.use('pdf')")
    if os.path.splitext(overlayfile)[1] != '.pdf':
        raise ValueError("overlayfile must end in .pdf: %s" % overlayfile)
    (cmap, mapping_d, mapper) = KyteDoolittleColorMapping()
    pts_per_inch = 72.0 # to convert between points and inches
    # some general properties of the plot
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('xtick', labelsize=9)
    matplotlib.rc('xtick', direction='out')
    matplotlib.rc('ytick', direction='out')
    matplotlib.rc('axes', linewidth=0.5)
    matplotlib.rc('ytick.major', size=3)
    matplotlib.rc('xtick.major', size=2.5)
    # define sizes (still in points)
    colorbar_bmargin = 20 # margin below color bars in points
    colorbar_tmargin = 15 # margin above color bars in points
    nlines = int(math.ceil(len(sites) / float(nperline)))
    lmargin = 25 # left margin in points
    barwidth = nperline * sitewidth
    figwidth = lmargin + rmargin + barwidth
    figheight = nlines * (logoheight + len(overlay) * (barheight + barspacing)) + barheight + colorbar_bmargin + colorbar_tmargin
    # set up the figure and axes
    fig = pylab.figure(figsize=(figwidth / pts_per_inch, figheight / pts_per_inch))
    # determine property types
    prop_types = {}
    for (prop_d, shortname, longname) in overlay:
        if all([isinstance(prop, str) for prop in prop_d.itervalues()]):
            proptype = 'discrete'
            propcategories = list(set(prop_d.values()))
            propcategories.sort()
            (vmin, vmax) = (0, len(propcategories) - 1)
        elif all ([isinstance(prop, (int, float)) for prop in prop_d.itervalues()]):
            proptype = 'continuous'
            propcategories = None
            (vmin, vmax) = (min(prop_d.values()), max(prop_d.values()))
            # If vmin is slightly greater than zero, set it to zero. This helps for RSA properties.
            if vmin >= 0 and vmin / float(vmax - vmin) < 0.05:
                vmin = 0.0
                # And if vmax is just a bit less than one, set it to that...
                if 0.9 <= vmax <= 1.0:
                    vmax = 1.0
        else:
            raise ValueError("Property %s is neither continuous or discrete. Values are:\n%s" % (shortname, str(prop_d.items())))
        prop_types[shortname] = (proptype, vmin, vmax, propcategories)
    assert len(prop_types) == len(overlay), "Not as many property types as overlays. Did you give the same name (shortname) to multiple properties in the overlay?"
    # loop over each line of the multi-lined plot
    prop_image = {}
    for iline in range(nlines):
        isites = sites[iline * nperline : min(len(sites), (iline + 1) * nperline)]
        xlength = len(isites) * sitewidth
        logo_ax = pylab.axes([lmargin / figwidth, ((nlines - iline - 1) * (logoheight + len(overlay) * (barspacing + barheight))) / figheight, xlength / figwidth, logoheight / figheight], frameon=False)
        logo_ax.yaxis.set_ticks_position('none')
        logo_ax.xaxis.set_ticks_position('none')
        pylab.yticks([])
        pylab.xlim(0.5, len(isites) + 0.5)
        pylab.xticks([])
        for (iprop, (prop_d, shortname, longname)) in enumerate(overlay):
            (proptype, vmin, vmax, propcategories) = prop_types[shortname]
            prop_ax = pylab.axes([lmargin / figwidth, ((nlines - iline - 1) * (logoheight + len(overlay) * (barspacing + barheight)) + logoheight + iprop * (barspacing + barheight)) / figheight, xlength / figwidth, barheight / figheight], frameon=True)
            prop_ax.xaxis.set_ticks_position('none')
            prop_ax.yaxis.set_ticks_position('left')
            pylab.xticks([])
            pylab.yticks([0], [shortname], size=8)
            pylab.xlim((0, len(isites)))
            pylab.ylim(-0.5, 0.5)
            propdata = pylab.zeros(shape=(1, len(isites)))
            propdata[ : ] = pylab.nan # set to nan for all entries
            for (isite, site) in enumerate(isites):
                if site in prop_d:
                    if proptype == 'continuous':
                        propdata[(0, isite)] = prop_d[site]
                    elif proptype == 'discrete':
                        propdata[(0, isite)] = propcategories.index(prop_d[site])
                    else:
                        raise ValueError('neither continuous nor discrete')
            prop_image[shortname] = pylab.imshow(propdata, interpolation='nearest', aspect='auto', extent=[0, len(isites), 0.5, -0.5], cmap=cmap, vmin=vmin, vmax=vmax)
            pylab.yticks([0], [shortname], size=8)
    # set up colorbar axes, then color bars
    if len(overlay) == 1:
        colorbarwidth = 0.4
        colorbarspacingwidth = 1.0 - colorbarwidth
    else:
        colorbarspacingfrac = 0.5 # space between color bars is this fraction of bar width
        colorbarwidth = 1.0 / (len(overlay) * (1.0 + colorbarspacingfrac)) # width of color bars in fraction of figure width
        colorbarspacingwidth = colorbarwidth * colorbarspacingfrac # width of color bar spacing in fraction of figure width
    propnames = {}
    for (icolorbar, (prop_d, shortname, longname)) in enumerate(overlay):
        (proptype, vmin, vmax, propcategories) = prop_types[shortname]
        if shortname == longname or not longname:
            propame = shortname
        else:
            propname = "%s (%s)" % (longname, shortname)
#        colorbar_ax = pylab.axes([colorbarspacingwidth * 0.5 + icolorbar * (colorbarwidth + colorbarspacingwidth), 1.0 - (colorbar_tmargin + barheight) / figwidth, colorbarwidth, barheight / figwidth], frameon=True)
        colorbar_ax = pylab.axes([colorbarspacingwidth * 0.5 + icolorbar * (colorbarwidth + colorbarspacingwidth), 1.0 - (colorbar_tmargin + barheight) / figheight, colorbarwidth, barheight / figheight], frameon=True)
        colorbar_ax.xaxis.set_ticks_position('bottom')
        colorbar_ax.yaxis.set_ticks_position('none')
        pylab.xticks([])
        pylab.yticks([])
        pylab.title(propname, size=9)
        if proptype == 'continuous':
            cb = pylab.colorbar(prop_image[shortname], cax=colorbar_ax, orientation='horizontal')
            # if range is close to zero to one, manually set tics to 0, 0.5, 1. This helps for RSA
            if -0.1 <= vmin <= 0 and 1.0 <= vmax <= 1.15:
                cb.set_ticks([0, 0.5, 1])
                cb.set_ticklabels(['0', '0.5', '1'])
        elif proptype == 'discrete':
            cb = pylab.colorbar(prop_image[shortname], cax=colorbar_ax, orientation='horizontal', boundaries=[i for i in range(len(propcategories) + 1)], values=[i for i in range(len(propcategories))])
            cb.set_ticks([i + 0.5 for i in range(len(propcategories))])
            cb.set_ticklabels(propcategories)
        else:
            raise ValueError("Invalid proptype")
    # save the plot
    pylab.savefig(overlayfile, transparent=True)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
