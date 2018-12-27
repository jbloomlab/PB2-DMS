import re
import os
import math

import natsort
import pandas
from pandas.api.types import CategoricalDtype
import numpy
import scipy.stats
import scipy.optimize
from statsmodels.sandbox.stats.multicomp import multipletests

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

from plotnine import *
# set ggplot theme
theme_set(theme_bw(base_size=12)) 

import seaborn
seaborn.set(context='talk',
            style='white',
            rc={
                'xtick.labelsize':15,
                'ytick.labelsize':15,
                'axes.labelsize':19,
                'font.family':'sans-serif',
                'font.sans-serif':['DejaVu Sans'],
                }
           )

from dms_tools2 import CODONS, AAS, AAS_WITHSTOP
import dms_tools2.utils

#: `color-blind safe palette <http://bconnelly.net/2013/10/creating-colorblind-friendly-figures/>`_
#: use by adding to your plots the following
#: `scale_fill_manual(COLOR_BLIND_PALETTE)` or
#: `scale_color_manual(COLOR_BLIND_PALETTE)`.
COLOR_BLIND_PALETTE = ["#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]


def latexSciNot(xlist):
    """Converts list of numbers to LaTex scientific notation.
    Useful for nice axis-tick formatting.
    Args:
        `xlist` (list)
            Numbers to format.
    Returns:
        List of latex scientific notation formatted strings.
    >>> latexSciNot([0, 3, 3120, -0.0000927])
    ['$0$', '$3$', '$3.1 \\\\times 10^{3}$', '$-9.3 \\\\times 10^{-5}$']
    >>> latexSciNot([0.001, 1, 1000, 1e6])
    ['$0.001$', '$1$', '$10^{3}$', '$10^{6}$']
    >>> latexSciNot([-0.002, 0.003, 0.000011])
    ['$-0.002$', '$0.003$', '$1.1 \\\\times 10^{-5}$']
    >>> latexSciNot([-0.1, 0.0, 0.1, 0.2])
    ['$-0.1$', '$0$', '$0.1$', '$0.2$']
    >>> latexSciNot([0, 1, 2])
    ['$0$', '$1$', '$2$']
    """
    formatlist = []
    for x in xlist:
        xf = "{0:.2g}".format(x)
        if xf[ : 2] == '1e':
            xf = "$10^{{{0}}}$".format(int(xf[2 : ]))
        elif xf[ : 3] == '-1e':
            xf = "$-10^{{{0}}}$".format(int(xf[3 : ]))
        elif 'e' in xf:
            (d, exp) = xf.split('e')
            xf = '${0} \\times 10^{{{1}}}$'.format(d, int(exp))
        else:
            xf = '${0}$'.format(xf)
        formatlist.append(xf)
    return formatlist

def annotateCodonCounts_accessible(counts, accessible):
    """Gets annotated `pandas.DataFrame` from codon counts.

    Some of the programs (e.g., `dms2_bcsubamplicons`) create 
    ``*_codoncounts.csv`` files when run with ``--chartype codon``.
    These CSV files have columns indicating the `site` and `wildtype`
    codon, as well as a column for each codon giving the counts for that 
    codon. This function reads that file (or a `pandas.DataFrame` read
    from it) to return a `pandas.DataFrame` where a variety of additional
    useful annotations have been added.

    Args:
        `counts` (str)
            Name of existing codon counts CSV file, or `pandas.DataFrame`
            holding counts.
        `accessible` (True/False)
            True : accessible by 1-nucleotide mutation
            False: not-accessible by 1-nucleotide mutation

    Returns:
        `df` (`pandas.DataFrame`)
            The DataFrame with the information in `counts` plus
            the following added columns for each site:

                `ncounts` : number of counts at site

                `mutfreq` : mutation frequency at site

                `nstop` : number of stop-codon mutations

                `nsyn` : number of synonymous mutations

                `nnonsyn` : number of nonsynonymous mutations

                `n1nt` : number of 1-nucleotide codon mutations

                `n2nt` : number of 2-nucleotide codon mutations
                
                `n3nt` : number of 3-nucleotide codon mutations

                `AtoC`, `AtoG`, etc : number of each nucleotide mutation
                type among codon mutations with **one** nucleotide change.

                `mutfreq1nt`, `mutfreq2nt`, `mutfreq3nt` : frequency
                of 1-, 2-, and 3-nucleotide codon mutations at site.

    """
    print('Annotating '+counts)
    if isinstance(counts, str):
        df = pandas.read_csv(counts)
    elif isinstance(counts, pandas.DataFrame):
        df = counts.copy()
    else:
        raise ValueError("invalid counts")
    assert set(dms_tools2.CODONS) <= set(df.columns), \
            "Did not find counts for all codons".format(counts)

    ## Added 20180424 SS to subset to accessible or non-accessible codons
    acccodons = {}
    for codon in dms_tools2.CODONS:
        acccodonlist = []
        for i, orignt in enumerate(codon):
            acccodonlist += [(codon[:i] + n + codon[i+1:]) 
                          for n in dms_tools2.NTS if n != orignt]
        acccodons[codon] = acccodonlist

    for (i, row) in df.iterrows():
        wt = row['wildtype']
        if accessible:
            for newcodon in dms_tools2.CODONS:
                if newcodon not in (acccodons[wt] + [wt]):
                    df.loc[df.index[i], newcodon] = 0
        else:
            for acccodon in acccodons[wt]:
                df.loc[df.index[i], acccodon] = 0
    ###
        
    df['ncounts'] = df[dms_tools2.CODONS].sum(axis=1)

    df['mutfreq'] = (((df['ncounts'] - df.lookup(df['wildtype'].index,
            df['wildtype'].values)) / df['ncounts'].astype('float'))
            .fillna(0))

    ntchanges = ['{0}to{1}'.format(nt1, nt2) for nt1 in dms_tools2.NTS
            for nt2 in dms_tools2.NTS if nt1 != nt2]

    nstoplist = []
    nsynlist = []
    nnonsynlist = []
    nXntlists = dict([(n + 1, []) for n in range(3)])
    nntchangeslists = dict([(ntchange, []) for ntchange in ntchanges])
    
    dfAAS = pandas.DataFrame(0, index=df.index, columns=(dms_tools2.AAS_WITHSTOP)) ### SS 20180424
    
    for (i, row) in df.iterrows():
        nstop = nsyn = nnonsyn = 0
        nXnt = dict([(n + 1, 0) for n in range(3)])
        nntchanges = dict([(ntchange, 0) for ntchange in ntchanges])
        wt = row['wildtype']
        wtaa = dms_tools2.CODON_TO_AA[wt]
        for c in dms_tools2.CODONS:
            if c == wt:
                continue
            aa = dms_tools2.CODON_TO_AA[c]
            dfAAS.loc[dfAAS.index[i], aa] += row[c] ### SS 20180424
            if aa == '*':
                nstop += row[c]
            elif aa == wtaa:
                nsyn += row[c]
            else:
                nnonsyn += row[c]
            ntdiffs = ['{0}to{1}'.format(nt1, nt2) for (nt1, nt2) 
                    in zip(wt, c) if nt1 != nt2]
            nXnt[len(ntdiffs)] += row[c]
            if len(ntdiffs) == 1:
                nntchanges[ntdiffs[0]] += row[c]
        nstoplist.append(nstop)
        nsynlist.append(nsyn)
        nnonsynlist.append(nnonsyn)
        for n in range(3):
            nXntlists[n + 1].append(nXnt[n + 1])
        for ntchange in ntchanges:
            nntchangeslists[ntchange].append(nntchanges[ntchange])
    df = df.assign(nstop=nstoplist, nsyn=nsynlist, nnonsyn=nnonsynlist)
    df = df.assign(n1nt=nXntlists[1], n2nt=nXntlists[2], n3nt=nXntlists[3])
    for ntchange in ntchanges:
        df[ntchange] = nntchangeslists[ntchange]

    for nnt in range(3):
        df['mutfreq{0}nt'.format(nnt + 1)] = (df['n{0}nt'.format(nnt + 1)]
            / df['ncounts'].astype('float')).fillna(0)
    
    df = pandas.concat([df, dfAAS], axis=1) ### SS 20180424

    return df

def plotCodonMutTypes_accessible(names, countsfiles, plotfile, accessible,
        classification='aachange', csvfile=None):
    """Plot average frequency codon mutation types.

    The averages are determined by summing counts for all sites.

    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `countsfiles` (list or series)
            ``*_codoncounts.csv`` files of type created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of created PDF plot file.
        `accessible` (True/False)
            True : accessible by 1-nucleotide mutation
            False: not-accessible by 1-nucleotide mutation
        `classification` (str)
            The method used to classify the mutation types. Can be:

                `aachange` : stop, synonymous, nonsynonymous

                `n_ntchanges` : number of nucleotide changes per codon

                `singlentchanges` : nucleotide change in 1-nt mutations
                
                `aaresult` : resultant amino acid after mutation

        `csvfile` (str or `None`)
            `None` or name of CSV file to which numerical data are written.
    """
    assert len(names) == len(countsfiles) == len(set(names))
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    counts = pandas.concat([annotateCodonCounts_accessible(f, accessible).assign(
            name=name) for (name, f) in zip(names, countsfiles)], 
            ignore_index=True)

    if classification == 'aachange':
        muttypes = {'stop':'nstop', 'synonymous':'nsyn', 
                'nonsynonymous':'nnonsyn'}
    elif classification == 'n_ntchanges':
        muttypes = dict([('{0} nucleotide'.format(n), 'n{0}nt'.format(n))
                for n in [1, 2, 3]])
    elif classification == 'singlentchanges':
        muttypes = dict([(ntchange, ntchange) for ntchange in [
                '{0}to{1}'.format(nt1, nt2) for nt1 in dms_tools2.NTS
                for nt2 in dms_tools2.NTS if nt1 != nt2]])
    ### SS 20180424
    elif classification == 'aaresult':
        muttypes = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                    '*Stop': '*'
                   }
    ###
    else:
        raise ValueError("Invalid classification {0}".format(classification))

    df = (counts[list(muttypes.values()) + ['ncounts', 'name']]
            .groupby('name', as_index=False)
            .sum(axis=1)
            .assign(ncounts=lambda x: x['ncounts'].astype('float'))
            )
    for (newcol, n) in muttypes.items():
        df[newcol] = (df[n] / df['ncounts']).fillna(0)

    if csvfile:
        df[['name'] + list(muttypes.keys())].to_csv(csvfile, index=False) 

    df = df.melt(id_vars='name', var_name='mutation type',
                value_vars=list(muttypes.keys()),
                value_name='per-codon frequency')

    ## SS 20180313: Make names a category so as to order bar plots in input order
    names_ordered = CategoricalDtype(categories=names, ordered=True)
    df['name_ordered'] = df['name'].astype(str).astype(names_ordered)
    print('Ordering names again!')
	##
    
    p = (ggplot(df)
            + geom_col(aes(x='name_ordered', y='per-codon frequency', ##SS 20180313: Replaced 'name' with 'name_ordered'
                fill='mutation type'), position='stack')
            + theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                    axis_title_x=element_blank())
            + scale_y_continuous(labels=latexSciNot)
            )
    if len(muttypes) <= len(COLOR_BLIND_PALETTE):
        p += scale_fill_manual(COLOR_BLIND_PALETTE)
    else:
        p += guides(fill=guide_legend(ncol=2))

    p.save(plotfile, height=2.7, width=(1.2 + 0.25 * len(names)),
            verbose=False)
    plt.close()

def plotMutFreq_customY(names, countsfiles, plotfile, ymax, maxcol=4):
    """Plot mutation frequency along primary sequence.
    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `countsfiles` (list or series)
            ``*_codoncounts.csv`` files of type created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of created PDF plot file.
        `maxcol` (int)
            Number of columns in faceted plot.
    SS 20180312: Added ymax variable
    """
    assert len(names) == len(countsfiles) == len(set(names))
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    counts = pandas.concat([dms_tools2.utils.annotateCodonCounts(f).assign(
            name=name) for (name, f) in zip(names, countsfiles)], 
            ignore_index=True).rename(columns={'mutfreq':'mutation frequency'})
    
    ncol = min(maxcol, len(names))
    nrow = math.ceil(len(names) / float(ncol))

    # make name a category to preserve order
    counts['name'] = counts['name'].astype('category', 
            categories=names)

    p = (ggplot(counts, aes(x='site', y='mutation frequency'))
            + geom_step(size=0.4)
#             + scale_y_continuous(labels=latexSciNot, 
#                     limits=(0, ymax)) ##SS 20180312: Replaced ylim with ymax variable
            + coord_cartesian(ylim=(0, ymax))
            + scale_x_continuous(limits=(counts['site'].min(),
                    counts['site'].max()))
            + facet_wrap('~name', ncol=ncol) 
            + theme(figure_size=(2.25 * (0.6 + ncol), 1.3 * (0.3 + nrow)))
            )
    p.save(plotfile, verbose=False)
    plt.close()
    
def plotCodonMutTypes_ordered(names, countsfiles, plotfile,
        classification='aachange', csvfile=None):
    """Plot average frequency codon mutation types.
    The averages are determined by summing counts for all sites.
    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `countsfiles` (list or series)
            ``*_codoncounts.csv`` files of type created by ``dms2_bcsubamp``.
        `plotfile` (str)
            Name of created PDF plot file.
        `classification` (str)
            The method used to classify the mutation types. Can be:
                `aachange` : stop, synonymous, nonsynonymous
                `n_ntchanges` : number of nucleotide changes per codon
                `singlentchanges` : nucleotide change in 1-nt mutations
        `csvfile` (str or `None`)
            `None` or name of CSV file to which numerical data are written.
    """
    assert len(names) == len(countsfiles) == len(set(names))
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    counts = pandas.concat([dms_tools2.utils.annotateCodonCounts(f).assign(
            name=name) for (name, f) in zip(names, countsfiles)], 
            ignore_index=True)

    if classification == 'aachange':
        muttypes = {'stop':'nstop', 'synonymous':'nsyn', 
                'nonsynonymous':'nnonsyn'}
    elif classification == 'n_ntchanges':
        muttypes = dict([('{0} nucleotide'.format(n), 'n{0}nt'.format(n))
                for n in [1, 2, 3]])
    elif classification == 'singlentchanges':
        muttypes = dict([(ntchange, ntchange) for ntchange in [
                '{0}to{1}'.format(nt1, nt2) for nt1 in dms_tools2.NTS
                for nt2 in dms_tools2.NTS if nt1 != nt2]])
    else:
        raise ValueError("Invalid classification {0}".format(classification))

    df = (counts[list(muttypes.values()) + ['ncounts', 'name']]
            .groupby('name', as_index=False)
            .sum(axis=1)
            .assign(ncounts=lambda x: x['ncounts'].astype('float'))
            )
    for (newcol, n) in muttypes.items():
        df[newcol] = (df[n] / df['ncounts']).fillna(0)

    if csvfile:
        df[['name'] + list(muttypes.keys())].to_csv(csvfile, index=False) 

    df = df.melt(id_vars='name', var_name='mutation type',
                value_vars=list(muttypes.keys()),
                value_name='per-codon frequency')
    
    ## SS 20180313: Make names a category so as to order bar plots in input order
    names_ordered = CategoricalDtype(categories=names, ordered=True)
    df['name_ordered'] = df['name'].astype(str).astype(names_ordered)
	##
	
    p = (ggplot(df)
            + geom_col(aes(x='name_ordered', y='per-codon frequency', ##SS 20180313: Replaced 'name' with 'name_ordered'
                fill='mutation type'), position='stack')
            + theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                    axis_title_x=element_blank())
            + scale_y_continuous(labels=latexSciNot)
            )
    if len(muttypes) <= len(COLOR_BLIND_PALETTE):
        p += scale_fill_manual(COLOR_BLIND_PALETTE)
    else:
        p += guides(fill=guide_legend(ncol=2))

    p.save(plotfile, height=2.7, width=(1.2 + 0.25 * len(names)),
            verbose=False)
    plt.close()

def plotCodonMutTypes_fromTable(df, muttypes, names, plotfile):
    """Plot average frequency codon mutation types.
    The averages are determined by summing counts for all sites.
    Args:
        `df` (codonmuttypes dataframe)
        `muttypes` (list)
        `names` (list, to order in)
        `plotfile` (str)
            Name of created PDF plot file.
    """

    df = df.melt(id_vars='name', var_name='mutation type',
                value_vars=muttypes,
                value_name='per-codon frequency')
    
    ## SS 20180313: Make names a category so as to order bar plots in input order
    names_ordered = CategoricalDtype(categories=names, ordered=True)
    df['name_ordered'] = df['name'].astype(str).astype(names_ordered)
	##
	
    p = (ggplot(df)
            + geom_col(aes(x='name_ordered', y='per-codon frequency', ##SS 20180313: Replaced 'name' with 'name_ordered'
                fill='mutation type'), position='stack')
            + theme(axis_text_x=element_text(angle=90, vjust=1, hjust=0.5),
                    axis_title_x=element_blank())
            + scale_y_continuous(labels=latexSciNot)
            )
    if len(muttypes) <= len(COLOR_BLIND_PALETTE):
        p += scale_fill_manual(COLOR_BLIND_PALETTE)
    else:
        p += guides(fill=guide_legend(ncol=2))

    p.save(plotfile, height=2.7, width=(1.2 + 0.25 * len(names)),
            verbose=False)
    plt.close()


def plotCorrMatrix2(names, infiles, plotfile, datatype,
        trim_unshared=True, title='', colors='black',
        contour=False, ncontours=10):
    """Plots correlations among replicates.

    Args:
        `names` (list or series)
            Names of samples for which we plot statistics.
        `infiles` (list or series)
            CSV files containing data. Format depends on `datatype`.
        `plotfile` (str)
            Name of created PDF plot file.
        `datatype` (str)
            Type of data for which we are plotting correlations:
                - `logprefs`: in format returned by ``dms2_prefs``, to be log transformed
                - `effect`: effect from ``dms_tools2.prefs.prefsToMutEffects``
                - `log2effect`: log2effect from ``dms_tools2.prefs.prefsToMutEffects``
        `trim_unshared` (bool)
            What if files in `infiles` don't have same sites / mutations?
            If `True`, trim unshared one and just analyze ones
            shared among all files. If `False`, raise an error.
        `title` (str)
            Title to place above plot.
        `colors` (str or list)
            Color(s) to color scatter points. If a string, should
            specify one color for all plots. Otherwise should be
            list of length `len(names) * (len(names) - 1) // 2`
            giving lists of colors for plots from top to bottom, 
            left to right.
        `contour` (bool)
            Show contour lines from KDE rather than points.
        `ncontours` (int)
            Number of contour lines if using `contour`.
    """
    assert len(names) == len(infiles) == len(set(names)) > 1
    assert os.path.splitext(plotfile)[1].lower() == '.pdf'

    if datatype == 'logprefs':
        # read prefs into dataframe, ensuring all have same characters
        prefs = [pandas.read_csv(f).assign(name=name) for (name, f) 
                in zip(names, infiles)]
        chars = set(prefs[0].columns)
        sites = set(prefs[0]['site'].values)
        for p in prefs:
            unsharedchars = chars.symmetric_difference(set(p.columns))
            if unsharedchars:
                raise ValueError("infiles don't have same characters: {0}"
                        .format(unsharedchars))
            unsharedsites = sites.symmetric_difference(set(p['site']))
            if trim_unshared:
                sites -= unsharedsites
            elif unsharedsites:
                raise ValueError("infiles don't have same sites: {0}".
                        format(unsharedsites))
        assert {'site', 'name'} < chars
        chars = list(chars - {'site', 'name'})
        # get measurements for each replicate in its own column
        df = (pandas.concat(prefs, ignore_index=True)
                    .query('site in @sites') # only keep shared sites
                    .melt(id_vars=['name', 'site'], var_name='char')
                    .pivot_table(index=['site', 'char'], columns='name')
                    )
        df.columns = df.columns.get_level_values(1)
        df = df.applymap(numpy.log2) ##SS 20180319: log2 preference values
    
    elif datatype in ['effect', 'log2effect']:
        mut_effect_df = [pandas.read_csv(f)
                         .assign(name=name)
                         .sort_values('mutation')
                         [['name', 'mutation', datatype]]
                    for (name, f) in zip(names, infiles)]
        muts = set(mut_effect_df[0]['mutation'].values)
        for m in mut_effect_df:
            unsharedmuts = muts.symmetric_difference(set(m['mutation']))
            if trim_unshared:
                muts -= unsharedmuts
            elif unsharedmuts:
                raise ValueError("infiles don't have same muts: {0}".
                        format(unsharedmuts))
        df = (pandas.concat(mut_effect_df, ignore_index=True)
                    .query('mutation in @muts') # only keep shared muts
                    .pivot_table(index='mutation', columns='name')
                    .dropna()
                    )
        df.columns = df.columns.get_level_values(1)

    else:
        raise ValueError("Invalid datatype {0}".format(datatype))

    ncolors = len(names) * (len(names) - 1) // 2
    if isinstance(colors, str):
        colors = [colors] * ncolors
    else:
        assert len(colors) == ncolors, "not {0} colors".format(ncolors)

    def corrfunc(x, y, contour, **kws):
        r, _ = scipy.stats.pearsonr(x, y)
        ax = plt.gca()
        ax.annotate('R = {0:.2f}'.format(r), xy=(0.05, 0.9), 
                xycoords=ax.transAxes, fontsize=19, 
                fontstyle='oblique')
        color = colors.pop(0)
        if contour:
            seaborn.kdeplot(x, y, shade=True, n_levels=ncontours)
        else:
#            plt.scatter(x, y, s=22, alpha=0.25, color=color,
#                    marker='o', edgecolor='none', rasterized=True)
            plt.scatter(x, y, s=10, alpha=0.5, marker='o', facecolors='none', edgecolors=color, rasterized=True)

    # map lower / upper / diagonal as here:
    # https://stackoverflow.com/a/30942817 for plot
    seaborn.set_style("ticks")
    p = seaborn.PairGrid(df)
    p.map_lower(corrfunc, colors=colors, contour=contour)
    if datatype == 'logprefs':
        (lowlim, highlim) = (df.values.min(), 0) ##SS 20180319: modified limits, ticks
        ticks = [0.001, 0.01, 0.1, 1]
        tickslog = [numpy.log2(i) for i in ticks]
        ticksstr = [str(i) for i in ticks]
        p.set(  xlim=(lowlim, highlim), 
                ylim=(lowlim, highlim), 
                xticks=tickslog,
                yticks=tickslog,
                xticklabels=ticksstr,
                yticklabels=ticksstr,
                )
        # hide upper, diag: https://stackoverflow.com/a/34091733
        for (i, j) in zip(*numpy.triu_indices_from(p.axes, 1)):
            p.axes[i, j].set_visible(False)
        for (i, j) in zip(*numpy.diag_indices_from(p.axes)):
            p.axes[i, j].set_visible(False)
            
    elif datatype in ['effect']:
        (lowlim, highlim) = (0, df.values.max())
        highlim += (highlim - lowlim) * 0.05
        lowlim -= (highlim - lowlim) * 0.05
        ticks = list(range(0,math.ceil(20.03),10))
        ticksstr = [str(i) for i in ticks]
        p.set(xlim=(lowlim, highlim), 
              ylim=(lowlim, highlim),
              xticks=ticks,
              yticks=ticks,
              xticklabels=ticksstr,
              yticklabels=ticksstr,
             )
        # hide upper, diag: https://stackoverflow.com/a/34091733
        for (i, j) in zip(*numpy.triu_indices_from(p.axes, 1)):
            p.axes[i, j].set_visible(False)
        for (i, j) in zip(*numpy.diag_indices_from(p.axes)):
            p.axes[i, j].set_visible(False)
            
            
    elif datatype in ['log2effect']:
        (lowlim, highlim) = (df.values.min(), df.values.max())
        highlim += (highlim - lowlim) * 0.05
        lowlim -= (highlim - lowlim) * 0.05
        p.set(xlim=(lowlim, highlim), ylim=(lowlim, highlim))
        # If I want to plot the non-log-transformed effect numbers
#         ticks = [0.001, 0.01, 0.1, 1, 10]
#         tickslog = [numpy.log2(i) for i in ticks]
#         ticksstr = [str(i) for i in ticks]
#         p.set(xlim=(lowlim, highlim), 
#               ylim=(lowlim, highlim),
#               xticks=tickslog,
#               yticks=tickslog,
#               xticklabels=ticksstr,
#               yticklabels=ticksstr,
#              )
        # hide upper, diag: https://stackoverflow.com/a/34091733
        for (i, j) in zip(*numpy.triu_indices_from(p.axes, 1)):
            p.axes[i, j].set_visible(False)
        for (i, j) in zip(*numpy.diag_indices_from(p.axes)):
            p.axes[i, j].set_visible(False)

    else:
        raise ValueError("invalid datatype")
        p.map_upper(seaborn.kdeplot, n_levels=20, cmap='Blues_d')

    if title:
        # following here: https://stackoverflow.com/a/29814281
        p.fig.suptitle(title)

    p.savefig(plotfile)
    plt.close()

    
    
def permtest(dataframe, metric, group, perm_no):
    """Performs permutation test, returns (p-value, eq)
    
    Takes as input:
    - 'dataframe': Dataframe which contains data
    - 'metric': column name for value being compared between two groups tested.
    e.g. the mutation differential selection value, 'mutdiffsel'.
    - 'group': column name for group being tested.
    e.g. whether of not a mutation has been experimentally verified, 'ExptVerified'.
    - 'perm_no': number of permutations to perform
    
    First, calculates test statistic for actual data.
    In this case, the test statistic is the difference in medians of the two groups.
    Then, shuffle column containing the values being compared.
    For each shuffle, recalculate test statistic.
    Variable 'count' tracks number of times test statistic >= actual data.
    If count > 0, p is count/perm_no, and eq is string '='.
    If count = 0, p 1/perm_no, and eq is string '<'.
    
    """   
    #calculate test statistic for actual data: difference of medians
    x = dataframe[dataframe[group] == 'Yes'][metric].median()
    y = dataframe[dataframe[group] == 'No'][metric].median()
    stat_actual = x - y
    #shuffle diffsel values
    count = 0
    for i in range(perm_no):
        df = dataframe.copy()
        df[metric] = numpy.random.permutation(df[metric])
        x = df[df[group] == 'Yes'][metric].median()
        y = df[df[group] == 'No'][metric].median()
        stat = x - y
#         if stat >= 0: # for testing
        if stat >= stat_actual:
            count +=1 
    if count > 0:
        p = count/perm_no
        eq = '='
    else:
        p = 1/perm_no
        eq = '<='
    return (p, eq)