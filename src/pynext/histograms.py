import numpy as np
import random
from dataclasses import dataclass

import matplotlib.pyplot as plt

from typing import Tuple, Optional
from typing      import TypeVar

from typing      import Tuple
from typing      import Dict
from typing      import List
from typing      import TypeVar
from typing      import Optional

from enum        import Enum

Number = TypeVar('Number', None, int, float)
Str   = TypeVar('Str', None, str)
Range = TypeVar('Range', None, Tuple[float, float])
Array = TypeVar('Array', List, np.array)
Int = TypeVar('Int', None, int)

NN = np.nan

from . stats import mean_and_std



@dataclass
class HistoPar:
    var    : np.array
    nbins  : int
    range  : Tuple[float]


@dataclass
class HistoPar2(HistoPar):
    var2    : np.array
    nbins2 : int
    range2 : Tuple[float]


@dataclass
class ProfilePar:
    x  : np.array
    y  : np.array
    xu : np.array
    yu : np.array


@dataclass
class PlotLabels:
    x     : str
    y     : str
    title : str


def labels(pl : PlotLabels):
    """
    Set x and y labels.
    """
    plt.xlabel(pl.x)
    plt.ylabel(pl.y)
    plt.title (pl.title)



def h1(x      : np.array,
       bins    : int,
       range   : Tuple[float],
       weights : Array = None,
       log     : bool  = False,
       normed  : bool  = False,
       color   : str   = 'black',
       width   : float = 1.5,
       style   : str   ='solid',
       stats   : bool  = True,
       lbl     : Optional[str]  = None):
    """
    histogram 1d with continuous steps and display of statsself.
    number of bins (bins) and range are compulsory.
    """

    mu, std = mean_and_std(x, range)

    if stats:
        entries  =  f'Entries = {len(x)}'
        mean     =  r'$\mu$ = {:7.2f}'.format(mu)
        sigma    =  r'$\sigma$ = {:7.2f}'.format(std)
        stat     =  f'{entries}\n{mean}\n{sigma}'
    else:
        stat     = ''

    if lbl == None:
        lab = ' '
    else:
        lab = lbl

    lab = stat + lab

    if color == None:
        n, b, p = plt.hist(x,
                       bins      = bins,
                       range     = range,
                       weights   = weights,
                       log       = log,
                       density   = normed,
                       histtype  = 'step',
                       linewidth = width,
                       linestyle = style,
                       label     = lab)

    else:

        n, b, p = plt.hist(x,
                       bins      = bins,
                       range     = range,
                       weights   = weights,
                       log       = log,
                       density   = normed,
                       histtype  = 'step',
                       edgecolor = color,
                       linewidth = width,
                       linestyle = style,
                       label     = lab)

    return n, b, mu, std


def plot_histo(pltLabels: PlotLabels, ax, legend= True,
               legendsize=10, legendloc='best', labelsize=11):

    if legend:
        ax.legend(fontsize= legendsize, loc=legendloc)
    ax.set_xlabel(pltLabels.x,fontsize = labelsize)
    ax.set_ylabel(pltLabels.y, fontsize = labelsize)
    if pltLabels.title:
        plt.title(pltLabels.title)


def h1d(x         : np.array,
        bins      : int,
        range     : Tuple[float],
        weights   : Array               = None,
        log       : bool                = False,
        normed    : bool                = False,
        color     : str                 = 'black',
        width     : float               = 1.5,
        style     : str                 ='solid',
        stats     : bool                = True,
        lbl       : Str                 = None,
        pltLabels : PlotLabels          =PlotLabels(x='x', y='y', title=None),
        legendloc : str                 ='best',
        figsize   : Tuple[float, float] =(6,6)):

    fig = plt.figure(figsize=figsize)
    ax      = fig.add_subplot(1, 1, 1)
    n, b, mu, std    = h1(x, bins=bins, range = range, stats = stats, lbl = lbl)
    plot_histo(pltLabels, ax, legendloc=legendloc)
    return n, b, mu, std


def h2(x         : np.array,
       y         : np.array,
       nbins_x   : int,
       nbins_y   : int,
       range_x   : Tuple[float],
       range_y   : Tuple[float],
       profile   : bool   = True):

    xbins  = np.linspace(*range_x, nbins_x + 1)
    ybins  = np.linspace(*range_y, nbins_y + 1)

    nevt, *_  = plt.hist2d(x, y, (xbins, ybins))
    plt.colorbar().set_label("Number of events")

    if profile:
        x, y, yu     = profile1d(x, y, nbins_x, range_x)
        plt.errorbar(x, y, yu, np.diff(x)[0]/2, fmt="kp", ms=7, lw=3)

    return nevt


def h2d(x         : np.array,
        y         : np.array,
        nbins_x   : int,
        nbins_y   : int,
        range_x   : Tuple[float],
        range_y   : Tuple[float],
        pltLabels : PlotLabels   = PlotLabels(x='x', y='y', title=None),
        profile  : bool          = False,
        figsize=(10,6)):

    fig = plt.figure(figsize=figsize)
    fig.add_subplot(1, 1, 1)

    nevt   = h2(x, y, nbins_x, nbins_y, range_x, range_y, profile)
    labels(pltLabels)
    return nevt
