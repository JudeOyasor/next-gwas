#!/usr/bin/env python
# coding: utf-8

# ## Reference(s):
# ---
# - https://github.com/khramts/assocplots

# In[319]:


import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.mstats import mquantiles
from scipy.stats import binom
from scipy.stats import norm
from scipy.stats import beta
from scipy.stats import linregress
from scipy.stats import chi2
import re
import argparse 


# In[311]:


def qqplot(data, labels, n_quantiles=100, alpha=0.95, error_type='theoretical', distribution = 'binomial', log10conv=True, color=['k', 'r', 'b'], fill_dens=[0.1, 0.1, 0.1], type = 'uniform', title='title'):
    '''
    Function for plotting Quantile Quantile (QQ) plots with confidence interval (CI)
    :param data: NumPy 1D array with data
    :param labels:
    :param type: type of the plot
    :param n_quantiles: number of quntiles to plot
    :param alpha: confidence interval
    :param distribution: beta/normal/binomial -- type of the error estimation. Most common in the literature is 'beta'.
    :param log10conv: conversion to -log10(p) for the figure
    :return: nothing
    '''
    xmax = 0
    ymax = 0
    np.seterr(divide = 'ignore') 
    if type == 'uniform':
        # we expect distribution from 0 to 1
        for j in range(len(data)):
            # define quantiles positions:
            q_pos = np.concatenate([np.arange(99.)/len(data[j]), np.logspace(-np.log10(len(data[j]))+2, 0, n_quantiles)])
            # define quantiles in data
            q_data = mquantiles(data[j], prob=q_pos, alphap=0, betap=1, limit=(0, 1)) # linear interpolation
            # define theoretical predictions
            q_th = q_pos.copy()
            # evaluate errors
            q_err = np.zeros([len(q_pos),2])
            if np.sum(alpha) > 0:
                for i in range(0, len(q_pos)):
                    if distribution == 'beta':
                        q_err[i, :] = beta.interval(alpha, len(data[j])*q_pos[i], len(data[j]) - len(data[j])*q_pos[i])
                    elif distribution == 'binomial':
                        q_err[i, :] = binom.interval(alpha=alpha, n=len(data[j]), p=q_pos[i])
                    elif distribution == 'normal':
                        q_err[i, :] = norm.interval(alpha, len(data[j])*q_pos[i], np.sqrt(len(data[j])*q_pos[i]*(1.-q_pos[i])))
                    else:
                        print('Distribution is not defined!')
                q_err[i, q_err[i, :] < 0] = 1e-15
                if (distribution == 'binomial') | (distribution == 'normal'):
                    q_err /= 1.0*len(data[j])
                    for i in range(0, 100):
                        q_err[i,:] += 1e-15
            # print(q_err[100:, :])
            slope, intercept, r_value, p_value, std_err = linregress(q_th, q_data)
            # print(labels[j], ' -- Slope: ', slope, " R-squared:", r_value**2)
            plt.plot(-np.log10(q_th[n_quantiles-1:]), -np.log10(q_data[n_quantiles-1:]), '-', color=color[j])
            plt.plot(-np.log10(q_th[:n_quantiles]), -np.log10(q_data[:n_quantiles]), '.', color=color[j], label=labels[j])
            xmax = np.max([xmax, - np.log10(q_th[1])])
            ymax = np.max([ymax, - np.log10(q_data[0])])
            # print(- np.log10(q_th[:]))
            if np.sum(alpha)>0:
                if error_type=='experimental':
                    plt.fill_between(-np.log10(q_th), -np.log10(q_data/q_th*q_err[:,0]), -np.log10(q_data/q_th*q_err[:,1]), color=color[j], alpha=fill_dens[j], label='%1.3f CI'%alpha)
        if np.sum(alpha)>0:
            if error_type=='theoretical':
                plt.fill_between(-np.log10(q_th), -np.log10(q_err[:,0]), -np.log10(q_err[:,1]), color=color[j], alpha=fill_dens[j], label='%1.3f CI'%alpha)
    plt.legend(loc=4)
    plt.xlabel('Theoretical -log10(p)')
    plt.ylabel('Experimental -log10(p)')
    plt.plot([0, 100], [0, 100],'--k')
    # print(xmax,ymax)
    plt.xlim([0, np.ceil(xmax)])
    plt.ylim([0, np.ceil(ymax*1.05)])
    plt.title(title)
    plt.tight_layout()
    np.seterr(divide = 'warn') 


# In[312]:


def sorted_nicely( l ):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


# In[313]:


def manhattan(p1, pos1, chr1, label1,
               p2=None, pos2=None, chr2=None, label2=None,
               plot_type='single',
               chrs_plot=None, chrs_names=None,
               cut = 2,
               colors = ['k', '0.5'],
               title='Title',
               xlabel='chromosome',
               ylabel='-log10(p-value)',
               top1 = 0,
               top2 = 0,
               lines = [10, 15],
               lines_colors = ['g', 'r'],
               zoom = None,
               scaling = '-log10'):
    '''
    Static Manhattan plot
    :param p1: p-values for the top panel
    :param pos1: positions
    :param chr1: chromosomes numbers
    :param label1: label
    :param p2: p-values for the bottom panel
    :param pos2: positions
    :param chr2: chromosomes numbers
    :param label2: label
    :param type: Can be 'single', 'double' or 'inverted'
    :param chrs_plot: list of chromosomes that should be plotted. If empty [] all chromosomes will be plotted
    :param cut: lower cut (default 2)
    :param colors: sequence of colors (default: black/gray)
    :param title: defines the title of the plot
    :param xlabel: defines the xlabel of the plot
    :param ylabel: defines the ylabel of the plot
    :param top: Defines the upper limit of the plot. If 0, it is detected automatically.
    :param lines: Horizontal lines to plot.
    :param lines_colors: Colors for the horizontal lines.
    :param zoom: [chromosome, position, range] Zooms into a region.
    :param scaling: '-log10' or 'none' (default -log10)
    :return:
    '''

    # Setting things up
    shift=np.array([0.0])
    plt.clf()

    # If chrs_plot is empty, we need to generate a list of chromosomes
    if chrs_plot is None:
        chrs_list = np.unique(chr1)
        if isinstance(chrs_list[0], str):
            chrs_list = sorted_nicely(chrs_list)
        else:
            chrs_list.sort()
    else:
        chrs_list = chrs_plot


    # If chrs_names is empty, we need to generate a list of names for chromosomes
    if chrs_names is None:
        chrs_names = [str(chrs_list[i]) for i in range(len(chrs_list))]

    plot_positions = False
    if len(chrs_list) == 1:
        plot_positions = True


    for ii, i in enumerate(chrs_list):
        if plot_type != 'single':
            ax1 = plt.subplot(2,1,1)
        else:
            plt.subplot(1,1,1)
        # print(i)
        filt = np.where(chr1==i)[0]
        x = shift[-1]+pos1[filt]
        if scaling=='-log10':
            y = -np.log10(p1[filt])
        elif scaling=='none':
            y = p1[filt]
        else:
            raise ValueError('Wrong "scaling" mode. Choose between "-log10" and "none"')
        plt.plot(x[y>cut], y[y>cut], '.', color=colors[ii % len(colors)])
        shift_f = np.max(x)

        if zoom is not None:
            if zoom[0] == i:
                zoom_shift = zoom[1] + shift[-1]

        if plot_type != 'single':
            plt.subplot(2,1,2)#, sharex=ax1)
            filt = np.where(chr2==i)[0]
            x = shift[-1]+pos2[filt]
            if scaling=='-log10':
                y = -np.log10(p2[filt])
            elif scaling=='none':
                y = p2[filt]
            else:
                raise ValueError('Wrong "scaling" mode. Choose between "-log10" and "none"')
            plt.plot(x[y>cut], y[y>cut], '.', color=colors[ii % len(colors)])
            shift_m = np.max(x)
        else:
            shift_m = 0

        shift = np.append(shift, np.max([shift_f, shift_m]))

        if plot_type != 'single':
            plt.subplot(2,1,1)
        else:
            plt.subplot(1,1,1)
        plt.plot([shift[-1], shift[-1]], [0, 1000], '-k', lw=0.5, color='lightgray')
        plt.xlim([0, shift[-1]])

        if plot_type != 'single':
            plt.subplot(2,1,2)
            plt.plot([shift[-1], shift[-1]], [0, 1000], '-k', lw=0.5, color='lightgray')
            plt.xlim([0, shift[-1]])
        # print(shift)

    # Defining top boundary of a plot
    if top1 == 0:
        if plot_type != 'single':
            if scaling == '-log10':
                top1 = np.ceil(np.max([np.max(-np.log10(p1)), np.max(-np.log10(p2))]))
            elif scaling == 'none':
                top1 = np.ceil(np.max([np.max(p1), np.max(p2)]))
            else:
                raise ValueError('Wrong "scaling" mode. Choose between "-log10" and "none"')
        else:
            if scaling == '-log10':
                top1 = np.ceil(np.max(-np.log10(p1)))
            elif scaling == 'none':
                top1 = np.ceil(np.max(p1))
            else:
                raise ValueError('Wrong "scaling" mode. Choose between "-log10" and "none"')


    if top2 == 0:
        if plot_type != 'single':
            top2 = top1

    # Setting up the position of labels:
    shift_label = shift[-1]
    shift = (shift[1:]+shift[:-1])/2.
    labels = chrs_names

    # Plotting horizontal lines
    for i, y in enumerate(lines):
        if plot_type != 'single':
            plt.subplot(2,1,1)
            plt.plot([0, shift_label], [y, y], color = lines_colors[i])
            plt.subplot(2,1,2)
            plt.plot([0, shift_label], [y, y], color = lines_colors[i])
        else:
            plt.subplot(1,1,1)
            plt.plot([0, shift_label], [y, y], color = lines_colors[i])

    if plot_type != 'single':
        plt.subplot(2,1,1)
        if not plot_positions:
            plt.xticks(shift, labels)
        plt.ylim([cut+0.05, top1])
    else:
        plt.subplot(1,1,1)
        plt.ylim([cut, top1])
    plt.title(title)
    if plot_type != 'single':
        plt.setp(plt.gca().get_xticklabels(), visible=False)
        if not plot_positions:
            plt.xticks(shift)
    else:
        if not plot_positions:
            plt.xticks(shift, labels)

    plt.text(shift_label*0.95,top1*0.95,label1,#bbox=dict(boxstyle="round", fc="1.0"),
            verticalalignment='top', horizontalalignment='right')

    if plot_type != 'single':
        plt.subplot(2,1,2)
        plt.ylim([cut, top2])
        if plot_type == 'inverted':
            plt.gca().invert_yaxis()
        if not plot_positions:
            plt.xticks(shift, labels)
        if plot_type == 'inverted':
            plt.text(shift_label*0.95,top2*0.95,label2,#bbox=dict(boxstyle="round", fc="1.0"),
                verticalalignment='bottom', horizontalalignment='right')
        else:
            plt.text(shift_label*0.95,top2*0.95,label2,#bbox=dict(boxstyle="round", fc="1.0"),
                verticalalignment='top', horizontalalignment='right')
        plt.ylabel(ylabel)
        plt.gca().yaxis.set_label_coords(-0.065,1.)
        plt.xlabel(xlabel)
        # plt.tight_layout(hspace=0.001)
        plt.subplots_adjust(hspace=0.00)
    else:
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

    if zoom is not None:
        if plot_type != 'single':
            plt.subplot(2,1,1)
            plt.xlim([zoom_shift-zoom[2], zoom_shift+zoom[2]])
            plt.subplot(2,1,2)
            plt.xlim([zoom_shift-zoom[2], zoom_shift+zoom[2]])
        else:
            plt.subplot(1,1,1)
            plt.xlim([zoom_shift-zoom[2], zoom_shift+zoom[2]])

    pass 


# ## Plot Data

# In[ ]:


if __name__ == "__main__" :
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--assoc', type=argparse.FileType('r', encoding='UTF-8'), required=True)
    parser.add_argument('--maxT', type=argparse.FileType('r', encoding='UTF-8'), required=True)
    args = parser.parse_args()
    
    assoc = pd.read_csv(filepath_or_buffer= args.assoc,delim_whitespace=  True)
    assoc.sort_values(by=['CHR','SNP','BP'],inplace=True)
    assoc = assoc.reset_index(drop =True)

    maxT = pd.read_csv(filepath_or_buffer= args.maxT, delim_whitespace=  True )
    maxT.sort_values(by=['CHR','SNP'],inplace=True)
    maxT = maxT.reset_index(drop =True)
    
    df = pd.merge(assoc, maxT, how ="inner",left_on=["CHR","SNP"], right_on=["CHR","SNP"] )
    
    mpl.rcParams['figure.dpi']=100
    mpl.rcParams['savefig.dpi']=100
    mpl.rcParams['figure.figsize']=5.375, 5.375

    qqplot( data = [df['P'].dropna()], 
            labels =[''], 
            color=['b'], 
            fill_dens=[0.2], 
            error_type='theoretical', 
            distribution='beta',
            title=''
          )
    plt.savefig('qq_plot.png', dpi=300)
    
    cmap = plt.get_cmap('seismic')
    colors = [cmap(i) for i in [0.0,0.33,0.67,0.90]]

    mpl.rcParams['figure.dpi']=100
    mpl.rcParams['savefig.dpi']=100
    mpl.rcParams['figure.figsize']=11, 5.375

    manhattan( p1= df['P'].dropna(), 
               pos1= df['BP'], 
               chr1= df['CHR'].astype(str), 
               label1= '',
               plot_type='single',
               chrs_plot= None,
               chrs_names=None,
               cut = 0,
               title='',
               xlabel='chromosome',
               ylabel='-log10(p)',
               colors = colors
             )

    plt.savefig('manhattan_plot.png', dpi=300)   

