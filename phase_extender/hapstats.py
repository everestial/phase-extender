import os
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import seaborn as sns

""" function that computes several stats from haplotype block before and after phase extension.
These data can be used to compare the improvements in phase extension, and plotting histogram.
This can include the size based on genome co-ordinate of the "POS" in each block.
This can also include the average number of variants per haplotype before/after phase extension.
and make a plot out of it - using pyPlot, MatlibPlot or R. """
plt.style.use('seaborn')


def compute_haplotype_stats(hap_data, soi, prefix, outputdir,logscale_x=None,logscale_y=None, show_plot=False):
    # compute stats from given hap_data
    stats_df_fname = prefix + "_haplotype_stats_" + soi + ".txt"
    stats_filepath = Path(outputdir, stats_df_fname)
    hap_stats = compute_stats_df(hap_data, soi, prefix, stats_filepath)
    plot_all_data(hap_stats, soi, prefix, outputdir, logscale_x=logscale_x, logscale_y=logscale_y, show_plot=False)
    # for plot total number of variants per-chromosome


def plot_all_data(
    hap_data_file, soi="", prefix="noprefix", outputdir=None,logscale_x=None, logscale_y=None,show_plot=False
):
    if outputdir is None:
        outputdir = Path(Path.cwd())
    if isinstance(hap_data_file, pd.DataFrame):
        hap_stats = hap_data_file
    else:
        print("computing stats from file")
        hap_stats = pd.read_csv(hap_data_file, sep="\t")

    var_path = Path(outputdir, "total_vars_" + soi + "_" + prefix + ".png")
    xlabel = "chromosomes"
    ylabel = "number of variants"
    suptitle = f"{prefix} number of variants for each chromosome"
    plot_bar_hapstats(
        hap_stats,
        xcol="CHROM",
        ycol="total_Vars",
        filepath=var_path,
        xlabel=xlabel,
        ylabel=ylabel,
        suptitle=suptitle,
        show_plot=show_plot,
    )
    # plots total number of haplotypes per-chromosome
    hap_path = Path(outputdir, "total_haps_" + soi + "_" + prefix + ".png")
    ylabel = "number of haplotypes"
    suptitle = f"{prefix} number of haplotypes for each chromosome"
    plot_bar_hapstats(
        hap_stats,
        xcol="CHROM",
        ycol="total_haplotypes",
        filepath=hap_path,
        xlabel=xlabel,
        ylabel=ylabel,
        suptitle=suptitle,
        show_plot=show_plot,
    )

    # plot histogram of haplotype size (by number of variants) per-chromosome
    # make different plots for each chromosome, but X-axis is shared.
    varsize_path = Path(outputdir, "hap_size_byVar_" + soi + "_" + prefix + ".png")
    genomic_path = Path(
        outputdir, "hap_size_byGenomicRange_" + soi + "_" + prefix + ".png"
    )

    if len(hap_stats) == 1:
        plot_hist_one_chr(
            hap_stats,
            hist_by="num_Vars_by_PI",
            filepath=varsize_path,
            prefix=prefix,
            logscale_x=logscale_x,
            logscale_y=logscale_y,
            show_plot=show_plot,
        )
        plot_hist_one_chr(
            hap_stats,
            hist_by="range_of_PI",
            filepath=genomic_path,
            prefix=prefix,
            logscale_x=logscale_x,
            logscale_y=logscale_y,
            show_plot=show_plot,
        )

    else:
        plot_hist_multi_chr(
            hap_stats,
            hist_by="num_Vars_by_PI",
            filepath=varsize_path,
            prefix=prefix,
            logscale_x=logscale_x,
            logscale_y=logscale_y,

            show_plot=show_plot,
        )
        plot_hist_multi_chr(
            hap_stats,
            hist_by="range_of_PI",
            filepath=genomic_path,
            prefix=prefix,
            logscale_x=logscale_x,

            logscale_y=logscale_y,
            show_plot=show_plot,
        )


def compute_stats_df(hap_data, soi, prefix, filepath):
    hap_stats = (
        hap_data.groupby(["CHROM", soi + ":PI"], sort=False)
        .size()
        .rename("num_Vars_by_PI")
        .reset_index(level=1)
        .astype(str)
        .groupby(level=0)
        .agg(",".join)
        .reset_index()
    )

    """ add other required columns """
    hap_stats["range_of_PI"] = (
        hap_data.groupby(["CHROM", soi + ":PI"], sort=False)["POS"]
        .apply(lambda g: g.max() - g.min())
        .rename("range")
        .reset_index(level=1)
        .astype(str)
        .groupby(level=0)
        .agg(",".join)
        .reset_index()["range"]
    )

    hap_stats["total_haplotypes"] = hap_stats.apply(
        lambda row: len(row[soi + ":PI"].split(",")), axis=1
    )

    hap_stats["total_Vars"] = hap_stats.apply(
        lambda row: sum([int(i) for i in row.num_Vars_by_PI.split(",")]), axis=1
    )

    """ write this concise statistics to a file. """
    hap_stats.to_csv(filepath, sep="\t", header=True, index=False)
    return hap_stats


def plot_bar_hapstats(
    hap_stats,
    xcol,
    ycol,
    filepath,
    xlabel=None,
    ylabel=None,
    suptitle=None,
    show_plot=False,
):
    ax = hap_stats.plot(x=xcol, y=ycol, kind="bar", title=suptitle)
    ax.set(xlabel=xlabel, ylabel=ylabel)
    if show_plot:
        plt.show()
    else:
        ax.figure.savefig(filepath)
        plt.close()
    return


def plot_hist_one_chr(hap_stats, hist_by, filepath, prefix,logscale_x=None,
            logscale_y=None, show_plot=False):
    data_i = hap_stats[hist_by].str.split(",")
    data_i = pd.Series([int(x) for x in data_i[0]])
    ax= sns.histplot(data_i, alpha=0.5, log_scale=True)

    # ax = data_i.plot(kind="hist", label=str(hap_stats["CHROM"]), alpha=0.5, log_scale=logscale_y)
    sub_title = (
        "number of variants" if hist_by == "num_Vars_by_PI" else "Genomic Distance"
    )
    title = f"{prefix} histogram of size of the haplotype (by {sub_title}) \n for each chromosome"
    ax.set_title(title)
    if logscale_x:
        sub_title += " (in log10 scale)"
    ax.set_xlabel(f"size of the haplotype (by {sub_title})")

    if logscale_y:

        y_label = "frequency of the haplotypes(in log10 scale) "
    else:
        y_label = "frequency of the haplotypes"
    
    ax.set_ylabel(y_label)

    # ax.text(.05, .5, 'frequency of the haplotypes', ha='center', va='center', rotation='vertical')
    if show_plot:
        plt.show()
    else:
        ax.figure.savefig(filepath)
        plt.close()
    return


def plot_hist_multi_chr(hap_stats, hist_by, filepath, prefix,logscale_x=None, logscale_y=None, show_plot=False):
    fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True)
    for i, data in hap_stats.iterrows():
        # first convert data to list of integers
        data_i = [int(x) for x in data[hist_by].split(",")]
        g= sns.histplot(data_i, alpha=0.5, ax= ax[i], log_scale=logscale_y)
        g.set(ylabel = None)
        ax[i].legend(str(data["CHROM"]))
        # data_i = [int(x) for x in data[hist_by].split(",")]
        # ax[i].hist(data_i, label=str(data["CHROM"]), alpha=0.5)
        # ax[i].legend()
    if logscale_y:
        y_label = "frequency of the haplotypes (in log10 scale)"
    else:
        y_label = "frequency of the haplotypes"


    fig.text(
        0.05,
        0.5,
        y_label,
        ha="center",
        va="center",
        rotation="vertical",
    )

    # add a common x-label
    sub_title = (
        "number of variants" if hist_by == "num_Vars_by_PI" else "Genomic Distance"
    )
    title = f"{prefix} histogram of size of the haplotype (by {sub_title}) \n for each chromosome"
    # plt.set_title(title)

    if logscale_x:
        x_label = f"size of the haplotype (by {sub_title}) (in log10 scale)"
    else:
        x_label = f"size of the haplotype (by {sub_title})"

    plt.xlabel(x_label)

    # plt.ylabel('frequency of the haplotypes')
    plt.suptitle(title)

    if show_plot:
        plt.show()
    else:
        plt.savefig(filepath)
        plt.close()
    return

def plot_stacked_haplotypes_variants(merged_df, outputdir= None, show_plot=True):
    if outputdir is None:
        outputdir = os.getcwd()

    fig, axs = plt.subplots(ncols=2, figsize=(12,6), gridspec_kw={'width_ratios': [3, 1]})
    # axs.autoscale(enable=True) 
    g =sns.barplot(y = 'total_Vars', data=merged,x='CHROM',color='#FFBD33', ax= axs[1].twinx())
    # g.set(yticklabels=[])
    # ax[1].set_y_axis('')
    sns.barplot(y = 'total_haplotypes', data=merged, x='CHROM', hue= 'stage', ax= axs[0])
    axs[0].set(title="Bar plot showing total number \n of haplotypes for each chromosomes")
    axs[0].set(ylabel='total number of haplotypes')
    axs[1].set(title="Bar plot showing total number \n  of variants for each chromosomes")
    axs[1].set(ylabel='total number of variants')

    plt.tight_layout()
    if show_plot:
        plt.show()
    fig.savefig(outputdir+"/total_haps_stacked_with_vars2.png")

def get_melted_df(df, col_name= 'range_of_PI'):
    splited = df[col_name].str.split(',', expand=True)
    splited['stage'] = ['initial', 'iteration01', 'iteration02']
    melted_df = splited.melt(id_vars= 'stage')
    melted_df.dropna(inplace=True)
    melted_df['value'] = pd.to_numeric(melted_df['value'])
    return melted_df


def plot_hist_one_chr_stacked(melted_df, hist_by=None, filepath=None, prefix= 'Stacked',logscale_x=True,
            logscale_y=True, show_plot=False):
    sub_title = (
            "number of variants" if hist_by == "num_Vars_by_PI" else "genomic distance"
        )
    title = f"{prefix} histogram of haplotype size by {sub_title} \n for each chromosome"
    ax = sns.displot(data=melted_df, x='value', hue='stage', log_scale=[logscale_x, logscale_y])
    # sns.move_legend(ax, "upper right", title='Stages')
    # ax.fig.suptitle(title)
    # if logscale_x:
    #     sub_title += " (in log2 scale)"
    # ax.set_xlabel(f"size of the haplotype (by {sub_title})")

    # if logscale_y:

    #     y_label = "frequency of the haplotypes(in log2scale) "
    # else:
    y_label = "frequency of the haplotypes"
    x_label = f"haplotype size by {sub_title}"
    
    # ax.set_ylabel(y_label)
    ax.set(title=title, xlabel =x_label, ylabel= y_label )
    plt.tight_layout()

    # ax.text(.05, .5, 'frequency of the haplotypes', ha='center', va='center', rotation='vertical')
    if show_plot:
        plt.show()
    if filepath:
        ax.figure.savefig(filepath)
        plt.close()

def get_merged_metrices(merged_df, outputdir=None):
    new_df = pd.DataFrame()
    new_df['CHROM'] = merged_df['CHROM']
    new_df['max_range_pi'] = merged_df['range_of_PI'].apply(lambda x: max((int(xs) for xs in x.split(','))))
    new_df['min_range_pi'] = merged_df['range_of_PI'].apply(lambda x: min((int(xs) for xs in x.split(','))))
    new_df['max_num_vars_pi']= merged_df['num_Vars_by_PI'].apply(lambda x: max((int(xs) for xs in x.split(','))))
    new_df['min_num_vars_pi'] = merged_df['num_Vars_by_PI'].apply(lambda x: min((int(xs) for xs in x.split(','))))
    new_df['stage'] = merged_df['stage']
    if outputdir:
        new_df.to_csv(outputdir+"/metrics.csv", index=False)
        return


    return new_df

# import matplotlib.pyplot as plt
# import pandas as pd

# ''' function that computes several stats from haplotype block before and after phase extension.
# These data can be used to compare the improvements in phase extension, and plotting histogram.
# This can include the size based on genome co-ordinate of the "POS" in each block.
# This can also include the average number of variants per haplotype before/after phase extension.
# and make a plot out of it - using pyPlot, MatlibPlot or R. '''
# def compute_haplotype_stats(hap_data, soi, prefix, outputdir) :

#     print()
#     ''' this function computes the stats of the haplotype file both before and after phase extension.
#         - we pass in "prefix" as "initial" vs. "final" for the haplotype block before and after
#           extension to compute statistics for respective phased file. '''

#     ''' Step 09 - A : write this haplotype information in a more concise format.
#         - in this step we are restructuring the data so proper statistics can be handled properly. '''
#     hap_stats = hap_data.groupby(['CHROM', soi +':PI'], sort=False)\
#         .size().rename('num_Vars_by_PI')\
#         .reset_index(level=1).astype(str)\
#         .groupby(level=0).agg(','.join).reset_index()

#     ''' add other required columns '''
#     hap_stats['range_of_PI'] = hap_data.groupby(['CHROM', soi +':PI'], sort=False)['POS']\
#         .apply(lambda g: g.max() - g.min()).rename('range').reset_index(level=1)\
#         .astype(str).groupby(level=0).agg(','.join).reset_index()['range']

#     hap_stats['total_haplotypes'] = hap_stats\
#         .apply(lambda row: len(row[soi +':PI'].split(',')), axis=1)

#     hap_stats['total_Vars'] = hap_stats \
#         .apply(lambda row: sum([int (i) for i in row.num_Vars_by_PI.split(',')]), axis=1)

#     ''' write this concise statistics to a file. '''
#     pd.DataFrame.to_csv(hap_stats, outputdir + '/' + prefix+'_haplotype_stats_' + soi + '.txt',
#                         sep='\t', header=True, index=False)


#     ## ** for future - plot haplotype STATs - may need improvement.
#     ''' Step 09 - B : Now, use the above dataframe to plot several statistical plots from here on.
#         - plot the haplotype statistics using matplotlib, pylab or ggplot (python).
#         - start with "with open() .... " to open a file, write the plot and then close. '''
#     # some good links for hist-plots: http://pandas.pydata.org/pandas-docs/version/0.18/visualization.html

#     # plots total number of variants per-chromosome
#     with open(outputdir + '/' + 'total_vars_'+ soi +'_'+ prefix+'.png', 'wb') as fig_initial:
#         hap_stats.plot(x='CHROM', y='total_Vars', kind='bar')
#         plt.xlabel('chromosomes')
#         plt.ylabel('number of variants')
#         plt.suptitle('number of variants for each chromosome')
#         plt.savefig(fig_initial)
#         # plt.show()  # optional: to show the plot to the console

#     # plots total number of haplotypes per-chromosome
#     with open(outputdir + '/' + 'total_haps_' + soi + '_' + prefix + '.png', 'wb') as fig_initial:
#         hap_stats.plot(x='CHROM', y='total_haplotypes', kind='bar')
#         plt.xlabel('chromosomes')
#         plt.ylabel('number of haplotypes')
#         plt.suptitle('number of haplotypes for each chromosome')
#         plt.savefig(fig_initial)


#     # plot histogram of haplotype size (by number of variants) per-chromosome
#     # make different plots for each chromosome, but X-axis is shared.
#     with open(outputdir + '/' + 'hap_size_byVar_' + soi + '_' + prefix + '.png', 'wb') as fig_initial:

#         ## making histogram plot.
#         # raising if-else when number of chromosome is 1 vs. more than 1. ** This can be optimize in the future.
#         fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True)

#         # when we have only one chromosome
#         if len(hap_stats) == 1:
#             # fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True, squeeze=False)
#             data_i = hap_stats['num_Vars_by_PI'].str.split(',')
#             data_i = pd.Series([int(x) for x in data_i[0]])
#             data_i.plot(kind='hist', label=str(hap_stats['CHROM']), alpha=0.5)
#             plt.ylabel('frequency of the haplotypes')

#         elif len(hap_stats) > 1:
#             for i, data in hap_stats.iterrows():
#                 # first convert data to list of integers
#                 data_i = [int(x) for x in data['num_Vars_by_PI'].split(',')]
#                 # print(data)
#                 # print(data_i)
#                 ax[i].hist(data_i, label=str(data['CHROM']), alpha=0.5)
#                 ax[i].legend()

#             fig.text(.05, .5, 'frequency of the haplotypes', ha='center', va='center', rotation='vertical')

#         # add a common x-label
#         plt.xlabel('size of the haplotype (by number of variants)')
#         # plt.ylabel('frequency of the haplotypes')
#         plt.suptitle(prefix + ' histogram of size of the haplotype (by number of variants) \n'
#                      'for each chromosome')
#         plt.savefig(fig_initial)


#     # plot histogram of haplotype size (by genomic ranges of haplotype) per-chromosome
#     # make different plots for each chromosome, but X-axis is shared.
#     with open(outputdir + '/' + 'hap_size_byGenomicRange_' + soi + '_' + prefix + '.png', 'wb') as fig_initial:

#         ## making histogram plot.
#         # raising if-else when number of chromosome is 1 vs. more than 1. ** This can be optimize in the future.
#         fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True)

#         # when we have only one chromosome
#         if len(hap_stats) == 1:
#             # fig, ax = plt.subplots(nrows=len(hap_stats), sharex=True, squeeze=False)
#             data_i = hap_stats['range_of_PI'].str.split(',')
#             data_i = pd.Series([int(x) for x in data_i[0]])
#             data_i.plot(kind='hist', label=str(hap_stats['CHROM']), alpha=0.5)
#             plt.ylabel('frequency of the haplotypes')

#         # when the dataframe has more than one chromosome
#         elif len(hap_stats) > 1:
#             for i, data in hap_stats.iterrows():
#                 # first convert data to list of integers
#                 data_i = [int(x) for x in data['range_of_PI'].split(',')]
#                 # print(data)
#                 # print(data_i)
#                 ax[i].hist(data_i, label=str(data['CHROM']), alpha=0.5)
#                 ax[i].legend()

#             # adding y-label separately if there are multiple chromosomes
#             fig.text(.05, .5, 'frequency of the haplotypes', ha='center', va='center', rotation='vertical')

#         # writing a common x-label
#         plt.xlabel('size of the haplotype (by Genomic Distance)')
#         #plt.ylabel('frequency of the haplotypes')
#         plt.suptitle(prefix + ' histogram of size of the haplotype (by Genomic Distance)\n'
#                      'for each chromosome')
#         plt.savefig(fig_initial)


#     #print('Global maximum memory usage: %.2f (mb)' % current_mem_usage())

#     # clear memory
#     del hap_data, hap_stats
