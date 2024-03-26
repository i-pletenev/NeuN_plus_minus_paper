import numpy as np
import pandas as pd
import cooler
import cooltools
import bioframe as bf

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib import ticker


##################################
# Neuronal dots-related functions
##################################

def group_expected(cvd, balance=True, smooth=False):
    """
    Convert expected values calculated per-chromosome to
    whole-genome ones
    
    Parameters
    ----------
    cvd : pd.DataFrame
        output of cooltools.expected_cis()
    balance : bool
        Use balanced or raw counts.
    smooth : bool
        If True, keep cvd['balanced.avg.smoothed.agg'] column
        in the output. Requires cvd['balanced.avg.smoothed.agg'] column
        
    Returns
    -------
    cvd_grp : pd.DataFrame
        Whole-genome expected values for each distance
    """
    col_prefix = 'balanced' if balance else 'count'
    sum_col = col_prefix + '.sum'
    agg_col = col_prefix + '.sum.agg'
    avg_col = col_prefix + '.avg.agg'
    
    # set rules for groupby
    agg_dict = {'n_valid': 'sum', 
                sum_col: 'sum'}
    if smooth == True:
        agg_dict['balanced.avg.smoothed.agg'] = 'first'
        
    # group by chromosome
    keep_cols = ['dist', 'n_valid', sum_col]
    if smooth == True:
        keep_cols.append('balanced.avg.smoothed.agg')
    cvd_grp = cvd[keep_cols]\
                .groupby('dist')\
                .agg(agg_dict)\
                .rename(columns={'n_valid': 'n_valid.agg',
                                 sum_col: agg_col})
    cvd_grp[avg_col] = cvd_grp[agg_col] / cvd_grp['n_valid.agg'] # whole-genome expected
       
    return cvd_grp


def expected_std(clr, cvd, exp_col='balanced.avg.agg', n_col='n_valid.agg', balance=True, sparse=False):
    """
    Compute standard deviation of expected
    
    Parameters
    ----------
    clr : cooler
        Cooler object to fetch data from
    cvd : pd.DataFrame
        Output of group_expected()
    exp_col : str
        'cvd' column name to extract expected value
    n_col : str
        'cvd' column name to extract number of valid pixels
    balance : bool
        Passed to 'clr.matrix()'
    sparse : bool
        Passed to 'clr.matrix()'. sparse=True slows down
        calculations but requires less RAM.
        
    Returns
    -------
    df['std'] : pd.Series
        Series with standard deviation for each distance
    """

    dist_out = []
    sum_sq_out = []
    chroms = clr.chromnames
        
    for chrom in chroms:
        mat = clr.matrix(balance=balance, sparse=sparse).fetch(chrom)
        n = mat.shape[0]
        dist_out.extend(list(range(n)))
        for i in range(n):
            exp = cvd.loc[i, exp_col]
            arr = mat.diagonal(i)
            arr = arr[(arr == arr) & (arr != 0)]
            
            sum_sq = np.nansum((arr - exp)**2)
            sum_sq_out.append(sum_sq)

    df = pd.DataFrame({'dist': dist_out,
                       'sum_sq': sum_sq_out})
    
    # aggregate by distance
    df = df.groupby('dist').sum()\
           .join(cvd[n_col])
    df['std'] = np.sqrt( df['sum_sq'] / (df[n_col] - 1) )
    
    return df['std']


def k27_hic_signal(k27_df, clr, cvd, view_df, min_sep=200000):
    """
    Get Hi-C values corresponding to h3k27me3 interactions
    1) assign peaks to corresponding hic bins. If more than one peak overlaps a bin, keep longest peak
    2) make all possible pairs out of k27 regions
    3) get corresponding hi-c obs/exp values
    
    k27_df - df with k27 peak coordinates.
    """
    k27_df = k27_df.copy(deep=True)
    if 'length' not in k27_df.columns:
        k27_df['length'] = k27_df['end'] - k27_df['start']
    
    # Convert peak coordinates to corresponding hic bins
    clr_bins = clr.bins()[:]\
                  .dropna(subset=['weight'])[['chrom', 'start', 'end']]
    k27_df_bins = bf.overlap(clr_bins, k27_df, how='inner', suffixes=('_bin', ''))\
                    .drop(columns=['chrom', 'start', 'end'])\
                    .rename(columns={'chrom_bin': 'chrom', 'start_bin': 'start', 'end_bin': 'end'})
    
    # Keep longest peak if multiple are in one bin
    k27_df_bins = k27_df_bins.sort_values(by='length', ascending=False)\
                             .drop_duplicates(subset=['chrom', 'start', 'end'], keep='first')
    
    pair_sites = bf.pair_by_distance(k27_df_bins, min_sep=min_sep, 
                                     max_sep=int(1e10), suffixes=('1', '2'))
    pup = cooltools.pileup(clr, pair_sites, view_df=view_df, 
                           expected_df=cvd, flank=0)
    pup_obs = cooltools.pileup(clr, pair_sites, view_df=view_df, 
                               expected_df=None, flank=0)
    pair_sites['obs_exp'] = pup.flatten()
    pair_sites['obs'] = pup_obs.flatten()
    
    return pair_sites


def bedpe_to_bed(df):
    """
    Extract unique bed regions from bedpe dataframe. Index is ignored.
    """
    df_reg1 = df[['chrom1', 'start1', 'end1']]
    df_reg2 = df[['chrom2', 'start2', 'end2']]

    df_reg1.columns = df_reg1.columns.str.rstrip('1')
    df_reg2.columns = df_reg2.columns.str.rstrip('2')

    df_bed = pd.concat([df_reg1, df_reg2], ignore_index=True)\
               .drop_duplicates(ignore_index=True)
    
    return df_bed   


def regs_to_bins(regs_df, binsize, mode='middle', return_index=True):
    """
    Assign bins to each genomic region
    
    Parameters
    ----------
    regs_df : pd.DataFrame
        Dataframe with genomic regions. 
        Required columns are: 'chrom', 'start', 'end'
    binsize : int
    mode : 'middle' or 'all'
        If 'middle': for each region, get a bin that overlap the middle of the region.
        If 'all': for each region, get all bins that overlap this region and write to separate rows 
    return_index : bool
        If True and mode='all', return index of input dataframe.
        
    Returns
    -------
    bins_df : pd.DataFrame
        Dataframe with bin coordinates written in
        'chrom', 'start', 'end' columns
    """
    
    if mode == 'middle':
        bins_df = regs_df.copy(deep=True)
        bins_df.loc[:, 'start'] = (bins_df['start'] + bins_df['end']) / 2 \
        - ((bins_df['start'] + bins_df['end']) / 2) % binsize
        bins_df.loc[:, 'end'] = bins_df['start'] + binsize
        bins_df.loc[:, ['start', 'end']] = bins_df.loc[:, ['start', 'end']].astype(int)
        
    elif mode == 'all':
        regs_cp = regs_df.copy(deep=True)
        regs_cp.loc[:, 'first_bin_start'] = regs_cp['start'] - (regs_cp['start'] % binsize)
        regs_cp.loc[:, 'last_bin_end'] = regs_cp['end'] - 1 - (regs_cp['end'] - 1) % binsize + binsize
        regs_cp.loc[:, 'nbins'] = (regs_cp['last_bin_end'] - regs_cp['first_bin_start']) // binsize
        regs_cp.loc[:, 'idx'] = regs_cp.index
        
        all_bins = np.empty((regs_cp['nbins'].sum(), 3), dtype=int)
        i = 0
        for reg_id in range(regs_cp.shape[0]):
            n, start, end, idx = regs_cp.loc[
                reg_id, ['nbins', 'first_bin_start', 'last_bin_end', 'idx']]
            all_bins[i:i+n, 0] = np.arange(start, end, binsize)
            all_bins[i:i+n, 1] = np.arange(start+binsize, end+binsize, binsize)
            all_bins[i:i+n, 2] = np.tile(idx, n)
            i+=n
            
        bins_df = pd.DataFrame(data=all_bins, columns=['start', 'end', 'idx'])
        regs_cp_drop_cols = ['start', 'end', 'first_bin_start', 'last_bin_end', 'nbins']
        bins_df = bins_df.merge(regs_cp.drop(columns=regs_cp_drop_cols), on='idx', how='left')
        bins_df.insert(0, 'chrom', bins_df.pop('chrom'))
        
        if not return_index:
            bins_df = bins_df.drop(columns='idx')
        else:
            bins_df.insert(bins_df.shape[1]-1, 'idx', bins_df.pop('idx'))
        
    else:
        raise ValueError()
        
    return bins_df


def obs_for_bin_pairs(bin_pairs, clr, obs_colname='obs', balance=True):
    """
    Get observed Hi-C values for pairs of bins
    
    Parameters
    ----------
    bin_pairs : pd.DataFrame
        Dataframe with pairs of genomic regions. 
        Required columns are: 
        'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'.
    clr : cooler.Cooler
        cooler to extract Hi-C counts
    obs_colname : str
        Name of column to store Hi-C values
    balance : bool
        Passed to 'clr.matrix()'
        
    Returns
    -------
    bin_pairs_obs : pd.DataFrame
        Input dataframe with new column 'obs'
    """
    bin_pairs_obs = bin_pairs
   
    mtx = clr.matrix(balance=balance)
    bin_pairs_obs[obs_colname] = bin_pairs_obs.apply(
        lambda crds: mtx.fetch([crds['chrom1'], crds['start1'], crds['end1']],
                               [crds['chrom2'], crds['start2'], crds['end2']])[0][0],
        axis=1)
    
    return bin_pairs_obs


def format_ticks_zoom_out(ax, base=2e7, x=True, y=True, rotate=True):
    bp_formatter = ticker.EngFormatter('b')
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base))
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base))
        ax.xaxis.tick_bottom()


def plot_hic_reg_zoom_out(clr, reg, ax, vmin, vmax, snips=None, title=None):
    """
    Plot region of a Hi-C heatmap. 
    Optionally, contour parts of a heatmap (snips) with rectangles.
    
    Parameters
    ----------
    clr : cooler
        Cooler object to fetch data from
    reg : list or tuple
        List or tuple with chromosome name, region start 
        and region end of the form (chrom, start, end)
    ax : matplotlib Axes
        Axes object to draw the plot onto
    vmin, vmax, cmap
        Passed on to `ax.matshow()`
    snips : list
        List of snips; each snip is a list itself
        conaining two region coordinates of the form 
        (chrom, start, end)
    title : str
        Title for the plot
    """
    # Create custom colormap
    cmap_colors = np.array([
        [255, 255, 255, 255],
        [245, 166, 35, 255],
        [208, 2, 27, 255],
        [0, 0, 0, 255]
    ]) / 255
    cmap_nodes = [0, 0.35, 0.6, 1]
    cust_cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(cmap_nodes, cmap_colors)))
    
    # Plot heatmap
    chrom, start, end = reg
    im = ax.matshow(
        clr.matrix(balance=True).fetch(reg),
        norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap=cust_cmap,
        extent=(start, end, end, start), 
        interpolation=None,
    );
    
    # Set title
    if title:
        ax.set_title(title)
    
    format_ticks_zoom_out(ax)
    
    # Draw rectangles
    if snips:
        for j, (reg1, reg2) in enumerate(snips):
            rect = patches.Rectangle(
                (reg2[1], reg1[1]), reg2[2]-reg2[1], reg1[2]-reg1[1], 
                edgecolor='k', facecolor='none', linewidth=1
            )
            ax.add_patch(rect)
    
    # Draw colorbar
    cb = plt.colorbar(im, ax=ax, shrink=0.5, location='right')
    cb.ax.set_ylabel('Contact probability', rotation=270, fontsize=8)
    
    return


def plot_hic_reg_zoom_in(clr, reg1, reg2, ax, vmin, vmax):
    start1, end1 = reg1[1], reg1[2]
    start2, end2 = reg2[1], reg2[2]
    
    # Create custom colormap
    cmap_colors = np.array([
        [255, 255, 255, 255],
        [245, 166, 35, 255],
        [208, 2, 27, 255],
        [0, 0, 0, 255]
    ]) / 255
    cmap_nodes = [0, 0.35, 0.6, 1]
    cust_cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(cmap_nodes, cmap_colors)))
    
    im = ax.matshow(
        clr.matrix(balance=True).fetch(reg1, reg2),
        norm=LogNorm(vmin, vmax),
        cmap=cust_cmap,
        extent=(start2, end2, end1, start1),
    );
    ax.set_xticks([])
    ax.set_yticks([])
    
    return im