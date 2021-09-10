#!/usr/bin/env python
# coding=utf-8
'''
Date         : 2021-08-26 10:29:53
LastEditTime : 2021-08-26 10:57:21
LastEditors  : windz
FilePath     : /pipelines/BS_seq_pipeline/script/metaplot_bs_seq.py
'''


import numpy as np
import pyBigWig
import pyranges as pr
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import scipy


# get last pa sites
# get single polya site gene
last_pa_bed = '/public/home/mowp/workspace/termination/cbRNA_pool/polya_sites/cbRNA.last_polya_cluster_summit.bed'
last_pa = pr.read_bed(last_pa_bed, as_df=True)
last_pa.columns = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'ratio']
last_pa.loc[:, 'Chromosome'] = last_pa.loc[:, 'Chromosome'].astype('str')

mask = last_pa.loc[:, 'Name'].str.contains('_1')
single_pa_site_gene = last_pa[mask].loc[:, 'Name'].map(lambda x: x.split('_')[0])

last_pa['Name'] = last_pa['Name'].map(lambda x: x.split('_')[0])
last_pa = last_pa.set_index(['Name'])


# get gene model
gene_model_bed = '/public/home/mowp/db/Arabidopsis_thaliana/isoform/araport11.gene.bed'
gene_model = pr.read_bed(gene_model_bed, as_df=True)
gene_model = gene_model.set_index(['Name'])


# load TSS
tss_bed = '/public/home/mowp/data/public_data/epigentics_data/cage_seq_2020/major_TSS.bed'
tss_bed = pr.read_bed(tss_bed, as_df=True)
tss_bed = tss_bed.set_index(['Name'])


# For bw file
def get_target_site(site_type: str, gene_id: str) -> int:
    if site_type == 'PAS':
        # polya site
        return last_pa.at[gene_id, 'End']
    elif site_type == 'TSS':
        try:
            values = tss_bed.loc[gene_id, :].values
            if values[4] == '+':
                return values[1]
            else:
                return values[2]
        except KeyError:
            values = gene_model.loc[gene_id, :].values
            if values[4] == '+':
                return values[1]
            else:
                return values[2]
    
    elif site_type == 'aTSS':
        # araport11注释的TSS
        values = gene_model.loc[gene_id, :].values
        if values[4] == '+':
            return values[1]
        else:
            return values[2]
        
    # elif site_type == 'TWE':
    #     return int(cb_pool_read_through_len.query('gene_id == @gene_id')['tts'])
            
    elif site_type == 'aTES':
        # araport11注释的TES
        # array(['1', 3629, 5899, '.', '+'], dtype=object)
        values = gene_model.loc[gene_id, :].values
        if values[4] == '+':
            return values[2]
        else:
            return values[1]
    else:
        raise KeyError


def get_bin_cov(methyratio: list, bins: int, threshold=4):
    methyratio = np.array(methyratio)
    results = []
    for i in range(0, len(methyratio), bins):
        bin_methyratio = methyratio[i:i + bins]
        if len(bin_methyratio[~np.isnan(bin_methyratio)]) >= threshold:
            # we selected for bins with at least 4 cytosines
            # that are each covered by at least 4 reads
            mean_methyratio = np.nanmean(bin_methyratio)
        else:
            mean_methyratio = np.nan
        results.append(mean_methyratio)

    return results


# site-point
def get_cov(infile: str, gene_id: str, site: str,
            before: int, after: int, bins: int,
            chrom_prefix: str):
    '''
    计算甲基化水平覆盖度
    '''
    chrom, start, end, _, strand = gene_model.loc[gene_id]

    if chrom in {'Pt', 'Mt', 'chrM', 'chrC'}:
        return None
    chrom = chrom_prefix + chrom

    target_site = get_target_site(site, gene_id)  # site1
    bwfile = pyBigWig.open(infile)
    try:
        if strand == '+':
            methyratio = bwfile.values(chrom, 
                                       target_site - before,
                                       target_site + after)
            cov = get_bin_cov(methyratio, bins)

        else:
            methyratio = bwfile.values(chrom, 
                                       target_site - after,
                                       target_site + before)[::-1]
            cov = get_bin_cov(methyratio, bins)

    except RuntimeError:
        return

    return cov, gene_id


def get_meta_site_result(infile,
                    gene_list,
                    site,
                    before=1000,
                    after=1000,
                    bins=100,
                    chrom_prefix='',
                    threads=64):
    '''
    Example:
    infile = '/public/home/mowp/test/BS_Seq/bw_files/WT_BS_GSM1242401.methratio.cg.bw'
    bins = 50
    b1 = 1000
    a1 = 1000
    tts_cov = get_meta_result(infile, all_gene, 'aTSS', bins=bins, before=b1, after=a1)

    b2 = 1000
    a2 = 1000
    pas_cov = get_meta_result(infile, all_gene, 'PAS', bins=bins, before=b2, after=a2)
    '''
    results = []
    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(gene_list) / threads)
        results = e.map(get_cov,
                        repeat(infile),
                        gene_list,
                        repeat(site),
                        repeat(before),
                        repeat(after),
                        repeat(bins),
                        repeat(chrom_prefix),
                        chunksize=chunksize)

    cov = []
    for res in results:
        if res is not None:
            cov_, gene_id = res
            cov.append(cov_)
    
    cov = np.nanmean(cov, axis=0)
    return cov


def set_ax(ax, b1, a1, b2, a2, bins, ylabel=None):
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].set_ylabel(ylabel)

    ax[1].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].yaxis.set_ticks_position('none')

    ax[0].set_xticks([0, b1//bins, (a1+b1)//bins])
    ax[0].set_xticklabels([f'-{b1//1000} kb', 'TSS', f'{a1//1000} kb'], rotation=90)

    ax[1].set_xticks([0, b2//bins, (a2+b2)//bins])
    ax[1].set_xticklabels([f'-{b2//1000} kb', 'PAS', f'{a2//1000} kb'], rotation=90)

    ax[0].axvline(b1//bins, ls='--', color='#555555')
    ax[1].axvline(b2//bins, ls='--', color='#555555')


# reference scale
def get_scale_cov(infile: str, gene_id: str, site1: str, site2: str,
                  before: int, after: int, regionbody: int, bins: int,
                  chrom_prefix: str):
    '''
    计算甲基化水平覆盖度
    '''
    chrom, start, end, _, strand = gene_model.loc[gene_id]

    if end - start < bins:
        return

    if chrom in {'Pt', 'Mt', 'chrM', 'chrC'}:
        return None
    chrom = chrom_prefix + chrom

    site1 = get_target_site(site1, gene_id)  # site1
    site2 = get_target_site(site2, gene_id)  # site1

    bwfile = pyBigWig.open(infile)
    try:
        if strand == '+':
            methyratio_5 = bwfile.values(chrom, site1 - before, site1)
            cov_5 = get_bin_cov(methyratio_5, bins)
            methyratio_3 = bwfile.values(chrom, site2, site2 + after)
            cov_3 = get_bin_cov(methyratio_3, bins)

            # gene_body_region
            methyratio_gb = bwfile.values(chrom, site1, site2)
            methyratio_gb = scipy.ndimage.zoom(methyratio_gb,
                                               regionbody / len(methyratio_gb),
                                               order=0,
                                               mode='nearest')
            cov_gb = get_bin_cov(methyratio_gb, bins)

        else:
            methyratio_5 = bwfile.values(chrom, site2 + before, site1)[::-1]
            cov_5 = get_bin_cov(methyratio_5, bins)
            methyratio_3 = bwfile.values(chrom, site2, site2 - after)[::-1]
            cov_3 = get_bin_cov(methyratio_3, bins)

            # gene_body_region
            methyratio_gb = bwfile.values(chrom, site1, site2)[::-1]
            methyratio_gb = scipy.ndimage.zoom(methyratio_gb,
                                               regionbody / len(methyratio_gb),
                                               order=0,
                                               mode='nearest')
            cov_gb = get_bin_cov(methyratio_gb, bins)

    except RuntimeError:
        return
    
    cov = np.concatenate([cov_5, cov_gb, cov_3])

    return cov, gene_id


def get_meta_scale_result(infile,
                          gene_list,
                          site1,
                          site2,
                          before=1000,
                          after=1000,
                          regionbody=1000,
                          bins=100,
                          chrom_prefix='',
                          threads=64):
    '''
    infile = '/public/home/mowp/test/BS_Seq/bw_files/WT_BS_GSM1242401.methratio.cg.bw'
    bins = 100
    b = 2000
    a = 2000
    m = 2000
    cov = get_meta_scale_result(infile, all_gene, 'aTSS', 'PAS', bins=bins, before=b, after=a, regionbody=m)
    '''
    results = []
    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(gene_list) / threads)
        results = e.map(get_scale_cov,
                        repeat(infile),
                        gene_list,
                        repeat(site1),
                        repeat(site2),
                        repeat(before),
                        repeat(after),
                        repeat(regionbody),
                        repeat(bins),
                        repeat(chrom_prefix),
                        chunksize=chunksize)

    cov = []
    for res in results:
        if res is not None:
            cov_, gene_id = res
            cov.append(cov_)

    cov = np.nanmean(cov, axis=0)
    return cov

