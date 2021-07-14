#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2020-07-06 18:15:02
LastEditTime : 2021-02-02 21:33:04
Description  : 
    Before using this script, 
    You should run bedtools intersect to get the input file
    eg. bedtools intersect -abam elongating_data.full_len.bam -b ~/db/Arabidopsis_thaliana/representative_gene_model/araport11.representative.gene_model.bed -wo -s -split -bed > elongating_data.full_len.bed.intersect
'''


from collections import defaultdict, namedtuple
import pickle
import pandas as pd
import numpy as np
import click

import logging
logging.basicConfig(level=logging.DEBUG,  
                    format='%(asctime)s %(filename)s: %(message)s',  
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    )

NAMES=[
    'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 
    'ThickStart', 'ThickEnd', 'ItemRGB', 'BlockCount', 'BlockSizes', 'BlockStarts', 
    'geneChromosome', 'geneStart', 'geneEnd', 'geneName', 'geneScore', 'geneStrand', 
    'geneThickStart', 'geneThickEnd', 'geneItemRGB', 'geneBlockCount', 'geneBlockSizes',   'geneBlockStarts', 'cov'
    ]

USECOLS = [
    'Chromosome', 'Start', 'End', 'Name', 'Strand', 'BlockSizes', 'BlockStarts',
    'geneStart', 'geneEnd', 'geneName', 'geneBlockCount', 'geneBlockSizes', 'geneBlockStarts', 
    ]

Splice_stats = namedtuple('Splice_stats', 'intron_total_count intron_count unsplice_count unsplice_intron')


def get_gene_info(item):
    '''
    获取基因内含子外显子信息
    '''
    start = item.geneStart
    strand = item.Strand
    exon = []
    intron = []
    n = 0
    for blockstart, blocksize in zip(item.geneBlockStarts, item.geneBlockSizes):
        n += 1
        if strand == '+':
            exon.append((start+blockstart, start+blockstart+blocksize))
            if n > 1:
                intron.append((exon[-2][-1], exon[-1][0]))
        else:
            exon.append(((start+blockstart+blocksize)*-1, (start+blockstart)*-1))
            if n > 1:
                intron.append((exon[-1][-1], exon[-2][0]))
                
    if strand == '-':
        exon = exon[::-1]
        intron = intron[::-1]
    return exon, intron


def get_read_info(item):
    '''
    获取read blocks信息
    '''
    start = item.Start
    strand = item.Strand
    blocks = []
    n = 0
    for blockstart, blocksize in zip(item.BlockStarts, item.BlockSizes):
        n += 1
        if strand == '+':
            blocks.append((start+blockstart, start+blockstart+blocksize))
        else:
            blocks.append(((start+blockstart+blocksize)*-1, (start+blockstart)*-1))
    if strand == '-':
        blocks = blocks[::-1]
    return blocks


def is_overlap(block, intron):
    minn = max(block[0], intron[0])
    maxn = min(block[1], intron[1])
    # if (maxn-minn)/(intron[1]-intron[0]) >= .5:
    if maxn-minn >= 10:
        return True
    else:
        return False
    
    
def get_splice_stats(blocks, introns):
    '''
    计算reads中没有剪切的intron个数
    返回:
        没有被剪切的intron个数
        跨过的总intron个数
    '''
    block_start_intron = 99999
    block_end_intron = -1
    for i in range(len(introns)):
        # read的第一个block
        if blocks[0][0]+10 <= introns[i][1] and block_start_intron > i:
            block_start_intron = i  # read从那个intron开始
        # read的最后一个block
        if blocks[-1][1]-10 >= introns[i][0] and block_end_intron < i:
            block_end_intron = i  # read从那个intron结束

    if block_start_intron == 99999 or block_end_intron == -1:  # read not overlap with intron
        return
    span_intron = block_end_intron-block_start_intron+1

    block_total_count = len(blocks)
    blocks = iter(blocks)
    introns = iter(introns)
    
    unsplice_intron = []  # 用于储存未剪切intron编号
    
    block_count, unsplice_count, intron_count, intron_number = 0, 0, 0, 0
    block = next(blocks)
    block_count += 1
    intron = next(introns)
    intron_number += 1  # 当前intron序号

    while True:
        '''
        忽略read没有跨过intron的情况
        read: _________________________
        gene: ________....._______..........._______
        ''' 
        if block[1] > intron[1] and intron_count < intron_number:
            intron_count = intron_number  # read跨过了多少个intron

        # if not(block_count == block_total_count and block[1] < intron[1]) and is_overlap(block, intron):
        if is_overlap(block, intron):
            unsplice_count += 1
            unsplice_intron.append(str(intron_number))
        try:
            if intron[0] >= block[1]:
                block = next(blocks)
                block_count += 1
            elif intron[0] <= block[1]:
                intron = next(introns)
                intron_number += 1
        except StopIteration:
            break
    
    if span_intron < unsplice_count:
        return

    return unsplice_count, span_intron, ','.join(unsplice_intron)


@click.command()
@click.option('-i', '--infile', required=True)
@click.option('-o', '--outfile', required=True)
def main(infile, outfile):
    logging.info('Start load bedtools intersect result')
    df = pd.read_csv(
        infile, 
        sep='\t', 
        names=NAMES,
        usecols=USECOLS,
        header=None
        )
    df['geneBlockSizes'] = df['geneBlockSizes'].map(lambda x: np.fromstring(x, sep=',', dtype='int'))
    df['geneBlockStarts'] = df['geneBlockStarts'].map(lambda x: np.fromstring(x, sep=',', dtype='int'))
    df['BlockSizes'] = df['BlockSizes'].map(lambda x: np.fromstring(x, sep=',', dtype='int'))
    df['BlockStarts'] = df['BlockStarts'].map(lambda x: np.fromstring(x, sep=',', dtype='int'))
    logging.info('Load bedtools intersect done!')

    logging.info('Main start')
    splice_stats_dict = defaultdict(lambda: {}) # 储存结果用的，pkl格式
    gene_info = {} # key=gene_id, value = (exon, intron)

    n = 0
    with open(outfile, 'w') as o:
        o.write('gene_id\tintron_total_count\tread_id\tintron_count\tunsplice_count\tunsplice_intron\n')
        for item in df.itertuples():
            
            n += 1
            if n % 100000 == 0:
                logging.info(f'Process {n} reads')

            gene_id = item.geneName.split('.')[0]
            intron_total_count = item.geneBlockCount-1

            # 跳过那些intronless的基因
            if intron_total_count == 0:
                continue
            # 获取基因剪切情况
            if gene_id not in gene_info:
                gene_info[gene_id] = get_gene_info(item)
            introns = gene_info[gene_id][1]
            # 获取read剪切情况
            blocks = get_read_info(item)

            splice_stats = get_splice_stats(blocks, introns)
            if splice_stats is not None:
                unsplice_count, intron_count, unsplice_intron = splice_stats
            else:
                continue
            '''
            item.Name: read_id
            unsplice_count: 未剪切intron个数
            intron_count: read跨了多少个intron
            intron_total_count: 代表性转录本的总intron个数
            unsplice_intron: 未剪切转录本编号
            item.geneName: 转录本号
            '''
            o.write(f'{item.geneName}\t{intron_total_count}\t{item.Name}\t{intron_count}\t{unsplice_count}\t{unsplice_intron}\n')
            # 下游分析时候注意intron_count=0的reads，一般都要去除
            splice_stats_dict[item.Name][gene_id] = Splice_stats(intron_total_count, intron_count, unsplice_count, unsplice_intron)
    
    with open(outfile+'.pkl', 'wb') as o:
        pickle.dump(dict(splice_stats_dict), o)


if __name__ == "__main__":
    main()