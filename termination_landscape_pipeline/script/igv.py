#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2021-01-13 10:40:39
LastEditTime : 2021-05-27 22:44:11
Description  : 
'''


import matplotlib.patches as mp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import pysam
import seaborn as sns


######
# plot gene model
######

def plot_gene_model(ax, gene_models, fig_start, fig_end, gene_color='k', y_space=1):
    """plot gene model in the axis

    Args:
    -----
        ax (matplotlib.axes): An axis object to plot.

        gene_models (pd.DataFrame): A bed12 like DataFrame:

            chrom     start       end      gene_id score strand thickStart  thickEnd  
                2  17922018  17924542  AT2G43110.1     0      +   17922066  17924291   

            rgb blockCount              blockSizes                blockStarts y_pos  
              0          6  361,196,73,50,114,372,  0,527,823,1802,1952,2152,     0  

        fig_start (int): figure xlim start.

        fig_end (int): figure xlim end.

        gene_color (str, optional): gene color. Defaults to 'k'.

        y_space (int, optional): the spaces between gene in y direction. Defaults to 1.
    """
    ylim = 0  # ax ylim的下限
    height = 3 # gene model 高度
    y_space = y_space+height*2
    for gene_model in gene_models.values:
        chrom, chromStart, chromEnd, gene_id, _, strand, thickStart, thickEnd, _, blockCount, blockSizes, blockStarts, y_pos = gene_model
        y_pos = -y_space*y_pos
        ylim = min(y_pos, ylim)
        
        # 数据类型转化
        chromStart = int(chromStart)
        chromEnd = int(chromEnd)
        thickStart = int(thickStart)
        thickEnd = int(thickEnd)
        blockSizes = np.fromstring(blockSizes, sep=',', dtype='int')
        blockStarts = np.fromstring(blockStarts, sep=',', dtype='int')+chromStart

        # 画转录起始位点及方向箭头
        small_relative = .06*(fig_end-fig_start)  # 箭头突出部分相对长度
        arrowprops = dict(arrowstyle="-|>", connectionstyle="angle", color=gene_color)
        if strand == '+':
            ax.annotate('', xy=(chromStart+small_relative, height*2+y_pos), xytext=(chromStart, y_pos), arrowprops=arrowprops)
        else:
            ax.annotate('', xy=(chromEnd-small_relative, height*2+y_pos), xytext=(chromEnd, y_pos), arrowprops=arrowprops)
        
        line = mp.Rectangle((chromStart, y_pos-height/8), chromEnd-chromStart, height/4, color=gene_color, linewidth=0) # 基因有多长这条线就有多长
        ax.add_patch(line)

        for exonstart, size in zip(blockStarts, blockSizes):
            if exonstart == chromStart and exonstart+size == chromEnd:
                utr_size = thickStart-chromStart
                utr = mp.Rectangle((exonstart, y_pos-height/2), utr_size, height, color=gene_color, linewidth=0)
                ax.add_patch(utr)
                utr_size = chromEnd-thickEnd
                utr = mp.Rectangle((thickEnd, y_pos-height/2), utr_size, height, color=gene_color, linewidth=0)
                ax.add_patch(utr)
                exon = mp.Rectangle((thickStart, y_pos-height), thickEnd-thickStart, height*2, color=gene_color, linewidth=0)
                ax.add_patch(exon)

            elif exonstart + size <= thickStart:
                # 只有5'/ 3'UTR
                utr = mp.Rectangle((exonstart, y_pos-height/2), size, height, color=gene_color, linewidth=0)
                ax.add_patch(utr)

            elif exonstart < thickStart and exonstart + size > thickStart:
                # 带有5' / 3' UTR的exon
                utr_size = thickStart-exonstart
                utr = mp.Rectangle((exonstart, y_pos-height/2), utr_size, height, color=gene_color, linewidth=0)
                exon = mp.Rectangle((exonstart+utr_size, y_pos-height), size-utr_size, height*2, color=gene_color, linewidth=0)
                ax.add_patch(utr)
                ax.add_patch(exon)

            elif exonstart >= thickStart and exonstart + size <= thickEnd:
                # 普通exon
                exon = mp.Rectangle((exonstart, y_pos-height), size, height*2, color=gene_color, linewidth=0)
                ax.add_patch(exon)

            elif exonstart < thickEnd and exonstart + size > thickEnd:
                # 带有3' / 5' UTR的exon
                utr_size = exonstart + size - thickEnd
                utr = mp.Rectangle((thickEnd, y_pos-height/2), utr_size, height, color=gene_color, linewidth=0)
                exon = mp.Rectangle((exonstart, y_pos-height), size-utr_size, height*2, color=gene_color, linewidth=0)
                ax.add_patch(utr)
                ax.add_patch(exon)

            elif exonstart >= thickEnd:
                # 只有3'/ 5'UTR
                utr = mp.Rectangle((exonstart, y_pos-height/2), size, height, color=gene_color, linewidth=0)
                ax.add_patch(utr)


        ax.annotate(gene_id, xy=((chromStart+chromEnd)/2, y_pos+height*1.5), ha='center')
    
    # set ax
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_ticks_position('none')
    
    ax.set_ylim(ylim-height, height+y_space)
    ax.set_xlim(fig_start, fig_end)


def get_gene_model(chrom, start, end, bed_path):
    """Get gene model information from the bed file.
        The bed file should be indexed by tabix.
        For details, see http://www.htslib.org/doc/tabix.html

    Args:
        chrom (str): chromosome id.

        start (int): the start of the region.

        end (int): the end of the region.

        bed_path (str): the PATH of the bed file.

    Returns:
        pd.DataFrame: A bed12 dataframe.
            Examples
            --------
            chrom     start       end      gene_id score strand thickStart  thickEnd  
                2  17922018  17924542  AT2G43110.1     0      +   17922066  17924291   

            rgb blockCount              blockSizes                blockStarts
              0          6  361,196,73,50,114,372,  0,527,823,1802,1952,2152,
    """    
    tbx = pysam.TabixFile(bed_path)
    gene_models = pd.DataFrame(
        [gene_model.split('\t') for gene_model in tbx.fetch(chrom, start, end)],
        columns = ['chrom', 'start', 'end', 'gene_id', 'score', 'strand', 'thickStart', 'thickEnd', 'rgb', 'blockCount', 'blockSizes', 'blockStarts']
    )
    gene_models.sort_values(['start', 'end'], inplace=True)
    
    return gene_models
    

######
# other function
######

def is_overlap(gene_a, gene_b, threshold=0):
    """To judge whether two region is overlap.

    Args:
        gene_a (tuple): (chrom, start, end)
        gene_b (tuple): (chrom, start, end)
        threshold (int, optional): the minimum space between two region. Defaults to 0.

    Returns:
        bool: if overlap True, else False 
    """
    minn = max(int(gene_a[1]), int(gene_b[1]))
    maxn = min(int(gene_a[2]), int(gene_b[2]))

    if maxn-minn >= -threshold:
        return True
    else:
        return False

    
def get_y_pos_discontinous(df, gene_list=None):
    """Get the y position of each region. Save the results in the df['y_pos'].

    Args:
        df (pd.DataFrame): a DataFrame must include first four columns: 
            chrom, start, end, gene_id.

            Examples
            --------
            chrom    start      end    gene_id strand  
                4  9105672  9106504  AT4G16100      +   
            
        gene_list (set, optional): a set contain which gene to plot. 
            When gene_list is None, plot all item in the df. Defaults to None.

    Returns:
        int: item counts in the y position.
    """

    # only plot reads/gene_model in the gene_list
    if gene_list is not None:
        df.query('gene_id in @gene_list', inplace=True)
    
    df.reset_index(drop=True, inplace=True)
    df['y_pos'] = df.index
        
    return len(df)


def filter_bam(
    df, strand=None, 
    start_before=None, start_after=None,
    end_before=None, end_after=None
    ):
    if strand is not None:
        df.query('strand == @strand', inplace=True)
    if start_before is not None:
        df.query('start <= @start_before', inplace=True)
    if start_after is not None:
        df.query('start >= @start_after', inplace=True)
    if end_before is not None:
        df.query('end <= @end_before', inplace=True)
    if end_after is not None:
        df.query('end >= @end_after', inplace=True)


def get_y_pos_continuous(df, gene_list=None, threshold=0):
    """Get the y position of each region. Save the results in the df['y_pos'].

    Args:
        df (pd.DataFrame): a DataFrame must include first four columns: 
            chrom, start, end, gene_id.

            Examples
            --------
            chrom    start      end    gene_id strand  
                4  9105672  9106504  AT4G16100      +   
            
        gene_list (set, optional): a set contain which gene to plot. 
            When gene_list is None, plot all item in the df. Defaults to None.
        
        threshold (int, optional): the minimum space between two region. Defaults to 0.

    Returns:
        int: item counts in the y position.
    """
    # initialization of y_pos columns
    df['y_pos'] = None

    read_list = []
    for index_id, item in enumerate(df.values):
        chrom, start, end, gene_id, *_ = item
        
        # only plot reads/gene_model in the gene_list
        if gene_list is not None and gene_id not in gene_list:
            df.drop(index_id, inplace=True)
            continue
        
        current_read = (chrom, start, end)
        is_add_read = False
        y_pos = -1
        for y_pos in range(len(read_list)):
            y_pos_read = read_list[y_pos]
            if not is_overlap(current_read, y_pos_read, threshold=threshold):
                read_list[y_pos] = current_read
                df.at[index_id, 'y_pos'] = y_pos
                is_add_read = True
                break
        if not is_add_read:
            read_list.append(current_read)
            y_pos += 1
            df.at[index_id, 'y_pos'] = y_pos
            
    return len(read_list)


######
# plot bam
######

def find_exon(read):
    """Find exon position in the pysam.read.

    Args:
        read (pysam.read): <pysam.libcalignedsegment.AlignedSegment>

    Returns:
        zip(exon_start, exon_size): the exon position in the read.
            exon_start: exon start coordinate, 0-based
            exon_size: exon length

        Examples
        --------
        [(9105672, 457), (9106400, 142)]
    """    
    BAM_CREF_SKIP = 3 #BAM_CREF_SKIP
    blockStart = []
    blockSize = []
    match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
    exon_start = read.reference_start
    length = 0
    for op, nt in read.cigartuples:
        if op in match_or_deletion:
            length += nt
        elif op == BAM_CREF_SKIP:
            blockStart.append(exon_start)
            blockSize.append(length)
            exon_start += length+nt
            length = 0
    blockStart.append(exon_start)
    blockSize.append(length)
    return zip(blockStart, blockSize)


def convert_bam(
    chrom, start, end, strand, 
    infile, subsample=None,
    gene_list=None,
    method='continuous',
    filter_strand=None,
    start_before=None, start_after=None,
    end_before=None, end_after=None
    ):
    """Conver pysam.reads into a DataFrame.

    Args:
    -----

        chrom (str): chromosome id

        start (int): start position of the region

        end (int): end position of the region

        strand (str): '+' or '-'

        infile (str): the PATH of the bam file

        subsample (int or float, optional): <number>|<frac> of items from axis to return.
            Defaults to None.
        
        gene_list (set, optional): a set contain which gene to plot. 
            When gene_list is None, plot all item in the df. Defaults to None.
        
        method ('continuous' | '3_end' | '5_end')

    Returns:
        DataFrame: A dataframe contain bam reads information
            Examples
            --------
            Return: pd.DataFrame
            chrom    start      end    gene_id strand  
                4  9105672  9106504  AT4G16100      +   

                                         read_id  gap  polya_len  
            37c91155-0ef5-47be-9e1b-d3a6753982f4  0.0          0   

                        exon     y_pos  
            [(9105672, 832)]      None  
    """    

    bam_data = []
    with pysam.AlignmentFile(infile, 'rb') as inbam:
        for read in inbam.fetch(chrom, start, end):
            gap = 0
            try:
                polya_len = read.get_tag('pa')
                gene_id = read.get_tag('gi')
                span_intron_count = read.get_tag('sn')  # span_intron_num
                unsplice_count = read.get_tag('rn')  # retention_intron_num
                unsplice_intron = read.get_tag('ri')  # retention_introns
                # only plot reads in gene_list
            except KeyError:
                # 如果bam文件没有加tag 则默认为0
                polya_len, gap, gene_id = 0, 0, 0
                span_intron_count, unsplice_count, unsplice_intron = None, None, None
            if read.is_supplementary or read.is_unmapped:
                continue
            
            read_strand = '-' if read.is_reverse else '+'
            exon = find_exon(read)

            bam_data.append(
                (read.reference_name, read.reference_start, read.reference_end,
                 gene_id, read_strand, read.query_name, gap, polya_len, list(exon), 
                 span_intron_count, unsplice_count, unsplice_intron)
            )
    bam_data = pd.DataFrame(
        bam_data, 
        columns=[
            'chrom', 'start', 'end', 'gene_id', 'strand', 'read_id', 'gap', 
            'polya_len', 'exon', 'span_intron_count', 'unsplice_count', 'unsplice_intron'
            ]
    )
    
    # 进行subsample
    if subsample is not None:
        if type(subsample) is int:
            bam_data = bam_data.sample(n=subsample, random_state=42)
        elif type(subsample) is float and subsample < 1:
            bam_data = bam_data.sample(frac=subsample, random_state=42)
    
    if method == 'continuous':
        # 类似于igv那样连续的排
        # add polya
        bam_data['polya_len'] = bam_data['polya_len'].map(lambda x: int(x) if x >= 15 else 0)
        bam_data = bam_data.assign(
            start = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['start'], bam_data['start']-bam_data['polya_len']], -1),
            end = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['end']+bam_data['polya_len'], bam_data['end']], -1),
        )
        if strand == '+':
            bam_data.sort_values(['start', 'end'], inplace=True)
        else:
            bam_data.sort_values(['end', 'start'], ascending=False, inplace=True)
        bam_data.reset_index(drop=True, inplace=True)
        get_y_pos_continuous(bam_data, gene_list=gene_list)
    
    elif method == '5_end':
        if strand == '+':
            bam_data.sort_values(['start', 'end'], inplace=True)
        else:
            bam_data.sort_values(['end', 'start'], ascending=False, inplace=True)
        
        # add polya
        bam_data['polya_len'] = bam_data['polya_len'].map(lambda x: int(x) if x >= 15 else 0)
        bam_data = bam_data.assign(
            start = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['start'], bam_data['start']-bam_data['polya_len']], -1),
            end = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['end']+bam_data['polya_len'], bam_data['end']], -1),
        )
        
        filter_bam(bam_data, strand=filter_strand, start_before=start_before, start_after=start_after, end_before=end_before, end_after=end_after)
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == '3_end':
        if strand == '+':
            bam_data.sort_values(['end', 'start'], inplace=True)
        else:
            bam_data.sort_values(['start', 'end'], ascending=False, inplace=True)

        # add polya
        bam_data['polya_len'] = bam_data['polya_len'].map(lambda x: int(x) if x >= 15 else 0)
        bam_data = bam_data.assign(
            start = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['start'], bam_data['start']-bam_data['polya_len']], -1),
            end = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['end']+bam_data['polya_len'], bam_data['end']], -1),
        )
        
        filter_bam(bam_data, strand=filter_strand, start_before=start_before, start_after=start_after, end_before=end_before, end_after=end_after)
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == 'gene_id':  # 按照gene_id 3' end排序
        if strand == '+':
            bam_data.sort_values(['gene_id', 'end', 'start'], inplace=True)
        else:
            bam_data.sort_values(['gene_id', 'start', 'end'], ascending=False, inplace=True)
        
        # add polya
        bam_data['polya_len'] = bam_data['polya_len'].map(lambda x: int(x) if x >= 15 else 0)
        bam_data = bam_data.assign(
            start = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['start'], bam_data['start']-bam_data['polya_len']], -1),
            end = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['end']+bam_data['polya_len'], bam_data['end']], -1),
        )
        get_y_pos_discontinous(bam_data, gene_list=gene_list)
    
    elif method == 'spliced':  # 只画完全剪切的reads
        bam_data.query('span_intron_count > 0 and unsplice_count == 0', inplace=True)
        if strand == '+':
            bam_data.sort_values(['end', 'start'], inplace=True)
        else:
            bam_data.sort_values(['start', 'end'], ascending=False, inplace=True)
        
        # add polya
        bam_data['polya_len'] = bam_data['polya_len'].map(lambda x: int(x) if x >= 15 else 0)
        bam_data = bam_data.assign(
            start = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['start'], bam_data['start']-bam_data['polya_len']], -1),
            end = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['end']+bam_data['polya_len'], bam_data['end']], -1),
        )
        filter_bam(bam_data, strand=filter_strand, start_before=start_before, start_after=start_after, end_before=end_before, end_after=end_after)
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == 'partially_spliced':  # 只画不完全剪切的reads
        bam_data.query('span_intron_count > 0 and unsplice_count > 0 and unsplice_count != span_intron_count ', inplace=True)
        if strand == '+':
            bam_data.sort_values(['end', 'start'], inplace=True)
        else:
            bam_data.sort_values(['start', 'end'], ascending=False, inplace=True)
        
        # add polya
        bam_data['polya_len'] = bam_data['polya_len'].map(lambda x: int(x) if x >= 15 else 0)
        bam_data = bam_data.assign(
            start = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['start'], bam_data['start']-bam_data['polya_len']], -1),
            end = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['end']+bam_data['polya_len'], bam_data['end']], -1),
        )
        filter_bam(bam_data, strand=filter_strand, start_before=start_before, start_after=start_after, end_before=end_before, end_after=end_after)
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == 'unspliced':  # 只画不完全剪切的reads
        bam_data.query('span_intron_count > 0 and unsplice_count == span_intron_count', inplace=True)
        if strand == '+':
            bam_data.sort_values(['end', 'start'], inplace=True)
        else:
            bam_data.sort_values(['start', 'end'], ascending=False, inplace=True)
        
        # add polya
        bam_data['polya_len'] = bam_data['polya_len'].map(lambda x: int(x) if x >= 15 else 0)
        bam_data = bam_data.assign(
            start = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['start'], bam_data['start']-bam_data['polya_len']], -1),
            end = np.select([bam_data['strand'] == '+', bam_data['strand'] == '-'], [bam_data['end']+bam_data['polya_len'], bam_data['end']], -1),
        )
        filter_bam(bam_data, strand=filter_strand, start_before=start_before, start_after=start_after, end_before=end_before, end_after=end_after)
        get_y_pos_discontinous(bam_data, gene_list=gene_list)


    return bam_data


def plot_bam(
    ax, bam_data, 
    fig_start, fig_end, 
    read_color='#5D93C4', 
    polya_color='lightcoral',
    y_space=1.5, # the space between reads in yaxis
    gene_list=None,
    pal=None
):
    """plot bam information to the ax

    Args:
    -----
        ax (matplotlib.axes): An axis object to plot
        
        bam_data (pd.DataFrame): A dataframe contain bam reads information
            Examples
            --------
            Return: pd.DataFrame
            chrom    start      end    gene_id strand  
                4  9105672  9106504  AT4G16100      +   

                                         read_id  gap  polya_len  
            37c91155-0ef5-47be-9e1b-d3a6753982f4  0.0          0   

                        exon     y_pos  
            [(9105672, 832)]         0  

        fig_start (int): figure xlim start.

        fig_end (int): figure xlim end.

        read_color (str, optional): read color. Defaults to '#5D93C4'.

        polya_color (str, optional): polya color. Defaults to 'lightcoral'.

        y_space (int, optional): the spaces between reads in y direction. Defaults to 1.5.

        pal (seaborn.palettes._ColorPalette, optional)
    """    
    if pal is None:
        pal = sns.color_palette("Paired")

    gene_color = {}
    gene_color_index = 0
    
    ylim = 0 # the start position of yaxis
    height = .5 # reads的高度
    for item in bam_data.values:
        chrom, start, end, gene_id, strand, read_id, gap, polya_len, exon, *_, ypos,  = item
        ypos = -y_space*ypos
        ylim = min(ypos, ylim)
        line = mp.Rectangle((start, ypos-height/4), end-start, height/2, color='#A6A6A6', linewidth=0)
        ax.add_patch(line)
        for block_start, block_size in exon:
            
            # set gene color
            # index为色板中的序列编号
            # 如果gene_list存在 则给gene_list里面都基因上色
            if gene_list is not None:
                if gene_id in gene_list:
                    if gene_id not in gene_color:
                        gene_color[gene_id] = gene_color_index
                        gene_color_index += 1
                    read_color_ = pal[gene_color[gene_id]]
                else:
                    # 不在gene_list里面的reads都设置成灰色
                    read_color_ = '#5D5D5D'
            # 如果没有则统一颜色
            else:
                read_color_ = read_color
                
                # TODO: 不同链的reads不同颜色
            exon = mp.Rectangle((block_start, ypos-height), block_size, height*2, color=read_color_, linewidth=0)
            ax.add_patch(exon)

        # plot polya
        if polya_len > 15:
            if strand == '+':
                polya_tail = mp.Rectangle((block_start+block_size, ypos-height), polya_len, height*2, color=polya_color, linewidth=0)
            else:
                polya_tail = mp.Rectangle((start, ypos-height), polya_len, height*2, color=polya_color, linewidth=0)        
            ax.add_patch(polya_tail)

    ax.set_ylim(ylim*1.1, height+y_space)
    #ax.set_xlim(fig_start, fig_end)


def set_ax(ax, plot_xaxis=False):
    """Set axes

    Args:
        ax (matplotlib.axes): An axis object to plot

        plot_xaxis (bool, optional): if True plot x_ticks and x_locator. 
            Defaults to False.
    """    
    ax.spines['right'].set_visible(False) # 去线
    ax.spines['left'].set_visible(False) # 去线
    ax.spines['top'].set_visible(False) # 去线
    ax.yaxis.set_major_locator(ticker.NullLocator()) # 去y数字

    if not plot_xaxis:
        ax.xaxis.set_major_locator(ticker.NullLocator()) # 去x数字
        ax.xaxis.set_ticks_position('none') # 去x刻度


######
# main Class
######

class IGV:
    """A IGV class

    Usage:
    ------
        chrom, start, end, strand = '5', 5371627, 5375616, '+'

        igv_plot = igv.IGV(chrom, start, end, strand=strand)

        araport11_isoform_path = '/public/home/mowp/db/Arabidopsis_thaliana/representative_gene_model/araport11.representative.gene_model.bed.gz'
        igv_plot.add_gene_model(
            araport11_isoform_path,
            gene_list = {'AT5G16440.1', 'AT5G16450.1'},
        )

        infile = '/public/home/mowp/workspace/termination/cbRNA_pool/elongating_data/elongating_data.bam'
        igv_plot.add_bam(
            infile,
            gene_list = {'AT5G16440', 'AT5G16450'},
        )

        igv_plot.plot(height=4, width=8)
    """    
    def __init__(self, chrom, start, end, strand='+'):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.bam_list = []
        
    
    def add_bam(
        self, bam_path, gene_list=None, subsample=None, 
        method='continuous',
        filter_strand=None,
        start_before=None, start_after=None,
        end_before=None, end_after=None
    ):
        bam_data = convert_bam(
            self.chrom, self.start, self.end, self.strand, 
            bam_path, subsample=subsample, gene_list=gene_list, 
            method=method, 
            filter_strand=filter_strand,
            start_before=start_before, start_after=start_after, 
            end_before=end_before, end_after=end_after
            )
        self.bam_list.append(bam_data)
        self.gene_list = gene_list
    

    def add_gene_model(self, anno, gene_list=None):
        self.anno = anno
        self.gene_models = get_gene_model(self.chrom, self.start, self.end, self.anno)
        self.gene_model_ylim = get_y_pos_continuous(self.gene_models, gene_list=gene_list, threshold=8)  # 不重叠y轴的基因数目
        
        
    def plot(
        self, 
        height=5, width=10, 
        gene_track_height=.5,
        bam_track_height=4,
        gene_color='k',
        extend_xlim_start = False,
        extend_xlim_end = False,
        polya_site=None,
    ):
                
        nrows = len(self.bam_list)+1  # bam_list track + gene model track
        height_ratios = [self.gene_model_ylim*gene_track_height]  # 设置基因track的高度
        height_ratios.extend([bam_track_height]*len(self.bam_list))  # 设置bam track的高度
        
        fig, ax = plt.subplots(
            nrows=nrows, 
            gridspec_kw={'height_ratios': height_ratios},
            figsize=(width, height),
            sharex = True,
        )

        # plot gene_model
        i = 0
        if len(self.bam_list) == 0:
            plot_gene_model(ax, self.gene_models, self.start, self.end, gene_color=gene_color, y_space=6)
        else:
            plot_gene_model(ax[0], self.gene_models, self.start, self.end, gene_color=gene_color, y_space=6)
        
        # plot bam files
        for i, bam_data in enumerate(self.bam_list, 1):
            plot_bam(ax[i], bam_data, self.start, self.end, gene_list=self.gene_list)
            plot_xaxis = i==nrows-1
            set_ax(ax[i], plot_xaxis=plot_xaxis)

            # add polya site
            if polya_site is not None:
                tbx = pysam.TabixFile(polya_site)
                for item in tbx.fetch(self.chrom, self.start, self.end):
                    chrom_, start_, end_, gene_id_, _, strand_, _ = item.split('\t') 
                    gene_id_ = gene_id_.split('_')[0]
                    if self.gene_list is not None:
                        if gene_id_ in self.gene_list:
                            ax[i].axvline(int(end_), ls='--', color='#555555')
                    else:
                        ax[i].axvline(int(end_), ls='--', color='#555555')

        
        # set last xaxis
        maxn = self.end
        minn = self.start
        for bam_data in self.bam_list:
            if len(bam_data) == 0:
                continue
            
            minn_ = min(bam_data['start'])
            maxn_ = max(bam_data['end'])

            if extend_xlim_start:
                if self.strand == '+':
                    minn = min(minn_, minn)
                else:
                    maxn = max(maxn_, maxn)

            if extend_xlim_end:
                if self.strand == '+':
                    maxn = max(maxn_, maxn)
                else:
                    minn = min(minn_, minn)

        if maxn-minn > 400:
            step = (maxn-minn)//400*100  # 坐标轴步长
        else:
            step = (maxn-minn)//40*10 
        if i == 0:
            ax_ = ax
        else:
            ax_ = ax[i]
        if self.strand == '+':
            xticks = np.arange(minn, maxn+step, step)
            ax_.set_xticks(xticks)
            ax_.set_xticklabels(xticks-minn)
            ax_.set_xlim(minn, maxn)
        else:
            xticks = np.arange(maxn, minn-step, -step)
            ax_.set_xticks(xticks)
            ax_.set_xticklabels(maxn-xticks)
            ax_.set_xlim(minn, maxn)
            ax_.invert_xaxis()
            
        return ax