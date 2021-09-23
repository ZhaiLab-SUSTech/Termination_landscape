#!/usr/bin/env python
# coding=utf-8
'''
Date         : 2021-08-30 10:10:24
LastEditTime : 2021-08-31 16:36:02
LastEditors  : windz
FilePath     : /tools/pyIGV/coverage_plot.py
'''


import os
import matplotlib.patches as mp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import pysam
import pyBigWig
from typing import Tuple, Optional


def plot_gene_model(
    ax,
    gene_models: pd.DataFrame,
    fig_start: int,
    fig_end: int,
    gene_color: str = "k",
    y_space: int = 1,
):
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
    height = 3  # gene model 高度
    y_space = y_space + height * 2
    for gene_model in gene_models.values:
        (
            chrom,
            chromStart,
            chromEnd,
            gene_id,
            _,
            strand,
            thickStart,
            thickEnd,
            _,
            blockCount,
            blockSizes,
            blockStarts,
            y_pos,
        ) = gene_model
        y_pos = -y_space * y_pos
        ylim = min(y_pos, ylim)

        # 数据类型转化
        chromStart = int(chromStart)
        chromEnd = int(chromEnd)
        thickStart = int(thickStart)
        thickEnd = int(thickEnd)
        blockSizes = np.fromstring(blockSizes, sep=",", dtype="int")
        blockStarts = np.fromstring(blockStarts, sep=",", dtype="int") + chromStart

        # 画转录起始位点及方向箭头
        small_relative = 0.06 * (fig_end - fig_start)  # 箭头突出部分相对长度
        arrowprops = dict(arrowstyle="-|>", connectionstyle="angle", color=gene_color)
        if strand == "+":
            ax.annotate(
                "",
                xy=(chromStart + small_relative, height * 2 + y_pos),
                xytext=(chromStart, y_pos),
                arrowprops=arrowprops,
            )
        else:
            ax.annotate(
                "",
                xy=(chromEnd - small_relative, height * 2 + y_pos),
                xytext=(chromEnd, y_pos),
                arrowprops=arrowprops,
            )

        line = mp.Rectangle(
            (chromStart, y_pos - height / 8),
            chromEnd - chromStart,
            height / 4,
            color=gene_color,
            linewidth=0,
        )  # 基因有多长这条线就有多长
        ax.add_patch(line)

        for exonstart, size in zip(blockStarts, blockSizes):
            if exonstart == chromStart and exonstart + size == chromEnd:
                utr_size = thickStart - chromStart
                utr = mp.Rectangle(
                    (exonstart, y_pos - height / 2),
                    utr_size,
                    height,
                    color=gene_color,
                    linewidth=0,
                )
                ax.add_patch(utr)
                utr_size = chromEnd - thickEnd
                utr = mp.Rectangle(
                    (thickEnd, y_pos - height / 2),
                    utr_size,
                    height,
                    color=gene_color,
                    linewidth=0,
                )
                ax.add_patch(utr)
                exon = mp.Rectangle(
                    (thickStart, y_pos - height),
                    thickEnd - thickStart,
                    height * 2,
                    color=gene_color,
                    linewidth=0,
                )
                ax.add_patch(exon)

            elif exonstart + size <= thickStart:
                # 只有5'/ 3'UTR
                utr = mp.Rectangle(
                    (exonstart, y_pos - height / 2),
                    size,
                    height,
                    color=gene_color,
                    linewidth=0,
                )
                ax.add_patch(utr)

            elif exonstart < thickStart and exonstart + size > thickStart:
                # 带有5' / 3' UTR的exon
                utr_size = thickStart - exonstart
                utr = mp.Rectangle(
                    (exonstart, y_pos - height / 2),
                    utr_size,
                    height,
                    color=gene_color,
                    linewidth=0,
                )
                exon = mp.Rectangle(
                    (exonstart + utr_size, y_pos - height),
                    size - utr_size,
                    height * 2,
                    color=gene_color,
                    linewidth=0,
                )
                ax.add_patch(utr)
                ax.add_patch(exon)

            elif exonstart >= thickStart and exonstart + size <= thickEnd:
                # 普通exon
                exon = mp.Rectangle(
                    (exonstart, y_pos - height),
                    size,
                    height * 2,
                    color=gene_color,
                    linewidth=0,
                )
                ax.add_patch(exon)

            elif exonstart < thickEnd and exonstart + size > thickEnd:
                # 带有3' / 5' UTR的exon
                utr_size = exonstart + size - thickEnd
                utr = mp.Rectangle(
                    (thickEnd, y_pos - height / 2),
                    utr_size,
                    height,
                    color=gene_color,
                    linewidth=0,
                )
                exon = mp.Rectangle(
                    (exonstart, y_pos - height),
                    size - utr_size,
                    height * 2,
                    color=gene_color,
                    linewidth=0,
                )
                ax.add_patch(utr)
                ax.add_patch(exon)

            elif exonstart >= thickEnd:
                # 只有3'/ 5'UTR
                utr = mp.Rectangle(
                    (exonstart, y_pos - height / 2),
                    size,
                    height,
                    color=gene_color,
                    linewidth=0,
                )
                ax.add_patch(utr)

        ax.annotate(
            gene_id, xy=((chromStart + chromEnd) / 2, y_pos + height * 1.5), ha="center"
        )

    # set ax
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_ticks_position("none")

    ax.set_ylim(ylim - height, height + y_space)
    ax.set_xlim(fig_start, fig_end)


def get_gene_model(chrom: str, start: int, end: int, bed_path: str, lambda_for_parse_geneId:Optional[str] = None) -> pd.DataFrame:
    """Get gene model information from the bed file.
        The bed file should be indexed by tabix.
        For details, see http://www.htslib.org/doc/tabix.html

    Args:
        chrom (str): chromosome id.

        start (int): the start of the region.

        end (int): the end of the region.

        bed_path (str): the PATH of the bed file.

        lambda_for_parse_geneId(str) : lambda expression for parse geneId. Defalut is x: x

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
        [gene_model.split("\t") for gene_model in tbx.fetch(chrom, start, end)],
        columns=[
            "chrom",
            "start",
            "end",
            "gene_id",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "rgb",
            "blockCount",
            "blockSizes",
            "blockStarts",
        ],
    )
    gene_models.sort_values(["start", "end"], inplace=True)
    if lambda_for_parse_geneId:
        lambda_for_parse_geneId_fc = eval(f"lambda {lambda_for_parse_geneId}")
        gene_models = gene_models.assign(gene_id = lambda df:df['gene_id'].map(lambda_for_parse_geneId_fc))
    return gene_models


def is_overlap(
    gene_a: Tuple[str, int, int], gene_b: Tuple[str, int, int], threshold: int = 0
) -> bool:
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

    if maxn - minn >= -threshold:
        return True
    else:
        return False


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
    df["y_pos"] = None

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
                df.at[index_id, "y_pos"] = y_pos
                is_add_read = True
                break
        if not is_add_read:
            read_list.append(current_read)
            y_pos += 1
            df.at[index_id, "y_pos"] = y_pos

    return len(read_list)


def set_ax(ax, plot_xaxis=False, data_range=None):
    """Set axes

    Args:
        ax (matplotlib.axes): An axis object to plot

        plot_xaxis (bool, optional): if True plot x_ticks and x_locator.
            Defaults to False.
    """
    ax.spines["right"].set_visible(False)  # 去线
    ax.spines["left"].set_visible(False)  # 去线
    ax.spines["top"].set_visible(False)  # 去线
    ax.yaxis.set_major_locator(ticker.NullLocator())  # 去y数字

    if not plot_xaxis:
        ax.xaxis.set_major_locator(ticker.NullLocator())  # 去x数字
        ax.xaxis.set_ticks_position("none")  # 去x刻度

    if data_range is not None:
        ax.set_ylim(data_range)


def plot_track(ax, track_data, start, end, track_color):
    ax.fill_between(np.linspace(start, end, end-start), track_data, color=track_color)


def find_exons(r):
    BAM_CREF_SKIP = 3 #BAM_CREF_SKIP
    res = []
    match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
    exon_start = r.reference_start
    length = 0
    for op, nt in r.cigartuples:
        if op in match_or_deletion:
            length += nt
        elif op == BAM_CREF_SKIP:
            res.append(np.array([exon_start, exon_start+length]))
            exon_start = res[-1][1]+nt
            length = 0
    res.append(np.array([exon_start, exon_start+length]))
    return res


def get_cov_from_bam(chrom, start, end, infile, gene_list=None, chrom_prefix=''):
    cov = np.zeros(end-start)
    chrom = chrom_prefix+chrom
    with pysam.AlignmentFile(infile, 'rb') as inbam:
        for read in inbam.fetch(chrom, start, end):
            if read.is_supplementary:
                continue

            if gene_list is not None:
                read_gene_id = read.get_tag('gi')
                if read_gene_id not in gene_list:
                    continue

            exons = find_exons(read)
            for exon in exons:
                exon = exon-start
                # 只取 self.start, self.end范围内的reads的coverage
                exon[np.where(exon > end-start)] = end-start
                exon[np.where(exon < 0)] = 0

                cov[exon[0]: exon[1]] += 1
        return cov
        

class COV:
    """A IGV class

    Usage:
    ------
        chrom, start, end, strand, force_tag_check = '5', 5371627, 5375616, '+', True

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
            read_color='r',
            bam_title='one bam'
        )

        igv_plot.plot(height=4, width=8)
    ------
    force_tag_check:
        if True, reads which don't have required tags will be ignored
    tag_required_dict:Mapping[str, str]:
        tag required. Key is the variable name, and value is corresponding tag.
        if force_tag_check is True, all tag should be stored in bam; if not, will assign values to variables in order, and when error is triggered, left variables will be ignored.
    """

    def __init__(
        self,
        chrom,
        start,
        end,
        strand="+",
    ):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.track_list = []
        self.color_list = []
        self.bamTitleList = []
        self.track_range = []


    def add_bam(
        self,
        bam_path,
        gene_list=None,
        color='#5D93C4',
        bam_title='',
        chrom_prefix='',
        data_range : list = None,
        log_scale = False,
    ):
        bam_data = get_cov_from_bam(
            self.chrom,
            self.start,
            self.end,
            bam_path,
            gene_list=gene_list,
            chrom_prefix=chrom_prefix,
        )
        if log_scale:
            bam_data = np.log2(bam_data+1)
        self.track_list.append(bam_data)
        self.color_list.append(color)
        self.bamTitleList.append(bam_title)
        self.track_range.append(data_range)


    def add_bw(
        self,
        bw_path,
        gene_list=None,
        color='#5D93C4',
        bam_title='',
        chrom_prefix='',
        data_range : list = None,
    ):
        bw = pyBigWig.open(bw_path)
        chrom = chrom_prefix+self.chrom
        bw_data = bw.values(chrom, self.start, self.end)
        bw_data = np.nan_to_num(bw_data)

        self.track_list.append(bw_data)
        self.color_list.append(color)
        self.bamTitleList.append(bam_title)
        self.track_range.append(data_range)


    def add_gene_model(self, anno:str, gene_list=None, auto_create_gz=False, lambda_for_parse_geneId:Optional[str]=None):
        """
        auto_create_gz: if provided bed is not compressed, the bgzip will be used to created the gz format file, corresponding index will also be built.
        lambda_for_parse_geneId: If not None, it will be used to parse the `gene_id` column in bed. None is equal to `x:x`
        """
        if anno.endswith('.gz'):
            self.anno = anno
        else:
            self.anno = anno + '.gz'
            if not os.path.exists(self.anno):
                assert auto_create_gz, "WRONG suffix, and <auto_create_gz> is False"
                print(f"WARNING: Try to create file {self.anno} and build index")
                os.system(f'sort -k1,1 -k2,2n {anno} > temp____ && bgzip -c  temp____ > {self.anno} && tabix -p bed {self.anno}')
        self.gene_models = get_gene_model(self.chrom, self.start, self.end, self.anno, lambda_for_parse_geneId=lambda_for_parse_geneId)
        self.gene_model_ylim = get_y_pos_continuous(
            self.gene_models, gene_list=gene_list, threshold=8
        )  # 不重叠y轴的基因数目


    def plot(
        self,
        height=5,
        width=10,
        gene_track_height=1,
        track_height=4,
        gene_color="k",
    ):

        nrows = len(self.track_list) + 1  # bam_list track + gene model track
        height_ratios = [gene_track_height]  # 设置基因track的高度
        height_ratios.extend([track_height] * len(self.track_list))  # 设置bam track的高度

        fig, ax = plt.subplots(
            nrows=nrows,
            gridspec_kw={"height_ratios": height_ratios},
            figsize=(width, height),
            sharex=True,
        )

        # plot gene_model
        i = 0
        if len(self.track_list) == 0:
            plot_gene_model(
                ax,
                self.gene_models,
                self.start,
                self.end,
                gene_color=gene_color,
                y_space=6,
            )
        else:
            plot_gene_model(
                ax[0],
                self.gene_models,
                self.start,
                self.end,
                gene_color=gene_color,
                y_space=6,
            )

        # plot bam files
        for i, (track_data, track_color, self.bamTitle, data_range) in enumerate(zip(self.track_list, self.color_list, self.bamTitleList, self.track_range), 1):
            plot_track(ax[i], track_data, self.start, self.end, track_color=track_color)
            plot_xaxis = i == nrows - 1
            set_ax(ax[i], plot_xaxis=plot_xaxis, data_range=data_range)


        # set last xaxis
        maxn = self.end
        minn = self.start
        for bam_data in self.track_list:
            if len(bam_data) == 0:
                continue

        if maxn - minn > 400:
            step = (maxn - minn) // 400 * 100  # 坐标轴步长
        else:
            step = (maxn - minn) // 40 * 10
        if i == 0:
            ax_ = ax
        else:
            ax_ = ax[i]
        if self.strand == "+":
            xticks = np.arange(minn, maxn + step, step)
            ax_.set_xticks(xticks)
            ax_.set_xticklabels(xticks - minn)
            ax_.set_xlim(minn, maxn)
        else:
            xticks = np.arange(maxn, minn - step, -step)
            ax_.set_xticks(xticks)
            ax_.set_xticklabels(maxn - xticks)
            ax_.set_xlim(minn, maxn)
            ax_.invert_xaxis()
        return ax
