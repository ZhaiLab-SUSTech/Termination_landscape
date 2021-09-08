#!/usr/bin/env python
# coding=utf-8
'''
Date         : 2021-04-19 09:45:44
LastEditTime : 2021-09-08 20:26:49
LastEditors  : windz
FilePath     : /public/home/mowp/test/nanopore_test/20210326_xrn2_test/script/add_tag_to_bam.py
'''


import pysam
import pandas as pd
import numpy as np
import click
import concurrent.futures


MIN_POLYA_LENGTH = 15
MAX_SPLICE_INTERMEDIATE_LEFT_LENGTH = -10
MAX_SPLICE_INTERMEDIATE_RIGHT_LENGTH = 10
READ_THROUGH_LENGTH = 10000
END3_ALIGNMENT_SCORE_LIMIT = (-5, 5)
PRIMER_TYPE = {'F-R', 'R-F', 'R-UF', 'UF-R'}


# get adapter_info
def extract_end3_score_value(df):
    df.eval('three_score = r_align_start-r_primer_start-genome_align_end', inplace=True)
    df.eval('five_score = genome_align_start-(f_align_end+f_primer_start)', inplace=True)
    is_correct_strand = df['rna_strand'] == df['read_align_strand']
    df['end3_values'] = df['five_score']
    df['end3_values'][is_correct_strand] = df['three_score'][is_correct_strand]


def get_adapter_info(infile):
    adapter_info = pd.read_csv(infile, sep='\t')
    adapter_info['read_id'] = adapter_info['read_core_id'].map(lambda x: x.split(',')[0])

    # remove multiple alignment reads
    adapter_info.drop_duplicates(['read_id'], inplace=True, keep=False)

    # only keep reads contained both adapters
    adapter_info.query('primer_type in @PRIMER_TYPE', inplace=True)

    # get 3' end score
    extract_end3_score_value(adapter_info)

    adapter_info_ = adapter_info[['read_id', 'polyA_type', 'end3_values']].copy()
    adapter_info_.set_index('read_id', inplace=True)
    adapter_info_ = adapter_info_.to_dict(orient='index')

    return adapter_info_

 
# get polya_info
def get_polya_info(infile):
    polya_info = pd.read_csv(infile, sep='\t')

    # remove multiple alignment reads
    polya_info['read_id'] = polya_info['read_core_id'].map(lambda x: x.split(',')[0])
    polya_info.drop_duplicates(['read_id'], inplace=True, keep=False)

    polya_info_ = polya_info[['read_id', 'polya_length']].copy()
    polya_info_.set_index('read_id', inplace=True)
    polya_info_ = polya_info_['polya_length'].to_dict()

    return polya_info_


# get read_info
def get_read_info(infile):
    read_info = pd.read_csv(infile, sep='\t', dtype={'retention_introns': str})
    # remove multiple alignment reads
    read_info['read_id'] = read_info['read_core_id'].map(lambda x: x.split(',')[0])
    read_info.drop_duplicates(['read_id'], inplace=True, keep=False)
    read_info['retention_introns'] = read_info['retention_introns'].fillna('')

    # find_splicing_intermediate
    read_info.eval(
        'end5ss_type = \
        (l_feature_type == "exon" and l_pos3 > @MAX_SPLICE_INTERMEDIATE_LEFT_LENGTH and l_feature_num != mRNA_intron_num+1) \
        or (l_feature_type != "exon" and l_feature_length+l_pos3 <= @MAX_SPLICE_INTERMEDIATE_RIGHT_LENGTH)',
        inplace=True
    )
    # get gene_id
    read_info['gene_id'] = read_info['mRNA'].map(lambda x: x.split('.')[0])
    read_info.set_index('read_id', inplace=True)

    read_info_ = read_info[['gene_id', 'end5ss_type', 'span_intron_num', 'retention_intron_num', 'retention_introns', 'mRNA_length']].copy()
    read_info_ = read_info_.to_dict(orient='index')

    return read_info_


def find_introns(r):
    BAM_CREF_SKIP = 3 #BAM_CREF_SKIP
    res = []
    match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
    base_position = r.reference_start
    for op, nt in r.cigartuples:
        if op in match_or_deletion:
            base_position += nt
        elif op == BAM_CREF_SKIP:
            junc_start = base_position
            base_position += nt
            res.append((junc_start, base_position))
    return res


def is_exceed_extend(read, gene_len):
    '''
    判断intron是否过长
    如果该intron比所在基因还长
    '''
    introns = np.array(find_introns(read))
    if len(introns) > 0:
        intron_len = introns[:, 1] - introns[:, 0]
        return (intron_len > gene_len).any()
    else:
        return False


@click.command()
@click.option('-i', '--infile', required=True)
@click.option('-o', '--outfile', required=True)
@click.option('--read_info', required=True)  
@click.option('--adapter_info', required=True)  
@click.option('--polya_info', required=True)
@click.option('-t', '--threads', default=10)
def main(infile, outfile, read_info, adapter_info, polya_info, threads):
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as e:
        read_info_ = e.submit(get_read_info, read_info)
        adapter_info_ = e.submit(get_adapter_info, adapter_info)
        polya_info_ = e.submit(get_polya_info, polya_info)
    
    read_info_ = read_info_.result()
    adapter_info_ = adapter_info_.result()
    polya_info_ = polya_info_.result()

    inbam = pysam.AlignmentFile(infile, 'rb') 
    outbam = pysam.AlignmentFile(outfile, 'wb', template=inbam)

    for read in inbam:
        if read.query_name in adapter_info_:
            polyA_type = adapter_info_[read.query_name]['polyA_type']
            end3_values = adapter_info_[read.query_name]['end3_values']
            polya_len = polya_info_[read.query_name]
            
            # get read_info
            if read.query_name in read_info_:
                gene_id = read_info_[read.query_name]['gene_id']
                is_splice_site = read_info_[read.query_name]['end5ss_type']
                span_intron_num = read_info_[read.query_name]['span_intron_num']
                retention_intron_num = read_info_[read.query_name]['retention_intron_num']
                retention_introns = read_info_[read.query_name]['retention_introns']
                mRNA_length = read_info_[read.query_name]['mRNA_length']
            else:
                # not overlap with annotated gene
                gene_id = 'None'
                is_splice_site = False
                span_intron_num = 0
                retention_intron_num = 0
                retention_introns = ''
                mRNA_length = 99999
                
            # remove read with low 3' end quality
            if polya_len < 15 and abs(end3_values) > 5:
                continue
            
            # remove splicing intermediate
            if is_splice_site and polya_len < 15:
                continue
            
            # 判断intron是否由于比对错误导致过长的情况
            if is_exceed_extend(read, mRNA_length):
                continue
            
            # adjust read strand
            if polyA_type == 'T':
                if read.is_reverse:
                    read.flag += -16
                else:
                    read.flag += 16
            
            read.set_tag('pa', polya_len)
            read.set_tag('gi', gene_id)
            read.set_tag('sn', span_intron_num)
            read.set_tag('rn', retention_intron_num)
            read.set_tag('ri', retention_introns)    
            
            outbam.write(read)


    inbam.close()
    outbam.close()

    # build bam index
    pysam.index(outfile, '-@ {threads}}')

if __name__ == "__main__":
    main()
    