#!/usr/bin/env python
# coding=utf-8
'''
Date         : 2021-08-24 21:05:33
LastEditTime : 2021-08-25 13:21:06
LastEditors  : windz
FilePath     : /pipelines/BS_seq_pipeline/script/get_bedgraph.py
'''

import gzip
import os
import click

@click.command()
@click.option('-i', '--infile', help='methratio.py output file.', required=True, type=click.Path(exists=True))
@click.option('-o', '--output_path', help='output filepath.', required=True)
@click.option('--chrom_prefix', help='chromosome name prefix', required=False, default='')
def main(infile: str, output_path: str, chrom_prefix: str):
    # infile = '/public/home/yuym/taiyi/data/BS_Seq/TAIR/public/WGBS_20210518/Methratio/SRX361939_methratio.txt.gz'
    # output_path = '/public/home/mowp/test/BS_Seq'

    filename = os.path.basename(infile).split('_methratio')[0]

    chh_out = open(f'{output_path}/{filename}.methratio.chh.bdg', 'w')
    chg_out = open(f'{output_path}/{filename}.methratio.chg.bdg', 'w')
    cg_out = open(f'{output_path}/{filename}.methratio.cg.bdg', 'w')

    outfile = {
        'CHH': chh_out,
        'CHG': chg_out,
        'CG': cg_out,
    }

    with gzip.open(infile, 'rt') as f:
        next(f)
        for line in f:
            line = line.rstrip().split('\t')
            chrom, pos, strand, context, ratio, eff_ct, c_counts, ct_count, *_ = line

            if chrom in {'chrM', 'chrC', 'Mt', 'Pt'}:
                continue

            pos = int(pos)
            ct_count = float(ct_count)
            c_counts = int(c_counts)

            if ct_count >= 4:
                ratio = round(c_counts/ct_count, 4)
                outfile[context].write(f'{chrom}\t{(pos-1)}\t{pos}\t{ratio}\n')

    chh_out.close()
    chg_out.close()
    cg_out.close()


if __name__ == "__main__":
    main()