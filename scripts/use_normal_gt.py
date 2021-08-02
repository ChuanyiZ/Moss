#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 chuanyi5 <chuanyi5@illinois.edu>
#
# Distributed under terms of the MIT license.

"""
Turn Strelka2, Mutect2 tumor VCFs into a single normal sample VCF
"""

import allel
import argparse
import numpy as np
import sys


def strelka2(inputs) -> dict:
    # chrom: pos: normal_genotype
    chrom_pos_gt = dict()
    cnt = 0
    for ifile in inputs:
        in_vcf = allel.read_vcf(ifile, fields='*')
        for chrom, pos, ref, alt, nt in zip(in_vcf["variants/CHROM"], in_vcf["variants/POS"], in_vcf["variants/REF"], in_vcf["variants/ALT"], in_vcf["variants/NT"]):
            chrom = str(chrom)
            pos = int(pos)
            alt = alt[0]

            if nt == 'ref':
                normal = ref + ref
            elif nt == 'het':
                normal = ref + alt
            elif nt == 'hom':
                normal = alt + alt
            else:
                continue

            if chrom in chrom_pos_gt:
                if pos in chrom_pos_gt[chrom]:
                    if chrom_pos_gt[chrom][pos][0] == chrom_pos_gt[chrom][pos][1] and normal[0] != normal[1]:
                        cnt += 1
                        chrom_pos_gt[chrom][pos] = normal
                else:
                    chrom_pos_gt[chrom][pos] = normal
            else:
                chrom_pos_gt[chrom] = {pos: normal}
    print(f"Disagreement on normal: {cnt} times.")
    return chrom_pos_gt


def mutect2(inputs, normal_name) -> dict:
    # chrom: pos: normal_genotype
    chrom_pos_gt = dict()
    cnt = 0
    cnt_het_hom = 0
    for ifile in inputs:
        # ["variants/CHROM", "variants/POS", "variants/REF", "variants/ALT", "calldata/GT"]
        in_vcf = allel.read_vcf(ifile, fields='*')
        idx_normal = np.argwhere(in_vcf["samples"] == normal_name)[0][0]
        zipped = zip(in_vcf["variants/CHROM"][in_vcf["variants/is_snp"]],
            in_vcf["variants/POS"][in_vcf["variants/is_snp"]],
            in_vcf["variants/REF"][in_vcf["variants/is_snp"]],
            in_vcf["variants/ALT"][in_vcf["variants/is_snp"]],
            in_vcf["calldata/GT"][in_vcf["variants/is_snp"]])
        for chrom, pos, ref, alt, gt in zipped:
            chrom = str(chrom)
            pos = int(pos)
            alt = alt[0]
            ref_alt = ref + alt
            normal = ref_alt[gt[idx_normal][0]] + ref_alt[gt[idx_normal][1]]
            if gt[idx_normal][0] != 0 or gt[idx_normal][1] != 0:
                cnt_het_hom += 1

            if chrom in chrom_pos_gt:
                if pos in chrom_pos_gt[chrom]:
                    if chrom_pos_gt[chrom][pos][0] == chrom_pos_gt[chrom][pos][1] and normal[0] != normal[1]:
                        cnt += 1
                        chrom_pos_gt[chrom][pos] = normal
                else:
                    chrom_pos_gt[chrom][pos] = normal
            else:
                chrom_pos_gt[chrom] = {pos: normal}
    print(f"Disagreement on normal: {cnt} times.")
    print(f"Not ref: {cnt_het_hom} times.")
    return chrom_pos_gt


def write_vcf(dict_chrom_pos_gt: dict, output):
    header = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\n"""
    with open(output, "w") as ofile:
        ofile.write(header)
        chunk = []
        cnt = 0
        for chrom, pos_gt in dict_chrom_pos_gt.items():
            for pos, gt in pos_gt.items():
                if gt[0] != gt[1]:
                    chunk.append(
                        f"{chrom}\t{pos}\t.\t{gt[0]}\t{gt[1]}\t.\t.\t.\tGT\t0|1\n")
                else:
                    chunk.append(
                        f"{chrom}\t{pos}\t.\t{gt[0]}\t.\t.\t.\t.\tGT\t0|0\n")
                cnt += 1
                if cnt % 1000 == 0:
                    ofile.writelines(chunk)
                    chunk = []
                    cnt = 0
        ofile.writelines(chunk)


def main(inputs, output, get_chrom_pos_gt):
    chrom_pos_gt = get_chrom_pos_gt(inputs)
    # ipdb.set_trace()
    write_vcf(chrom_pos_gt, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input tumor result VCF file", action="append")
    parser.add_argument("--normal-name", help="name of the normal sample in the VCF file, only used for Mutect")
    parser.add_argument("-f", "--input-files", help="input tumor result VCF file list")
    parser.add_argument("-t", "--tool", help="[M|m|Mutect] or [S|s|Strelka]", required=True)
    parser.add_argument("-o", "--output", help="output normal genotype VCF file")
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])



    if args.input_files is not None:
        with open(args.input_files) as ifile:
            if args.input is None:
                args.input = []
            for line in ifile:
                args.input.append(line.strip())
    
    
    if args.tool[0].lower() == 'm':
        if args.normal_name is None:
            exit()
        chrom_pos_gt = mutect2(args.input, args.normal_name)
    elif args.tool[0].lower() == 's':
        chrom_pos_gt = strelka2(args.input)
    write_vcf(chrom_pos_gt, args.output)
