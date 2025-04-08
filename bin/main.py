#!/usr/bin/env python
from maf import filter_maf
from nr import process_noredundant
from pre import process_snp_dists
from pmrcal import *
from wtdenoise import *
import argparse 
import os
import polars as pl
parser = argparse.ArgumentParser(description='Calculate recombination along chromsome')
parser.add_argument('-m', '--maf', type=float, metavar='', help='max major allele frequency')
parser.add_argument('-n', '--nonredundant', type=float, metavar='', help='criterion for nonredundant strains')
parser.add_argument('-cc', '--cold_criterion', type=float, metavar='', help='criterion for recombination cold defult -0.4')
parser.add_argument('-hc', '--hot_criterion', type=float, metavar='', help='criterion for recombination hot defult 0.6')
parser.add_argument('-i','--input', metavar='', type=str,required=True, help='input path snp matrix ')
parser.add_argument('-o','--output', metavar='', type=str,required=True, help='output path')
args = parser.parse_args()
i=0
pl.enable_string_cache()
if args.cold_criterion:
    cc=args.cold_criterion
else:
    cc=-0.4
breakpoint()
if args.hot_criterion:
    hc=args.hot_criterion
else:
    hc=0.6
if args.nonredundant:
    path=process_noredundant(args.input, args.output,args.nonredundant)
    i=1
if args.maf:
    if i==0:
        path=filter_maf(args.input, args.maf, args.output)
        i=1
    else:
        path=filter_maf(path, args.maf, path)
if i==0:
    process_snp_dists(args.input,args.output)
    process_snp_data(args.input,args.output)
    plot_and_save(get_crange, get_hrange,wavelet_transform, args.output,hc=hc,lc=cc)
else:    
    process_snp_dists(path,path)
    process_snp_data(path,path)
    plot_and_save(get_crange, get_hrange,  wavelet_transform, path,hc=hc,lc=cc)
