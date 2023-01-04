#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 18:57:46 2021

@author: april
"""

import re
import numpy as np
import codecs
import argparse
import os, sys

parser = argparse.ArgumentParser(description='取同一测序深度，计算测到的基因数——取同一测序深度，对Aligned.out.sam文件进行过滤： \
                                              python2 DepthGene_v2.py \
                                              -b Aligned.out.sam \
                                              -d depth \
                                              -p prefix')


parser.add_argument('-b', '--bam', type=str, required=True, help='file Aligned.out.bam')
parser.add_argument('-d', '--depth', type=int, required=True, help='depth')
parser.add_argument('-p', '--prefix', type=str, required=True, help='prefix')

args = parser.parse_args()

os.system("samtools view -h %s > %s.sam " % (args.bam, args.prefix))
sam=str(args.prefix)+".sam"
ref_dict={}
count=0
with open(sam, 'r') as ref:
    for line in ref:
        if line.startswith("@"):
            continue
        count += 1
        info = line.split("\t")
        ID = info[0]
        ref_dict[ID] = 1
print("total reads number (paired reads count for 1): "+str(len(ref_dict)))
print("total reads number: "+str(count))


d = int(args.depth/2)
select_ID = np.random.choice(list(ref_dict.keys()), d, replace = False)
wri_id = {}
for i in select_ID:
    wri_id[i] = 1
print("select read ID: " + str(len(wri_id)))

n = 0
out = str(args.prefix) + '_' + str(args.depth/1000000) + '_Depth.sam'
with open(sam, "r") as sam:
    with open(out, 'w') as OUT:
        for line in sam:
            line = line.rstrip()
            if line.startswith("@"):
                OUT.write(line+"\n")
            else:
                tmp = line.split("\t")
                sam_id = tmp[0]
                flag = wri_id.get(sam_id, 0)
                if flag == 1:
                    n += 1
                    OUT.write(line+"\n")
                else:
                    continue
print("Output sam file number: "+str(n))

bam2 = str(args.prefix) + '_' + str(args.depth/1000000) + '_Depth.bam'
os.system("samtools view -b -S %s > %s " % (out, bam2))
