# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 14:44:58 2017

@author: Justin
"""
#This script acted as one validation step for our analytical pipeline to check for false negatives - 
#varients that were excluded in our filtering process that should not have been.

#Find the number of usable bases common to two samples - 
#	arg 1 is first sample, arg 2 is second (parent)

#Sample files are generated from a VCF file using VCF tools.
#	Example: vcftools --vcf your.annotated.vcf --recode --out sample1 --indv sampleName
#Filter out any non-heterozygotes:
#	Example: grep "0/1" sample1.recode.vcf > sample1_heteros 
		
#ex: python depth_adj.py sample1_heteros sample2_heteros

import sys

sampleFile = sys.argv[1]
parentFile = sys.argv[2]
sample_list = open(sampleFile, 'r')
parent_list = open(parentFile, 'r')

samplecount = 0
parentcount = 0

#step 1) create sample set
sample_set = set()
for line in sample_list:
    samplecount = samplecount + 1
    line = line.strip().split()
    sample_set.add(str(line[0])+str(line[1]))
print("sample list compiled")

#step 2) create parent set
parent_set = set()
for line in parent_list:
    parentcount = parentcount + 1
    line = line.strip().split()
    parent_set.add(str(line[0])+str(line[1]))
print("parent list compiled")

#step 3) check for bp that meet the standeard in both sample and parent
heteros = len(sample_set & parent_set)
percent = heteros/parentcount
print(str(sampleFile) + "\t" + "Heterozygotes found: " + str(heteros))
print("Sample: " + str(samplecount))
print("Parent: " + str(parentcount))
print("Percent found: " + str(percent))
