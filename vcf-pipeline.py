# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 09:46:34 2017

@author: Justin
"""
#Arguments required: 1)vcf file to filter 2)name of parent 3)outFile name
#	1) vcf file to filter: This is the file containing the varients. Must be in VCF format. For more details copy link below:
#		http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/
#	2) Name of parent: This argument lets you indicate which of your samples is the parent sample that the
#		others decsend from.
#	3) outFile Name: indicate where you would like to store results. 
#		This script will make a filtered vcf file with the passed name: outFileEx.vcf
#		and a text file with tab delimited data easily imported into Excel or R: outFile.txt

#Ex: python vcf_filter_pipeline.py myvcf.vcf 2460_V out

#count the number of samples with variant
#In: a list of the sample data for a given varient
#Out: a list containing the number of each geneotype from the sample list for a given varient
#	heterogeneous, homogeneous to the reference geneotype, homogeneous but different from the reference, unreadable 

import sys

def count_variant_samples(sample_list):
    hetero = 0
    homo_ref = 0
    homo_alt = 0
    unreadable = 0
    for sample in sample_list:
        sample_data = sample.split(":")
        #count the number of 0/1, 0/0, 1/1, oand ./.
        if sample_data[0] == "0/1" or sample_data[0] == "0|1":
            hetero = hetero + 1
        elif sample_data[0] == "0/0" or sample_data[0] == "0|0":
            homo_ref = homo_ref + 1
        elif sample_data[0] == "1/1" or sample_data[0] == "1|1":
            homo_alt = homo_alt + 1
        elif sample_data[0] == "./." or sample_data[0] == ".|.":
            unreadable = unreadable + 1
        #Error message
		else:
            print("GT not recognized")
    count = [hetero, homo_ref, homo_alt, unreadable]
    return count

#Return the genotype of a given sample
#In: sample to get genotype of, list of sample for a varient
#Out: genotype of sample
def return_GT(sample_index, sample_list):
    target_sample = sample_list[sample_index]
    return target_sample[0:3]

#determine what to do with a given variant based on the hetero/homo counts
#Case 0)More than half of the samples are unreadable or the variant doesn't meet uniquiness criteria
#Case 1)Parent is homo_ref & there is only one variant present (homo or hetero)
#Case 2)Parent is homo_alt & all samples are homo_alt except 1
#Case 3)2 Parents are hetero & all samples are hetero except 1 (homo_ref, or homo_alt)
def choose_case(variant_counts, parent_index, sample_list):
    sample_count = len(sample_list)
    readable_samples = sample_count - variant_counts[3]
    req_sample_number = readable_samples - 1
    parent_GT = return_GT(parent_index, sample_list)
    #Case 0 test
    if float(variant_counts[3])/float(len(sample_list)) > .5:
        return 0
    #Case 1 test
    elif parent_GT == "0/0" or parent_GT == "0|0":
        if variant_counts[0] == 1 or variant_counts[2] == 1:
            return 1
        else:
           return 0
    #Case 2 test   
    elif parent_GT == "1/1" or parent_GT == "1|1":
       if variant_counts[2] == (req_sample_number):
            return 2
       else:
            return 0
    #Case 3 test   
    elif parent_GT == "0/1" or parent_GT == "0|1":
        if variant_counts[0] == (req_sample_number):
            return 3
        else:
            return 0
    #All remainders
    else:
        return 0

#Return index number of AD in the format list
def find_ADindex(formatlist):
    names = formatlist.split(":")
    return(names.index("AD"))

#Return index number of the sample containing the variant
def find_sample(sample_list):
    index_counter = 0
    for sample in sample_list:
        sample_data = sample.split(":")
        if sample_data[0] == "0/1" or sample_data[0] == "1/1" or sample_data[0] == "0|1" or sample_data[0] == "1|1":
            sample_index = index_counter + 9
            return sample_index
        else:
            index_counter = index_counter + 1
    #if no no variant with sufficient data is found 
    failed = -1
    return failed

#check AD 
#both the reference and alternate must be greater than the cutoff 
# (use 2 for a defaul cutoff AD value)
def check_AD(AD_index, sample_data, cut_off):
    segmented_data = sample_data.split(":")
    AD = segmented_data[AD_index]
    splitAD = AD.split(",")
    if int(splitAD[0]) >= cut_off and int(splitAD[1]) >= cut_off:
        return True
    else:
        return False

#Ensure there are enough samples with sufficiant read depth for quality results
def count_pass(sample_list):
    counter = 0
    unreadable = 0
    for sample in sample_list:
        counter = counter + 1
        segmented_data = sample.split(":")
        GT = segmented_data[0]
        if GT == "./.":
            unreadable = unreadable + 1
    usable_samples = counter - unreadable
    return usable_samples


#MAIN
#open vcf file to read - You can input the name directly below or
#	default read file passed as first argument
#InFileName = "vcfsample3.txt"
InFileName = sys.argv[1]
InFile = open(InFileName, 'r')

#open vcf file to write
OutFileName = "testout"
OutVCF = open(OutFileName + ".vcf", 'w')
OutTable = open(OutFileName + ".txt", 'w')

#initialize variables to defaults
AD_index_assigned = False
#parent name can be passed as an argument or added below
#parent = "2460_V"
parent = sys.argv[2]
raw_parent_index = -1
AD_cutoff = 2
usable_sample_min = 5
pre_sample_col = 9

for line in InFile:
    #first identify header lines and write to vcf file directly
    if line[0] == "#":
            OutVCF.write(line)
            #If the header line contains the column names, extract the names
            if line[0] + line[1] == "#C":
                column_headers = line.strip().split()
                sample_count = len(column_headers) - pre_sample_col
        
                #Find the index number of the progenetor sample for later comparison
                if parent in column_headers:
                    raw_parent_index = column_headers.index(parent)
                else:
                    print("Error: could not determine progenetor sample")
                    break
    
    #pass data lines to be parsed appropriately            
    else:
        tab_delim = line.strip().split()
        adj_parent_index = raw_parent_index - pre_sample_col
        sample_list = tab_delim[pre_sample_col:]
        variant_counts = count_variant_samples(sample_list)
        case = choose_case(variant_counts, adj_parent_index, sample_list)
        #First check if you know the index number of the AD values and assign it if it is the first loop
        if AD_index_assigned == False:
            AD_index = find_ADindex(tab_delim[8])
            AD_index_assigned = True
        
       
        
        #Once you have the AD index go straight to the next step for each line
        #1) Find the sample with the variant (hetero or homozygote) -> if it is the progenetor remove variant. Ensure line has at least one valid variant sample
        #2) Find AD values for ref and alt -> if either is below 2 remove variant
        #3) if variant passes both tests write to outVCF and outTable
        if find_sample(tab_delim[pre_sample_col:]) >= 0:
            sample_data = tab_delim[find_sample(tab_delim[9:])]
            variant_index = tab_delim.index(sample_data)
            #print("pass find_sample")
            if variant_index != adj_parent_index:
               # print("pass parent check")
                if check_AD(AD_index, sample_data, AD_cutoff) == True:
                    #print("pass check_AD")
                    if count_pass(tab_delim[pre_sample_col:]) >= usable_sample_min:       
                       # print("pass read min")
                        OutVCF.write(line)
                        OutTable.write(tab_delim[0] + "\t" + tab_delim[1] + "\t" + column_headers[variant_index] + "\n")
            
InFile.close()
OutVCF.close()
OutTable.close()
print("done")
