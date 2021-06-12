# vcf-pipeline
Scripts used in publication: "Low genetic variation is associated with low mutation rate in the giant duckweed" 18 Mar 2019
  https://www.nature.com/articles/s41467-019-09235-5#data-availability
  See especially: Methods -> Mutation rate estimation and false-negative calculations
 
 This pipeline filters a list of varients in vcf format and represents several steps of data quality control.
 Vcf files contain a list of varients and information about those varients including the genotype of all samples at that point. 
 More information about vcf file formats can be found here:
  http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/

The files in this repository are not representitive of every step of data analysis for this publication, but rather are selected
to highlight some of the projects I put together while working on this publication.

Information about included files:

1) pipeline.py
Arguments required: 1)vcf file to filter 2)name of parent 3)outFile name
	1) vcf file to filter: This is the file containing the varients. Must be in VCF format. 
	2) Name of parent: This argument lets you indicate which of your samples is the parent sample that the
		others decsend from.
	3) outFile Name: indicate where you would like to store results. 
		This script will make a filtered vcf file with the passed name: outFileEx.vcf
		and a text file with tab delimited data easily imported into Excel or R: outFile.txt

Ex: python vcf_filter_pipeline.py myvcf.vcf 2460_V out

count: the number of samples with variant
In: a list of the sample data for a given varient
Out: a list containing the number of each geneotype from the sample list for a given varient
  heterogeneous, homogeneous to the reference geneotype, homogeneous but different from the reference, unreadable 

2)false_neg.py
  This script acted as one validation step for our analytical pipeline to check for false negatives - 
  varients that were excluded in our filtering process that should not have been.

  Find the number of usable bases common to two samples - 
	  arg 1 is first sample, arg 2 is second (parent)

  Sample files are generated from a VCF file using VCF tools.
  	Example: vcftools --vcf your.annotated.vcf --recode --out sample1 --indv sampleName
  Filter out any non-heterozygotes:
  	Example: grep "0/1" sample1.recode.vcf > sample1_heteros 
		
  ex: python depth_adj.py sample1_heteros sample2_heteros
