#!/usr/bin/python
import subprocess
import os

def bash(cmd):
	process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
	output = process.communicate()[0]
	return output

def copy_to_s3(files, local_dir, s3_dir):
	if not isinstance(files, list):
		files = [files]
	for f in files:
		aws_loc = 's3://rna-seq/%s/'%s3_dir + f
		bash('aws s3 cp %s/'%local_dir+f+' '+aws_loc)

i='DRR033936-len-q-t'

copy_to_s3([i+'.bam', i+'.bam.bai'], '.', 'bam-files')
copy_to_s3(['%s_genes.fpkm_tracking'%i, '%s_isoforms.fpkm_tracking'%i,'%s_transcripts.gtf'%i,'%s_skipped.gtf'%i],\
				     	'.','cufflinks-out')

