#!/usr/bin/python

import os
import subprocess

import RNAseq_analysis
reload(RNAseq_analysis)

current_obj = RNAseq_analysis.RNAseq(['1585_sequence_short.fastq',\
										'3578-E_sequence_short.fastq', '3578_sequence_short.fastq',\
										'9399-E_sequence_short.fastq', '9399_sequence_short.fastq'])
current_obj.clean_split_seqs()
current_obj.map_reads(genome='NCBI')
current_obj.map_reads(genome='miRNA')


def bash(cmd):
	process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
	output = process.communicate()[0]
	return output

def make_temp_seqs():
	if not os.path.exists('temp_seqs/'):
		os.makedirs('temp_seqs/')
	else:
		os.system('rm temp_seqs/*')


def copy_to_s3(files, local_dir, s3_dir):
	if not isinstance(files, list):
		files = [files]
	for f in files:
		aws_loc = 's3://rna-seq/%s/'%s3_dir + f
		bash('aws s3 cp %s/'%local_dir+f+' '+aws_loc)

def copy_from_s3(files,local_dir, s3_dir):
	if not isinstance(files, list):
		files = [files]
	for f in files:
		aws_loc = 's3://rna-seq/%s/'%s3_dir + f
		bash('aws s3 cp ' + aws_loc + ' %s/'%local_dir + f)

def write_short(seqs, short_seqs):
	
	with open('temp_seqs/%s'%seqs) as f:
		with open('temp_seqs/%s'%short_seqs, 'w') as of:
			counter=0
			for line in f:
				if counter<40000000:
					of.write(line)
					counter+=1
				else:
					break



seqs = ['3578-E_sequence.fastq.gz','3578_sequence.fastq.gz']




for s in seqs:
	make_temp_seqs()

	print "copying %s from s3"%s
	copy_from_s3([s],'temp_seqs','raw-seq-files')
	os.system('gunzip temp_seqs/%s'%s)
	big_file = s[:-3]
	small_file = big_file.split('.')[0]+'_short.fastq'
	write_short(big_file, small_file)
	copy_to_s3([small_file], 'temp_seqs', 'raw-seq-files')
	os.system('rm -r temp_seqs')





