#!/usr/bin/python
import os
import subprocess

import RNAseq_analysis
reload(RNAseq_analysis)


# spit = RNAseq_analysis.RNAseq(['DRR033935.fastq', 'DRR033936.fastq'])
# spit.clean_split_seqs(min_len=30)

def bash(cmd):
	process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
	output = process.communicate()[0]
	return output

def copy_to_s3_1(files, local_dir, s3_dir):
	if not isinstance(files, list):
		files = [files]
	print files
	for f in files:
		aws_loc = 's3://rna-seq/%s/'%s3_dir + f
		bash('aws s3 cp %s/'%local_dir+f+' '+aws_loc)

def copy_from_s3_1(files,local_dir, s3_dir):
	if not isinstance(files, list):
		files = [files]
	for f in files:
		aws_loc = 's3://rna-seq/%s/'%s3_dir + f
		bash('aws s3 cp ' + aws_loc + ' %s/'%local_dir + f)

def make_down_seqs():
	if not os.path.exists('down_seqs/'):
		os.makedirs('down_seqs/')
	else:
		os.system('rm down_seqs/*')


def fastq_dump(f, clean=True):
	print 'fastq dumping'
	os.system('./sratoolkit.2.5.0-1-ubuntu64/bin/fastq-dump --split-3 -O ./down_seqs %s'%f)
	print 'copying to S3'
	copy_to_s3_1([f+'_1.fastq',f+'_2.fastq'], 'down_seqs', 'clean-seq-files')
	copy_to_s3_1([f+'_1.fastq',f+'_2.fastq'], 'down_seqs', 'raw-seq-files')
	# if clean:
	# 	SeqIO.convert(f[:-4]+'.fastq', 'fastq', f[:-4]+'.fasta', 'fasta')
	# 	copy_to_s3(f[:-4]+'.fastq', '.', 'clean-repertoire-data')


with open('RNA_seq_ids_1.txt', 'r') as f:
	f.readline()
	errors=[]
	for line in f:

		ids = line.split('\n')[0]
		if ids == 'SRR1550982':
			# make_down_seqs()
			# print "downloading %s from SRA"%line
			# # try:
			# # 	os.system('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/%s/%s/%s.sra -O /home/ubuntu/down_seqs/%s.sra'%(ids[:6], ids, ids, ids))
			# # except:
			# # 	errors.append('error downloading %s'%ids)

			# try:
			# 	fastq_dump(ids)
			# except:
			# 	errors.append('error fastq dumping %s'%ids)
			# os.system('rm -r down_seqs')

			try:
				current_obj = RNAseq_analysis.RNAseq(['%s.fastq'%ids])

			except:
				errors.append('error initializing %s'%ids)

			# try:
			# 	current_obj.clean_split_seqs(min_len=30)
			# except:
			# 	errors.append('error cleaning %s'%ids)

			# try:
			# 	current_obj.map_reads(genome='NCBI')
			# except:
			# 	errors.append('error mapping %s'%ids)

			# try:
			current_obj.run_HTseq(genome='NCBI')
			# except: 
			# 	errors.append('error HTseq %s'%ids)

			# try:
			# 	current_obj.run_cufflinks(genome='NCBI')
			# except:
			# 	errors.append('error fastq dumping %s'%ids)
			# print errors


			










