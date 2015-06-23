#!/usr/bin/python

import RNAseq_analysis
reload(RNAseq_analysis)


# spit = RNAseq_analysis.RNAseq(['DRR033935.fastq', 'DRR033936.fastq'])
# spit.clean_split_seqs(min_len=30)

def make_down_seqs():
	if not os.path.exists('down_seqs/'):
		os.makedirs('down_seqs/')
	else:
		os.system('rm down_seqs/*')


def fastq_dump(f, clean=True):
	os.system('./sratoolkit.2.5.0-1-ubuntu64/bin/fastq-dump --split-3 -O ./down_seqs ./%s'%f)
	copy_to_s3(f[:-4]+'.fastq', '.', 'patient_repertoire_data')
	if clean:
		SeqIO.convert(f[:-4]+'.fastq', 'fastq', f[:-4]+'.fasta', 'fasta')
		copy_to_s3(f[:-4]+'.fastq', '.', 'clean-repertoire-data')


with open('RNA_seq_ids.txt', 'r') as f:

	for line in f:
		make_down_seqs()
		print "downloading %s from SRA"%ids
		os.system('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/%s/%s/%s.sra -O /home/ubuntu/down_seqs/%s.sra'%(line[:6], line, line, line))
		fastq_dump(ids+'.sra')

		current_obj = RNAseq_analysis.RNAseq(['down_seqs/%s.fastq'%line])
		current_obj.clean_split_seqs(min_len=30)
		current_obj.map_reads(genome='NCBI')
		current_obj.run_HTseq(genome='NCBI')
		current_obj.run_cufflinks(genome='NCBI')

		os.system('rm -r down_seqs')