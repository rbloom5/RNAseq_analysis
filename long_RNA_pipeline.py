 #!/usr/bin/python
import os
import subprocess
import webbrowser
import multiprocessing
import HTSeq
import collections
import sys

# object for RNA seq cleaning and analysis

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

def copy_from_s3(files,local_dir, s3_dir):
	if not isinstance(files, list):
		files = [files]
	for f in files:
		aws_loc = 's3://rna-seq/%s/'%s3_dir + f
		bash('aws s3 cp ' + aws_loc + ' %s/'%local_dir + f)


def make_temp_seqs():
	if not os.path.exists('temp_seqs/'):
		os.makedirs('temp_seqs/')
	else:
		os.system('rm temp_seqs/*')


def trim_barcodes(barcode_file, keep_files, output_prefix, data_dir):
	with open(barcode_file) as bf:
		bc_lens={}
		for line in bf:
			bc_lens[line.split('\t')[0]]=len(line.split('\t')[1])

	f_keep_files = []
	for kf in keep_files:
		for key in bc_lens:
			if output_prefix+key == kf:
				print "trimming barcode from "+key
				bash('fastx_trimmer -f %s -i %s -o %s'%(bc_lens[key]+1,data_dir+'/'+kf, data_dir+'/'+kf+'-len-q-t'))
				f_keep_files.append(kf+'-len-q-t')

	return f_keep_files



def split_barcodes(i, barcode_file, output_prefix,file_dir):
	print "splitting sequences by barcode"
	split_cmd = 'cat temp_seqs/%s-len-q.fastq | fastx_barcode_splitter.pl --bcfile %s \
			--bol --prefix temp_seqs/%s --suffix .fastq'%(i,barcode_file,output_prefix)

	# gives each barcode it's own file and makes a list of files with more than 1000 sequences to keep
	keep_files = []
	for line in subprocess.Popen(split_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[1:-3]:
		if int(line.split('\t')[1])>1000:
			keep_files.append(output_prefix+line.split('\t')[0]+'.fastq')

	keep_files = trim_barcodes(barcode_file, keep_files, output_prefix, file_dir)




def get_genome_files(genome):
	if genome == 'Ensembl':
		genome_loc = '/home/ubuntu/genomes/Ensembl/GRCh37/Sequence/Bowtie2Index/genome'
		gtf_loc = '/home/ubuntu/genomes/Ensembl/GRCh37/Annotation/Genes/genes.gtf'
		gtf_index = '/home/ubuntu/genomes/current_annotations/Ensemble/Ensemble_annotations'

	elif genome == 'UCSC':
		genome_loc = '/home/ubuntu/genomes/UCSC/hg19/Sequence/Bowtie2Index/genome'
		gtf_loc = '/home/ubuntu/genomes/UCSC/hg19/Annotation/Genes/genes.gtf'
		gtf_index = '/home/ubuntu/genomes/current_annotations/UCSC/UCSC_annotations'

	elif genome == 'NCBI':
		genome_loc = '/home/ubuntu/genomes/NCBI/build37.2/Sequence/Bowtie2Index/genome'
		gtf_loc = '/home/ubuntu/genomes/NCBI/build37.2/Annotation/Genes/genes.gtf'
		gtf_index = '/home/ubuntu/genomes/current_annotations/NCBI/NCBI_annotations'

	elif genome == 'miRNA':
		genome_loc = '/home/ubuntu/genomes/Ensembl/GRCh37/Sequence/Bowtie2Index/genome'
		gtf_loc = '/home/ubuntu/genomes/Ensembl/GRCh37/Annotation/SmallRNA/mature.fa'
		gtf_index = '/home/ubuntu/genomes/current_annotations/miRNA/miRNA_annotations'

	return genome_loc, gtf_loc, gtf_index





def HTseq_count(bam_file, gtf_file, out_dir, identifier, parallel = True ):
	gtf_file = HTSeq.GFF_Reader(gtf_file)
	features = HTSeq.GenomicArrayOfSets( "auto", stranded=True )

	print "extracting features from gtf file"
	for feature in gtf_file:
		# if feature.type == "exon":
		features[feature.iv] += feature.attr[identifier]

	counts = collections.Counter( )

	almnt_file = HTSeq.SAM_Reader(bam_file)
	counts = collections.Counter( )
	for bundle in HTSeq.pair_SAM_alignments( almnt_file, bundle=True ):
		if len(bundle) != 1:
			continue  # Skip multiple alignments
		first_almnt, second_almnt = bundle[0]  # extract pair
		if not first_almnt.aligned and second_almnt.aligned:
			count[ "_unmapped" ] += 1
			continue
		gene_ids = set()
		for iv, val in features[ left_almnt.iv ].steps():
			gene_ids |= val
		for iv, val in features[ right_almnt.iv ].steps():
			gene_ids |= val
		if len(gene_ids) == 1:
			gene_id = list(gene_ids)[0]
			counts[ gene_id ] += 1
		elif len(gene_ids) == 0:
			counts[ "_no_feature" ] += 1
		else:
			counts[ "_ambiguous" ] += 1

	for gene_id in counts:
		print gene_id, counts[ gene_id ]





def clean_seqs(files, output_prefix='', min_len=25, qual_score=20):
	# gets rid of low quality sequences and converts to FASTA, puts into clean-rna-seq
	# a list of the file ids that comes out of the splitting
	#assume files are in directory 'temp_seqs'
	if not isinstance(files, list):
		files = [files]

	# loop through files, copy to ec2, process based on input parameters, and copy output back to s3
	new_file_ids = []
	for f in files:
		
		i = f.split('/')[-1].split('.')[0] # i stands for id - the file name without the .fastq extension
		f_tail = f.split('/')[-1] # f_tail is just the filename, no directory
		#unzip if needed
		if f_tail[-3:]=='.gz':
			os.system('gunzip temp_seqs/%s'%f_tail)
			f_tail = f_tail[:-3]
		
		#remove sequences with quality threshold below the specified input parameters
		print "removing sequences that are too short"
		bash('fastx_clipper -i temp_seqs/%s -o temp_seqs/%s-len.fastq -l %s'%(f_tail,i,min_len))
		print "removing low quality sequences"
		bash('fastq_quality_filter -q %s -p 80 -i temp_seqs/%s-len.fastq -o temp_seqs/%s-len-q.fastq'%(qual_score,i,i))


		# generate fastqc reports on new cleaned sequence files:
		print 'generating fastQC report'
		bash("fastqc temp_seqs/" + i + "-len-q.fastq -outdir=./temp_seqs")
		copy_to_s3(['%s_fastqc.html'%(i+'-len-q')], 'temp_seqs', 'fastqc-reports')


		# copy back to s3
		print "copying %s files back to s3"%i
		copy_to_s3(i+'-len-q.fastq', 'temp_seqs', 'clean-seq-files')

		#get rid of temp files
		os.system('rm -r temp_seqs')

		# append files without extensions as ids to save in object for all the other methods
		new_file_ids.append(i+"-len-q")

	return new_file_ids




def map_reads(files, genome='UCSC', paired_end=False, fusions=False, isoforms=False):
	genome_loc, gtf_loc, gtf_index = get_genome_files(genome)

	new_ids = []	
	for i in files:
		make_temp_seqs()
		print "copying %s from s3"%(i+'.fastq')
		if paired_end:
			copy_from_s3([i+'_1.fastq',i+'_2.fastq'],'temp_seqs','clean-seq-files')
		else:
			copy_from_s3([i+'.fastq'],'temp_seqs','clean-seq-files')


		num_cores = multiprocessing.cpu_count()
		print "running tophat on %s cores"%num_cores
		tophat_base = 'tophat2 -p %s -o /home/ubuntu --no-coverage-search -G %s --transcriptome-index %s'%(num_cores,gtf_loc,gtf_index)
		
		if fusions:
			tophat_base += ' --fusion-search --bowtie1'
		if not isoforms:
			tophat_base += ' --no-novel-juncs'

		if paired_end:
			tophat_string = tophat_base+' %s %s %s'%\
							(genome_loc,'temp_seqs/'+i+'_1.fastq','temp_seqs/'+i+'_2.fastq')
		else:
			tophat_string = tophat_base+' %s %s '%\
							(genome_loc,'temp_seqs/'+i+'.fastq')				
		print tophat_string
		bash(tophat_string)
		os.system('rm -r logs')


		new_i = i+'_' + genome
		if fusions:
			fusion_string = 'tophat-fusion-post --num-threads %s %s'%(num_cores, genome_loc)
			bash(fusion_string)
			os.system('mv tophatfusion_out/result.txt tophatfusion_out/%s'%(new_i+'_fusions.txt'))
			copy_to_s3(new_i+'_fusions.txt', 'tophatfusion_out', 'fusions')
			os.system('rm -r tophatfusion_out')

		#rename and move output from tophat to s3
		tophat_out = ['accepted_hits.bam',  'align_summary.txt',  'deletions.bed',  'fusions.out', \
						 'insertions.bed',  'junctions.bed',  'prep_reads.info',  'unmapped.bam']

		for out in tophat_out:
			new_out = new_i+'_'+out
			os.system('mv %s %s'%(out, new_out))
			copy_to_s3(new_out, '.', 'bam-files')
			os.system('rm %s'%new_out)

		new_ids.append(new_i)

		sys.stdout.flush()

	os.system('rm -r temp_seqs')
	print 'mapping completed!'
	return new_ids





def run_cufflinks(files, genome='UCSC'):
	genome_loc, gtf_loc, gtf_index = get_genome_files(genome)

	for i in files:
		make_temp_seqs()
		print "copying %s from s3"%(i+'_accepted_hits.bam')
		copy_from_s3(i+'_accepted_hits.bam','temp_seqs','bam-files')
		num_cores = multiprocessing.cpu_count()

		# run cufflinks to count reads and estimate gene expression
		bash('cufflinks -p %s -N -G %s -u --compatible-hits-norm temp_seqs/%s'\
			%(num_cores, gtf_loc, i+'_accepted_hits.bam'))

		
		cufflinks_out = ['genes.fpkm_tracking', 'isoforms.fpkm_tracking','transcripts.gtf','skipped.gtf']
		for out in cufflinks_out: 
			new_out = i+'_'+out
			os.system('mv %s %s'%(out,new_out))
			copy_to_s3(new_out, '.', 'cufflinks-out')
			os.system('rm %s'%new_out)

		sys.stdout.flush()

	os.system('rm -r temp_seqs')
	print "cufflinks complete!"
	return




if __name__ == "__main__":
	seq_file = sys.argv[1]
	new_id = clean_seqs(seq_file)
	new_id = map_reads(new_id, fusions=True, genome='NCBI')
	run_cufflinks(new_id, genome='NCBI')













