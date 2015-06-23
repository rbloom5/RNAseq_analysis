 #!/usr/bin/python
import os
import subprocess
import webbrowser
import multiprocessing
import HTSeq
import collections

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


def sam_index(f, outdir=''):
	print "sam index"
	os.system('samtools index %s'%f)
	# if outdir:
	# 	os.system('mv %s.bai %s/'%(f.split('/')[-1], outdir))
	os.system('rm %s.bai'%f)



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
		if len(bundle) != 1
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




class RNAseq:
	def __init__(self, files, s3_folder='raw-seq-files', aws=True):
		# Initialize the object with the files you want to analyze
		# files is a list of files (including the extension) of the sequencing data

		# if the files are local (i.e. not in s3), then give the full path to the 
		# files on your computer and set aws=False
		if not isinstance(files, list):
			files = [files]
		if aws:
			prefix = 's3://rna-seq/' + s3_folder +'/'
		else:
			prefix = ''

		self.file_locations = [prefix+f for f in files]
		self.ids = [f.split('.')[0] for f in files]

	def clean_split_seqs(self, output_prefix='', barcode_file=None, min_len=15, qual_score=20, aws=True):
		# gets rid of low quality sequences and converts to FASTA, puts into clean-rna-seq

		# a list of the file ids that comes out of the splitting
		self.ids=[]

		# loop through files, copy to ec2, process based on input parameters, and copy output back to s3
		for f in self.file_locations:
			make_temp_seqs()
			# i stands for id - the file name without the .fastq extension
			# f_tail is just the filename, no directory
			i = f.split('/')[-1].split('.')[0]
			f_tail = f.split('/')[-1]

			if aws:
				#if aws, need to copy files from s3 to ec2 first
				print "copying %s from s3"%f
				os.system('aws s3 cp %s temp_seqs/'%f)
				
				print "removing sequences that are too short"
				bash('fastx_clipper -i temp_seqs/%s -o temp_seqs/%s-len.fastq -l %s'%(f_tail,i,min_len))

			else:
				print "removing sequences that are too short"
				bash('fastx_clipper -i %s -o temp_seqs/%s-len.fastq -l %s'%(f,i,min_len))


			#remove sequences with quality threshold below the specified input parameters
			print "removing low quality sequences"
			bash('fastq_quality_filter -q %s -p 80 -i temp_seqs/%s-len.fastq -o temp_seqs/%s-len-q.fastq'%(qual_score,i,i))

			#split up files by barcode
			if barcode_file:
				keep_files = split_barcodes(i, barcode_file, output_prefix, 'temp_seqs')

			else:
				# if no barcodes, just copy back the 
				keep_files = [f for f in os.listdir('temp_seqs') if f.endswith('-len-q.fastq')]
				for kf in keep_files:
					os.system('fastx_trimmer -f %s -i %s -o %s'%(13,'temp_seqs/'+kf, 'temp_seqs/'+kf[:-6]+'-t.fastq'))

				keep_files = [f for f in os.listdir('temp_seqs') if f.endswith('-len-q-t.fastq')]



			# generate fastqc reports on new cleaned sequence files:
			print 'generating fastQC report'
			for kf in keep_files:
				bash("fastqc temp_seqs/" + kf + " -outdir=./temp_seqs")
				copy_to_s3(['%s_fastqc.html'%kf.split('.')[0]], 'temp_seqs', 'fastqc-reports')


			# copy back to s3
			print "copying %s files back to s3"%i
			copy_to_s3(keep_files, 'temp_seqs', 'clean-seq-files')

			#get rid of temp files
			os.system('rm -r temp_seqs')

			# append files without extensions as ids to save in object for all the other methods
			for f in keep_files:
				self.ids.append(f.split('.')[0])


	def map_reads(self, genome='hg19', coverage_search=False):
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
			genome_loc = '/home/ubuntu/genomes/Homo_sapiens/NCBI/build37.2/Annotation/SmallRNA/hairpin.fa'
			gtf_loc = '/home/ubuntu/genomes/Homo_sapiens/NCBI/build37.2/Annotation/SmallRNA/mirs.gff'
			gtf_index = '/home/ubuntu/genomes/current_annotations/miRNA/miRNA_hairpin_annotations'

			

		for i in self.ids:
			make_temp_seqs()
			print "copying %s from s3"%(i+'.fastq')
			copy_from_s3([i+'_1.fastq',i+'_2.fastq'],'temp_seqs','clean-seq-files')
			num_cores = multiprocessing.cpu_count()

			print "running tophat on %s cores"%num_cores

			tophat_string = 'tophat2 -p %s --no-coverage-search -G %s --transcriptome-index %s --no-novel-juncs %s %s %s'%\
								(num_cores,gtf_loc,gtf_index,genome_loc,'temp_seqs/'+i+'_1.fastq','temp_seqs/'+i+'_2.fastq')

			print tophat_string
			bash(tophat_string)

			#rename and move output from tophat to s3
			os.system('mv tophat_out/accepted_hits.bam ~/tophat_out/%s'%(i+'.bam'))
			os.system('mv tophat_out/align_summary.txt ~/tophat_out/%s'%(i+'_align_summary.txt'))

			#index the bam file for IGV viewing
			sam_index('~/tophat_out/%s'%(i+'.bam'), outdir='tophat_out')

			#copy to s3
			copy_to_s3([i+'.bam',i+'_align_summary.txt', i+'.bam.bai'], 'tophat_out', 'bam-files')

			os.system('rm -r tophat_out/*')

		os.system('rm -r temp_seqs')


	def run_cufflinks(self, genome='Ensembl'):
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
			genome_loc = '/home/ubuntu/genomes/Homo_sapiens/NCBI/build37.2/Annotation/SmallRNA/hairpin.fa'
			gtf_loc = '/home/ubuntu/genomes/Homo_sapiens/NCBI/build37.2/Annotation/SmallRNA/mirs.gff'
			gtf_index = '/home/ubuntu/genomes/current_annotations/miRNA/miRNA_hairpin_annotations'

		for i in self.ids:
			make_temp_seqs()
			print "copying %s from s3"%(i+'.bam')
			copy_from_s3(i+'.bam','temp_seqs','bam-files')
			num_cores = multiprocessing.cpu_count()

			# run cufflinks to count reads and estimate gene expression
			bash('cufflinks -p %s -N -G %s -u --compatible-hits-norm temp_seqs/%s'\
				%(num_cores, gtf_loc, i+'.bam'))
			# bash('cufflinks -p %s temp_seqs/%s'\
			# 	%(num_cores, i+'.bam'))

			os.system('mv genes.fpkm_tracking %s_genes.fpkm_tracking'%i)
			os.system('mv isoforms.fpkm_tracking %s_isoforms.fpkm_tracking'%i)
			os.system('mv transcripts.gtf %s_transcripts.gtf'%i)
			os.system('mv skipped.gtf %s_skipped.gtf'%i)

			print "copying output to s3"
			copy_to_s3(['%s_genes.fpkm_tracking'%i, '%s_isoforms.fpkm_tracking'%i,'%s_transcripts.gtf'%i,'%s_skipped.gtf'%i],\
				     	'.','cufflinks-out')

			for i in ['%s_genes.fpkm_tracking'%i, '%s_isoforms.fpkm_tracking'%i,'%s_transcripts.gtf'%i,'%s_skipped.gtf'%i]:
				os.system('rm %s'%i)
		os.system('rm -r temp_seqs')



	def run_HTseq(self, genome='Ensembl'):
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
			genome_loc = '/home/ubuntu/genomes/Homo_sapiens/NCBI/build37.2/Annotation/SmallRNA/hairpin.fa'
			gtf_loc = '/home/ubuntu/genomes/Homo_sapiens/NCBI/build37.2/Annotation/SmallRNA/mirs.gff'
			gtf_index = '/home/ubuntu/genomes/current_annotations/miRNA/miRNA_hairpin_annotations'

		for i in self.ids:
			make_temp_seqs()
			print "copying %s from s3"%(i+'.bam')
			copy_from_s3(i+'.bam','temp_seqs','bam-files')

			counts = HTseq_count('temp_seqs/'+i+'.bam', gtf_loc, 'temp_seqs', 'transcript_id')

			with open('%s_HTseq_counts.txt'%i,'w') as f:
				for c in counts:
					f.write('%s\t%s\n'%(c, counts[c]))

			copy_to_s3('%s_HTseq_counts.txt'%i,'.','htseq-out')
			os.system('rm %s_HTseq_counts.txt'%i)
			os.system('rm temp_seqs/%s.bam'%i)












		




