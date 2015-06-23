#!/usr/bin/python

import collections
import HTSeq

gtf_file = '/home/ubuntu/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf'


gtf_file = HTSeq.GFF_Reader(gtf_file)
features = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
feature_dict = collections.Counter()
print "extracting features from gtf file"
for feature in gtf_file:
	# if feature.type == "exon":
	features[feature.iv] += feature.attr["gene_id"]
	feature_dict[feature.type] +=1

print feature_dict


gtf_file = '/home/ubuntu/genomes/Homo_sapiens/NCBI/build37.2/Annotation/Genes/genes.gtf'

gtf_file = HTSeq.GFF_Reader(gtf_file)
features = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
feature_dict = collections.Counter()
print "extracting features from gtf file"
for feature in gtf_file:
	# if feature.type == "exon":
	features[feature.iv] += feature.attr["gene_id"]
	feature_dict[feature.type] +=1

print feature_dict
