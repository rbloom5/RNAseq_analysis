#!/usr/bin/python

import urllib
import os

urllib.urlretrieve('ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz', 'genomes/Ensembl_GRCh37.tar.gz')
urllib.urlretrieve('ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/NCBI/build37.2/Homo_sapiens_NCBI_build37.2.tar.gz', 'genomes/NCBI_build37.2.tar.gz')
urllib.urlretrieve('ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz', 'genomes/UCSC_hg19.tar.gz')

os.system('tar -zxvf genomes/Ensembl_GRCh37.tar.gz')
os.system('tar -zxvf genomes/NCBI_build37.2.tar.gz')
os.system('tar -zxvf genomes/UCSC_hg19.tar.gz')