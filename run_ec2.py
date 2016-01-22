# Program requires a a system argument with the path to the key to use on the instance

import boto.ec2
import time
import sys


"""modify for current data and pipeline"""
data_locs=['test/test_seqs.fastq'] #can be a list of all the S3 locs of files to be processed
script_loc = 'scripts/RNAseq.py'




"""  MAIN SCRIPT """

key = sys.argv[1]
ami = sys.argv[2]
script_name = script_file.split('/')[-1]
# Enhanced creation now with the addition of 'user_data'
for data_loc in data_locs:
	data_name = data_file.split('/')[-1]

	conn = boto.ec2.connect_to_region('us-west-1')

	#user_data_script requires that the first sys argv to script_file is the data_file to be processed
	user_data_script = """#!/bin/bash 
	mkdir temp_seqs
	aws s3 cp s3://rna-seq/%s /home/ubuntu/temp_seqs/%s --region us-west-1
	aws s3 cp s3://rna-seq/%s /home/ubuntu/%s --region us-west-1
	python %s temp_seqs/%s """%(data_loc,data_name,script_loc, script_name, script_name, data_name)



	# Red Hat Enterprise Linux 6.4 (ami-7d0c6314)
	new_reservation = conn.run_instances(
	                        ami,
	                        key_name=key,
	                        instance_type='m4.2xlarge',
	                        security_groups=['default'],
	                        user_data=user_data_script,
	                        instance_profile_name='myinstanceprofile_rnaseq')
	print "New instance created."

	# Add a Name to the instance, then loop to wait for it to be running.
	instance = new_reservation.instances[0]
	conn.create_tags([instance.id], {"Name":"In Out Instance Test"})
	while instance.state == u'pending':
	    print "Instance state: %s" % instance.state
	    time.sleep(10)
	    instance.update()

	print "Instance state: %s" % instance.state
	print "Public dns: %s" % instance.public_dns_name

