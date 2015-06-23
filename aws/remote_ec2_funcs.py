#!/usr/bin/python
from optparse import OptionParser
import logging
import time
import subprocess
from boto.ec2.connection import EC2Connection
import boto


# sra_file = ''

# sra_ids = []
# with open(sra_file, 'r') as f:
# 	for line in f:
# 		sra_ids.append(line.split('\n')[0])

# n=18
# num = float(len(sra_ids))/n
# chunks = [ sra_ids [i:i + int(num)] for i in range(0, (n-1)*int(num), int(num))]
# #l is a list of lists - each sub list containing the sraids for that chunk
# chunks.append(sra_ids[(n-1)*int(num):])

# print chunks

def terminate_instances(ec2, instance_ip):
    for reservation in ec2.get_all_instances():
        for instance in reservation.instances:
            if instance.state not in ['shutting-down', 'terminated'] and instance.ip_address == instance_ip:
                logging.info('terminating instance %s (%s)', instance.id, instance.state)
                ec2.terminate_instances([instance.id])

def print_instances(ec2):
    for reservation in ec2.get_all_instances():
        for instance in reservation.instances:
            print instance.id, instance.state, instance.ip_address

# Returns True if everything went OK 
#
def execute_command(cmd):
    p = subprocess.Popen(cmd,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    (output, error) = p.communicate()
    logging.info('%s exited; returncode: %s; stdout: %s; stderr: %s', cmd, p.returncode, repr(output), repr(error))
    return p.returncode == 0


# Create an instance and wait for it to be ssh-able
#
def run_instance(ec2, image_id, instance_type, key_name, key_file, timeout=300):
	start_time = time.time()
	reservation = ec2.run_instances(image_id=image_id, instance_type=instance_type, key_name=key_name)
	instance = reservation.instances[0]
	logging.info('started instance %s', instance.id)
	while time.time() - start_time < timeout:
	    if instance.update() == 'running':
	        instance = ec2.get_all_instances([instance.id])[0].instances[0]
	        logging.info('instance %s now running at %s', instance.id, instance.ip_address)
	        while time.time() - start_time < timeout:
	            cmd = ['ssh', '-o', 'StrictHostKeyChecking=no', '-i', key_file, 'ubuntu@' + instance.ip_address, 'echo']
	            if execute_command(cmd):
	                logging.info('instance %s responding %s', instance.id, instance.ip_address)
	                return instance.ip_address
	            time.sleep(2)
	        pass
	    time.sleep(2)

