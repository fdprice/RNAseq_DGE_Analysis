from __future__ import print_function
import os
import os.path
import sys
import gzip
import itertools
import SbatchController
import subprocess

sbatch_memreq = 8

# Parse command line parameters
if len(sys.argv) != 3:
    print("Usage: python " + sys.argv[0] + " <sample map file> <qc dir>",file=sys.stderr)
    sys.exit()

sample_map_filename, qc_dir = sys.argv[1:]

# make QC dir
if not os.path.exists(qc_dir):
    print ("creating " + qc_dir + " QC directory\n")
    os.makedirs(qc_dir)

# make FastQC cmds
def qc_cmd(r_path, qc_dir):
	return " ".join(["fastqc -f fastq -o", qc_dir, r_path])

cmd_list = list()
with open(sample_map_filename, "rU") as sample_map:
    for line in sample_map:
        sample_id, subsample_id, r1_filename, r2_filename = line.split()
        for r_filename in [r1_filename, r2_filename]:
        	r_path = os.path.join(os.path.dirname(sample_map_filename), r_filename)
        	cmd_list.append(qc_cmd(r_path, qc_dir))
       
#print(cmd_list)

controller = SbatchController.SbatchController(cmd_list, memory=sbatch_memreq, cmds_per_node=1, mount_test=qc_dir)
controller.run_slurm_submission()

failed_cmds = controller.get_failed_cmds()
if failed_cmds:
	for failed in failed_cmds:
		print(failed[0] + ' in job ' + str(failed[1]) + ' failed with ret ' + str(failed[2]) + ', see ' + str(failed[3]) + '.stdout and ' + str(failed[3]) + '.stderr')
	raise Exception('Ending process due to failed QC commands.')	
else:
	print('No failed QC commands, finished.')

controller.clean_logs()
