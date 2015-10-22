from __future__ import print_function
import os
import os.path
import sys
import SbatchController
import argparse
import time

print(time.strftime("%c"))
print(" ".join(sys.argv))

parser = argparse.ArgumentParser()

parser.add_argument('sample_map', help='location of sample map file', type=str)
parser.add_argument('reference', help='reference genome: Human|Mouse|Rat|Chicken', type=str)
parser.add_argument('barcodes', help='barcode plate: P1|P2||P3|P1P2|Trugrade_384_set1|Trugrade_96_set1|Trugrade_96_set2|Trugrade_96_set3|Trugrade_96_set4|Truegrade_96_1234', type=str)
parser.add_argument('alignment_dir', help='directory to process alignments', type=str)
parser.add_argument('analysis_dir', help='directory to calculate gene expression', type=str)
parser.add_argument('--short_slurm_queue', help='slurm_queue to run short jobs. default=short', type=str)
parser.add_argument('--long_slurm_queue', help='slurm_queue to run long jobs. default=long', type=str)
parser.add_argument('--loose_barcodes', help='allows well barcodes to have 1 mismatch', action='store_true')
parser.add_argument('--cleanup', help='removes resulting fq and sam files upon successful completion', action='store_true') 

args = parser.parse_args()

sample_map_filename = args.sample_map
species = args.reference
barcode_plate = args.barcodes
alignment_dir = args.alignment_dir
dge_dir = args.analysis_dir
loose_barcodes = str(args.loose_barcodes)
cleanup = args.cleanup

pythoncmd = 'python';

sbatch_queue = "long"
short_queue = "short"

if args.short_slurm_queue:
	short_queue = args.short_slurm_queue
if args.long_slurm_queue:
	sbatch_queue = args.long_slurm_queue

sbatch_memreq = 32		################# previously 4
debug_flag = 1          # whether or not to show debug messages

bindir = os.path.abspath(os.path.dirname(sys.argv[0]))		# locate scripts directory
reference_dir =  os.path.join(bindir, "../Reference")		# locate "Reference" dir relative to scripts dir 

reference_prefix_map = {"Human": os.path.join(reference_dir, "Human_RefSeq", "refMrna_ERCC_polyAstrip.hg19.fa"),  
                        "Mouse": os.path.join(reference_dir, "Mouse_RefSeq", "refMrna_ERCC_polyAstrip.mm10.fa"),
                        "Chicken": os.path.join(reference_dir, "Chicken_RefSeq", "refMrna_ERCC_polyAstrip.gg4.fa"), 
                        "Rat": os.path.join(reference_dir, "Rat_RefSeq", "refMrna_ERCC_polyAstrip.rn5.fa")}     
			
sym2ref = {"Human": os.path.join(reference_dir, "Human_RefSeq", "refGene.hg19.sym2ref.dat"),
           "Mouse": os.path.join(reference_dir, "Mouse_RefSeq", "refGene.mm10.sym2ref.dat"),
           "Chicken": os.path.join(reference_dir, "Chicken_RefSeq", "refGene.gg4.sym2ref.dat"),
	   "Rat": os.path.join(reference_dir, "Rat_RefSeq", "refGene.rn5.sym2ref.dat")}

barcodes = {"P1": os.path.join(reference_dir, "barcodes_plate1.dat"),
            "P2": os.path.join(reference_dir, "barcodes_plate2.dat"),
            "P3": os.path.join(reference_dir, "barcodes_plate3.dat"),
            "P1P2": os.path.join(reference_dir, "barcodes_plates12.dat"),
			"Trugrade_384_set1": os.path.join(reference_dir, "barcodes_trugrade_384_set1.dat"),
			"Trugrade_96_set1": os.path.join(reference_dir, "barcodes_trugrade_96_set1.dat"),
			"Trugrade_96_set2": os.path.join(reference_dir, "barcodes_trugrade_96_set2.dat"),
			"Trugrade_96_set3": os.path.join(reference_dir, "barcodes_trugrade_96_set3.dat"),
			"Trugrade_96_set4": os.path.join(reference_dir, "barcodes_trugrade_96_set4.dat"),
			"Trugrade_96_set1234": os.path.join(reference_dir, "barcodes_trugrade_96_sets1234_combined.dat")}

ercc_fasta = os.path.join(reference_dir, "ERCC92.fa")

# Error handling:  Check for existence of Reference directory
if not os.path.exists(reference_dir):
    raise Exception("\nFatal error: " + reference_dir + 
        " not found.\n\"Reference\" directory must exist within" + bindir + "\n")

# Error handling:  Check for existence of alignment directory.  Create directory if not found.
if not os.path.exists(alignment_dir):
    print ("creating " + alignment_dir + " alignment directory\n")
    os.makedirs(alignment_dir)

# Error handling:  Check for existence of dge_dir directory.  Create directory if not found.
if not os.path.exists(dge_dir):
    print ("creating " + dge_dir + " dge directory\n")
    os.makedirs(dge_dir)

#  Error handling:  Check for existence of sample_map_filename
if not ( os.path.isfile(sample_map_filename) ):                         # check for sample_map_filename
    raise Exception("\n\nFatal error: " + sample_map_filename + " sample map not found.\n")

# Error handling:  Check for species matching either [Mm]ouse or [Hh]uman or [Rr]at or [Cc]hicken
species = species.capitalize()
if species in ("Mouse", "Human", "Rat", "Chicken"):
    print ("Species: " + species + "\n")
else:
    raise Exception("\n\nFatal error: " + species + " not a recognized species.\n")

reference_prefix = reference_prefix_map[species]
sym2ref = sym2ref[species]
barcodes = barcodes[barcode_plate]

#
## run split and align
#

# Construct command for split_and_align job
def sa_cmd(sample_id, subsample_id, r1_path, r2_path, alignment_dir, reference_prefix, short_queue):
	split_call = "".join([bindir, "/split_and_align2.py"])
	return " ".join([pythoncmd, split_call, sample_id, subsample_id, r1_path, 
	    r2_path, alignment_dir, reference_prefix, short_queue])

# Read though sample map and launch split_and_align jobs for each entry
cmd_list = list()
with open(sample_map_filename, "rU") as sample_map:
    for line in sample_map:
        sample_id, subsample_id, r1_filename, r2_filename = line.split()
	### Error handling for fastq files
	if r1_filename == r2_filename:
		raise Exception("\nFatal error: Read 1 filename shhould not equal read 2 filename\n")
	if not ( os.path.isfile(r1_filename) ):
		raise Exception("\nFatal error: " + r1_filename + " not found. Check " + sample_map_filename + "\n")
	elif not r1_filename.endswith(".fastq") and not r1_filename.endswith(".fastq.gz"):
		raise Exception("\nFatal error: " + r1_filename + " must end with filename extension of .fastq or .fastq.gz\n")
	else:
	    print (subsample_id + " " + r1_filename + " input fastq file found")   
	if not ( os.path.isfile(r2_filename) ):
	    raise Exception("\nFatal error: " + r2_filename + " not found. Check " + sample_map_filename + "\n")
	elif not r2_filename.endswith(".fastq") and not r2_filename.endswith(".fastq.gz"):
	    raise Exception("\nFatal error: " + r2_filename + " must end with filename extension of .fastq or .fastq.gz\n")
	else :
	    print (subsample_id + " " + r2_filename + " input fastq file found")
        ### 

        r1_path = os.path.join(os.path.dirname(sample_map_filename), r1_filename)
        r2_path = os.path.join(os.path.dirname(sample_map_filename), r2_filename)
        cmd_list.append(sa_cmd(sample_id, subsample_id, r1_path, r2_path, alignment_dir, 
            reference_prefix, short_queue))

# print out status and debug messages
print('#### Split and Align ####')
print("%s command(s) to run" % len(cmd_list))
if debug_flag:
    for index, item in enumerate(cmd_list):
        print(index, item)

controller = SbatchController.SbatchController(cmd_list, queue=sbatch_queue, 
    memory=sbatch_memreq, cmds_per_node=1, see=True, debug_flag=debug_flag) # , mount_test=alignment_dir) ### removed mount_test as argument
controller.run_slurm_submission()

failed_cmds = controller.get_failed_cmds()
if failed_cmds:
	for failed in failed_cmds:
		print(failed[0] + ' in job ' + str(failed[1]) + ' failed with ret ' 
		    + str(failed[2]) + ', see ' + str(failed[3]) + '.stdout and ' + str(failed[3]) + '.stderr')
	raise Exception('Ending process due to failed alignment commands.')	
else:
	print('No failed aligment commands, continuing.')

if not debug_flag:
    controller.clean_logs()
    
#
## run merge and count
#

sbatch_memreq = 80   ################# previously 8


## Helpers

# Contruct slurm sbatch command for merge_and_count job
def mc_cmd(sample_id, sym2ref, ercc_fasta, barcodes, alignment_dir, dge_dir, loose_barcodes):
	merge_call = "".join([bindir, "/merge_and_count.py"])
	return " ".join([pythoncmd, merge_call, sample_id, sym2ref, ercc_fasta, barcodes, alignment_dir, dge_dir, loose_barcodes])

# Read though sample map and launch merge_and_count jobs for each sample_id
merge_cmd_list = list()
with open(sample_map_filename, "rU") as sm_file:
    for sample_id in set([line.split()[0] for line in sm_file]):
        merge_cmd_list.append(mc_cmd(sample_id, sym2ref, ercc_fasta, barcodes, alignment_dir, dge_dir, loose_barcodes))

# print out status and potential debug messages
print('#### Merge and Count ####')
print("%s command(s) to run" % len(merge_cmd_list))
if debug_flag:
    for index, item in enumerate(merge_cmd_list):
        print(index, item)


controller2 = SbatchController.SbatchController(merge_cmd_list, queue=sbatch_queue, memory=sbatch_memreq, cmds_per_node=1, debug_flag=debug_flag) # , mount_test=dge_dir)  ### removed mount_test as argument
controller2.run_slurm_submission()


failed_cmds = controller2.get_failed_cmds()
success = False
if failed_cmds:
	for failed in failed_cmds:
		print(failed[0] + ' in job ' + str(failed[1]) + ' failed with ret ' + str(failed[2]) + ', see ' + str(failed[3]) + '.stdout and ' + str(failed[3]) + '.stderr')
else:
	print('No failed merge and count commands, finished.')
	success = True

if success and not debug_flag:
	controller2.clean_logs()
	
	if cleanup and not debug_flag:
		for align_file in [f for f in os.listdir(alignment_dir) if (f.endswith(".fastq") or f.endswith(".sam"))]:
			aligned_path = os.path.join(alignment_dir, align_file)
			os.remove(aligned_path)
		print('Alignment files removed.')


