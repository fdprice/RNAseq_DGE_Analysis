from __future__ import print_function
import os
import os.path
import sys
import gzip
import itertools
import SbatchController

## Parameters


sbatch_memreq = 8 ####
bwa_aln_flags = "-l 24"
batch_size = 5000000
r1_length = 16
#r2_length = 33
debug_flag = 1

## Helpers

# Open .fastq or .fastq.gz files for reading
def open_fastq_or_gz(filename):
    if filename.endswith(".fastq") and os.access(filename, os.F_OK):
        return open(filename, "rU")
    elif filename.endswith(".fastq.gz") and os.access(filename, os.F_OK):
        return gzip.open(filename, "rb")
    elif filename.endswith(".fastq") and os.access(filename + ".gz", os.F_OK):
        return gzip.open(filename + ".gz", "rb")
    elif filename.endswith(".fastq.gz") and os.access(filename[:-3], os.F_OK):
        return open(filename[:-3], "rU")
    raise IOError, "Unknown file: " + filename

# Contruct SLURM sbatch command for BWA alignment job
# def sbatch_bwa_cmd(fastq_path, reference_prefix, alignment_dir):
#     return " ".join(["sbatch", "-q", sbatch_queue, "-o",  os.path.join(alignment_dir, "bwa.sbatch"), "-J bwa", \
#                     "-R", sbatch_memreq, "-R", sbatch_ioreq, \
#                     "\"bwa aln", bwa_aln_flags, reference_prefix, fastq_path, "|", "bwa samse", \
#                     reference_prefix, "-", fastq_path, ">", ".".join([fastq_path, "sam"]) + "\""])
                    
def bwa_cmd(fastq_path, reference_prefix, alignment_dir):  ## separated bwa from sbatch
    return " ".join(["bwa aln", bwa_aln_flags, reference_prefix, fastq_path, "|", "bwa samse -n 20", \
                    reference_prefix, "-", fastq_path, ">", ".".join([fastq_path, "sam"])])


# Write fastq batch and submit BWA alignment job to SLURM
# def write_fastq_and_align(alignment_dir, sample_id, subsample_id, read_count, name_seq_qual_list, reference_prefix):
#     fastq_path = os.path.join(alignment_dir, ".".join([sample_id, subsample_id, str(read_count), "fastq"]))
#     with open(fastq_path, "w") as out:
#         for name, seq, qual in name_seq_qual_list:
#             print("\n".join(["@" + name, seq, "+", qual]), file=out)
#     os.system(sbatch_bwa_cmd(fastq_path, reference_prefix, alignment_dir))

def write_fastq_and_return_align_cmd(alignment_dir, sample_id, subsample_id, read_count, name_seq_qual_list, reference_prefix):
    fastq_path = os.path.join(alignment_dir, ".".join([sample_id, subsample_id, str(read_count), "fastq"]))
    with open(fastq_path, "w") as out:
        for name, seq, qual in name_seq_qual_list:
            print("\n".join(["@" + name, seq, "+", qual]), file=out)
    #os.system(sbatch_bwa_cmd(fastq_path, reference_prefix, alignment_dir))   
    return bwa_cmd(fastq_path, reference_prefix, alignment_dir)

# Mask sequence by quality score
def mask(seq, qual, min_qual=10):
    return "".join((b if (ord(q) - 33) >= min_qual else "N") for b, q in itertools.izip(seq, qual))

## Main

# Parse command line parameters
if len(sys.argv) != 8:
    print("Usage: python " + sys.argv[0] + " <sample id> <subsample id> <read1 .fastq or .fastq.gz> " +\
            "<read2 fastq or .fastq.gz> <alignment dir> <reference prefix> <sbatch_queue>",file=sys.stderr)
    sys.exit()

sample_id, subsample_id, r1_fastq, r2_fastq, alignment_dir, reference_prefix, sbatch_queue = sys.argv[1:]

# Read through R1 and R2 fastq files in parallel, add R1 barcode to R2 name, and launch batched
# BWA alignments via SLURM with output written to alignment_dir
print('### Reading and re-writing FASTQ data ###')
with open_fastq_or_gz(r1_fastq) as r1_file, open_fastq_or_gz(r2_fastq) as r2_file:
    read_count = 0
    buf = list()
    cmd_list = list() ### added

    r1_r2 = itertools.izip(r1_file, r2_file)
    for header1, header2 in r1_r2:
        seq1, seq2 = r1_r2.next()
        plus1, plus2 = r1_r2.next()
        qual1, qual2 = r1_r2.next()

        read_name1, read_name2 = header1.split()[0][1:], header2.split()[0][1:]
        assert read_name1 == read_name2
        seq2, qual2 = seq2.rstrip(), qual2.rstrip()
        barcode, seq, qual = mask(seq1[0:6], qual1[0:6], min_qual=10) + mask(seq1[6:r1_length], qual1[6:r1_length]), seq2, qual2
        barcoded_name = ":".join([read_name2, barcode])

        buf.append((barcoded_name, seq, qual))
        read_count += 1
        if read_count % batch_size == 0:
            cmd_to_add = write_fastq_and_return_align_cmd(alignment_dir, sample_id, subsample_id, read_count, buf, reference_prefix)
            cmd_list.append(cmd_to_add)
            buf = list()

    if len(buf) > 0:
    	cmd_to_add = write_fastq_and_return_align_cmd(alignment_dir, sample_id, subsample_id, read_count, buf, reference_prefix)
    	cmd_list.append(cmd_to_add)
	
	# print out status and debug messages
    print('#### FASTQs written. Starting Alignments Calls ####')
    print("%s command(s) to run" % len(cmd_list))
    if debug_flag:
        for index, item in enumerate(cmd_list):
            print(index, item)

	
    controller = SbatchController.SbatchController(cmd_list, queue=sbatch_queue, memory=sbatch_memreq, cmds_per_node=1, see = True, debug_flag=debug_flag ) # , mount_test=alignment_dir) ### removed mount_test as argument

    controller.run_slurm_submission()

    failed_cmds = controller.get_failed_cmds()
    if failed_cmds:
        for failed in failed_cmds:
            print(failed[0] + ' in job ' + str(failed[1]) + ' failed with ret ' + str(failed[2]) + ', see ' + str(failed[3]) + '.stdout and ' + str(failed[3]) + '.stderr')
        raise Exception('Ending process due to failed commands.')
    else:
        print('No failed commands.')

    if not debug_flag:
        controller.clean_logs()
