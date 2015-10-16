#!/bin/bash
#
#SBATCH -p general                              # partition (queue)
#SBATCH -N 1                                    # number of nodes
#SBATCH -n 1                                    # number of cores
#SBATCH --mem 1000                              # memory pool for all cores
#SBATCH -t 2-0:00                               # time (D-HH:MM)
#SBATCH -o DGE_run.%N.%j.out                    # STDOUT
#SBATCH -e DGE_run.%N.%j.err                    # STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=feodor_price@harvard.edu    # send-to address

# guarantee our consistent defaults
source new-modules.sh
module load python/2.7.6-fasrc01 
module load bwa/0.7.9a-fasrc01
module load fastqc/0.11.3-fasrc01
module load samtools/1.2-fasrc01

# load our Anaconda environment
source activate PYTHON_DGE

# template script. Assumes this script will be in run directory
python $HOME/git/RNAseq_DGE_Analysis/run_DGE_analysis.py \
   --short_slurm_queue "short" --long_slurm_queue "long" \
   --loose_barcodes --cleanup \
   sample_map.txt \
   Human \
   Trugrade_384_set1 \
   /n/regal/rubin_lab/fprice/DGE_processing/Alignment \
   DGE_OCT15

