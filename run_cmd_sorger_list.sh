
# guarantee our consisten defaults
source new-modules.sh
module load python/2.7.6-fasrc01 
module load bwa/0.7.9a-fasrc02
module load fastqc/0.11.3-fasrc01
module load samtools/1.2-fasrc01

# load our Anaconda environment
source activate PYTHON_DGE


python ./dge-prod/Scripts/run_DGE_analysis.py --short_slurm_queue "short" --long_slurm_queue "long" --loose_barcodes --cleanup sample_map.txt Human Trugrade_384_set1 /n/regal/rubin_lab/fprice/DGE_processing/Alignment DGE_July15


