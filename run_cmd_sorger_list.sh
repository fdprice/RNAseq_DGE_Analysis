

module add seq/bwa/0.7.8

python ./dge-prod/Scripts/run_DGE_analysis.py --short_lsf_queue "sorger_12h" --long_lsf_queue "sorger_unlimited" --loose_barcodes --cleanup sample_map.txt Human Trugrade_384_set1 /groups/sorger/Marc/DGE_processing/Alignment DGE_July15
