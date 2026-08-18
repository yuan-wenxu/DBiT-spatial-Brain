#!/bin/bash

# RCTD reference
rctd_reference_dir=
rctd_reference_cache=
rctd_reference_barcode_column=V1
rctd_reference_numi_column=nCount_RNA
rctd_cell_type_column=C66_named
rctd_reference_dataset_column=Dataset
rctd_reference_gene_column=2
rctd_max_cells_per_type=200
rctd_cores=8
rctd_mode=full
rctd_ref_min_umi=100
rctd_spa_min_umi=100

# BANKSY expression integration
banksy_expression_hvg_count=3000
banksy_expression_components=20
banksy_expression_weight=0.3

# Slurm resources (used when execution_mode=hpc)
sbatch_domain_cpus=${rctd_cores}
sbatch_domain_partition=
sbatch_domain_mem=64G
sbatch_domain_time=24:00:00
