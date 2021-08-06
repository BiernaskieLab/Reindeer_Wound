#!/bin/bash
#
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --partition=bigmem
#SBATCH --time=24:00:00

export PATH=/work/biernaskie_lab/sarthak_sinha/cellrangeratacv120transfer/cellranger-atac_v1.2.0/cellranger-atac-1.2.0:$PATH
export PATH=/home/sarthak.sinha1/bcl2fastq2-v2.20/bin:$PATH



cellranger-atac aggr --id=Day_7_Ectopic_Aggr \
                 --csv=ATAC_aggr.csv \
                 --normalize=depth \
		 --reference=/work/biernaskie_lab/sarthak_sinha/cellrangeratacv120transfer/cellranger-atac_v1.2.0/refdata-cellranger-atac-GRCh38-1.2.0