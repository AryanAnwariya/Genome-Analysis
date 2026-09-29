#!/bin/bash
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
NAME=$1; ORDER=$2; SAMPLES=$3
export APPTAINER_CACHEDIR=/scratch/adas/apptainer/cache APPTAINER_TMPDIR=/scratch/adas/apptainer/tmp
W=/scratch/adas/genome-analysis/results/container_demo_$NAME
mkdir -p $W
cd $W
ln -sfn /sci/LINE LINE
ln -sfn /sci/chromosome_sizes chromosome_sizes
apptainer exec --bind $W:/work /scratch/adas/apptainer/sci_container.sif bash -c "
cd /work
export LD_LIBRARY_PATH=/sci/gsl/lib:\${LD_LIBRARY_PATH:-}
. /root/miniconda3/etc/profile.d/conda.sh
conda activate sci
echo start: \$(date)
python -m sci.sci -n demo -f /sci-data/demo_data.txt -r 100000 -g chromosome_sizes/hg19.chrom.sizes -o $ORDER -s $SAMPLES -k 5
echo end: \$(date)
ls -l /work
"
