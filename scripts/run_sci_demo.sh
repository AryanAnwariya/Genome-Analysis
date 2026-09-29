#!/bin/bash
#SBATCH --job-name=sci_demo
#SBATCH --partition=compute
#SBATCH --qos=cpulimit
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/adas/genome-analysis/logs/sci_demo_%j.log
export APPTAINER_CACHEDIR=/scratch/adas/apptainer/cache APPTAINER_TMPDIR=/scratch/adas/apptainer/tmp
W=/scratch/adas/genome-analysis/results/container_demo
mkdir -p $W
cd $W
ln -sfn /sci/LINE LINE
ln -sfn /sci/chromosome_sizes chromosome_sizes
apptainer exec --bind $W:/work /scratch/adas/apptainer/sci_container.sif bash -c '
cd /work
export LD_LIBRARY_PATH=/sci/gsl/lib:${LD_LIBRARY_PATH:-}
for c in /opt/conda /opt/miniconda3 /opt/anaconda3 /root/miniconda3 /root/anaconda3 /usr/local/conda /miniconda /conda; do
  if [ -f $c/etc/profile.d/conda.sh ]; then . $c/etc/profile.d/conda.sh; echo "conda: $c"; break; fi
done
conda activate sci 2>&1 || echo "conda activate failed"
which python; python --version 2>&1
python -c "import sci; print(sci.__file__)" 2>&1
echo "=== run start: $(date)"
python -m sci.sci -n demo -f /sci-data/demo_data.txt -r 100000 -g chromosome_sizes/hg19.chrom.sizes -o both -s 25 -k 5
echo "=== run end: $(date)"
ls -l /work
'
