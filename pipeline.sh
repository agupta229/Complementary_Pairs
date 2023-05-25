#!/bin/bash
set -euo pipefail
if [ $# -lt 4 ]
then
  echo "Usage: $0 <input.bam> <forward.bg> <reverse.bg> <output_prefix>"
  exit
fi

IN=$1
FILE=$(basename $IN)
NAME=${FILE%%.*}
echo "#!/bin/bash -l
#SBATCH -J cmp_${NAME}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=32G

cd /homes9/akansha/meyersonlab/complementary_count

module load samtools
module load gcc/7.3.0 bedtools/2.29.0
./scripts/get_complementary_count_multiple.sh $1 $2 $3 $4
" > job_cmp_${NAME}.slurm
