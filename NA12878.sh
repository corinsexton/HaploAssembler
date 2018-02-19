#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=384G # memory per CPU core
#SBATCH -J "NA12878 full vcf"

module load samtools
#time ./FastHaploCaller.py -b ~/compute/NA12878/chr21_NA12878_S1.bam -v ~/compute/NA12878/chr21_S1_total.vcf -c chr21 > chr21_newvcf.total
time ./FastHaploCaller.py -b ~/compute/NA12878/chr21_NA12878_S1.bam -v ~/compute/NA12878/NA12878_GIAB.vcf -c chr21:9411318-9500318 > x
time ./FastHaploCaller.py -b ~/compute/NA12878/chr21_NA12878_S1.bam -v ~/compute/NA12878/chr21_NA12878_S1_noO.vcf -c chr21:9411318-9500318 > y
#time ./HaploCaller.py -b ~/compute/NA12878/chr21_NA12878_S1.bam -v ~/compute/NA12878/chr21_NA12878_S1_noO.vcf -c chr21:9411318-9503055 > nofast
#time ./HaploCaller.py -b ~/compute/NA12878/chr21_NA12878_S1.bam -v ~/compute/NA12878/chr21_NA12878_S1_noO.vcf -c chr21 > total.out

