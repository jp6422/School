#!/bin/bash
#SBATCH --job-name=salmon # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jp6422@nyu.edu # Where to send mail
#SBATCH --ntasks=4 # Run on a single CPU
#SBATCH --mem=64gb # Job memory request
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH -p cpu_medium
#SBATCH --output=/gpfs/scratch/jp6422/final/salmon.log # Standard output and error log

module load salmon/1.4.0

salmon index -i sindex -t /gpfs/data/courses/bminga3004/2022/Practicum5/Homo_sapiens.GRCh38.cdna.all.fa


salmon quant -i /gpfs/scratch/jp6422/final/sindex -o /gpfs/scratch/jp6422/final/kallisto/pramehi-ctrl1/ -l A -1 /gpfs/scratch/jp6422/final/PRAMEhi-ctrl1_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMEhi-ctrl1_R2_trimmed.fq 
salmon quant -i /gpfs/scratch/jp6422/final/sindex -o /gpfs/scratch/jp6422/final/kallisto/pramehi-ctrl2/ -l A -1 /gpfs/scratch/jp6422/final/PRAMEhi-ctrl2_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMEhi-ctrl2_R2_trimmed.fq 
salmon quant -i /gpfs/scratch/jp6422/final/sindex -o /gpfs/scratch/jp6422/final/kallisto/pramehi-ctrl3/ -l A -1 /gpfs/scratch/jp6422/final/PRAMEhi-ctrl3_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMEhi-ctrl3_R2_trimmed.fq 
salmon quant -i /gpfs/scratch/jp6422/final/sindex -o /gpfs/scratch/jp6422/final/kallisto/pramehi-sirna1 -l A -1 /gpfs/scratch/jp6422/final/PRAMEhi-siRNA1_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMEhi-siRNA1_R2_trimmed.fq
salmon quant -i /gpfs/scratch/jp6422/final/sindex -o /gpfs/scratch/jp6422/final/kallisto/pramehi-sirna2 -l A -1 /gpfs/scratch/jp6422/final/PRAMEhi-siRNA2_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMEhi-siRNA2_R2_trimmed.fq
salmon quant -i /gpfs/scratch/jp6422/final/sindex -o /gpfs/scratch/jp6422/final/kallisto/pramehi-sirna3 -l A -1 /gpfs/scratch/jp6422/final/PRAMEhi-siRNA3_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMEhi-siRNA3_R2_trimmed.fq
salmon quant -i /gpfs/scratch/jp6422/final/sindex  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-ctrl1 -l A -1 /gpfs/scratch/jp6422/final/PRAMElo-ctrl1_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMElo-ctrl1_R2_trimmed.fq
salmon quant -i /gpfs/scratch/jp6422/final/sindex  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-ctrl2 -l A -1 /gpfs/scratch/jp6422/final/PRAMElo-ctrl2_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMElo-ctrl2_R2_trimmed.fq
salmon quant -i /gpfs/scratch/jp6422/final/sindex  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-ctrl3 -l A -1 /gpfs/scratch/jp6422/final/PRAMElo-ctrl3_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMElo-ctrl3_R2_trimmed.fq
salmon quant -i /gpfs/scratch/jp6422/final/sindex  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-oe1  -l A -1 /gpfs/scratch/jp6422/final/PRAMElo-oe1_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMElo-oe1_R2_trimmed.fq
salmon quant -i /gpfs/scratch/jp6422/final/sindex  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-oe2 -l A -1 /gpfs/scratch/jp6422/final/PRAMElo-oe2_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMElo-oe2_R2_trimmed.fq
salmon quant -i /gpfs/scratch/jp6422/final/sindex  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-oe3 -l A -1 /gpfs/scratch/jp6422/final/PRAMElo-oe3_R1_trimmed.fq -2 /gpfs/scratch/jp6422/final/PRAMElo-oe3_R2_trimmed.fq



