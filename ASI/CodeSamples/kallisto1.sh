#!/bin/bash 
#SBATCH --job-name=Kallisto # Job name 
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=jp6422@nyu.edu # Where to send mail 
#SBATCH --ntasks=4 # Run on a single CPU 
#SBATCH --mem=64gb # Job memory request 
#SBATCH --time=24:00:00 # Time limit hrs:min:sec 
#SBATCH -p cpu_medium
#SBATCH --output=/gpfs/scratch/jp6422/final/kallisto.log # Standard output and error log



#make sure you cd into your directory and change all paths to match your own paths
cd /gpfs/scratch/jp6422/final/


module load kallisto/0.44.0

### Indexing needed for Kallisto, point to the one in course directory or create your own 
kallisto index -i index.idx  /gpfs/data/courses/bminga3004/2022/Practicum5/Homo_sapiens.GRCh38.cdna.all.fa

#run kallisto

kallisto quant -i /gpfs/scratch/jp6422/final/index.idx -o /gpfs/scratch/jp6422/final/kallisto/pramehi-ctrl1/ -b 100/ --bias  /gpfs/scratch/jp6422/final/PRAMEhi-ctrl1_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMEhi-ctrl1_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx -o /gpfs/scratch/jp6422/final/kallisto/pramehi-ctrl2/ -b 100 --bias /gpfs/scratch/jp6422/final/PRAMEhi-ctrl2_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMEhi-ctrl2_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx -o /gpfs/scratch/jp6422/final/kallisto/pramehi-ctrl3/ -b 100 --bias /gpfs/scratch/jp6422/final/PRAMEhi-ctrl3_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMEhi-ctrl3_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx -o /gpfs/scratch/jp6422/final/kallisto/pramehi-sirna1 -b 100 --bias /gpfs/scratch/jp6422/final/PRAMEhi-siRNA1_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMEhi-siRNA1_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx -o /gpfs/scratch/jp6422/final/kallisto/pramehi-sirna2 -b 100 --bias /gpfs/scratch/jp6422/final/PRAMEhi-siRNA2_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMEhi-siRNA2_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx -o /gpfs/scratch/jp6422/final/kallisto/pramehi-sirna3 -b 100 --bias /gpfs/scratch/jp6422/final/PRAMEhi-siRNA3_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMEhi-siRNA3_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-ctrl1 -b 100 --bias /gpfs/scratch/jp6422/final/PRAMElo-ctrl1_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMElo-ctrl1_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-ctrl2 -b 100 --bias /gpfs/scratch/jp6422/final/PRAMElo-ctrl2_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMElo-ctrl2_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-ctrl3 -b 100 --bias /gpfs/scratch/jp6422/final/PRAMElo-ctrl3_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMElo-ctrl3_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-oe1 -b 100 --bias /gpfs/scratch/jp6422/final/PRAMElo-oe1_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMElo-oe1_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-oe2 -b 100 --bias /gpfs/scratch/jp6422/final/PRAMElo-oe2_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMElo-oe2_R2_trimmed.fq
kallisto quant -i /gpfs/scratch/jp6422/final/index.idx  -o /gpfs/scratch/jp6422/final/kallisto/pramelo-oe3 -b 100 --bias /gpfs/scratch/jp6422/final/PRAMElo-oe3_R1_trimmed.fq /gpfs/scratch/jp6422/final/PRAMElo-oe3_R2_trimmed.fq



 





