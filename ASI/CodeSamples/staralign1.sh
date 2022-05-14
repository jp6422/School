#!/bin/bash
#SBATCH --job-name=aligner2 # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jonah.poczobutt@nyulangone.org # Where to send mail
#SBATCH --ntasks=1 # Run on a single nodes
#SBATCH --mem=64gb # Job memory request
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=/gpfs/scratch/jp6422/final/align2.log # Standard output and error log
#SBATCH -p cpu_medium

module load gcc/10.2.0
module load star/2.7.7a 

cd /gpfs/scratch/jp6422/final

#trim here
module load trimgalore/0.5.0

module load python/cpu/2.7.15-ES

trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-ctrl1_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-ctrl1_R2.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-ctrl2_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-ctrl2_R2.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-ctrl3_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-ctrl3_R2.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-siRNA1_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-siRNA1_R2.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-siRNA2_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-siRNA2_R2.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-siRNA3_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMEhi-siRNA3_R2.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMElo-ctrl1_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMElo-ctrl1_R2.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMElo-ctrl2_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMElo-ctrl2_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMElo-ctrl3_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMElo-ctrl3_R2.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/FinalProjects/Project2/PRAMElo-oe1_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/Finalprojects/Project2/PRAMElo-oe1_R2.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/Finalprojects/Project2/PRAMElo-oe2_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/Finalprojects/Project2/PRAMElo-oe2_R2.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/Finalprojects/Project2/PRAMElo-oe3_R1.fastq.gz  –q 30
trim_galore /gpfs/data/courses/bminga3004/2022/Finalprojects/Project2/PRAMElo-oe3_R2.fastq.gz  –q 30

#unzip files
gunzip PRAMEhi-ctrl1_R1_trimmed.fastq.gz
gunzip PRAMEhi-ctrl1_R2_trimmed.fastq.gz
gunzip PRAMEhi-ctrl2_R1_trimmed.fastq.gz
gunzip PRAMEhi-ctrl2_R2_trimmed.fastq.gz
gunzip PRAMEhi-ctrl3_R1_trimmed.fastq.gz
gunzip PRAMEhi-ctrl3_R1_trimmed.fastq.gz
gunzip PRAMEhi-siRNA1_R1_trimmed.fastq.gz
gunzip PRAMEhi-siRNA1_R2_trimmed.fastq.gz
gunzip PRAMEhi-siRNA2_R1_trimmed.fastq.gz
gunzip PRAMEhi-siRNA2_R2_trimmed.fastq.gz
gunzip PRAMEhi-siRNA3_R1_trimmed.fastq.gz
gunzip PRAMEhi-siRNA3_R2_trimmed.fastq.gz
gunzip PRAMElo-ctrl1_R1_trimmed.fastq.gz
gunzip PRAMElo-ctrl1_R2_trimmed.fastq.gz
gunzip PRAMElo-ctrl2_R1_trimmed.fastq.gz
gunzip PRAMElo-ctrl2_R2_trimmed.fastq.gz
gunzip PRAMElo-ctrl3_R1_trimmed.fastq.gz
gunzip PRAMElo-ctrl3_R2_trimmed.fastq.gz
gunzip PRAMElo-oe1_R1_trimmed.fastq.gz 
gunzip PRAMElo-oe1_R2_trimmed.fastq.gz
gunzip PRAMElo-oe2_R1_trimmed.fastq.gz
gunzip PRAMElo-oe2_R2_trimmed.fastq.gz
gunzip PRAMElo-oe3_R1_trimmed.fastq.gz
gunzip PRAMElo-oe3_R2_trimmed.fastq.gz




    

#generate index

STAR --runMode genomeGenerate --genomeFastaFiles /gpfs/scratch/jp6422/Practicum9/hg38.fa --genomeDir GenomeDir


#align unzipped fq.gz files for all samples

STAR --genomeDir GenomeDir --readFilesIn PRAMEhi-ctrl1_R1_trimmed.fq PRAMEhi-ctrl1_R2_trimmed.fq --outFileNamePrefix PRAMEhi-ctrl1
STAR --genomeDir GenomeDir --readFilesIn PRAMEhi-ctrl2_R1_trimmed.fq PRAMEhi-ctrl2_R2_trimmed.fq --outFileNamePrefix PRAMEhi-ctrl2
STAR --genomeDir GenomeDir --readFilesIn PRAMEhi-ctrl3_R1_trimmed.fq PRAMEhi-ctrl3_R2_trimmed.fq --outFileNamePrefix PRAMEhi-ctrl3
STAR --genomeDir GenomeDir --readFilesIn PRAMEhi-siRNA1_R1_trimmed.fq PRAMEhi-siRNA1_R2_trimmed.fq --outFileNamePrefix PRAMEhi-siRNA1
STAR --genomeDir GenomeDir --readFilesIn PRAMEhi-siRNA2_R1_trimmed.fq PRAMEhi-siRNA2_R2_trimmed.fq --outFileNamePrefix PRAMEhi-siRNA2
STAR --genomeDir GenomeDir --readFilesIn PRAMEhi-siRNA3_R1_trimmed.fq PRAMEhi-siRNA3_R2_trimmed.fq --outFileNamePrefix PRAMEhi-siRNA3
STAR --genomeDir GenomeDir --readFilesIn PRAMElo-ctrl1_R1_trimmed.fq PRAMElo-ctrl1_R2_trimmed.fq --outFileNamePrefix PRAMElo-ctrl1
STAR --genomeDir GenomeDir --readFilesIn PRAMElo-ctrl2_R1_trimmed.fq PRAMElo-ctrl2_R2_trimmed.fq --outFileNamePrefix PRAMElo-ctrl2
STAR --genomeDir GenomeDir --readFilesIn PRAMElo-ctrl3_R1_trimmed.fq PRAMElo-ctrl3_R2_trimmed.fq --outFileNamePrefix PRAMElo-ctrl3
STAR --genomeDir GenomeDir --readFilesIn PRAMElo-oe1_R1_trimmed.fq PRAMElo-oe1_R2_trimmed.fq --outFileNamePrefix PRAMElo-oe1
STAR --genomeDir GenomeDir --readFilesIn PRAMElo-oe2_R1_trimmed.fq PRAMElo-oe2_R2_trimmed.fq --outFileNamePrefix PRAMElo-oe2
STAR --genomeDir GenomeDir --readFilesIn PRAMElo-oe3_R1_trimmed.fq PRAMElo-oe3_R2_trimmed.fq --outFileNamePrefix PRAMElo-oe3


