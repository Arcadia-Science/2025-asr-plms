# ESM2 Scoring

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

## Purpose

Code for calculating ESM2 pseudo-perplexity scores for protein sequence for different ESM2 models.

## Setup

ESM2 pseudo-perplexity scores were calculated using [esm2_pppl_calculator.py](esm2_pppl_calculator.py) on a GPU-enable AWS EC2 instance with the following specifications:

AMI Name: **Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 22.04)**  
Instance Type: **g4dn.2xlarge (NVIDIA T4 GPU, 32 GiB memory)**  
Root Volume Size: **100 GiB**  
Storage Type: **EBS (General Purpose SSD, gp3)**

Once launched, set up the Python environment:

Install Miniconda:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

Create the Conda environment from the provided [esm2_env.yml](esm2_env.yml) file:

```bash
conda env create -n esm_env -f esm2_env.yml
conda activate esm_env
```

## Usage

Calculate pseudo-perplexity scores using ESM2:

```bash
python esm2_pppl_calculator.py --input ../ASR/ADA1_ASR/outputs/consensus_ADA1_ancestors.fa --output consensus_ADA1_ancestors_scores.csv --model esm2_t33_650M_UR50D 
```

## Reproduce

To reproduce the ESM2 pseudo-perplexity scores for the ADA1 and isomaltase sequences, run the following commands:

```bash
python esm2_pppl_calculator.py --input ADA1_native_and_anc.fasta --output ADA1_all_esm2_scores_8M.csv --model esm2_t6_8M_UR50D 
python esm2_pppl_calculator.py --input ADA1_native_and_anc.fasta --output ADA1_all_esm2_scores_35M.csv --model esm2_t12_35M_UR50D 
python esm2_pppl_calculator.py --input ADA1_native_and_anc.fasta --output ADA1_all_esm2_scores_150M.csv --model esm2_t30_150M_UR50D 
python esm2_pppl_calculator.py --input ADA1_native_and_anc.fasta --output ADA1_all_esm2_scores_650M.csv --model esm2_t33_650M_UR50D 
python esm2_pppl_calculator.py --input ADA1_native_and_anc.fasta --output ADA1_all_esm2_scores_3B.csv --model esm2_t36_3B_UR50D 

python esm2_pppl_calculator.py --input consensus_ancestors/consensus_ADA1_ancestors.fa --output consensus_ADA1_ancestors_esm2_scores_650M.csv --model esm2_t33_650M_UR50D 
python esm2_pppl_calculator.py --input consensus_ancestors/consensus_isomaltase_ancestors.fa --output consensus_isomaltase_ancestors_esm2_scores_650M.csv --model esm2_t33_650M_UR50D

python esm2_pppl_calculator.py --input isomaltase_native_and_anc.fasta --output isomaltase_all_esm2_scores_8M --model esm2_t6_8M_UR50D
python esm2_pppl_calculator.py --input isomaltase_native_and_anc.fasta --output isomaltase_all_esm2_scores_35M --model esm2_t12_35M_UR50D
python esm2_pppl_calculator.py --input isomaltase_native_and_anc.fasta --output isomaltase_all_esm2_scores_150M --model esm2_t30_150M_UR50D
python esm2_pppl_calculator.py --input isomaltase_native_and_anc.fasta --output isomaltase_all_esm2_scores_3B --model esm2_t36_3B_UR50D
python esm2_pppl_calculator.py --input isomaltase_native_and_anc.fasta --output isomaltase_all_esm2_scores_650M --model esm2_t33_650M_UR50D

python esm2_pppl_calculator.py --input ADA1_remaining.fasta --output ADA1_remaining_esm2_scores_650M --model esm2_t33_650M_UR50D
python esm2_pppl_calculator.py --input isomaltase_remaining.fasta --output isomaltase_remaining_esm2_scores_650M --model esm2_t33_650M_UR50D
```
