# Notebook Publication Template

This repo is a template for Jupyter notebook publications. The template produces a publication rendered and hosted by Quarto, which can be viewed at [this demo URL](https://arcadia-science.github.io/notebook-pub-template/).

## Template Documentation

All the learning resources for this template can be found in `developer-docs/`.

- [Quickstart Guide](developer-docs/QUICKSTART.md) - **The most efficient way to get started** is to follow this guide (the rest can wait)
- [Environment Setup Guide](developer-docs/ENVIRONMENT_SETUP.md) - How to set up your development environment
- [Publishing Guide](developer-docs/PUBLISHING_GUIDE.md) - How to publish your notebook publication
- [Template Architecture](developer-docs/TEMPLATE_ARCHITECTURE.md) - Understanding the template's structure

---

**NOTE: When ready to publish, fill in the information below, then delete this line and everything above it.**

# Exploring ESM2 Pseudo-perplexity of Maximum Likelihood Ancestral Protein Sequences

This code repository contains or points to all materials required for creating and hosting the publication entitled, *"[PUB-TITLE]"*.

The publication is hosted at [this URL]([PUB-URL]).

## Data Description

**Ancestral Sequence Generation**

Ancestral sequences were reconstructed using the workflow in [ASR/ASR_notebook.ipynb](ASR/ASR_notebook.ipynb).  All sequence inputs and outputs are in the [ASR/](ASR/) directory.

**ESM2 Pseudo-perplexity Calculation**

ESM2 pseudo-perplexity scores were calculated using the [esm2_pppl_calculator.py](esm2_pppl_calculator.py) on a GPU-enable AWS EC2 instance with the following specifications to ensure full reproducibility:

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

Create the Conda environment from the provided [esm2_env.yml] file:
```bash
conda env create -n esm_env -f esm2_env.yml
conda activate esm_env
```

Calculate pseudo-perplexity scores using ESM2:
```bash
python esm2_pppl_calculator.py --input consensus_ADA1_ancestors.fa --output consensus_ADA1_ancestors_scores.csv --model esm2_t33_650M_UR50D 
```

## Reproduce

Please see [SETUP.qmd](SETUP.qmd).

## Contribute

Please see [CONTRIBUTING.qmd](CONTRIBUTING.qmd).
