# Ancestral Sequence Reconstruction (ASR)

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

## Purpose

This README details how we calculated the ASR sequences in the pub. Alignment was with [MAFFT](https://mafft.cbrc.jp/alignment/software/), tree building with [IQTree](https://iqtree.github.io/), reconstruction with codeml package from [PAML](https://github.com/abacus-gene/paml), gap inference with DOWNPASS algorithm from [PASTML](https://github.com/evolbioinfo/pastml).

## Setup

This repository uses conda to manage the computational and build environment. If you don’t have it installed (check with `conda --version`), you can find operating system-specific instructions for installing miniconda [here](https://www.anaconda.com/docs/getting-started/miniconda/main). After installing, run the following commands to create and activate the environment.

Create the environment from the provided [asr_env.yml](./asr_env.yml) file:

```bash
conda env create -n asr_env -f asr_env.yml
conda activate asr_env
```

## Reproduce

This publication uses ASR on both ADA1 and Isomaltase. The inputs/ouputs for each are found in `ADA1_ASR/` and `Isomaltase_ASR/`, respectively.

To reproduce the ASR analysis for ADA1, run:

```bash
papermill ASR_notebook.ipynb ASR_notebook_ADA1.ipynb -p input_fasta ADA1_ASR/inputs/ADA1_curated_022525_under420.fa -p outgroup_path ADA1_ASR/inputs/outgroup_names.txt
```

The output notebook will be written to `ASR_notebook_ADA1.ipynb`.

To reproduce the ASR analysis for Isomaltase, run:

```bash
papermill ASR_notebook.ipynb ASR_notebook_Isomaltase.ipynb -p input_fasta Isomaltase_ASR/inputs/isomaltase_dereplicated_final.fa -p outgroup_path Isomaltase_ASR/inputs/outgroup_names.txt
```

The output notebook will be written to `ASR_notebook_Isomaltase.ipynb`.

**Note** When finalizing the publication, we realized that the ASR workflow is non-deterministic because IQ-TREE uses a stochastic hill-climbing algorithm for tree search. The results presented in the paper were run with an unknown seed, so can't be reproduced exactly. If you are interested in fixing a seed, you can add a `-seed 1234` parameter to each `iqtree2` subprocess call in `ASR_notebook.ipynb`.

## Data Description

Key output files:

* **ancestor_tree.txt**: phylogenetic tree of input sequences with ancestral nodes labelled. You can visualize in FigTree or other tree viewing software to get ancestral node numbers of interest, or can visualize in the notebooks directly.
* **ML_ancestors.fa**: final ML sequences for all ancestral nodes
* **posterior_probabilities.json**: probabilities of all 20 amino acids at every position for all ancestors

Other output files:

* **recoding_dict.txt**: key for recoding of sequences to meet PAML 10 character seq id limit
* **ancestor_cladogram.txt**: cladogram of input sequences with ancestral nodes labelled
* **ancestor_recoded_cladogram.txt**: cladogram of input sequences recoded for paml
ancestor_recoded_tree.txt: phylogen* **ancestor_recoded_cladogram.txt**etic tree of input sequences recoded for paml
* **codeml.ctl**: config file for codeml, written by notebook from example file
* **ML_ancestors_with_gaps.fa**: ML sequences for all ancestral nodes before removing gaps with DOWNPASS
* **posterior_probabilities_no_gaps.json**: probabilities of all 20 amino acids at every position for all ancestors before gaps removed
* **rst**: raw paml output file parsed by the notebook to generate all other output files

## Compute Specifications

This was run locally on macOS (osx64) with 36 GM RAM.
