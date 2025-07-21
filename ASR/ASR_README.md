# 2025-asr

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)


## Purpose

Method for Ancestral Sequence Reconstruction (ASR) of protein sequences. Alignment with MAFFT, tree building with IQTree, reconstruction with codeml package from PAML, gap inference with DOWMPASS algorithm from PASTML.

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
TODO: Replace <NAME> with the name of your environment
mamba env create -n <NAME> --file envs/dev.yaml
conda activate <NAME>
```

<details><summary>Developer Notes (click to expand/collapse)</summary>

1. Install your pre-commit hooks:

    ```{bash}
    pre-commit install
    ```

    This installs the pre-commit hooks defined in your config (`./.pre-commit-config.yaml`).

2. Export your conda environment before sharing:

    As your project develops, the number of dependencies in your environment may increase. Whenever you install new dependencies (using either `pip install` or `mamba install`), you should update the environment file using the following command.

    ```{bash}
    conda env export --no-builds > envs/dev.yml
    ```

    `--no-builds` removes build specification from the exported packages to increase portability between different platforms.
</details>

## Data

Inputs: 
(1) Fasta file of input sequences including known outgroup sequence(s)
(2) Text file with seq ids for outgroup sequences

## Overview

### Description of the folder structure

### Methods

Update second cell in ASR_notebook.ipynb with path to input sequences and outgroup name files.  Run all to run the whole notebook.

Visualize ancestor_tree in FigTree or other tree viewing software to get ancestral node numbers of interest, or can visualize in the notebook but a bit clunky.


Key output files:
ancestor_tree.txt: phylogenetic tree of input sequences with ancestral nodes labelled
ML_ancestors.fa: final ML sequences for all ancestral nodes
posterior_probabilities.json: probabilities of all 20 amino acids at every position for all ancestors


Other output files:
recoding_dict.txt: key for recoding of sequences to meet PAML 10 character seq id limit
ancestor_cladogram.txt: cladogram of input sequences with ancestral nodes labelled
ancestor_recoded_cladogram.txt: cladogram of input sequences recoded for paml
ancestor_recoded_tree.txt: phylogenetic tree of input sequences recoded for paml
codeml.ctl: config file for codeml, written by notebook from example file
ML_ancestors_with_gaps.fa: ML sequences for all ancestral nodes before removing gaps with DOWNPASS
posterior_probabilities_no_gaps.json: probabilities of all 20 amino acids at every position for all ancestors before gaps removed
rst: raw paml output file parsed by the notebook to generate all other output files


### Compute Specifications

Run locally on macOS (osx64) with 36 GM RAM and no paralellization. 

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).

---
## For Developers

This section contains information for developers who are working off of this template. Please adjust or edit this section as appropriate when you're ready to share your repo.

### GitHub templates
This template uses GitHub templates to provide checklists when making new pull requests. These templates are stored in the [.github/](./.github/) directory.

### VSCode
This template includes recommendations to VSCode users for extensions, particularly the `ruff` linter. These recommendations are stored in `.vscode/extensions.json`. When you open the repository in VSCode, you should see a prompt to install the recommended extensions.

### `.gitignore`
This template uses a `.gitignore` file to prevent certain files from being committed to the repository.

### `pyproject.toml`
`pyproject.toml` is a configuration file to specify your project's metadata and to set the behavior of other tools such as linters, type checkers etc. You can learn more [here](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/)

### Linting
This template automates linting and formatting using GitHub Actions and the `ruff` linter. When you push changes to your repository, GitHub will automatically run the linter and report any errors, blocking merges until they are resolved.
