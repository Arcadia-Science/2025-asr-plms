
# Do protein language models understand evolution? Mixed evidence from ancestral sequences and ESM2.

This code repository contains or points to all materials required for creating and hosting the publication entitled, *"[Do protein language models understand evolution?  Mixed evidence from ancestral sequences and ESM2.]"*.

The publication is hosted at [this URL]([https://arcadia-science.github.io/2025-asr-plms/]).

## Data Description

**Ancestral Sequence Generation**

Ancestral sequences were reconstructed using the workflow in [ASR/ASR_notebook.ipynb](ASR/ASR_notebook.ipynb).  All sequence inputs and outputs are in the [ASR/](ASR/) directory.

**ESM2 Pseudo-perplexity Calculation**

ESM2 pseudo-perplexity scores were calculated using [ESM2_scoring/esm2_pppl_calculator.py](ESM2_scoring/esm2_pppl_calculator.py) on a GPU-enable AWS EC2 instance. All inputs and outputs are in the [ESM2_scoring/](ESM2_scoring/) directory.

## Reproduce

Please see [SETUP.qmd](SETUP.qmd).

## Contribute

Please see [CONTRIBUTING.qmd](CONTRIBUTING.qmd).
