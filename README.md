# Human Genetics GWAS Analysis

A Python data-analysis project exploring genome-wide association study
(GWAS) records across eight genes: BRCA1, BRCA2, HBB, CFTR, APOE,
LDLR, TCF7L2 and F5.

## Objectives

- Inspect dataset structure, data types and missing values.
- Summarise association records by gene.
- Explore genome-wide significance labels.
- Parse numerical beta values and effect directions.
- Reconstruct log10 p-values from mantissas and exponents.
- Check reported odds ratios against confidence intervals.
- Export derived datasets, summary tables and plots.

## Tools

- Python 3
- pandas
- NumPy
- Matplotlib
- openpyxl

## Project structure

- `analyse_gwas.py` — initial summaries and visualisation.
- `inspect_gwas.py` — dataset inspection and quality checks.
- `advanced_analysis.py` — beta parsing, p-value analysis and OR/CI checks.
- `results/` — exported CSV reports and PNG plots.

## Dataset

The working dataset contains 731 association records across eight genes.
The raw Excel workbook is excluded from version control.

Obtain the dataset from its original source and verify its licence,
attribution requirements and redistribution permissions before sharing it.

The scripts expect the workbook to be named:

`human_genetics_gwas_associations copy.xlsx`

Place the workbook in the project directory before running the scripts.

## Setup

Install the required packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install pandas numpy matplotlib openpyxl
```

Run the scripts from the project directory:

```bash
python3 inspect_gwas.py
python3 analyse_gwas.py
python3 advanced_analysis.py
```

## Current findings

- 731 records and 731 unique association IDs were observed.
- 661 beta fields were parsed into numerical values.
- Three extremely small p-values required log-scale reconstruction
  because their floating-point values underflowed to zero.
- 49 records had an odds ratio and both confidence limits available.
  All 49 estimates fell within their reported limits.
- The collected records contained different numbers of associations
  for different genes.

These findings describe this dataset snapshot, not all published
research on these genes.

## Scientific limitations

- The dataset may be filtered or incomplete and may not represent all
  published genetic associations.
- Association does not establish causation.
- Beta values may use different effect measures and units. They should
  not be pooled or compared without checking their definitions.
- Missing values and placeholders can affect summaries.
- Reconstructed p-values should be checked against the original source
  when precision is important.
- Passing an odds-ratio/confidence-interval check does not independently
  validate a study or its findings.
- Interpret results alongside original publications and documentation.

## Reproducibility

The scripts and generated reports are intended to make the analysis
inspectable and repeatable. Reproducibility depends on obtaining the
appropriate source dataset and compatible package versions.

## Author

Personal bioinformatics and scientific Python portfolio project.
