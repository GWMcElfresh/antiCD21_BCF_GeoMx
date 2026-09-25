# BCF anti-CD21 transcriptional modules

This repository contains the manuscript companion script for the B-cell follicle
figures in Chan CN, Busman-Sahay K, Nekorchuk M, et al.,
"Anti-CD21/CR2 antibody treatment disrupts SIVmac239 deposition on follicular
dendritic cells during antiretroviral therapy."

Authors of the companion script: GW McElfresh, Daniel Newhouse, Leanne Whitmore.

The plotted values are row z-scores of voom log2 expression from the BCF
limma–voom analysis, and the mean of those z-scores for each module.
Figure 3B plots log2 fold change with anti-CD21 in the positive direction.
Points are colored when Benjamini–Hochberg adjusted *P* < 0.05 and the
absolute log2 fold change exceeds 0.58.

## Reproduce the analysis

The lockfile records R 4.4.1. Restore the locked renv library from the repository root:

```sh
Rscript -e 'renv::restore()'
```

Run the script:

```sh
Rscript bcf_cd21_module_figures.R
```

The script writes these tracked manuscript artifacts under `figures/`:

- `figure3A_bcf_heatmap.png` and `.pdf` — BCF heatmap
- `figure3B_bcf_volcano.png` and `.pdf` — BCF volcano, anti-CD21 positive
- `figure3D_immediate_early.png` and `.pdf` — immediate-early heatmap and score
- `figure3E_tfh.png` and `.pdf` — Tfh heatmap and score
- `figureS11B_nfkb.png` and `.pdf` — Supplemental Figure 11B, NF-kB
- `figureS11C_ifn_isg.png` and `.pdf` — Supplemental Figure 11C, IFN/ISG
- `data/Source_data.xlsx`

## Source data

`data/Source_data.xlsx` is the source-data file. Heatmap sheets are the row
z-scores of voom log2 expression used to color each heatmap (genes by BCF ROI).
Score sheets are the mean of those z-scores, one value per ROI, which is the
quantity drawn in the score plots. The Figure 3B sheet is every gene in the
BCF contrast: log2 fold change, −log10(Benjamini–Hochberg adjusted *P*), and
the color class drawn on the plot. Positive log2 fold change is higher in
anti-CD21.

The script reads these inputs from `data/`:

- `voom_log2_bcf.tsv` — voom log2 expression for the plotted genes and BCF ROIs
- `sample_metadata.tsv` — BCF ROI, animal, and treatment
- `module_genes.tsv` — gene symbols in each module
- `bcf_limma_topTable.tsv` — BCF anti-CD21 versus no-treatment statistics
- `module_statistics.tsv` — limma t, raw *P*, and Holm-adjusted *P* drawn on the score plots

Figure 3A contains genes with Benjamini–Hochberg adjusted *P* < 0.05 and
absolute log2 fold change of at least 0.58. Holm adjustment on the score plots
is across the immediate-early, Tfh, NF-kB, and IFN/ISG modules.
