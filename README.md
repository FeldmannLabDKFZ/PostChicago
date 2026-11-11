
# PostChicago

<!-- badges: start -->
<!-- badges: end -->

The goal of PostChicago is to integrate, visualize and normalize Capture-*C data with cross-sample significant interactions and other genomic intervals, such as epigenetic peaks: 

Functionally linking enhancers to gene promoters at high resolution often entails analysis of 3C datasets, including viewpoint-based methods such as Capture-(Hi)C. Quantitative comparison across samples and experiments, together with visualization and integration of promoter interactions with other datasets are of critical importance. While existing statistical tools identify significant interactions and quantify underlying signal, systematically comparing the shape and extent of this signal across large sets of intervals remains challenging. Moreover, experimental artefacts, such as inefficient crosslinking, can affect downstream interpretation, but are difficult to detect systematically. Here we introduce PostChicago, an R/Bioconductor package that fills this gap. Building upon the widely used interaction-calling package Chicago, PostChicago provides a complete workflow for quality control, visualization, integration and quantitative analysis of promoter interactions between multiple samples, including annotation and quantification of genomic regions based on external datasets. PostChicago extends visual comparisons beyond individual viewpoints by computing normalized metaprofiles around large sets of intervals which display the aggregate interaction signal, enabling comparison between narrow and broad interaction patterns. We further demonstrate how PostChicago can identify experimental crosslinking artefacts. Together, PostChicago enables integrative, shape-sensitive analysis and quality control of promoter Capture-(Hi)C datasets at a scale that would otherwise require extensive custom scripting.


## Installation

You can install the development version of PostChicago like so:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## The following initializes the development version of Bioconductor
BiocManager::install(version = "devel")

BiocManager::install("PostChicago")
```

## Documentation

The vignette and reference manual are available on the
[Bioconductor package page](https://bioconductor.org/packages/PostChicago).
After installation, open the vignette with:

```r
browseVignettes("PostChicago")
```


## Citation

If you use PostChicago, please cite: Krall S, Blum B, Smialkovska V, Ozcan E, Dachtler S, Oldfield G, Feldmann A (2026). “Capture-(Hi)C analysis with PostChicago: Integrative metaprofiling and quantitative comparison of chromatin interactions.” bioRxiv. doi:10.64898/2026.10.08.757510. Preprint.