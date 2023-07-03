SpaceANOVA
================
Souvik Seal

This is an R package implementing the proposed method from the paper,
“SpaceANOVA: Spatial co-occurrence analysis of cell types in multiplex
imaging data using point process and functional ANOVA”. It can be used
to study differential spatial co-occurrence of pairs of cell types
across multiple groups, such as case/control.

## Loading the package

We install and load the developmental version of SpaceANOVA from GitHub.

``` r

suppressMessages(devtools::install_github('sealx017/SpaceANOVA'))
# Warning in readLines(old_path): incomplete final line found on
# '/Users/sealso/.R/Makevars'
# ── R CMD build ─────────────────────────────────────────────────────────────────
#      checking for file ‘/private/var/folders/8k/tj8qsmdj7_74drqc_slz0rkh0000gp/T/Rtmp18oVYA/remotes13519366d9af3/sealx017-SpaceANOVA-1f7d397/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/8k/tj8qsmdj7_74drqc_slz0rkh0000gp/T/Rtmp18oVYA/remotes13519366d9af3/sealx017-SpaceANOVA-1f7d397/DESCRIPTION’
#   ─  preparing ‘SpaceANOVA’:
#    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#   ─  checking for LF line-endings in source and make files and shell scripts
#   ─  checking for empty or unneeded directories
#        NB: this package now depends on R (>= 3.5.0)
#        WARNING: Added dependency on R >= 3.5.0 because serialized objects in
#      serialize/load version 3 cannot be read in older versions of R.
#      File(s) containing such objects:
#        ‘SpaceANOVA/Data/IMC_T1DM.rda’
#   ─  building ‘SpaceANOVA_0.1.0.tar.gz’
#      
# 
require(spatstat)
require(tidyr)
require(dplyr)
require(fda.usc)
require(ggplot2)
require(cowplot)
require(SpaceANOVA)
```

## Loading the example dataset

Next, we import the type 1 diabetis melitus (T1DM) IMC dataset \[1\].
The input data should have a simple form: a) rows should correspond to
individual cells, and b) six columns should correspond to the group,
subject, and image indices of the cell, along with it’s x,y co-ordinates
and cell type.

``` r

data("IMC_T1DM")
knitr::kable(head(IMC_T1DM), format="markdown")
```

| Group | ID   | imageID |   x |   y | cellType |
|:------|:-----|:--------|----:|----:|:---------|
| Onset | 6362 | A01     |  53 |   0 | Others   |
| Onset | 6362 | A01     | 128 |   0 | Others   |
| Onset | 6362 | A01     | 135 |   0 | Others   |
| Onset | 6362 | A01     | 450 |   0 | Others   |
| Onset | 6362 | A01     | 458 |   0 | Others   |
| Onset | 6362 | A01     | 551 |   0 | Others   |

## Compute the summary functions and perform association analysis

This function can be used to perform the entire analysis at one go,
starting from summary function estimation to the univariate and
multivariate FANOVA-based association analyses. The function takes
several input parameters out of which the most important ones are
displayed below. “fixed_r” corresponds to the grid values of radius r.
“Summary_function” corresponds to which summary function to be used: g,
K, or L. “Hard_ths” corresponds to the lowest number of a particular
cell type, say A, that an image can have (the image will be dropped if
not for any pair of cell types involving A). “perm” is used to employ
the permutation envelope-based adjustment to the summary functions in
presence of holes in images and “nPerm” denotes the number of
permutations to be used. “cores” is the number of cores to be used, if
more than 1, mclapply (parallel package) is called.

``` r

Final_result = All_in_one(data = IMC_T1DM, fixed_r = seq(0, 100, by = 5), Summary_function = "g", Hard_ths = 10, perm = TRUE, nPerm = 20, cores = 8)
# [1] 2
# [1] 3
# [1] 4
# [1] 5
# [1] 6
# [1] 7
# [1] 8
# [1] 9
```

## Extracting p-values

Extract the p-values of the univariate and multivariate tests:
SpaceANOVA Univ. and SpaceANOVA Mult. and print.

``` r
p_res = p_extract(Final_result)
Univ_p = p_res[[1]]
Mult_p = p_res[[2]]
print(Univ_p)
#              Others           Tc           Th        alpha        delta
# Others 7.006061e-01 1.562276e-01 6.379778e-01 6.553735e-15 5.859720e-05
# Tc     1.579690e-01 1.211840e-05 6.589536e-01 7.197470e-10 2.211829e-01
# Th     6.605541e-01 6.698687e-01 2.171020e-01 2.750797e-03 6.357745e-05
# alpha  5.399850e-15 4.792938e-10 4.860598e-03 1.848498e-12 1.636848e-04
# delta  7.142692e-05 1.990115e-01 2.744462e-05 1.535029e-04 4.588658e-05
# beta   7.643171e-02 1.759561e-01 3.301321e-02 5.786298e-02 2.033933e-04
#                beta
# Others 7.182447e-02
# Tc     1.784679e-01
# Th     2.087199e-02
# alpha  3.074615e-02
# delta  7.044844e-06
# beta   2.849975e-01
```

## Heatmap of -log10 of p-values

Display the -log10 of the p-values using heatmap.

``` r

Plot.heatmap(Univ_p, main = "SpaceANOVA Univ.")
```

<img src="README_files/figure-gfm/P-value visualization-1.png" width="100%" />

``` r
Plot.heatmap(Mult_p, main = "SpaceANOVA Mult.")
```

<img src="README_files/figure-gfm/P-value visualization-2.png" width="100%" />

## Visualizing summary functions

Next, we can check the summary functions corresponding to those pairs of
cell types whose p-values turned out to be significant. In this case, we
visualize the g-functions of pairs: (alpha, alpha) and (beta, Th). We
also display the point-wsie F-values which might help to understand
which values of the radius r were most influential.

``` r

Pointwise_Fvals = Final_result[[1]][[2]]
Functional_results = Final_result[[2]]

Plot.functions(Functional_results, pair = c("alpha", "alpha"), Fvals = Pointwise_Fvals)
# Warning: The `fun.y` argument of `stat_summary()` is deprecated as of ggplot2 3.3.0.
# ℹ Please use the `fun` argument instead.
# ℹ The deprecated feature was likely used in the SpaceANOVA package.
#   Please report the issue to the authors.
# This warning is displayed once every 8 hours.
# Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
# generated.
# [[1]]
```

<img src="README_files/figure-gfm/plotting summary functions-1.png" width="100%" />

    # 
    # [[2]]

<img src="README_files/figure-gfm/plotting summary functions-2.png" width="100%" />

## Visualizing cellular organization

We can also visualize cellular organization in different images of every
subject.

## References

1\) Damond, N., Engler, S., Zanotelli, V. R., Schapiro, D., Wasserfall,
C. H., Kusmartseva, I., Nick, H. S., Thorel, F., Herrera, P. L.,
Atkinson, M. A., et al. (2019). A map of human type 1 diabetes
progression by imaging mass cytometry. Cell metabolism, 29(3), 755–768.

2\) Baddeley, A., Rubak, E., and Turner, R. (2015). Spatial point
patterns: methodology and applications with R. CRC press.

3\) Zhang, J.-T. (2013). Analysis of variance for functional data. CRC
press.
