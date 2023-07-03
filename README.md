SpaceANOVA
================
Souvik Seal

This is an R package implementing the proposed method from the paper,
“SpaceANOVA: Spatial co-occurrence analysis of cell types in multiplex
imaging data using point process and functional ANOVA”. It can be used
to study differential spatial co-occurrence of pairs of cell types
across multiple groups, such as case/control, in multiplex imaging
datasets.

## Loading the package

We install and load the developmental version of SpaceANOVA from GitHub.

``` r

suppressWarnings(devtools::install_github('sealx017/SpaceANOVA'))
# 
# ── R CMD build ─────────────────────────────────────────────────────────────────
#      checking for file ‘/private/var/folders/8k/tj8qsmdj7_74drqc_slz0rkh0000gp/T/RtmptMKohV/remotes50f274ef85d7/sealx017-SpaceANOVA-d3ff624/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/8k/tj8qsmdj7_74drqc_slz0rkh0000gp/T/RtmptMKohV/remotes50f274ef85d7/sealx017-SpaceANOVA-d3ff624/DESCRIPTION’
#   ─  preparing ‘SpaceANOVA’:
#      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
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
require(SpaceANOVA)

require(spatstat)
require(tidyr)
require(dplyr)
require(fda.usc)
require(ggplot2)
require(cowplot)
require(gridExtra)
```

## Loading the example dataset

Next, we import the type 1 diabetis melitus (T1DM) IMC dataset \[1\].
The input data should have a simple form: a) rows should correspond to
individual cells, and b) six columns should correspond to the group,
subject, and image indices of the cell, along with it’s x,y co-ordinates
and cell type. In the original dataset, there were three groups of
subjects: “Non-diabetic”, “Onset”, and “Long-duration” (4 in each
group). We truncated the dataset to have only the first two groups for
simplicity. It should be highlighted that our method can handle any
number of groups (more than 2).

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

``` r
IMC_T1DM %>% group_by(Group, ID) %>% reframe(n = n()) # check group and subject IDs, and respective cell counts (spread across multiple images)
# # A tibble: 8 × 3
#   Group        ID         n
#   <fct>        <fct>  <int>
# 1 Non-diabetic 6126  158693
# 2 Non-diabetic 6134  181609
# 3 Non-diabetic 6386  127350
# 4 Non-diabetic 6278  185399
# 5 Onset        6362  202272
# 6 Onset        6414  166052
# 7 Onset        6228  142106
# 8 Onset        6380  129456
```

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
more than 1, mclapply (parallel package) is called to construct the
permutation envelope.

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
# Others 7.030910e-01 1.572046e-01 6.421546e-01 5.140401e-15 6.434097e-05
# Tc     1.589949e-01 1.143203e-05 6.626439e-01 6.437610e-10 2.164065e-01
# Th     6.596069e-01 6.778576e-01 1.819694e-01 2.989368e-03 1.031986e-07
# alpha  4.078163e-15 4.141564e-10 5.548258e-03 1.492405e-12 1.638608e-04
# delta  7.833031e-05 1.901258e-01 1.451958e-08 1.532716e-04 4.315377e-05
# beta   8.732236e-02 2.218902e-01 1.569129e-02 4.735376e-02 1.247981e-04
#                beta
# Others 8.273657e-02
# Tc     2.228232e-01
# Th     1.006248e-02
# alpha  2.653912e-02
# delta  3.232963e-06
# beta   2.983041e-01
```

## Heatmap of -log10 of p-values

Display the -log10 of the p-values using heatmap.

``` r
Plot.heatmap(Univ_p, main = "SpaceANOVA Univ.")
```

<img src="README_files/figure-gfm/P-value visualization-1.png" width="40%" />

``` r
Plot.heatmap(Mult_p, main = "SpaceANOVA Mult.")
```

<img src="README_files/figure-gfm/P-value visualization-2.png" width="40%" />

## Visualizing summary functions

Next, we can check the summary functions corresponding to those pairs of
cell types whose p-values turned out to be significant. Here, we
visualize the g-functions of pairs: (alpha, alpha) and (beta, Th). For
every pair, the first panel shows subject-level mean functions (averaged
across the images of a subject), while the second panel shows
image-level summary functions. We also display the point-wise F-values
which might help to understand which particular values of the radius r
are most influential, i.e., where the maximum difference between the
group-level summary functions are observed.

### Pair: (alpha, alpha)

We notice that the spatial co-occurrence of alpha cells is positive in
both the groups but decreases in the “Onset” group.

``` r
Pointwise_Fvals = Final_result[[1]][[2]]
Functional_results = Final_result[[2]]

# Pair: (alpha, alpha)
Plot.functions(Functional_results, pair = c("alpha", "alpha"), Fvals = Pointwise_Fvals)
# [[1]]
```

<img src="README_files/figure-gfm/plotting summary functions 1-1.png" width="90%" />

    # 
    # [[2]]

<img src="README_files/figure-gfm/plotting summary functions 1-2.png" width="90%" />

### Pair: (beta, Th)

We notice that the spatial co-occurrence of (beta, Th) is negative in
both the groups meaning that the cell types avoid each other. But the
degree of avoidance decreases in the “Onset” group.

``` r
# Pair: (beta, Tc)
Plot.functions(Functional_results, pair = c("beta", "Tc"), Fvals = Pointwise_Fvals)
# [[1]]
```

<img src="README_files/figure-gfm/plotting summary functions 2-1.png" width="90%" />

    # 
    # [[2]]

<img src="README_files/figure-gfm/plotting summary functions 2-2.png" width="90%" />

## Visualizing cellular organization

We canvisualize cellular organization in different images of every
subject. Here, 4 images each are shown for two subjects, “6126” from
group “Non-diabetic” and “6414” from group “Onset”.

``` r
table(IMC_T1DM$cellType)
# 
#   alpha    beta   delta  Others      Tc      Th 
#   73368   56901   13623 1132887   11949    4209
palette = c("darkorchid1","red", "cyan", "grey", "blue", "green") #assign colors to cell types 

Plot.cellTypes(data = IMC_T1DM, ID = "6126", palette = palette)
```

<img src="README_files/figure-gfm/plotting celluar organization-1.png" width="90%" />

    # TableGrob (2 x 2) "arrange": 4 grobs
    #   z     cells    name           grob
    # 1 1 (1-1,1-1) arrange gtable[layout]
    # 2 2 (1-1,2-2) arrange gtable[layout]
    # 3 3 (2-2,1-1) arrange gtable[layout]
    # 4 4 (2-2,2-2) arrange gtable[layout]
    Plot.cellTypes(data = IMC_T1DM, ID = "6414", palette = palette)

<img src="README_files/figure-gfm/plotting celluar organization-2.png" width="90%" />

    # TableGrob (2 x 2) "arrange": 4 grobs
    #   z     cells    name           grob
    # 1 1 (1-1,1-1) arrange gtable[layout]
    # 2 2 (1-1,2-2) arrange gtable[layout]
    # 3 3 (2-2,1-1) arrange gtable[layout]
    # 4 4 (2-2,2-2) arrange gtable[layout]

## References

1\) Damond, N., Engler, S., Zanotelli, V. R., Schapiro, D., Wasserfall,
C. H., Kusmartseva, I., Nick, H. S., Thorel, F., Herrera, P. L.,
Atkinson, M. A., et al. (2019). A map of human type 1 diabetes
progression by imaging mass cytometry. Cell metabolism, 29(3), 755–768.

2\) Baddeley, A., Rubak, E., and Turner, R. (2015). Spatial point
patterns: methodology and applications with R. CRC press.

3\) Zhang, J.-T. (2013). Analysis of variance for functional data. CRC
press.
