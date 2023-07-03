SpaceANOVA
================
Souvik Seal


This is an R package implementing the proposed method from the paper, "SpaceANOVA: Spatial co-occurrence analysis of cell types in multiplex imaging data using point
process and functional ANOVA". It can be used to study differential spatial co-occurrence of pairs of cell types across multiple groups, such as case/control, in multiplex imaging datasets.

## Loading the package

We install and load the developmental version of SpaceANOVA from GitHub.

```r

suppressWarnings(devtools::install_github('sealx017/SpaceANOVA'))
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

Next, we import the type 1 diabetis melitus (T1DM) IMC dataset [1]. The input data should have a simple form: a) rows should correspond to individual cells, and b) six columns should correspond to the group, subject, and image indices of the cell, along with it's x,y co-ordinates and cell type. In the original dataset, there were three groups of subjects: "Non-diabetic", "Onset", and "Long-duration" (4 in each group). We truncated the dataset to have only the first two groups for simplicity. It should be highlighted that our method can handle any number of groups (more than 2). 

```r 
data("IMC_T1DM")
knitr::kable(head(IMC_T1DM), format="markdown")
IMC_T1DM %>% group_by(Group, ID) %>% reframe(n = n()) # check group and subject IDs, and respective cell counts (spread across multiple images)

```

## Compute the summary functions and perform association analysis

This function can be used to perform the entire analysis at one go, starting from summary function estimation to the univariate and multivariate FANOVA-based association analyses. 
The function takes several input parameters out of which the most important ones are displayed below. "fixed_r" corresponds to the grid values of radius r. "Summary_function" corresponds to which summary function to be used: g, K, or L.  "Hard_ths" corresponds to the lowest number of a particular cell type, say A, that an image can have (the image will be dropped if not for any pair of cell types involving A). "perm" is used to employ the permutation envelope-based adjustment to the summary functions in presence of holes in images and "nPerm" denotes the number of permutations to be used. "cores" is the number of cores to be used, if more than 1, mclapply (parallel package) is called to construct the permutation envelope.


```r 

Final_result = All_in_one(data = IMC_T1DM, fixed_r = seq(0, 100, by = 5), Summary_function = "g", Hard_ths = 10, perm = TRUE, nPerm = 20, cores = 8)

```

## Extracting p-values 

Extract the p-values of the univariate and multivariate tests: SpaceANOVA Univ. and SpaceANOVA Mult. and print.

```r
p_res = p_extract(Final_result)
Univ_p = p_res[[1]]
Mult_p = p_res[[2]]
print(Univ_p)

```

## Heatmap of -log10 of p-values 


Display the -log10 of the p-values using heatmap.

```r
Plot.heatmap(Univ_p, main = "SpaceANOVA Univ.")
Plot.heatmap(Mult_p, main = "SpaceANOVA Mult.")

```

## Visualizing summary functions

Next, we can check the summary functions corresponding to those pairs of cell types whose p-values turned out to be significant. Here, we visualize the g-functions of pairs: (alpha, alpha) and (beta, Th). For every pair, the first panel shows subject-level mean functions (averaged across the images of a subject), while the second panel shows image-level summary functions. We also display the point-wise F-values which might help to understand which particular values of the radius r are most influential, i.e., where the maximum difference between the group-level summary functions are observed. 

### Pair: (alpha, alpha)

We notice that the spatial co-occurrence of alpha cells is positive in both the groups but
decreases in the "Onset" group.

```r
Pointwise_Fvals = Final_result[[1]][[2]]
Functional_results = Final_result[[2]]

# Pair: (alpha, alpha)
Plot.functions(Functional_results, pair = c("alpha", "alpha"), Fvals = Pointwise_Fvals)

```

### Pair: (beta, Th)

We notice that the spatial co-occurrence of (beta, Th) is negative in both the groups meaning that the cell types avoid each other. But the degree of avoidance decreases in the "Onset" group.

```r
# Pair: (beta, Tc)
Plot.functions(Functional_results, pair = c("beta", "Tc"), Fvals = Pointwise_Fvals)

```

## Visualizing cellular organization

We canvisualize cellular organization in different images of every subject. Here, 4 images each are shown for two subjects, "6126" from group "Non-diabetic" and "6414" from group "Onset". 

```r
table(IMC_T1DM$cellType)
palette = c("darkorchid1","red", "cyan", "grey", "blue", "green") #assign colors to cell types 

Plot.cellTypes(data = IMC_T1DM, ID = "6126", palette = palette)
Plot.cellTypes(data = IMC_T1DM, ID = "6414", palette = palette)

```


## References

1\) Damond, N., Engler, S., Zanotelli, V. R., Schapiro, D., Wasserfall, C. H.,
Kusmartseva, I., Nick, H. S., Thorel, F., Herrera, P. L., Atkinson, M. A., et al.
(2019). A map of human type 1 diabetes progression by imaging mass cytometry.
Cell metabolism, 29(3), 755–768.

2\) Baddeley, A., Rubak, E., and Turner, R. (2015). Spatial point patterns: methodology and applications with R. CRC press.

3\) Zhang, J.-T. (2013). Analysis of variance for functional data. CRC press.



