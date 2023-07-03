#' @title Generates heatmaps of the -log10 of p-values
#'
#' @param p_matrix is a dtaframe with p-values corresponding to every pair of cell types.
#'
#' @return a heatmap object from the pheatmap package.
#'
#' @export
#'

Plot.heatmap <- function(p_matrix, main = "SpaceANOVA Univ."){
  breaks <- c(0, 10, 0.5)
  breaks <- seq(from = breaks[1], to = breaks[2], by = breaks[3])
  colours = c("#4575B4", "white", "#D73027")
  pal <- grDevices::colorRampPalette(colours)(length(breaks))

  return(pheatmap::pheatmap(-log10(p_matrix), cluster_rows = F, cluster_cols = F,
                     color=pal, main=main, fontsize = 8))
}

#' @title Generates plots of summary functions for a pair of cell types
#'
#' @param Functional_results is a list of functional data generated using All_in_one function.
#' @param pair is the pair of cell types to be visualized.
#' @param Fvals is NULL or a list of point-wise F-values corresponding to SpaceANOVA Univ. and
#' SpaceANOVA Mult. obtained from All_in_one function.
#' @param range_lim is NULL or a vector of radius values on which the functions would be displayed.
#' Must be in increasing order and a subset of fixed_r.
#' @return a ggplot object with plots of mean and individual summary functions across groups.
#'
#' @export
#'

Plot.functions<- function(Functional_results = Functional_results,
                        pair = pair, Fvals = NULL, range_lim = NULL)
{
  Mean_frame = Functional_results[[1]]
  Wmean_frame = Functional_results[[2]]
  Depth_mean_frame = Functional_results[[3]]
  Full_dat_frame = Functional_results[[4]]
  celltypes = Functional_results[[5]]
  ID_groups = Functional_results[[6]]


  m1 = which(celltypes == pair[1])
  m2 = which(celltypes == pair[2])
  mean_dat = as.data.frame(Mean_frame[ , m1, m2, ])
  good_IDs = ID_groups[, 1, m1, m2][which(rowSums(abs(mean_dat)) > 0)]
  good_IDs_group = ID_groups[, 2, m1, m2][which(rowSums(abs(mean_dat)) > 0)]
  good_IDs_group = factor(good_IDs_group, levels = unique(good_IDs_group))
  if(is.null(Fvals) == FALSE){
    Fvals_dat = Fvals[ , m1, m2]
  }

  mean_dat = mean_dat[which(rowSums(abs(mean_dat)) > 0), ]
  mean_dat_ID = cbind.data.frame(mean_dat, as.factor(good_IDs), as.factor(good_IDs_group))


  All_means_long = gather(mean_dat_ID, range, Mean, (1:dim(mean_dat)[2]))
  colnames(All_means_long)[1:2] = c("ID", "Group")
  All_means_long = All_means_long[order(All_means_long$ID), ]
  All_means_long$range = factor(as.numeric(gsub("V", "", All_means_long$range)))


  idc = 1
  All_func_dat = NULL
  for(i in ID_groups[, 1, m1, m2])
  {
    func_dat = Full_dat_frame[idc, m1, m2][[1]]
    if(length(func_dat) != 1){
      good_rows = apply(func_dat, 1, function(x) is.character(x))
      if(length(good_rows)>0){
        func_dat = func_dat[good_rows, ]
        All_func_dat = rbind(All_func_dat, data.frame(ID = func_dat[,2], Group = func_dat[,3],
                                                      imageID = func_dat[,4],
                                                      range = as.numeric(func_dat[,1]), func = as.numeric(func_dat[,5])))
      }
    }
    idc = idc + 1
  }

  All_func_dat$ID = factor(All_func_dat$ID, levels = unique(All_func_dat$ID))
  All_func_dat$imageID = factor(All_func_dat$imageID, levels = unique(All_func_dat$imageID))
  All_func_dat$Group = factor(All_func_dat$Group, levels = unique(All_func_dat$Group))

  levels(All_means_long$range) = unique(All_func_dat$range)
  All_means_long$range =   as.numeric(levels(All_means_long$range))[All_means_long$range]
  All_means_long$Group = factor(All_means_long$Group, levels = c(levels(All_func_dat$Group)))

  if(is.null(range_lim) == FALSE){
    All_means_long = All_means_long[All_means_long$range < range_lim, ]
    All_func_dat = All_func_dat[All_func_dat$range < range_lim, ]
  }

  simple_mean = ggplot(All_means_long, aes(x = range, y = Mean, group = ID))+
    geom_line(linewidth = 0.5) + facet_wrap(~Group, labeller = label_both) +
    ggtitle(paste0('Mean of centered summary functions of pair:  \n (', paste0(pair, collapse = ", ") , ") ",
                   "across subjects")) +
    stat_summary(aes(y = Mean, group=1), fun = mean, colour = "red",
                 geom = "line", group=1, size = 1)+
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          legend.title = element_text(size=10), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "bottom",
          panel.background = element_blank()) + labs(x = "r", y = "g(r)") +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    geom_hline(yintercept = 0, linetype="dotted", color = "blue", linewidth = 0.5)

  #print(simple_mean)


  mult_mean = ggplot(All_func_dat, aes(x = range, y = func, group = imageID))+
    geom_line(linewidth=0.5) + facet_wrap(~Group, labeller = label_both) +
    stat_summary(aes(y = func, group=1),
                 fun = mean, colour = "red", size = 1, geom = "line", group=1) +
    ggtitle(paste0('Centered summary functions of pair: \n  (', paste0(pair, collapse = ", ") , ")",
                   " across images")) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          legend.title = element_text(size=10), #change legend title font size
          legend.text = element_text(size=10),
          legend.position = "bottom",
          panel.background = element_blank())+ labs(x = "r", y = "g(r)") +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    geom_hline(yintercept = 0, linetype="dotted", color = "blue", linewidth = 0.5)

  if(is.null(Fvals) == TRUE) {return(list(simple_mean, mult_mean))
  }else{
    unidat = data.frame(Fval = c(0, unlist(Fvals_dat[[1]])[1:(length(unique(All_func_dat$range))-1)]),
                        r = unique(All_func_dat$range))
    UniF = ggplot(unidat, aes(x = r, y = Fval))+
      geom_point(size=0.5) +
      ggtitle(paste0('Accompanying point-wise F values')) +
      theme(plot.title = element_text(hjust = 0.5, size = 10),
            legend.title = element_text(size=10), #change legend title font size
            legend.text = element_text(size=10),
            legend.position = "bottom",
            panel.background = element_blank())+ labs(x = "r", y = "F(r)") +
      guides(colour = guide_legend(override.aes = list(size=5)))

    Muldat = data.frame(Fval = c(0, unlist(Fvals_dat[[2]])[1:(length(unique(All_func_dat$range))-1)]),
                        r = unique(All_func_dat$range))
    MulF = ggplot(Muldat, aes(x = r, y = Fval))+
      geom_point(size=0.5) +
      ggtitle(paste0('Accompanying point-wise F values')) +
      theme(plot.title = element_text(hjust = 0.5, size = 10),
            legend.title = element_text(size=10), #change legend title font size
            legend.text = element_text(size=10),
            legend.position = "bottom",
            panel.background = element_blank())+ labs(x = "r", y = "F(r)")

    simple_mean = cowplot::plot_grid(simple_mean, UniF, align = "h",  axis = "b",
                                     nrow = 1, rel_widths = c(2/3, 1/3))

    mult_mean = cowplot::plot_grid(mult_mean, MulF, align = "h",  axis = "b",
                                   nrow = 1, rel_widths = c(2/3, 1/3))


    return(list(simple_mean, mult_mean))
  }
}


#' @title Generates plots of cellular organization of images of a particular subject
#'
#' @param data is the original data frame with cell level data.
#' @param ID is a character denoting the ID of the subject to consider.
#' @return a ggplot object of cellular organization (cells at their
#' locations with colors according to their types) in different images of a subject.
#'
#' @export
#'



Plot.cellTypes <- function(data = data, ID = ID, palette = NULL)
{
  one_subject = data[data$ID == ID, ]
  subject_images = unique(one_subject$imageID)
  check_proportions_full = table(one_subject$cellType)
  good_image_counter = matrix(0, nrow = length(subject_images), ncol = 1)

  image_info = suppressMessages(one_subject %>%
                                  group_by(imageID, cellType) %>% summarise(n = n()))

  image_info_wide = reshape(as.data.frame(image_info), idvar = "imageID",
                            timevar = "cellType", direction = "wide")

  r = 1
  for(image_id in image_info_wide[1:3, 1]){
    one_subject_one_image = one_subject[one_subject$imageID == image_id, ]
    if(is.null(palette) == FALSE){assign(paste0("simple_mean_", r), one_subject_one_image  %>%
             ggplot( aes(x = x, y = y, color = cellType)) + geom_point(size = 1) +
             theme(legend.position = "bottom",
             plot.title = element_text(hjust = 0.5)) +
             scale_color_manual(values = palette) +
             guides(colour = guide_legend(override.aes = list(size = 2))) +
             ggtitle(paste0("Image ",  image_id, " from subject ", ID))
           )}
    else{assign(paste0("simple_mean_", r), one_subject_one_image  %>%
             ggplot( aes(x = x, y = y, color = cellType)) + geom_point(size = 1) +
             theme(legend.position = "bottom",
                   plot.title = element_text(hjust = 0.5)) +
             guides(colour = guide_legend(override.aes = list(size = 2))) +
             ggtitle(paste0("Image ",  image_id, " from subject ", ID)))}

    r = r+1
  }

  three_plots = grid.arrange(grobs = list(simple_mean_1, simple_mean_2, simple_mean_3),
                             ncol = 3, main = paste0("Celltypes: ", paste0(pairs, collapse = ", ")))
  return(three_plots)
}


