#' @title Sets up the estimated functions as appropriate functional data object
#'
#' @param Final_result is a list output from the function "Spat.summary".
#' @param Summary_function is a character denoting which of the three summary
#' functions to use: "g" (default), "K", and "L" in the final test.
#' @param perm is TRUE (or, FALSE) if permutation-based adjustement is used (or, not).
#' @param subset is either NULL (default) or a vector radius values on which the analysis
#'  will be restricted.
#' @param print is TRUE (or, FALSE) if progress is to be displayed (or, not).
#' @param trimming_cutoff is a fraction between 0 and 1, used only in computing depth-based
#' mean of the multiple functions of a subject.
#' @return It returns a list with appropriate functional data frames to be used by
#' downstream functions.
#'
#' @export


Functional.objects <- function(Final_result = Final_result,
                               Summary_function = "g",
                               perm = TRUE, subset = NULL,
                               print = F, trimming_cutoff = 0.1)
{

  if(Summary_function == "g"){
    Sum_function = Final_result[[3]]
  }else if(Summary_function == "K"){
    Sum_function = Final_result[[1]]
  }else if(Summary_function == "L" ){
    Sum_function = Final_result[[2]]
  }
  Image_counts = Final_result[[4]]
  fixed_r = Final_result[[5]]
  IDs = Final_result[[6]]
  celltypes = Final_result[[7]]
  Hard_ths = Final_result[[8]]
  n_celltypes = length(celltypes)
  R = length(fixed_r)

  Mean_frame = Wmean_frame = Depth_mean_frame = array(0, dim = c(length(IDs),  n_celltypes,
                                                                 n_celltypes, R-1))
  if(is.null(subset) == FALSE){
    Mean_frame = Wmean_frame = Depth_mean_frame = array(0, dim = c(length(IDs),  n_celltypes,
                                                                   n_celltypes, length(subset)))
  }
  ID_groups = array(0, dim = c(length(IDs), 2,  n_celltypes,  n_celltypes))
  image_dat_frame =  array(list(), dim = c(length(IDs),  n_celltypes,  n_celltypes))

  groups = list()

  idc = 1
  for(i in IDs){
    image_info = do.call(rbind, Image_counts[[idc]])
    dat = Sum_function[[idc]]
    for(cell1 in celltypes){
      m1 = which(celltypes == cell1)
      for(cell2 in celltypes){
        m2 = which(celltypes == cell2)
        good_images =  which(image_info[ , as.character(cell1)] >= Hard_ths &
                               image_info[ , as.character(cell2)] >= Hard_ths)
        idc_spec = NULL
        if(length(good_images) > 0){
          if(perm == TRUE){
            for(x in good_images){
              idc_spec = rbind(idc_spec, dat[[x]][[1]][m1, m2, ] - dat[[x]][[2]][m1, m2, ]) #
            }
          }else{
            for(x in good_images){
              idc_spec = rbind(idc_spec, dat[[x]][[1]][m1, m2, ] - dat[[x]][[2]]) #
            }
          }


          dat_sub_wide = cbind(fixed_r[2:R], t(idc_spec))
          if(is.null(subset) == FALSE){
            dat_sub_wide = dat_sub_wide[subset, ]
          }
          weights = rowSums(image_info[, c(as.character(cell1),
                                           as.character(cell2))])[good_images]
          arg_range = as.vector(dat_sub_wide[,1])
          image_funcs = as.matrix(dat_sub_wide[, -1])
          image_funcs[is.na(image_funcs)] = 0
          colnames(image_funcs) =   image_info$imageID[good_images]


          usc_obj = fda.usc::fdata(rangeval = arg_range, mdata = t(image_funcs))
          depth_based_mean = suppressMessages(t(func.trim.mode(usc_obj, trim = trimming_cutoff)$data))

          Mean_frame[ idc, m1, m2, ] = rowMeans(image_funcs, na.rm = T)
          Wmean_frame[ idc, m1, m2, ] = WM(image_funcs, weights)
          Depth_mean_frame[ idc, m1, m2, ] = depth_based_mean

          wide_dat = data.frame(image_funcs, range = t(t(arg_range)),
                                ID = i, Group = unique(image_info$Group))

          long_dat = as.matrix(gather(wide_dat, im, func, 1:dim(image_funcs)[2]))

          image_dat_frame[idc, m1, m2] = list(long_dat)

          ID_groups[idc, , m1, m2] = c(as.character(i), as.character(unique(image_info$Group)))
        }else{
          Mean_frame[ idc, m1, m2, ] = Wmean_frame[ idc, m1, m2, ] = Depth_mean_frame[ idc, m1, m2, ] = 0
          image_dat_frame[idc, m1, m2] = list(0)
          ID_groups[idc, , m1, m2] = c( as.character(i), as.character(unique(image_info$Group)))
        }
      }
    }
    idc = idc + 1
    #print(idc)
  }
  Functional_results = list( Mean_frame, Wmean_frame, Depth_mean_frame,
                             image_dat_frame, celltypes, ID_groups)
  if(print == T) {print("Functions set up for FANOVA!")}
  return(Functional_results)
}




#' @title Performs pairwise FANOVA both Univ. and Mult.
#'
#' @param Functional_results is a list output from the function "Functional.objects".
#' @param which_mean is NULL (default) or "weighted" or "depth" depending on which type of
#' mean to be used. NULL corresponds to simple mean, while "weighted" and "depth" correspond to
#' weighted and functional depth based means, respectively.
#' @param pairs is NULL (default) or a specific pair of cell types for which the differential
#' analysis is to be performed.
#' @param print is TRUE or FALSE based on whether progression details are to be shown once the algorithm starts till completeion.
#'
#' @return It returns a list with the estimated summary functions and other input parameters to be passed on to
#' the downstream functions.
#'
#' @export


Pairwise.FANOVA.short <- function(Functional_results = Functional_results,
                                  which_mean = NULL,
                                  pairs = NULL, print = F)
{
  if(which_mean == "weighted"){
  Mean_frame = Functional_results[[2]]
  }else if(which_mean == "depth"){
  Mean_frame = Functional_results[[3]]
  }else{
  Mean_frame = Functional_results[[1]]
  }
  Full_dat_frame = Functional_results[[4]]
  celltypes = Functional_results[[5]]
  ID_groups = Functional_results[[6]]

  if(is.null(pairs) == T){
    res = array(NA, dim = c(3, length(celltypes), length(celltypes)))
    Fvals = array(list(), dim = c(2, length(celltypes), length(celltypes)))
    for(cell1 in celltypes){
      m1 = which(celltypes == cell1)
      for(cell2 in celltypes){
        m2 = which(celltypes == cell2)
        mean_dat = as.data.frame(Mean_frame[ , m1, m2, ])
        mean_dat = mean_dat[, colSums(abs(mean_dat)) > 0]
        good_IDs = ID_groups[, 1, m1, m2][which(rowSums(abs(mean_dat)) > 0)]
        good_IDs_group = ID_groups[, 2, m1, m2][which(rowSums(abs(mean_dat)) > 0)]
        good_IDs_group = factor(good_IDs_group, levels = unique(good_IDs_group))

        if(any(table(good_IDs_group) == 1) | length(unique(good_IDs_group)) < 2){
          res[ , m1, m2] = rep(NA, 3)
        }else{
          data = cbind(good_IDs, good_IDs_group,
                       mean_dat[which(rowSums(abs(mean_dat)) > 0), ])
          colnames(data)[1:2] = c("ID", "Group")
          FuncANOVA1 = UniFANOVA(data)

          idc = 1
          All_func_dat = NULL
          for(i in ID_groups[, 1, m1, m2])
          {
            func_dat = Full_dat_frame[idc, m1, m2][[1]]
            if(length(func_dat) != 1){
              good_rows = apply(func_dat, 1, function(x) is.character(x))
              if(length(good_rows)>0){
                func_dat = func_dat[good_rows, ]
                func_dat = data.frame(ID = func_dat[,2], Group = func_dat[,3],
                                      imageID = func_dat[,4],
                                      range = as.numeric(func_dat[,1]),
                                      func = as.numeric(func_dat[,5]))

                All_func_dat = rbind(All_func_dat, func_dat)

              }

            }
            idc = idc + 1
          }

          All_func_dat$ID = factor(All_func_dat$ID, levels = unique(All_func_dat$ID))
          All_func_dat$imageID = factor(All_func_dat$imageID,
                                        levels = unique(All_func_dat$imageID))
          All_func_dat$Group = factor(All_func_dat$Group, levels = unique(All_func_dat$Group))
          MultiFANOVA = MFanova.tests(data = All_func_dat[All_func_dat$range > 0,])

          res[ , m1, m2] = c(Mean_FB = FuncANOVA1$pval,
                             MFANOVA_GPF_int = MultiFANOVA$GPF_int$pval,
                             MFANOVA_GPF =  MultiFANOVA$GPF$pval)
          Fvals[, m1, m2] = list(FuncANOVA1$Fvals, MultiFANOVA$Fvals)
        }}}
  }else{
    m1 = which(celltypes == pairs[1])
    m2 = which(celltypes == pairs[2])
    mean_dat = as.data.frame(Mean_frame[ , m1, m2, ])
    mean_dat = mean_dat[, colSums(abs(mean_dat)) > 0]
    good_IDs = ID_groups[, 1, m1, m2][which(rowSums(abs(mean_dat)) > 0)]
    good_IDs_group = ID_groups[, 2, m1, m2][which(rowSums(abs(mean_dat)) > 0)]
    good_IDs_group = factor(good_IDs_group, levels = unique(good_IDs_group))

    if(any(table(good_IDs_group) == 1) | length(unique(good_IDs_group)) < 2){
      res = rep(NA, 3)
    }else{
      data = cbind(good_IDs, good_IDs_group,
                   mean_dat[which(rowSums(abs(mean_dat)) > 0), ])
      colnames(data)[1:2] = c("ID", "Group")
      FuncANOVA1 = UniFANOVA(data)

      idc = 1
      All_func_dat = NULL
      for(i in ID_groups[, 1, m1, m2])
      {
        func_dat = Full_dat_frame[idc, m1, m2][[1]]
        if(length(func_dat) != 1){
          good_rows = apply(func_dat, 1, function(x) is.character(x))
          if(length(good_rows)>0){
            func_dat = func_dat[good_rows, ]
            func_dat = data.frame(ID = func_dat[,2], Group = func_dat[,3],
                                  imageID = func_dat[,4],
                                  range = as.numeric(func_dat[,1]),
                                  func = as.numeric(func_dat[,5]))
            All_func_dat = rbind(All_func_dat, func_dat)

          }

        }
        idc = idc + 1
      }

      All_func_dat$ID = factor(All_func_dat$ID, levels = unique(All_func_dat$ID))
      All_func_dat$imageID = factor(All_func_dat$imageID,
                                    levels = unique(All_func_dat$imageID))
      All_func_dat$Group = factor(All_func_dat$Group, levels = unique(All_func_dat$Group))



      MultiFANOVA = MFanova.tests(data = All_func_dat[All_func_dat$range > 0,])

      res = c(Mean_GPF = FuncANOVA1$pval,
              MFANOVA_GPF_int = MultiFANOVA$GPF_int$pval,
              MFANOVA_GPF =  MultiFANOVA$GPF$pval)
      Fvals = list(FuncANOVA1$Fvals, MultiFANOVA$Fvals)
    }}
  if(print == T) {print("FANOVA complete!")}
  return(list(res, Fvals, celltypes))
}




#' @title Performs the entire analysis from computing summary functions to FANOVA in a
#' streamlined way
#'
#' @param data is the data matrix with columns named, "Group", "ID",  "imageID","cellType", "x", "y",
#' and rows representing the cells.
#' @param fixed_r is the grid of radius r on which the summary functions will be estimated. Should be
#' increasing in order and start with 0.
#' @param Summary_function is a character denoting which of the three summary
#' functions to use: "g" (default), "K", and "L" in the final test.
#' @param ID_subset is a vector of subject IDs on which the analysis will be restricted. If NULL (default), all
#' of the subjects are used.
#' @param celltypes is a vector of cell types on which the analysis will be restricted. If NULL (defaukt),
#' all unique cell types (co-occurrence of pairs of all unique cell types) are considered.
#' @param Hard_ths is a constant denoting the lowest number of a particular cell type that an image can
#' have. For example, Hard_ths = 10 means that if an image has less than 10
#' cells of a particular cell type A, in all pairwise comparisons involving A, such as (A, A), (A, B), (A, C),.., the image
#' would be dropped from the analysis (B, C,.. are other cell types).
#' @param perm is TRUE or FALSE denoting if permutation based adjustment will be performed or not, respectively.
#' @param nPerm is an integer denoting the umber of permutations to be used. Only used if perm = TRUE o
#' @param print is TRUE or FALSE based on whether progression details are to be shown once the algorithm starts till completeion.
#' @param parallel is an integer denoting the number of cores to be used.
#'
#' @return It returns a list with the estimated summary functions and other input parameters to be passed on to
#' the downstream functions.
#'
#' @export

All_in_one <- function(data, fixed_r = seq(0, 100, by = 1), Summary_function = "g",
                       ID_subset = NULL, celltypes = NULL, Hard_ths = 10,
                       perm = TRUE, nPerm = 50,  print = F, cores = 8){
  Final_result = Spat.summary(data = data, fixed_r, ID_subset,
                              celltypes, Hard_ths, perm = perm, nPerm = nPerm,
                              print = print, cores = cores)
  Functional_results = Functional.objects(Final_result = Final_result,
                                          Summary_function = Summary_function,
                                          perm = perm, print = print)
  res = Pairwise.FANOVA.short(Functional_results = Functional_results, print = print)

  return(res)
}


WM<-function(Func_mat, weights){
  Func_mat[is.na(Func_mat)] = 0;  weights[is.na(weights)] = 0
  return(Func_mat %*% weights/sum(weights))
}


