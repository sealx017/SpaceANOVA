
#' @title Compute the spatial summary functions: K, L, and g, for each image of every subject in all groups and
#' adjust the functions using permutation-based envelope (optional).
#'
#' @param data is the data matrix with columns named, "Group", "ID",  "imageID","cellType", "x", "y",
#' and rows representing the cells.
#' @param fixed_r is the grid of radius r on which the summary functions will be estimated. Should be
#' increasing in order and start with 0.
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
#'
#' @return It returns a list with the estimated summary functions and other input parameters to be passed on to
#' the downstream functions.
#'
#' @export

Spat.summary <- function(data = data, fixed_r = seq(0, 100, by = 1), ID_subset = NULL,
                         celltypes = NULL, Hard_ths = 10, perm = TRUE, nPerm = 19, print = F)
{
  L = K = g = Image_counts = list()
  if(is.null(ID_subset) == T){IDs = unique(data$ID)
  }else{IDs = ID_subset}

  if(is.null(celltypes) == T){celltypes = unique(data$cellType)
  }else{celltypes = celltypes}

  n_celltypes = length(celltypes)
  n_celltype_pairs = n_celltypes^2
  R = length(fixed_r)
  idc = 1
  for(i in IDs){
    one_subject = data[data$ID == i,]
    subject_images = unique(one_subject$imageID)
    check_proportions_full = table(one_subject$cellType)
    good_image_counter = matrix(0, nrow = length(subject_images), ncol = 1)
    rowsum_tracker = matrix(0, nrow = length(subject_images), ncol = 4)
    K_all = L_all = g_all = Image_all = list()

    r = 1
    for(image_id in subject_images){
      one_subject_one_image = one_subject[one_subject$imageID == image_id, ]
      check_proportions = table(one_subject_one_image$cellType)
      cp = matrix(0, 1, length(celltypes)); colnames(cp) = celltypes
      cp[, names(check_proportions)] = check_proportions

      if(all(check_proportions < Hard_ths)){
        K_all[[r]] = L_all[[r]] =  g_all[[r]] = NULL
      }else{
        good_phenotypes = names(which(check_proportions >= Hard_ths))
        one_subject_one_image = one_subject_one_image[one_subject_one_image$cellType %in% good_phenotypes , ]
        x = one_subject_one_image$x; xrange = range(x)
        y = one_subject_one_image$y; yrange = range(y)
        PP_obj = ppp(x, y, xrange, yrange)
        m = marks(PP_obj) = factor(one_subject_one_image$cellType, good_phenotypes)

        Kall <- alltypes(PP_obj, Kcross, r = fixed_r, correction = "isotropic")
        gall <- pcf(Kall, spar=1, method="c", divisor ="d")
        K_theo <- Kall$fns[[1]]$theo[2:R]; L_theo <- fixed_r[2:R]; g_theo = 1
        K_in = L_in = g_in = array(0, dim = c(n_celltypes, n_celltypes, R-1))
        subset = which(celltypes %in% good_phenotypes)

        if(perm == "TRUE"){
          Perm = Perm_spat(PP_obj, n_celltypes, subset, fixed_r, R, nPerm = nPerm)
          K_theo = Perm[[1]]; L_theo = Perm[[2]]; g_theo = Perm[[3]]
        }

        s = 1
        for(r1 in subset){
          for(r2 in subset){
            K_in[ r1, r2, ] = Kall$fns[[s]]$iso[2:R]
            L_in[ r1, r2, ] =  sqrt(K_in[ r1, r2, ]/pi)
            gfunc <- gall$fns[[s]]$pcf[2:R]
            gfunc[is.na(gfunc)] <- 1
            g_in[ r1, r2, ] = gfunc
            s = s+1
          }
        }

        K_all[[r]] = list(K_in, K_theo)
        L_all[[r]] = list(L_in, L_theo)
        g_all[[r]] = list(g_in, g_theo)
      }
      Image_all[[r]] = cbind(data.frame(unique(one_subject_one_image$Group),
                                        unique(one_subject_one_image$ID),
                                        unique(one_subject_one_image$imageID)), cp)
      colnames(Image_all[[r]])[1:3] = c("Group", "ID", "imageID")
      r = r+1
    }

    K[[idc]] = K_all; L[[idc]] = L_all; g[[idc]] = g_all
    Image_counts[[idc]] = Image_all

    idc = idc + 1
    print(idc)
  }

  Final_result = list(K, L, g, Image_counts, fixed_r, IDs, celltypes, Hard_ths)
  if(print == T) {print("Summary functions computed!")}
  return(Final_result)
}
