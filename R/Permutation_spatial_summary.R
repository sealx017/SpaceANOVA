#' @title Computes permutation-based mean of spatial summary functions: K, L, and g, for a single image
#'
#' @param PP_obj is a spatstat object imported from the Spatial_summary function.
#' @param n_celltypes is the total number of cell types.
#' @param subset is the subset of cell types considered.
#' @param fixed_r is the vector of grid values of radius r.
#' @param R is the maximum of the grid values.
#' @param nPerm is an integer denoting the umber of permutations to be used. Only used if perm = TRUE o
#'
#' @return It returns a list with the permutation-mean of summary functions
#'
#' @export


Perm_spat <-function(PP_obj, n_celltypes, subset, fixed_r, R, nPerm = 19)
{
  perm_func <- function(i){
    PP_perm <- PP_obj
    PP_perm$marks <- sample(PP_obj$marks)
    Kall <- alltypes(PP_perm, Kcross, r = fixed_r, correction = "isotropic")
    gall <- pcf(Kall, spar = 1, method="c", divisor ="d")
    return(list(Kall, gall))
  }
  #res <-  parallel::mclapply(X = as.list(1:nPerm), perm_func, mc.cores = 8)
  res <-  lapply(X = as.list(1:nPerm), perm_func)

  K_perm = L_perm = g_perm = array(0, dim = c(n_celltypes, n_celltypes, R-1))
  s = 1
  for(r1 in subset){
    for(r2 in subset){
      mean_Kall = mean_gall = 0
      for (i in 1:nPerm) {
        mean_Kall = mean_Kall + res[[i]][[1]]$fns[[s]]$iso[2:R]/nPerm
        mean_gall = mean_gall + res[[i]][[2]]$fns[[s]]$pcf[2:R]/nPerm
      }
      K_perm[ r1, r2, ] = mean_Kall #Kall$fns[[s]]$border Kall$fns[[s]]$iso
      L_perm[ r1, r2, ] =  sqrt(mean_Kall/pi) #sqrt(Kall$fns[[s]]$iso/pi)
      gfunc <- mean_gall
      gfunc[is.na(gfunc)] <- 1
      g_perm[ r1, r2, ] = gfunc
      s = s+1
    }
  }

  return(list(K_perm, L_perm, g_perm))
}
