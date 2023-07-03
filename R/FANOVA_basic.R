UniFANOVA = function(data = data, contrast = NULL){
  group.label0 = unique(data$Group)
  subject.label0 = unique(data$ID)

  G = length(group.label0) #Number of groups
  sum_m_g = length(subject.label0)  #Total number of subjects
  p = (dim(data)[2] - 2) #The number of radii used

  group_sub_table = suppressMessages(data %>% group_by(Group, ID) %>% summarise(n = n()))
  group_table_subject_count = group_sub_table %>% group_by(Group) %>% summarise(sub_count = n())

  #Model matrix W_gr based on group ID (same for every r)--------
  W_gr = model.matrix(~Group - 1, data = data[, c(1:2)])
  WTW_gr_inv = solve(t(W_gr)%*%W_gr)

  if(is.null(contrast) == T){
    C_matrix = matrix(0, nrow = G - 1, ncol = G)
    C_matrix[, 1] = 1
    for(gr in 1:(G-1)){
      C_matrix[gr, (gr + 1)] = -1
    }
  }else{
    if(nrow(contrast) > (G - 1)){print("Contrast has more rows than groups!")
    }else{
      C_matrix = contrast
    }
  }

  contrast_rows = (dim(C_matrix)[1])
  z = matrix(0, nrow = sum_m_g, ncol = p)
  sumse = sumg = sumt = frat  = 0

  r = 1
  for(rval in 1:p)
  {
    y = data[, (rval + 2)]
    mu_hat = WTW_gr_inv%*%t(W_gr)%*%y
    test = t(C_matrix %*% mu_hat)%*%solve(C_matrix%*%WTW_gr_inv%*%t(C_matrix))%*%(C_matrix %*% mu_hat)

    sumg = c(sumg, test)
    z[,r] = y - W_gr %*% mu_hat

    se = sum(z[,r]^2)
    tot = sum((y - mean(y))^2)

    sumse = c(sumse, se)
    sumt = c(sumt, tot)

    frat = c(frat, test/se*(sum_m_g - G)/contrast_rows)
    r = r+1
  }

  sumg = sum(sumg)
  sumse = sum(sumse)
  sumt = sum(sumt)

  z = z/(as.matrix(rep(1, sum_m_g)) %*% sqrt(colSums(z^2)))
  if(sum_m_g >= p){z = t(z) %*% z}else{z = z %*% t(z)}
  statGPF = mean(frat[-1])


  trace_gammasq = sum(z * z)
  betaGPF = (trace_gammasq/p^2/contrast_rows)/((sum_m_g - G)/(sum_m_g-G-2))
  dGPF = ((sum_m_g - G)/(sum_m_g-G-2))^2/((trace_gammasq)/p^2/contrast_rows)

  pvalueGPF = pchisq(statGPF/betaGPF,dGPF, lower.tail = FALSE)

  res = list(pval = pvalueGPF, Fvals = list(frat[-1]))
}

MFanova.tests = function(data = data, contrast = NULL, parallel = FALSE, nslaves = NULL)
{
  group.label0 = unique(data$Group)
  subject.label0 = unique(data$ID)
  image.label0 = unique(data$imageID)

  G = length(group.label0) #Number of groups
  sum_m_g = length(subject.label0)  #Total number of subjects
  N = length(image.label0) #Total number of images
  p = length(unique(data$range)) #The number of radii used

  group_sub_im_table = suppressMessages(data %>% group_by(Group, ID, imageID) %>% summarise(n = n()))
  group_sub_table = suppressMessages(group_sub_im_table %>% group_by(Group, ID) %>% summarise(n = n()))
  group_table_image_count = group_sub_table %>% group_by(Group) %>% summarise(image_count = sum(n))
  group_table_subject_count = group_sub_table %>% group_by(Group) %>% summarise(sub_count = n())

  #Cellmeans model matrix W based on subject ID (same for every r)--------
  r_init = unique(data$range)[1]
  W = model.matrix(~ID - 1, data = data[data$range == r_init, ])
  W = W[,paste0("ID", subject.label0)] #Reordering based on original ID order
  WTW_inv = solve(t(W)%*%W)

  #Model matrix W_gr based on group ID (same for every r)--------
  W_gr = model.matrix(~Group - 1, data = data[data$range == r_init, ])
  WTW_gr_inv = solve(t(W_gr)%*%W_gr)
  #gs = diag(t(W_gr)%*%W_gr) #Number of images per group

  small_hi = group_sub_table$n/rep(group_table_image_count$image_count, group_table_subject_count$sub_count)
  #Number of images per subject
  #scaled by the # images per group
  if(is.null(contrast) == T){
    C_matrix = matrix(0, nrow = G - 1, ncol = sum_m_g)
    C_matrix[, 1:group_table_subject_count$sub_count[1]] = 1
    for(gr in 1:(G-1)){
      C_matrix[gr, (cumsum(group_table_subject_count$sub_count)[gr] + 1):(cumsum(group_table_subject_count$sub_count)[gr + 1])] = -1
      C_matrix[gr, ] = C_matrix[gr, ]*small_hi
    }
  }else{
    if(nrow(contrast) > (G - 1)){print("Contrast has more rows than groups!")
    }else{
      C_matrix = matrix(0, nrow = nrow(contrast), ncol = sum_m_g)
      for(gr in 1:nrow(contrast)){
        C_matrix[gr, ] = contrast[gr, ]*small_hi
      }
    }
  }

  contrast_rows = (dim(C_matrix)[1])
  z = z_noint = matrix(0, nrow = N, ncol = p)
  frat = f_int_rat = f_ind_rat = 0

  r = 1
  for(rval in unique(data$range))
  {
    dat = data[data$range == rval, ]
    y = dat$func
    mu_hat = WTW_inv%*%t(W)%*%y
    test = t(C_matrix %*% mu_hat)%*%solve(C_matrix%*%WTW_inv%*%t(C_matrix))%*%(C_matrix %*% mu_hat)

    z[,r] = y - W %*% mu_hat
    z_noint[,r] = y - W_gr %*% WTW_gr_inv%*%t(W_gr)%*%y

    se = sum(z[,r]^2)
    tot = sum((y - mean(y))^2)
    int = tot - test - se

    frat = c(frat, test/se*(N - sum_m_g)/contrast_rows)
    f_int_rat = c(f_int_rat, int/se*(N - sum_m_g)/(sum_m_g - G))
    f_ind_rat = c(f_ind_rat, test/(int + se)*(N - contrast_rows)/contrast_rows)
    r = r+1
  }

  z = z/(as.matrix(rep(1, N)) %*% sqrt(colSums(z^2)))
  if(N >= p){z = t(z) %*% z}else{z = z %*% t(z)}
  statGPF = mean(frat[-1])
  intGPF = mean(f_int_rat[-1])
  indGPF = mean(f_ind_rat[-1])


  trace_gammasq = sum(z * z)
  betaGPF_int = (trace_gammasq/p^2/(sum_m_g - G))/((N-sum_m_g)/(N-sum_m_g-2))
  dGPF_int = ((N-sum_m_g)/(N-sum_m_g-2))^2/((trace_gammasq)/p^2/(sum_m_g - G))

  betaGPF = (trace_gammasq/p^2/contrast_rows)/((N-sum_m_g)/(N-sum_m_g-2))
  dGPF = ((N-sum_m_g)/(N-sum_m_g-2))^2/((trace_gammasq)/p^2/contrast_rows)

  z_noint = z_noint/(as.matrix(rep(1, N)) %*% sqrt(colSums(z_noint^2)))
  if(N >= p){z_noint = t(z_noint) %*% z_noint}else{z_noint = z_noint %*% t(z_noint)}

  trace_gammasq_noint = sum(z_noint * z_noint)
  betaGPF_noint = (trace_gammasq_noint/p^2/contrast_rows)/((N-G)/(N-G-2))
  dGPF_noint = ((N-G)/(N-G-2))^2/((trace_gammasq_noint)/p^2/contrast_rows)

  pvalueGPF_int = pchisq(intGPF/betaGPF_int, dGPF_int, lower.tail = FALSE)
  pvalueGPF = ifelse(pvalueGPF_int < 0.05, pchisq(indGPF/betaGPF_noint,dGPF_noint,
                                                  lower.tail = FALSE),
                     pchisq(statGPF/betaGPF,dGPF, lower.tail = FALSE))

  res = list(GPF_int = list(pval = pvalueGPF_int),
             GPF = list(pval = pvalueGPF),
             Fvals = list(frat[-1]))
  return(res)
}
