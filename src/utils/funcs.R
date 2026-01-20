calculate_residuals = function(par, data, dt=0.01) {
  lambda = par[['lambda']]
  t = seq(0,max(data$day),dt)
  R = interpolate_rna(data$day, data$rna, t)
  P_pred = Pt(R, t, lambda)[match(data$day, t)]
  residuals = data$protein - P_pred
  return (residuals)
}

calculate_residuals_mixture = function(par, data, dt=0.01, K=2) {
  lambdas = par[paste0('lambda',1:K)]
  weights_raw = par[paste0('w',1:(K-1))]
  weights = c(weights_raw, 1-sum(weights_raw))
  if(any(weights < 0) | any(weights > 1)) {
    return(rep(1e10, length(data$day)))
  }
  t = seq(0,max(data$day),dt)
  R = interpolate_rna(data$day, data$rna, t)
  P_pred = Pt_mixture(R, t, lambdas, weights)$total[match(data$day, t)]
  residuals = data$protein - P_pred
  return (residuals)
}

Pt_mixture = function(R, t, lambdas, weights) {
  K = length(lambdas)
  P_components = matrix(0, nrow=length(t), ncol=K)
  P_components[1,] = weights
  P = numeric(length(t))
  P[1] = 1
  for (i in 2:length(t)) {
    dt_val = (t[i] - t[i-1])
    for (k in 1:K) {
      dP = lambdas[k] * dt_val * R[i-1] * weights[k] - lambdas[k] * dt_val * P_components[i-1, k]
      P_components[i, k] = P_components[i-1, k] + dP
    }
    P[i] = sum(P_components[i,])
  }
  return(list(total=P, components=P_components))
}


proportion_labeled_mixture = function(halflives, proportions, t, avails) {
  K = length(halflives)
  lambdas = log(2) / halflives

  # Track each subpopulation separately
  prop_heavy_sub = matrix(0, nrow=length(t), ncol=K)

  for (k in 1:K) {
    lambda_k = lambdas[k]
    for (i in 2:length(t)) {
      dt_val = t[i] - t[i-1]
      protein_turned_over = lambda_k * dt_val
      original = (1 - protein_turned_over) * prop_heavy_sub[i-1, k]
      nascent = protein_turned_over * avails(t[i])
      prop_heavy_sub[i, k] = original + nascent
    }
  }

  # Weight by proportions
  total_prop_heavy = prop_heavy_sub %*% proportions

  return(as.vector(total_prop_heavy))
}



runs_test = function(x) {
  x = x[!is.na(x) & x != 0]
  n = length(x)
  n_pos = sum(x > 0)
  n_neg = sum(x < 0)
  runs = rle(sign(x))
  n_runs = length(runs$lengths)
  expected_runs = 1 + (2 * n_pos * n_neg) / n
  var_runs = (2 * n_pos * n_neg * (2 * n_pos * n_neg - n)) / (n^2 * (n - 1))
  z_stat = (n_runs - expected_runs) / sqrt(var_runs)
  p_value = 2 * pnorm(-abs(z_stat))
  list(statistic = z_stat, p.value = p_value, n_runs = n_runs, n_pos = n_pos, n_neg = n_neg)
}

upper = function(x, ci=0.95) {
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x, na.rm=T) + sds*sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
}

lower = function(x, ci=0.95) {
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x, na.rm=T) - sds*sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
}

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x > 0 & signed, '+', ''),formatC(100*x,format='f',digits=digits),'%'))

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}



parxpdt = function(expr) {
  par(xpd=T)
  expr
  par(xpd=F)
}

rbind_files = function(path, grepstring) {
  all_files = list.files(path, full.names=T)
  these_files = all_files[grepl(grepstring,all_files)]
  if (exists('rbound_table')) rm('rbound_table')
  for (this_file in these_files) {
    this_tbl = read_delim(this_file, col_types=cols()) %>% clean_names()
    this_tbl$file = gsub('.*\\/','',gsub('\\.[tc]sv','',this_file))
    if (exists('rbound_table')) {
      rbound_table = rbind(rbound_table, this_tbl)
    } else {
      rbound_table = this_tbl
    }
  }
  return (rbound_table)
}

free_lysine = function(t) {
  pmax(0,1-0.503*(exp(-t*log(2)*0.799))-0.503*exp(-t*log(2)/39.423))
}

always_all = function(t) {
  1
}

proportion_labeled = function(thalf, t, avails) {
  lambda = log(2) / thalf
  proportion_heavy = numeric(length(t))
  proportion_heavy[1] = 0
  for (i in 2:length(t)) {
    protein_turned_over = lambda * dt
    original = (1 - protein_turned_over) * proportion_heavy[i - 1]
    nascent = protein_turned_over * avails(t[i])
    proportion_heavy[i] = original + nascent
  }
  return (proportion_heavy)
}


proportion_labeled_t_unary = function(thalf, t, which_t, avails) {
  return (proportion_labeled(thalf, t, avails)[t==which_t])
}

proportion_labeled_t = function(thalf, t, which_t, avails) {
  result = numeric(length(thalf))
  for (i in 1:length(result)) {
    result[i] = proportion_labeled(thalf[i], t, avails)[t==which_t]
  }
  return (result)
}


function_for_uniroot = function(thalf, prop_labeled, t, which_t, avails) {
  proportion_labeled_t(thalf, t, which_t, avails) - prop_labeled
}
find_thalf_unary = function(prop_labeled, t=seq(0,8,.01), which_t=8, avails=free_lysine) {
  max_possible = proportion_labeled_t(thalf=0.01, t=t, which_t=which_t, avails=avails)
  if (prop_labeled > max_possible) {
    return (as.numeric(NA))
  }
  uniroot(function_for_uniroot, lower=0.1, upper=100, tol=0.01, 
          prop_labeled = prop_labeled, t = t, which_t = which_t, avails=avails)$root
}

# vectorized version
find_thalf = function(prop_labeled, t=seq(0,8,.01), which_t=8, avails=free_lysine) {
  result = numeric(length(prop_labeled))
  for (i in 1:length(prop_labeled)) {
    result[i] = find_thalf_unary(prop_labeled[i], t, which_t, avails)
  }
  return (result)
}

interpolate_rna = function(day, rna, xout) {
  tibble(day, rna) %>%
    group_by(day) %>%
    summarize(.groups='keep', rna=mean(rna)) %>%
    ungroup() %>%
    add_row(day=0, rna=1, .before=1) -> actual_data_points
  approx(x=actual_data_points$day, y=actual_data_points$rna, xout=xout)$y
}

Pt = function(R, t, lambda) {
  P = numeric(length(t))
  dP = numeric(length(t))
  P[1] = 1
  for (i in 2:length(t)) {
    dt = (t[i]- t[i-1])
    dP[i] = lambda * dt * R[i-1] - lambda * dt * P[i-1]
    P[i] = P[i-1] + dP[i]
  }
  return (P)
}

calculate_residuals = function(par, data, dt=0.01) {
  lambda = par[['lambda']]
  t = seq(0,max(data$day),dt)
  R = interpolate_rna(data$day, data$rna, t)
  P_pred = Pt(R, t, lambda)[match(data$day, t)]
  residuals = data$protein - P_pred
  return (residuals)
}



label_difference = function(prop_labeled, chow_days, genotype) {
  # consolidate into 1 vs. 1 comparison of FFI vs. both control groups
  grp = case_when(genotype=='129(TT-3F4-FFI)HOZ' ~ 'FFI',
                  genotype!='129(TT-3F4-FFI)HOZ' ~ 'controls')
  # halflife calculation
  test_grp_thalf = fit_isotopic_thalf(chow_days[grp=='FFI'], prop_labeled[grp=='FFI'])
  ctrl_grp_thalf = fit_isotopic_thalf(chow_days[grp=='controls'], prop_labeled[grp=='controls'])
  thalf_ratio = test_grp_thalf/ctrl_grp_thalf
  # betareg model
  prop_labeled = pmin(pmax(prop_labeled,1e-6),1-1e-6) # betareg rejects 0 values so allow 1 ppm tolerance
  br_obj = betareg(prop_labeled ~ chow_days * grp)
  br_ffi_pval = summary(br_obj)$coefficients$mean['grpFFI','Pr(>|z|)']
  br_interaction_pval = summary(br_obj)$coefficients$mean['chow_days:grpFFI','Pr(>|z|)']
  # edit the return statement to decide which to use
  return (tibble(thalf_ratio=thalf_ratio, br_ffi_pval=br_ffi_pval, br_interaction_pval=br_interaction_pval))
}

