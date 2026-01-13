
## Figure 4 #### 
tell_user('done.\nCreating Figure 4...')
resx=300
png('display_items/figure-4.png',width=6.5*resx,height=4.5*resx,res=resx)

layout_matrix = matrix(c(1,1,1,1,
                         2,3,4,5,
                         6,6,6,6,
                         7,8,9,10), nrow=4, byrow=T)
layout(layout_matrix, heights=c(0.07, 1, 0.07, 1))

panel = 1

ffi_all = read_tsv('data/iqp/ffi_turnover.tsv', col_types=cols()) %>%
  mutate(total = light + heavy) %>%
  mutate(prop_labeled = heavy / (total))

leg = tibble(genotype=c("129(TT-3F4-FFI)HOZ","129(TT-3F4WT)","B6/N"),
             disp = c('ki-3F4-FFI','ki-3F4-WT','C57BL/6N WT'),
             color=c('#D95F02','#22127A','#77127A'),
             xgeno = c(1,2,3),
             ttest_grouping = c('test','control','control'))

# Read pre-computed genotypic differences from Julia analysis for thalf_effective
genotypic_diffs_julia = read_tsv('data/fig4/genotypic_diffs.tsv', col_types=cols()) %>%
  select(protein, peptide, age, ratio_thalf_julia = ratio_thalf)

ffi_all %>% 
  inner_join(leg, by='genotype') %>%
  group_by(protein, peptide) %>%
  mutate(n_ages_detected = length(unique(age))) %>%
  filter(n_ages_detected == 2) %>%
  ungroup() %>%
  select(-n_ages_detected) %>%
  group_by(protein, peptide, age) %>%
  summarize(.groups='keep',
            n = n(),
            pval_total = t.test(total[ttest_grouping=='test'], total[ttest_grouping=='control'])$p.value,
            ratio_total = mean(total[ttest_grouping=='test']) / mean(total[ttest_grouping=='control']),
            label_obj = label_difference(prop_labeled, chow_days, genotype),
            pval_label_ffi = label_obj$br_ffi_pval, 
            pval_label_interaction = label_obj$br_interaction_pval) %>%
  ungroup() %>%
  select(-label_obj) %>%
  left_join(genotypic_diffs_julia, by = c("protein", "peptide", "age")) %>%
  mutate(ratio_thalf = ratio_thalf_julia) %>% # Use Julia's ratio
  mutate(log2_ratio_total = log2(ratio_total),
         log2_ratio_thalf = log2(ratio_thalf)) %>%
  mutate(total_color = case_when(log2_ratio_total < 0 ~ '#FF0000',
                                 log2_ratio_total > 0 ~ '#0000FF',
                                 log2_ratio_total ==0 ~ '#000000')) %>%
  mutate(label_ffi_color = case_when(log2_ratio_thalf < 0 ~ '#FF0000',
                                     log2_ratio_thalf > 0 ~ '#0000FF',
                                     log2_ratio_thalf ==0 ~ '#000000'))%>%
  mutate(pepnick = substr(peptide, 1, 4)) -> genotypic_diffs_by_age

# Read mixture model parameters
mixture_params = read_tsv('data/fig4/table-mixture-model-results.tsv', col_types=cols())

# Define mixture model prediction function
proportion_labeled_mixture = function(thalves, proportions, t, avails_func=free_lysine) {
  # Simulate each component separately
  K = length(thalves)
  prop_heavy_components = matrix(0, nrow=length(t), ncol=K)
  
  for (k in 1:K) {
    lambda_k = log(2) / thalves[k]
    for (i in 2:length(t)) {
      protein_turned_over = lambda_k * dt
      prop_heavy_components[i, k] = (1 - protein_turned_over) * prop_heavy_components[i-1, k] +
        (protein_turned_over * avails_func(t[i]))
    }
  }
  
  # Weight by proportions
  total_prop_heavy = prop_heavy_components %*% proportions
  return(as.vector(total_prop_heavy))
}


for (this_age in c('young','aged')) {
  
  ### age header row ####
  
  par(mar=c(0,0,0,0))
  plot(NA, NA, xlim=0:1, ylim=0:1, axes=F, ann=F, xaxs='i', yaxs='i')
  if (this_age=='young') {
    mtext(side=3, line=-2, text='young (9 weeks)', cex=1)
  } else if (this_age=='aged') {
    mtext(side=3, line=-2, text='aged (63 weeks)', cex=1)
  }
  
  ### A/E abundance volcano ####
  
  genotypic_diffs_by_age %>%
    filter(age==this_age) -> ffi
  
  ffi$total_color = alpha(ffi$total_color, case_when(ffi$pval_total > 0.05/nrow(ffi) ~ 0.1,
                                                     TRUE ~ 1))
  ffi$label_ffi_color = alpha(ffi$label_ffi_color, case_when(ffi$pval_label_ffi > 0.05/nrow(ffi) ~ 0.1,
                                                             TRUE ~ 1))  
  
  par(mar=c(3,2.5,3,0.5))
  xlims = c(-3, 3)
  xbigs = -3:3
  xats = seq(-3, 3, .5)
  ylims = c(0, 10)
  ybigs = 0:10
  ybiglabs = as.character(0:10) %>% gsub('10','≥10', .)
  yats = seq(0, 10, .5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xbigs, labels=NA, tck=-0.05)
  axis(side=1, at=xbigs, lwd=0, line=-0.25)
  axis(side=1, at=xats, labels=NA, tck=-0.02)
  mtext(side=1, line=1.6, text='L2FC abundance', cex=0.6)
  axis(side=2, at=ybigs, labels=NA, tck=-0.05)
  axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.4)
  axis(side=2, at=yats, labels=NA, tck=-0.02)
  mtext(side=2, line=1.4, text='-log10(P)', cex=0.6)
  abline(h=-log10(0.05 / nrow(ffi)))
  par(xpd=T)
  points(x=ffi$log2_ratio_total, y=pmin(-log10(ffi$pval_total), max(ylims)), pch=20, col=ffi$total_color)
  ffi %>% 
    filter(protein=='PRNP') -> to_label
  points(x=to_label$log2_ratio_total, y=pmin(-log10(to_label$pval_total), max(ylims)), pch=1)
  for (i in 1:nrow(to_label)) {
    text(x=to_label$log2_ratio_total[i], y=pmin(-log10(to_label$pval_total[i]), max(ylims)), pos=4, labels=bquote(italic(.(to_label$protein[i])) * ' ' * .(to_label$pepnick[i])), cex=0.6)
  }
  par(xpd=F)
  mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1
  
  ### B/F half-life volcano ####
  
  par(mar=c(3,0.5,3,2.5))
  xlims = c(-3, 3)
  xbigs = -3:3
  xats = seq(-3, 3, .5)
  ylims = c(0, 10)
  ybigs = 0:10
  ybiglabs = as.character(0:10) %>% gsub('10','≥10', .)
  yats = seq(0, 10, .5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xbigs, labels=NA, tck=-0.05)
  axis(side=1, at=xbigs, lwd=0, line=-0.25)
  axis(side=1, at=xats, labels=NA, tck=-0.02)
  mtext(side=1, line=1.6, text='L2FC half-life', cex=0.6)
  axis(side=2, at=yats, labels=NA, tck=-0.02)
  abline(h=-log10(0.05 / nrow(ffi)))
  par(xpd=T)
  points(x=ffi$log2_ratio_thalf, y=pmin(-log10(ffi$pval_label_ffi), max(ylims)), pch=20, col=ffi$label_ffi_color)
  ffi %>% 
    filter(protein=='PRNP') -> to_label
  # amazingly, bquote is not vectorized? the only way to do this is in a loop??
  # text(x=to_label$log2_ratio_total, y=pmin(-log10(to_label$pval_total), max(ylims)), pos=4, labels=bquote(italic(.(to_label$protein)) * ' ' * .(to_label$pepnick)), cex=0.6)
  points(x=to_label$log2_ratio_thalf, y=pmin(-log10(to_label$pval_label_ffi), max(ylims)), pch=1)
  for (i in 1:nrow(to_label)) {
    is_vveq_f_panel = to_label$pepnick[i] == 'VVEQ' & panel == 6
    text(x=to_label$log2_ratio_thalf[i], 
      y=pmin(-log10(to_label$pval_label_ffi[i]), max(ylims)), 
      pos=ifelse(is_vveq_f_panel, 3, 4), # Keep original left/right positioning
      adj=if (is_vveq_f_panel) c(1.5, -0.5) else NULL, # Nudge VVEQ only
      labels=bquote(italic(.(to_label$protein[i])) * ' ' * .(to_label$pepnick[i])), 
      cex=0.6)
  }
  par(xpd=F)
  mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1
  
  ### C/D/G/H label accumulation by day ####
  
  use_peptides = c('GENFTETDVK','VVEQMCVTQYQK')
  
  for (this_peptide in use_peptides) {
    
    first_panel = this_peptide == use_peptides[1]
    if (first_panel) par(mar=c(3,3,3,0.25)) # make the first panel narrow on right
    else par(mar=c(3,0.25,3,3)) # make the second panel narrow on left
    
    xlims = c(-0.5, 8.5)
    xbigs = c(0, 2, 4, 6, 8)
    xats = 0:8
    ylims = c(0, .50)
    ybigs = 0:10/10
    ybiglabs = percent(ybigs)
    yats = 0:20/20
    plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
    axis(side=1, at=xbigs, labels=NA, tck=-0.05)
    axis(side=1, at=xbigs, lwd=0, line=-0.25)
    axis(side=1, at=xats, labels=NA, tck=-0.02)
    mtext(side=1, line=1.6, text='day', cex=0.6)
    axis(side=2, at=yats, labels=NA, tck=-0.02)
    if (first_panel) {
      axis(side=2, at=ybigs, labels=NA, tck=-0.05)
      axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.5)
      mtext(side=2, line=2.5, text='proportion labeled', cex=0.6)
    }
    mtext(side=3, text=substr(this_peptide,1,4), line=0, cex=0.6)
    
    for (this_genotype in unique(leg$genotype)) {
      
      subs = ffi_all %>%
        filter(genotype==this_genotype & peptide==this_peptide & age==this_age) %>%
        inner_join(leg, by='genotype')
      subs %>%
        group_by(protein, peptide, genotype, disp, chow_days, color) %>%
        summarize(.groups='keep',
                  mean = mean(prop_labeled),
                  l95 = lower(prop_labeled),
                  u95 = upper(prop_labeled)) %>%
        ungroup() -> this_smry
      
      # Plot individual data points
      points(x=subs$chow_days, y=subs$prop_labeled, col=subs$color, pch=1, bg='#FFFFFF')
      
      # Plot mixture model fit from t=0 to 8
      model_t = seq(0, 8, dt)
      params = mixture_params %>%
        filter(peptide==this_peptide & genotype==this_genotype & age==this_age)
      
      if (nrow(params) > 0 && !is.na(params$thalf_fast) && !is.na(params$thalf_slow)) {
        model_pred = proportion_labeled_mixture(
          thalves = c(params$thalf_fast, params$thalf_slow),
          proportions = c(params$prop_fast, params$prop_slow),
          t = model_t,
          avails_func = free_lysine
        )
        points(x=model_t, y=model_pred, type='l', col=this_smry$color[1], lty=1, lwd=2)
      }
      
      legend('topleft', leg$disp, col=leg$color, lwd=1, text.col=leg$color, bty='n', cex=0.7)
      
    }
    mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1
    
  }
  
}

silence_is_golden = dev.off()
### end Fig 4 ####
