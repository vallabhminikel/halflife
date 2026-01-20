

## Figure S3 ####
tell_user('done.\nCreating Figure S3...')
resx=300
png('display_items/figure-s3.png',width=13*resx,height=16*resx,res=resx)

layout_matrix = matrix(1:40, nrow=10, byrow=F)
layout(layout_matrix, heights=rep(c(0.3,1),5))

loq = read_tsv('data/loq.tsv', col_types=cols())

reference_dilution = 12

loq %>%
  mutate(total_area_fragment = replace_na(total_area_fragment,0)) %>%
  group_by(peptide, protein, dilution, attomoles_peptide_analyzed) %>%
  summarize(.groups='keep',
            mean_area = mean(total_area_fragment),
            cv_area = sd(total_area_fragment)/mean(total_area_fragment)) %>%
  ungroup() %>%
  group_by(peptide, protein) %>%
  mutate(expected_area = (0.5^(12-dilution))*mean_area[dilution==12]) %>%
  mutate(relative_error = (mean_area - expected_area)/expected_area) %>%
  ungroup() -> loq_stats

cumcv_threshold = 0.20

# determine CV
loq_stats %>%
  arrange(protein, peptide, desc(dilution)) %>%
  group_by(protein, peptide) %>%
  mutate(cumcv = cummean(cv_area)) %>% 
  mutate(lloq = suppressWarnings(pmin(12,min(dilution[cumcv < cumcv_threshold], na.rm=T)))) %>%
  ungroup() %>% 
  mutate(is_lloq = dilution==lloq) -> lloq_determination

lloq_determination %>%
  filter(is_lloq) %>%
  select(protein, peptide, mean_area) %>%
  mutate(gene = gsub('_.*','',protein)) -> lloqs

wt_means = pivot_tab %>% 
  inner_join(iqp_meta, by=c('identifier'='animal_id')) %>%
  filter(genotype=='C57BL/6N', labeled=='yes') %>%
  group_by(tissue, peptide) %>%
  summarize(.groups='keep', mean_heavy=mean(heavy)) %>%
  ungroup()

wt_means %>%
  inner_join(lloqs, by=c('peptide')) %>%
  rename(wt_mean_heavy = mean_heavy, lloq_area=mean_area) %>%
  mutate(lloq_ratio = lloq_area / wt_mean_heavy) %>%
  select(tissue, peptide, lloq_ratio) -> lloq_vs_wt_stats

write_supp_table(lloq_vs_wt_stats, 'Ratio of LLOQ to mean heavy peptide area in WT animals labeled for 8 days.')

for (pep in unique(loq_stats$peptide)) {
  
  this_gene = lloqs$gene[lloqs$peptide==pep]
  
  subs = loq_stats %>% filter(peptide==pep)
  xlims = c(10,50000)
  label_x = 3e1
  par(mar=c(0.5, 5, 3, 5))
  ylims = c(0,1)
  plot(NA, NA, xlim=xlims, ylim=ylims, log='x', axes=F, ann=F, xaxs='i', yaxs='i')
  points(x=subs$attomoles_peptide_analyzed, y=subs$cv_area, type='h', col='gray', lwd=10, lend=1)
  axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
  axis(side=2, at=0:2/2, labels=percent(0:2/2), las=2)
  mtext(side=2, line=3.5, text='%CV')
  abline(h=.20, col='red', lty=3, lwd=0.5)
  mtext(side=3, text=this_gene, font=3, line=1, cex=0.9)
  mtext(side=3, text=pep, line=0, cex=0.7)
  
  par(mar=c(4,5,0.5,5))
  ylims = c(1e3, 1e9)
  plot(NA, NA, xlim=xlims, ylim=ylims, log='xy', axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=subs$attomoles_peptide_analyzed, labels=formatC(subs$attomoles_peptide_analyzed,format='e',digits=0))
  mtext(side=1, line=2.5, text='attomoles of peptide')
  axis(side=2, las=2)
  mtext(side=2, line=3.5, text='mean area')
  points(subs$attomoles_peptide_analyzed, subs$mean_area, pch=20)
  points(subs$attomoles_peptide_analyzed, subs$expected_area, col='gray', lty=3, type='l')
  wt_means %>%
    filter(peptide==pep) -> subs
  abline(h=subs$mean_heavy, col='blue', lty=3)
  mtext(side=4, at=subs$mean_heavy, col='blue', text=paste0(subs$tissue, ''), las=2, cex=0.8)
  lloqs %>%
    filter(peptide==pep) %>%
    select(mean_area) %>%
    pull() -> this_lloq
  abline(h=this_lloq, col='red', lty=3)
  mtext(side=4, at=this_lloq, col='red', text='LLQ', las=2, cex=0.8)
  
  
}
silence_is_golden = dev.off()
### end Figure S3 ####
