tell_user('done.\nCreating Figure 2...')


crl_prot = read_tsv('data/crl_proteomics.tsv', col_types=cols()) %>%
  mutate(heavy = replace_na(heavy, bql)) %>%
  mutate(total = light + heavy) %>%
  mutate(proportion_heavy = heavy / total) %>%
  inner_join(tissue_meta, by ='tissue')


pivot_tab %>%
  inner_join(iqp_meta, by=c('identifier'='animal_id')) %>%
  mutate(total = heavy + light) %>%
  mutate(prop_labeled = heavy / total) %>%
  inner_join(lloq_stats, by=c('peptide','protein'='gene')) %>%
  mutate(above_lloq = pmin(light,heavy) >= lloq) -> iqp_data


genotype_meta = tibble(genotype = c("C57BL/6N","Tga20 het","ZH3/+","Tg25109 het; ZH3/ZH3"),
                       explanation = c('WT','Tga20','het KO','HuPrP'),
                       xgeno = c(1,2,3,4))

prnp_pep_meta = tibble(peptide = c("VVEQMCVTQYQK", "GENFTETDVK", "QHTVTTTTK"),
                       short = c('VVEQ','GENF','QHTV'),
                       pep_color = c("#e5b6b6", "#849cb5", "#a2bca2"),
                       xpep = c(2,1,3))

tissue_meta = tibble(tissue=c('brain','colon'),
                     color = c('#CCCF00','#EEBB77'),
                     xtiss=c(1,2))

adjusted_free_lysine = function(x) {
  return (function(t) { pmin(1,x * free_lysine(t)) })
}


paired_samples = c(12, 23, 24)

tissue_colors = c(brain='#CCCF00', colon='#EEBB77')


read_tsv('data/iqp/quantified_lh_hh.tsv', col_types=cols()) %>%
  clean_names() %>%
  mutate(tissue = tolower(tissue)) %>%
  filter(iqp_id %in% paired_samples) %>%
  select(tissue, iqp_id, pg_genes, pg_group_label, peptide, charge, protein_abundance, 
         light_light, light_heavy, heavy_heavy, log_min_signal_binned) %>%
  mutate(ria = 2/(2 + (light_heavy / heavy_heavy))) -> lh_hh

write_supp_table(lh_hh, 'LH/HH data used to calculate relative isotope abundance (RIA).')


lh_hh %>%
  inner_join(tissue_meta, by='tissue') %>%
  rename(x = xtiss) %>%
  group_by(tissue, iqp_id, x, color) %>% # summarize across all peptides
  summarize(.groups='keep',
            n    = n(),
            mean_ria = mean(ria, na.rm=T),
            l95 = lower(ria),
            u95 = upper(ria)) %>%
  ungroup() -> individual_lh_hh

individual_lh_hh %>%
  group_by(tissue, x, color) %>%
  summarize(.groups='keep',
            n    = n(),
            mean = mean(mean_ria),
            l95 = lower(mean_ria),
            u95 = upper(mean_ria)) %>%
  ungroup() -> lh_hh_tissue_smry 

colon_ria_multiplier = lh_hh_tissue_smry$mean[lh_hh_tissue_smry$tissue=='colon'] / lh_hh_tissue_smry$mean[lh_hh_tissue_smry$tissue=='brain']


resx=300
png('display_items/figure-2.png',width=6.5*resx,height=7.0*resx,res=resx)

layout_matrix = matrix(c(1,1,1,1,1,1,1,1,1,
                         2,2,2,3,3,3,4,4,4,
                         5,5,5,5,6,7,7,8,8), nrow=3, byrow=T)
layout(layout_matrix, 
       heights = c(1, 1, 1),
       widths =  c(1,1,1,0.5,1.5,1,1,1,1))

panel = 1

### peptide diagram ####
par(mar=c(0.25,0,1,0.25))
peptide_panel = image_convert(image_read('data/misc/peptide-diagram.png'),'png')
plot(as.raster(peptide_panel))
mtext(side=3, adj=0.05, text=LETTERS[panel], line=-0.75); panel = panel + 1

### A CRL percent labeled ####



par(mar=c(3,4,3,3))
xlims = c(-1, 9)
xbigs = c(0, 2, 4, 6, 8)
xats = 0:8
ylims = c(0, .70)
ybigs = 0:10/10
ybiglabs = percent(ybigs)
yats = 0:20/20
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, lwd=0, line=-0.25)
axis(side=1, at=xats, labels=NA, tck=-0.02)
mtext(side=1, line=1.6, text='day')
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.5)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=2.5, text='proportion labeled', cex=0.8)
points(crl_prot$day, crl_prot$proportion_heavy, col=crl_prot$color, pch=21, cex=0.5, bg='#FFFFFF')
for (this_tissue in unique(tissue_meta$tissue)) {
  subs = crl_prot %>% filter(tissue==this_tissue)
  subs %>%
    filter(heavy==bql & day==0) %>%
    summarize(effective_bql = max(proportion_heavy)) %>%
    pull(effective_bql) -> ef_bql
  abline(h=ef_bql, lty=3)
  mtext(side=4, at=ef_bql, text=paste0(this_tissue,' LLQ'), cex=0.6, las=2)
  subs %>%
    group_by(day, color) %>%
    summarize(.groups='keep',
              mean = mean(proportion_heavy),
              l95 = lower(proportion_heavy),
              u95 = upper(proportion_heavy)) -> this_smry
  points(x=this_smry$day, y=this_smry$mean, type='l', lwd=2, col=this_smry$color)
  polygon(x=c(this_smry$day, rev(this_smry$day)), y=c(this_smry$l95, rev(this_smry$u95)), col=alpha(this_smry$color, ci_alpha), border=NA)
  mtext(side=4, at=this_smry$mean[this_smry$day==8], line=-0.25, text=this_tissue, las=2, col=this_smry$color, cex=0.8)
}
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1


crl_prot %>% 
  group_by(tissue, day, color) %>%
  summarize(.groups='keep',
            mean = mean(proportion_heavy),
            l95 = lower(proportion_heavy),
            u95 = upper(proportion_heavy)) %>%
  ungroup()  -> crl_smry

write_supp_table(crl_smry, 'Labeling of VVEQ peptide in Charles River proteomics data.')

### B free lysine model ####
par(mar=c(3,4,3,2))
xlims = c(-0.5, 8)
xbigs = c(0, 2, 4, 6, 8)
xats = 0:8
ylims = c(0, 0.7)
ybigs = 0:10/10
ybiglabs = percent(ybigs)
yats = 0:20/20
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, lwd=0, line=-0.25, cex.axis=0.8)
axis(side=1, at=xats, labels=NA, tck=-0.02)
mtext(side=1, line=1.6, text='day')
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.5)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=2.5, text='expected proportion labeled', cex=0.8)

lysine_meta = tibble(model = c('100%','plasma free'),
                     color = c("#429898", "#A90101"))

points(x=t, y=free_lysine(t), type='l', lwd=2, col=lysine_meta$color[lysine_meta$model=='plasma free'])
text(x=4, y=.55, col=lysine_meta$color[lysine_meta$model=='plasma free'], labels='plasma free lysine', cex=0.7)
par(xpd=T)
text(x=max(xlims)*1.15, y=.40, srt=270, col='black', las=2, 'protein half-life')
par(xpd=F)

t = seq(0, 8, dt)
thalf_values = c(1:5, 8, 10, 15, 20)

alt_hl_color = '#A9A9A9'
for (thalf in thalf_values) {
  lines(t, proportion_labeled(thalf, t, free_lysine), lwd = 1,col=alt_hl_color)
  par(xpd=T)
  text(t[length(t)], proportion_labeled_t(thalf, t, t[length(t)], free_lysine), labels = thalf, pos = 4,cex=0.7)
  par(xpd=F)
}
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1



### C ratios to 5-day ####

par(mar=c(3,4,3,2))
xlims = c(-0.5, 16)
xbigs = 0:8*2
xats = 0:16
ylims = c(0, 2)
ybigs = 0:4/2
ybiglabs = ybigs
yats = 0:8/4
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, lwd=0, line=-0.25, cex.axis=0.8)
axis(side=1, at=xats, labels=NA, tck=-0.02)
mtext(side=1, line=1.6, text='day')
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.25, cex.axis=0.8)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=2.25, text='ratio vs. 5-day half-life', cex=0.8)
abline(h=c(0.5, 1.5), lty=3,col ="red")
longer_t = seq(0, 16, dt)
answer_for_5 = proportion_labeled(5, longer_t, free_lysine)
for (thalf in thalf_values) {
  if (thalf==5) {
    this_color = '#000000'
  } else {
    this_color = alt_hl_color
  }
  answer_for_this = proportion_labeled(thalf, longer_t, free_lysine)
  ratio_this_to_5 = answer_for_this / answer_for_5
  lines(longer_t, ratio_this_to_5, lwd = 1,col=this_color)
  par(xpd=T)
  text(longer_t[length(longer_t)], ratio_this_to_5[length(longer_t)], labels = thalf, pos = 4, cex=0.5)
  par(xpd=F)
}
par(xpd=T)
text(x=max(xlims)*1.15, y=1.0, srt=270, col='black', las=2, 'protein half-life')
par(xpd=F)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

### D by PRNP genotype ####
iqp_data %>%
  filter(labeled == 'yes') %>%
  filter(protein == 'PRNP') %>%
  filter(!(genotype=='Tg25109 het; ZH3/ZH3' & peptide=='VVEQMCVTQYQK')) %>%
  filter(above_lloq) %>%
  inner_join(prnp_pep_meta, by='peptide') %>%
  inner_join(tissue_meta, by='tissue') %>%
  inner_join(genotype_meta, by='genotype') %>%
  mutate(x = (xtiss-1)*8 + (xgeno-1)*2.2 + (xpep-1)) -> all_prnp_plab


all_prnp_plab %>%
  group_by(genotype, explanation, tissue, peptide, x, color, pep_color, short) %>%
  summarize(.groups='keep',
            n = n(),
            mean = mean(prop_labeled),
            l95 = lower(prop_labeled),
            u95 = upper(prop_labeled)) %>%
  ungroup() %>%
  mutate(estimated_thalf_unadjusted = find_thalf(mean, avails=free_lysine),
    estimated_thalf_adjusted = find_thalf(mean, avails=adjusted_free_lysine(colon_ria_multiplier))) %>%
  mutate(estimated_thalf = case_when(tissue=='brain' ~ estimated_thalf_unadjusted,
                                     tissue=='colon' ~ estimated_thalf_adjusted)) %>%
  select(-estimated_thalf_unadjusted, -estimated_thalf_adjusted) -> smry

smry$pval = as.numeric(NA)
for (i in 1:nrow(smry)) {
  this_tissue = smry$tissue[i]
  this_peptide = smry$peptide[i]
  this_genotype = smry$genotype[i]
  all_prnp_plab %>%
    filter(tissue==this_tissue, peptide==this_peptide) %>%
    filter(genotype=='C57BL/6N') %>%
    pull(prop_labeled) -> controls
  all_prnp_plab %>%
    filter(tissue==this_tissue, peptide==this_peptide) %>%
    filter(genotype==this_genotype) %>%
    pull(prop_labeled) -> test_group
  tobj = t.test(controls, test_group)
  smry$pval[i] = tobj$p.value
}
testing_burden = sum(smry$genotype != 'C57BL/6N')
smry$pbonf = pmin(1, smry$pval*testing_burden)

write_supp_table(smry, 'PrP peptide labeling and estimated half-life by genotype in IQ Proteomics data.')

smry %>% 
  group_by(tissue) %>% 
  summarize(.groups='keep',
            mean_labeled=mean(mean)) %>%
  ungroup() -> prp_smry_by_tissue

write_supp_table(prp_smry_by_tissue, 'Mean PrP peptide labeling by tissue in IQ Proteomics data.')


prnp_iqp_smry = smry # save to use in next panel

par(mar=c(5,4,3,1))
xlims = range(smry$x) + c(-0.7, 0.7)
ylims = c(0, 0.7)
ybigs = 0:10/10
ybiglabs = percent(ybigs)
yats = 0:20/20
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, line=0.2, at=smry$x, text=smry$short, cex=0.35, las=2)
smry %>%
  group_by(genotype, explanation, tissue) %>%
  summarize(.groups='keep',
            minx = min(x),
            maxx = max(x),
            midx = mean(x)) %>%
  ungroup() -> tranches
overhang = 0.5
tranche_line = 2.0
for (i in 1:nrow(tranches)) {
  axis(side=1, line=tranche_line, at=c(tranches$minx[i]-overhang,tranches$maxx[i]+overhang), tck=0.03, labels=NA)
  mtext(side=1, line=tranche_line-0.4, padj=1, at=c(tranches$midx[i]), text=gsub(' ', '\n', tranches$explanation[i]), cex=0.4)
}
smry %>%
  group_by(tissue) %>%
  summarize(.groups='keep',
            minx = min(x),
            maxx = max(x),
            midx = mean(x)) %>%
  ungroup() -> tranches
tranche_line = 0.3
for (i in 1:nrow(tranches)) {
  axis(side=3, line=tranche_line, at=c(tranches$minx[i]-overhang,tranches$maxx[i]+overhang), tck=0.03, labels=NA)
  mtext(side=3, line=tranche_line+0.2, padj=0, at=c(tranches$midx[i]), text=gsub(' ','\n',tranches$tissue[i]), cex=0.7)
}
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.5, cex.axis=0.8)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=2.5, text='proportion labeled', cex=0.8)
barwidth=0.9
rect(xleft=smry$x-barwidth/2, xright=smry$x+barwidth/2, ybottom=rep(0,nrow(smry)), ytop=smry$mean, col=alpha(smry$color,ci_alpha), lwd=1.5, border=NA)
par(xpd=T)
set.seed(1)
points(jitter(all_prnp_plab$x,amount=0.1), all_prnp_plab$prop_labeled, col=all_prnp_plab$color, pch=21, bg='#FFFFFF')
arrows(x0=smry$x, y0=smry$l95, y1=smry$u95, code=3, angle=90, length=0.05, col='#000000', lwd=1.5)
par(xpd=F)
mtext(side=3, adj=0, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(individual_lh_hh, 'Summary of free lysine calculated from LH/HH ratios by individual animal and tissue.')

write_supp_table(lh_hh_tissue_smry, 'Summary of free lysine calculated from LH/HH ratios by tissue.')

par(mar=c(3,3,3,1))
xlims = c(0.5, 2.5)
ylims = c(0, 1)
ybigs = 0:4/4
yats  = 0:10/10
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=percent(ybigs), lwd=0, las=2, line=-0.5)
axis(side=2, at=yats,  labels=NA, tck=-0.02)
mtext(side=2, line=2.5, text='RIA', cex=0.8)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, line=-0.1, at=tissue_meta$xtiss, text=tissue_meta$tissue, cex=0.65)
points(x=individual_lh_hh$x, y=individual_lh_hh$mean_ria,
       pch=21, col=individual_lh_hh$color, bg='#FFFFFF', cex=1)
arrows(x0=lh_hh_tissue_smry$x, y0=lh_hh_tissue_smry$l95, y1=lh_hh_tissue_smry$u95,
       code=3, angle=90, length=0.05, lwd=2, col=lh_hh_tissue_smry$color)
barwidth = 0.4
segments(x0=lh_hh_tissue_smry$x-barwidth, x1=lh_hh_tissue_smry$x+barwidth, y0=lh_hh_tissue_smry$mean, lwd=1.5, col=lh_hh_tissue_smry$color)
mtext(side=3, adj=0, text=LETTERS[panel], line=0.5); panel = panel + 1


### G & H theoretical vs. diff prots ####

reported = lit_half %>%
  mutate(reported_halflife = case_when(is.na(cerebellum) ~ cortex,
                                       is.na(cortex) ~ cerebellum,
                                       TRUE ~ (cortex + cerebellum)/2)) %>%
  inner_join(name_map, by=c('gene'='fornasiero_name')) %>%
  filter(iqp_name %in% iqp_data$protein)

# check none are missing
# length(unique(iqp_data$protein))
# setdiff(unique(iqp_data$protein), reported$iqp_name)

for (this_tissue in c('brain','colon')) {
  
  if (this_tissue == 'brain') {
    par(mar = c(3, 3, 3, 0.5))
  } else {
    par(mar = c(3, 0.5, 3, 3))
  }
  
  xlims = c(0, 12)
  ylims = c(0, 1)
  
  plot(NA, NA, xlim = xlims, ylim = ylims, xaxs = 'i', yaxs = 'i', axes = F, ann = F)
  axis(1, at = xlims, labels = NA, tck = 0)
  axis(1, at = c(0:10 * 2), labels = c(0:10 * 2), cex.axis=0.7, lwd=0, line=-0.5)
  axis(1, at = c(0:10 * 2), tck = -0.05, labels=NA)
  axis(1, at = c(0:12 * 1), tck = -0.025, labels = NA)
  mtext(side = 1, line = 1.6, text = 'protein half life', cex = 0.7)
  
  if (this_tissue=='brain') {
    axis(side = 2, at = 0:5 / 4, labels = scales::percent(0:5 / 4), las = 2, lwd = 0, line=-0.5, cex.axis=0.8)
    mtext(side = 2, line = 2.25, text = 'proportion labeled', cex = 0.7)
  }
  axis(2, at = ylims, labels = NA)
  axis(2, at = c(0:4 * 0.25), labels = NA, cex = 0.8, tck = -0.05)
  axis(2, at = c(0:8 * 0.125), labels = NA, cex = 0.8, tck = -0.025)
  mtext(side = 3, line = 0, text = this_tissue, cex = 0.8)
  
  prnp_iqp_smry %>%
    filter(tissue==this_tissue) %>%
    filter(genotype=='C57BL/6N') %>%
    arrange(desc(mean)) -> subs
  
  if (this_tissue=='colon') {
    avails=adjusted_free_lysine(colon_ria_multiplier)
    subs -> colon_subs
  } else if (this_tissue=='brain') {
    avails=free_lysine
    subs -> brain_subs
  }
  
  thalf_values = seq(0, 12, .1)
  points(thalf_values, proportion_labeled_t(thalf_values, t=t, which_t = 8, avails=avails), type = 'l', col = lysine_meta$color[2], lwd = 2)
  
  prnp_color = '#0001CD'
  
  label_x = 1.2
  segments(x0=rep(0, nrow(subs)), x1=find_thalf(subs$mean,avails=avails), y0=subs$mean, col=prnp_color, lty=3)
  text(x=rep(label_x,2), y=subs$mean+c(-0.02,0.02), pos=c(3,1), labels=subs$short, col="#0001CD", cex=0.5)
  segments(x0=find_thalf(subs$mean,avails=avails), y0=rep(0,nrow(subs)), y1=subs$mean, col=prnp_color, lty=3)
  
  lloq_meta = tibble(above_lloq = c(1,0,0),
                     color = c('#000000','#C9C9C9','#BFB7EF'),
                     disp = c('OK','below LLOQ','controls'))
  
  iqp_data %>%
    filter(genotype=='C57BL/6N',
           tissue == this_tissue) %>%
    filter(protein != 'PRNP') %>%
    group_by(protein, peptide, labeled) %>%
    summarize(.groups='keep',
              n = n(),
              mean = mean(prop_labeled, na.rm=T),
              prop_above_lloq = sum(above_lloq)/n) %>%
    ungroup() %>%
    mutate(is_prnp = protein=='PRNP') %>%
    mutate(disp = case_when(prop_above_lloq >= 0.5 ~ 'OK',
                            labeled == 'no' ~ 'controls',
                            TRUE ~ 'below LLOQ')) %>%
    inner_join(lloq_meta, by=c('disp')) %>%
    inner_join(reported, by=c('protein'='iqp_name')) -> smry
  
  write_supp_table(smry, paste0('Labeling versus reported half-life in ',this_tissue,'.'))
  
  points(smry$reported_halflife, smry$mean, pch=1, cex = 1.2, col =smry$color)
  
  legend("topright", legend = lloq_meta$disp, col =lloq_meta$color, pch = 1, cex = 0.5, bty='n')
  
  mtext(side=3, adj=0, text=LETTERS[panel], line=0.5); panel = panel + 1
}
silence_is_golden = dev.off()
### end Figure 2 ####






resx=300
png('display_items/figure-s4.png',width=6.5*resx,height=4.0*resx,res=resx)

layout_matrix = matrix(c(1,2), nrow=1, byrow=T)
layout(layout_matrix)

panel = 1

tissue_colors = c(brain='#CCCF00', colon='#EEBB77')

x_map = tibble(tissue = rep(c('brain','colon'), each=length(paired_samples)),
               iqp_id = rep(paired_samples, 2),
               x      = c(1:3, 5:7)) %>%
  mutate(color = tissue_colors[tissue],
         label = paste('animal', match(iqp_id, paired_samples)))

plot_dat = lh_hh %>%
  filter(iqp_id %in% paired_samples) %>%
  mutate(ria = 2/(2 + (light_heavy / heavy_heavy))) %>%
  inner_join(x_map, by=c('tissue','iqp_id'))


plot_dat %>%
  distinct(peptide, charge) %>%
  summarize(n()) %>% pull() -> unique_peptides_for_lh_hh

smry = plot_dat %>%
  group_by(tissue, iqp_id, x, color) %>%
  summarize(.groups='keep',
            n    = n(),
            mean = mean(ria, na.rm=T),
            l95 = lower(ria),
            u95 = upper(ria)) %>% 
  ungroup()

smry -> lh_hh_smry_by_animal_and_tissue


  
par(mar=c(3,4,3,1))
xlims = c(0.3, 7.7)
ylims = c(0, 1)
ybigs = 0:4/4
yats  = 0:10/10
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=ybigs, labels=NA, tck=-0.04)
axis(side=2, at=ybigs, labels=percent(ybigs), lwd=0, las=2, line=-0.25)
axis(side=2, at=yats,  labels=NA, tck=-0.02)
mtext(side=2, line=2.5, text='RIA', cex=0.8)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, line=0.2, at=x_map$x, text=gsub(' ','\n',x_map$label), padj=0, cex=0.4)
for (this_tissue in c('brain','colon')) {
  xs = x_map$x[x_map$tissue==this_tissue]
  axis(side=1, at=c(min(xs)-0.6, max(xs)+0.6), line=1.2, tck=0.03, labels=NA)
  mtext(side=1, line=1.5, at=mean(xs), text=this_tissue, cex=0.8)
}
set.seed(1)
points(x=jitter(plot_dat$x, amount=0.15), y=plot_dat$ria,
       pch=21, col=alpha(plot_dat$color, 0.2), bg=alpha('#FFFFFF',0.2), cex=0.5)
arrows(x0=smry$x, y0=smry$l95, y1=smry$u95,
       code=3, angle=90, length=0.05, lwd=2, col=smry$color)
barwidth = 0.4
segments(x0=smry$x-barwidth, x1=smry$x+barwidth, y0=smry$mean, lwd=1.5, col=smry$color)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

bin_map = lh_hh %>%
  filter(iqp_id %in% paired_samples) %>%
  distinct(tissue, log_min_signal_binned) %>%
  arrange(tissue, log_min_signal_binned) %>%
  mutate(x     = ifelse(tissue=='brain', row_number(), row_number() + 1),
         color = tissue_colors[tissue],
         label = sprintf('%.2f', log_min_signal_binned))

plot_dat2 = lh_hh %>%
  filter(iqp_id %in% paired_samples) %>%
  mutate(ria = 2/(2 + (light_heavy / heavy_heavy))) %>%
  inner_join(bin_map, by=c('tissue','log_min_signal_binned'))

smry2 = plot_dat2 %>%
  group_by(tissue, log_min_signal_binned, x, color) %>%
  summarize(.groups='keep',
            n_peptides    = n(),
            mean = mean(ria, na.rm=T),
            l95 = lower(ria),
            u95 = upper(ria)) %>%
  ungroup() 

smry2 -> lh_hh_smry_by_bin_and_tissue
write_supp_table(lh_hh_smry_by_bin_and_tissue, 'Summary of free lysine calculated from LH/HH ratios by intensity bin and tissue.')

plot_dat2 %>%
  distinct(peptide, charge) %>%
  summarize(n()) %>% pull() -> unique_peptides_for_lh_hh_by_bin


plot_dat2 %>%
  distinct(pg_genes) %>%
  summarize(n()) %>% pull() -> unique_genes

par(mar=c(3,4,3,1))
xlims = c(0.3, max(bin_map$x) + 0.7)
ylims = c(0, 1)
ybigs = 0:4/4
yats  = 0:10/10
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=ybigs, labels=NA, tck=-0.04)
axis(side=2, at=ybigs, labels=percent(ybigs), lwd=0, las=2, line=-0.25)
axis(side=2, at=yats,  labels=NA, tck=-0.02)
mtext(side=2, line=2.5, text='RIA', cex=0.8)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, line=0.2, at=bin_map$x, text=paste0('bin\n',bin_map$label), padj=0, cex=0.4)
for (this_tissue in c('brain','colon')) {
  xs = bin_map$x[bin_map$tissue==this_tissue]
  axis(side=1, at=c(min(xs)-0.6, max(xs)+0.6), line=1.2, tck=0.03, labels=NA)
  mtext(side=1, line=1.5, at=mean(xs), text=this_tissue, cex=0.8)
}
set.seed(1)
points(x=jitter(plot_dat2$x, amount=0.15), y=plot_dat2$ria,
       pch=21, col=alpha(plot_dat2$color, 0.2), bg=alpha('#FFFFFF', 0.2), cex=0.5)
arrows(x0=smry2$x, y0=smry2$l95, y1=smry2$u95,
       code=3, angle=90, length=0.05, lwd=2, col=smry2$color)
barwidth = 0.4
segments(x0=smry2$x-barwidth, x1=smry2$x+barwidth, y0=smry2$mean, lwd=1.5, col=smry2$color)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

dev.off()




