tell_user('done.\nCreating Figure 2...')
resx=300
png('display_items/figure-2.png',width=6.5*resx,height=4.5*resx,res=resx)

layout_matrix = matrix(c(1,1,1,2,2,2,3,3,3,
                         4,4,4,4,4,5,5,6,6), nrow=2, byrow=T)
layout(layout_matrix)

panel = 1

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
  mutate(estimated_thalf = find_thalf(mean),
         estimated_thalf_if_100 = find_thalf(mean, avails=always_all)) -> smry

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

write_supp_table(smry, 'PrP peptide labeling by genotype in IQ Proteomics data.')

prnp_iqp_smry = smry # save to use in next panel

par(mar=c(3,4,3,1))
xlims = range(smry$x) + c(-0.7, 0.7)
ylims = c(0, 0.7)
ybigs = 0:10/10
ybiglabs = percent(ybigs)
yats = 0:20/20
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, line=-0.1, at=smry$x, text=smry$short, cex=0.35)
smry %>%
  group_by(genotype, explanation, tissue) %>%
  summarize(.groups='keep',
            minx = min(x),
            maxx = max(x),
            midx = mean(x)) %>%
  ungroup() -> tranches
overhang = 0.5
tranche_line = 1.0
for (i in 1:nrow(tranches)) {
  axis(side=1, line=tranche_line, at=c(tranches$minx[i]-overhang,tranches$maxx[i]+overhang), tck=0.03, labels=NA)
  mtext(side=1, line=tranche_line-0.4, padj=1, at=c(tranches$midx[i]), text=tranches$explanation[i], cex=0.6)
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




### E/F theoretical vs. diff prots ####

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
  
  thalf_values = seq(0, 12, .1)
  points(thalf_values, proportion_labeled_t(thalf_values, t=t, which_t = 8, avails=free_lysine), type = 'l', col = lysine_meta$color[2], lwd = 2)
  if (this_tissue=='colon') {
    points(thalf_values, proportion_labeled_t(thalf_values, t=t, which_t = 8, avails=always_all), type = 'l', col = lysine_meta$color[1], lwd = 2)
  }
  
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
  
  prnp_iqp_smry %>%
    filter(tissue==this_tissue) %>%
    filter(genotype=='C57BL/6N') %>%
    arrange(desc(mean)) -> subs
  prnp_color = '#0001CD'
  
  if (this_tissue == 'brain') {
    avails = free_lysine
    label_x = 1.5
  } else if (this_tissue == 'colon') {
    avails = always_all
    label_x = 1.5
  }
  segments(x0=rep(0, nrow(subs)), x1=find_thalf(subs$mean,avails=avails), y0=subs$mean, col=prnp_color, lty=3)
  text(x=rep(label_x,2), y=subs$mean+c(-0.02,0.02), pos=c(3,1), labels=subs$short, col="#0001CD", cex=0.5)
  segments(x0=find_thalf(subs$mean,avails=avails), y0=rep(0,nrow(subs)), y1=subs$mean, col=prnp_color, lty=3)
  
  legend("topright", legend = lloq_meta$disp, col =lloq_meta$color, pch = 1, cex = 0.5, bty='n')
  
  mtext(side=3, adj=0, text=LETTERS[panel], line=0.5); panel = panel + 1
}

silence_is_golden = dev.off()
### end Figure 2 ####