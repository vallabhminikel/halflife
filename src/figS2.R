
## Figure S2 ####
tell_user('done.\nCreating Figure S2...')
resx=300
png('display_items/figure-s2.png',width=6.5*resx,height=3.0*resx,res=resx)

layout_matrix = matrix(c(1,2,2), nrow=1, byrow=T)
layout(layout_matrix)

panel = 1

### A CRL abundance ####

par(mar=c(3,4,3,1))
pep_stds = read_tsv('data/crl_pepstds.tsv', col_types=cols())
bql = pep_stds$back_calc_ng_ml[1] / 10 # TO DO - check this with Abdul Basit
tissue_meta = tibble(tissue=c('brain','colon'),
                     color = c('#CCCF00','#EEBB77'),
                     x=c(1,2))

crl_prot = read_tsv('data/crl_proteomics.tsv', col_types=cols()) %>%
  mutate(heavy = replace_na(heavy, bql)) %>%
  mutate(total = light + heavy) %>%
  mutate(proportion_heavy = heavy / total) %>%
  inner_join(tissue_meta, by ='tissue')

xlims = c(0.5, 2.5)
ylims = c(0, 8)
ybigs = 0:10
ybiglabs = ybigs
yats = 0:100/10
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, at=tissue_meta$x, text=tissue_meta$tissue, cex=0.8)
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.25)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=1.6, text='total peptide\nconcentration (ng/mL)', cex=0.8)

crl_prot %>%
  group_by(x, color) %>%
  summarize(.groups='keep',
            mean = mean(total),
            l95 = lower(total),
            u95 = upper(total)) -> smry
barwidth=0.8
rect(xleft=smry$x-barwidth/2, xright=smry$x+barwidth/2, ybottom=rep(0,nrow(smry)), ytop=smry$mean, col=alpha(smry$color,ci_alpha), lwd=1.5, border=NA)
set.seed(1)
points(jitter(crl_prot$x,amount=0.25), crl_prot$total, col=crl_prot$color, pch=21, bg='#FFFFFF')
arrows(x0=smry$x, y0=smry$l95, y1=smry$u95, code=3, angle=90, length=0.05, col='#000000', lwd=1.5)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(smry, 'Abundance of PrP peptide in Charles River proteomics data.')

### B IQ Proteomics ####

par(mar=c(3,4,3,1))
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


iqp_data %>%
  filter(labeled == 'yes') %>%
  filter(protein == 'PRNP') %>%
  filter(!(genotype=='Tg25109 het; ZH3/ZH3' & peptide=='VVEQMCVTQYQK')) %>%
  filter(above_lloq) %>%
  inner_join(prnp_pep_meta, by='peptide') %>%
  inner_join(tissue_meta, by='tissue') %>%
  inner_join(genotype_meta, by='genotype') %>%
  mutate(x = (xtiss-1)*8 + (xpep-1)*4 + (xgeno-1)) -> all_prnp_plab

all_prnp_plab %>%
  group_by(genotype, explanation, tissue, peptide, x, color, pep_color, short) %>%
  summarize(.groups='keep',
            mean = mean(total),
            l95 = lower(total),
            u95 = upper(total)) %>%
  ungroup() -> smry

smry %>%
  filter(tissue=='brain') %>%
  filter(explanation %in% c('WT','Tga20')) %>%
  select(explanation, peptide, mean) %>%
  pivot_wider(names_from = explanation, values_from=mean) %>%
  mutate(tga20_rel = `Tga20`/`WT`) -> tga20_expression_level

smry %>%
  select(explanation, tissue, peptide, mean) %>%
  pivot_wider(names_from=tissue, values_from=mean) %>%
  mutate(brain_colon_ratio = brain/colon) -> iqp_bc_ratio

xlims = range(smry$x) + c(-0.7, 0.7)
ylims = c(3e4, 4e8)
ybigs = 10^(4:8)
ybiglabs = gsub('\\+0','',formatC(10^(4:8),format='e',digits=0))
yats = rep(1:9, 5) * rep(10^(4:8),each=9)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, line=0.0, at=smry$x, text=smry$explanation, cex=0.35)
smry %>%
  group_by(short, tissue) %>%
  summarize(.groups='keep',
            minx = min(x),
            maxx = max(x),
            midx = mean(x)) %>%
  ungroup() -> tranches
overhang = 0.45
tranche_line = 1.0
for (i in 1:nrow(tranches)) {
  axis(side=1, line=tranche_line, at=c(tranches$minx[i]-overhang,tranches$maxx[i]+overhang), tck=0.015, labels=NA)
  mtext(side=1, line=tranche_line-0.4, padj=1, at=c(tranches$midx[i]), text=tranches$short[i], cex=0.6)
}
smry %>%
  group_by(tissue) %>%
  summarize(.groups='keep',
            minx = min(x),
            maxx = max(x),
            midx = mean(x)) %>%
  ungroup() -> tranches
tranche_line = 0.5
for (i in 1:nrow(tranches)) {
  axis(side=3, line=tranche_line, at=c(tranches$minx[i]-overhang,tranches$maxx[i]+overhang), tck=0.015, labels=NA)
  mtext(side=3, line=tranche_line+0.2, padj=0, at=c(tranches$midx[i]), text=gsub(' ','\n',tranches$tissue[i]), cex=0.7)
}
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=0, cex.axis=0.8)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=2.5, text='total peptide intensity', cex=0.8)
barwidth=0.9
rect(xleft=smry$x-barwidth/2, xright=smry$x+barwidth/2, ybottom=rep(1,nrow(smry)), ytop=smry$mean, col=alpha(smry$color,ci_alpha), lwd=1.5, border=NA)
set.seed(1)
points(jitter(all_prnp_plab$x,amount=0.1), all_prnp_plab$total, col=all_prnp_plab$color, pch=21, bg='#FFFFFF')
arrows(x0=smry$x, y0=pmax(smry$l95,1), y1=smry$u95, code=3, angle=90, length=0.05, col='#000000', lwd=1.5)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(smry, 'Abundance of PrP peptides in IQ Proteomics data.')

silence_is_golden = dev.off()
### end Figure S2 ####