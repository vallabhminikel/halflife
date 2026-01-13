
## Figure S1 ####
tell_user('done.\nCreating Figure S1...')
resx=300
png('display_items/figure-s1.png',width=6.5*resx,height=3.0*resx,res=resx)

layout_matrix = matrix(1:3, nrow=1, byrow=T)
layout(layout_matrix)
panel = 1

par(mar = c(3,3,3,2))

### A effect of tube size ####

rbind_files('data/','18[56]_summary.tsv') %>%
  mutate(plate = as.integer(substr(file,1,3))) %>%
  mutate(spin = grepl('spin',sample)) %>%
  mutate(sample_id = gsub('VL23-01 ','',gsub(' \\+ spin','',sample))) %>%
  filter(grepl('VL23-01',sample)) -> elisa
cohort = read_tsv('data/colon_qc_cohort.tsv',col_types=cols()) %>%
  clean_names() %>%
  mutate(sample_id = as.character(sample_id)) %>%
  mutate(tube_size = as.integer(substr(homogenization,1,1))) %>%
  mutate(genotype = gsub('PrP ','',genotype))
gt_meta = tibble(genotype = c('KO','WT'),
                 offset = c(.2, -.2),
                 color = c('#FBC74A','#00A7CD'))
tube_meta = tibble(tube_x=1:2,
                   tube_size=c(2,7))
cohort %>%
  inner_join(elisa, by=c('sample_id')) %>%
  filter(plate==185) %>%
  select(animal_id, sample_id, tube_size, genotype, ngml_av) %>%
  inner_join(gt_meta, by='genotype') %>%
  inner_join(tube_meta, by='tube_size') %>%
  mutate(x = tube_x + offset) -> main
main %>%
  group_by(x, tube_x, offset, color, tube_size, genotype) %>%
  summarize(.groups='keep',
            mean = mean(ngml_av),
            l95 = lower(ngml_av),
            u95 = upper(ngml_av)) %>%
  ungroup() -> smry

xlims = c(0.5, 2.5)
ylims = c(0, 52)
yats = 0:60
ybigs = 0:6*10
ybiglabs = ybigs
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, at=smry$x, text=smry$genotype, cex=0.6, line=0.25)
mtext(side=1, at=tube_meta$tube_x, text=paste0(tube_meta$tube_size, ' mL'), cex=0.8, line=1.5)
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.25)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=1.6, text='PrP (ng/mL)', cex=0.8)
llq = 0.05 * 100
abline(h=llq, lty=3, col='black')
mtext(side=4, at=llq, text='LLQ', cex=0.6, line=0.1, las=2)
barwidth=0.4
rect(xleft=smry$x-barwidth/2, xright=smry$x+barwidth/2, ybottom=rep(0,nrow(smry)), ytop=smry$mean, col=alpha(smry$color,ci_alpha), lwd=1.5, border=NA)
set.seed(1)
points(jitter(main$x,amount=0.15), main$ngml_av, col=main$color, pch=21, bg='#FFFFFF')
arrows(x0=smry$x, y0=smry$l95, y1=smry$u95, code=3, angle=90, length=0.025, col='#000000', lwd=1)
mtext(side=3, adj=0.0, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(smry, 'Effect of homogenization tube size on colon PrP quantification.')

### B effect of spin ####

gt_meta = tibble(genotype = c('KO','WT'),
                 gt_x = c(0, 5),
                 color = c('#FBC74A','#00A7CD'))
tube_meta = tibble(tube_x=c(-1,1),
                   tube_size=c(2,7))
spin_meta = tibble(spin_x = c(-0.5, 0.5),
                   spin = c(F,T),
                   spin_disp = c('',' + spin'))
cohort %>%
  inner_join(elisa, by=c('sample_id')) %>%
  filter(plate==186) %>%
  select(animal_id, sample_id, tube_size, genotype, ngml_av, spin) %>%
  inner_join(gt_meta, by='genotype') %>%
  inner_join(tube_meta, by='tube_size') %>%
  inner_join(spin_meta, by='spin') %>%
  mutate(x = gt_x + tube_x + spin_x) %>%
  mutate(disp = paste0(tube_size, ' mL',spin_disp)) -> main
main %>%
  group_by(x, tube_x, gt_x, spin_x, disp, color, tube_size, genotype) %>%
  summarize(.groups='keep',
            mean = mean(ngml_av),
            l95 = lower(ngml_av),
            u95 = upper(ngml_av)) %>%
  ungroup() -> smry


xlims = range(smry$x) + c(-0.5, 0.5)
ylims = c(0, 52)
yats = 0:60
ybigs = 0:6*10
ybiglabs = ybigs
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
smry %>%
  group_by(gt_x, genotype) %>%
  summarize(.groups='keep',
            minx = min(x),
            maxx = max(x),
            midx = mean(x)) %>%
  ungroup() -> tranches
tranche_line = 0.5
overhang = 0.4
for (i in 1:nrow(tranches)) {
  axis(side=3, line=tranche_line, at=c(tranches$minx[i]-overhang,tranches$maxx[i]+overhang), tck=0.015, labels=NA)
  mtext(side=3, line=tranche_line+0.2, padj=0, at=c(tranches$midx[i]), text=gsub(' ','\n',tranches$genotype[i]), cex=0.7)
}
par(xpd=T)
text(x=smry$x, y=0, adj=1, srt=45, labels=paste0(smry$disp,'  '), cex=0.8)
par(xpd=F)
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.25)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=1.6, text='PrP (ng/mL)', cex=0.8)
llq = 0.05 * 100
abline(h=llq, lty=3, col='black')
mtext(side=4, at=llq, text='LLQ', cex=0.6, line=0.1, las=2)
barwidth=0.4
rect(xleft=smry$x-barwidth/2, xright=smry$x+barwidth/2, ybottom=rep(0,nrow(smry)), ytop=smry$mean, col=alpha(smry$color,ci_alpha), lwd=1.5, border=NA)
set.seed(1)
points(jitter(main$x,amount=0.15), main$ngml_av, col=main$color, pch=21, bg='#FFFFFF')
suppressWarnings(arrows(x0=smry$x, y0=smry$l95, y1=smry$u95, code=3, angle=90, length=0.025, col='#000000', lwd=1))
mtext(side=3, adj=0.0, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(smry, 'Effect of centrifugation on colon PrP quantification.')

### C stability study ####
stab_raw = read_tsv('data/223_summary.tsv', col_types=cols()) %>%
  filter(grepl('colon',sample))
smry = tribble(
  ~x, ~name1, ~name2, ~disp,
  1, '', 'normal', 'normal',
  3, 'thaw twice', 'thaw twice', 'thawed twice',
  2, 'thaw once', 'thaw once', 'thawed once',
  4, '4C O/N', '4°C', '4°C O/N',
  5, 'RT O/N', 'RT', 'RT O/N'
)
samps = tribble(
  ~color, ~qclevel, ~qcdisp,
  '#FF0000', 'Hi', 'High QC (WT)',
  '#FF8800', 'Mid', 'Mid QC (Het KO)',
  '#00FF00', 'Lo', 'Low QC (10% WT, 90% hom KO)',
  '#0000FF', 'Neg', 'Neg QC (Hom KO)'
)
stab_raw %>%
  mutate(name1 = trimws(gsub('Mo Pos (Hi|Mid|Lo|Neg) QC colon','',sample))) %>%
  inner_join(smry, by='name1') %>%
  mutate(qclevel = gsub(' .*','',gsub('Mo Pos ','',sample))) %>%
  inner_join(samps, by='qclevel') %>%
  select(x, disp, qclevel, qcdisp, color, ngml_av) -> stab


xlims = range(stab$x) + c(-0.5, 0.5)
ylims = c(0, 20)
yats = 0:60
ybigs = 0:6*10
ybiglabs = ybigs
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
par(xpd=T)
text(x=smry$x, y=0, adj=1, srt=45, labels=paste0(smry$disp, '  '), cex=0.7)
par(xpd=F)
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.25)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=1.6, text='PrP (ng/mL)', cex=0.8)
llq = 0.02048 * 100
abline(h=llq, lty=3, col='black')
mtext(side=4, at=llq, text='LLQ', cex=0.6, line=0.1, las=2)
points(stab$x, stab$ngml_av, col=stab$color, pch=20)
parxpdt(legend(x=2.5,y=23,samps$qcdisp, col=samps$color, pch=20, cex=0.6))
mtext(side=3, adj=0.0, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(stab, 'Stability of colon PrP in ELISA.')

silence_is_golden = dev.off() ### end Figure S1 ####
