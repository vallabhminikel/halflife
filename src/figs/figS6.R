
## Figure S6 Tga20 expression level #### 
tell_user('done.\nCreating Figure S6...')
resx=300
png('display_items/figure-s6.png',width=resx*3.25,height=3.5*resx, res=resx)

par(mar=c(1,3,1,1))

rbind_files('data/','21[46]_summary.tsv') %>%
  mutate(plate = as.integer(substr(file,1,3))) -> elisa
cohort = read_tsv('data/tg_expression_cohort.tsv',col_types=cols()) %>%
  mutate(animal = as.character(animal))
meta = tibble(genotype = c('WT','Tga20'),
              x = c(1,2),
              color = c('#545454','#78AB46'))
cohort %>%
  inner_join(elisa, by=c('animal'='sample', 'plate')) %>%
  select(animal, plate, genotype, ngml_av) %>%
  inner_join(meta, by='genotype') %>%
  group_by(plate) %>%
  mutate(rel = ngml_av / mean(ngml_av[genotype=='WT'])) %>%
  select(-ngml_av) -> tgexp
tgexp %>%
  group_by(x, color, genotype) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            l95 = lower(rel),
            u95 = upper(rel)) %>%
  ungroup() -> smry

write_supp_table(smry, 'PrP expression in Tga20 vs. WT mice.')

xlims = c(0.5, 2.5)
ylims = c(0, 3)
ybigs = 0:10
ybiglabs = ybigs
yats = 0:100/10
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
mtext(side=1, at=meta$x, text=meta$genotype, cex=0.8)
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.25)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=1.6, text='PrP (fold WT)', cex=0.8)
abline(h=1, lty=3)
barwidth=0.8
rect(xleft=smry$x-barwidth/2, xright=smry$x+barwidth/2, ybottom=rep(0,nrow(smry)), ytop=smry$mean, col=alpha(smry$color,ci_alpha), lwd=1.5, border=NA)
set.seed(1)
points(jitter(tgexp$x,amount=0.25), tgexp$rel, col=tgexp$color, pch=21, bg='#FFFFFF')
arrows(x0=smry$x, y0=smry$l95, y1=smry$u95, code=3, angle=90, length=0.05, col='#000000', lwd=1.5)


silence_is_golden = dev.off() ### end Fig S6 Tga20 #### 
