
## Figure S5 Ki817 mouse ####
tell_user('done.\nCreating Figure S5...')
resx=300
png('display_items/figure-s5.png',width=resx*6.5,height=5*resx, res=resx)

layout_matrix = matrix(1:4, nrow=4, byrow=T)
layout(layout_matrix, heights=c(0.8, 0.2, 0.8, 0.2))

grch38_coding_offset = 4699605 - 4680251
depth = read_tsv('data/TACONICM10_1_gene_depth.tsv',col_types=cols()) %>%
  mutate(pos = pos + grch38_coding_offset)
all_possible_pos = tibble(pos=seq(4570000,4740000),by=1)

depth %>%
  select(-chrom) %>%
  right_join(all_possible_pos, by='pos') %>%
  mutate(depth=replace_na(depth,0)) %>%
  arrange(pos) -> depthtbl

depthtbl %>%
  mutate(pos100 = floor(pos/100)*100) %>%
  group_by(pos100) %>%
  summarize(.groups='keep',
            p30 = mean(depth >= 30),
            p10 = mean(depth >= 10),
            mn = mean(depth)) %>%
  ungroup() -> depth100

write_supp_table(depth100, 'Sequencing depth of ki817 mouse against GRCh38 human genome reference.')

xlims = c(4570000, 4740000)

zooms = tibble(zoomed = c('out','in'),
               xmin = c(4570000,4686456-1000),
               xmax = c(4740000, 4701588+1000))

for (i in 1:nrow(zooms)) {
  xlims = c(zooms$xmin[i], zooms$xmax[i])
  
  
  pseudozero = 1
  ylims = c(pseudozero,1e5)
  ybigs = c(1, 10, 100, 1000, 10000)
  yats = rep(1:9, 6) * 10^(rep(-1:4, each=9))
  ybiglabs = c('≤1', '10', '100', '1K','10K')
  xbigs = seq(min(xlims), max(xlims), 10000)
  xats = seq(min(xlims), max(xlims), 1000)
  
  par(mar=c(1,4,1,2))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
  axis(side=1, at=xbigs, labels=NA,tck=-0.05)
  axis(side=1, at=xats, labels=NA,tck=-0.02)
  axis(side=1, at=xbigs, line=-0.5, labels=paste0(formatC(xbigs/1e6,digits=2,format='f'),'M'), lwd=0)
  axis(side=2, at=ybigs, labels=NA,tck=-0.04)
  axis(side=2, at=yats, labels=NA,tck=-0.02)
  axis(side=2, at=ybigs, line=-0.25, las=2, labels=ybiglabs, lwd=0)
  mtext(side=2, line=2.5, text='sequencing depth', cex=0.8)
  mtext(side=1, line=1.6, adj=0, text='chr20 position', cex=0.8)
  polygon(c(depth100$pos100,max(depth100$pos100),rev(depth100$pos100)), c(pmax(depth100$mn,pseudozero),pseudozero,rep(0,nrow(depth100))), lwd=3, col='#CEAB1277', border='#CEAB12')
  
  
  exon1_start = c(4686456,4721909 )# transcription start site
  exon1_end = c(4686512,4721969)
  exon2_start =c(4699211,4724541)
  exon2_end = c(4701588,4728460) # transcription end site
  cds_start = c(4699221,4724552)
  cds_end = c(4699982,4725079)
  
  landmarks = tibble(gene = c('PRNP','PRND'),
                     exon1_start = exon1_start,
                     exon2_start = exon2_start,
                     exon1_end   = exon1_end,
                     exon2_end   = exon2_end,
                     cds_start   = cds_start,
                     cds_end     = cds_end,
                     exon2_fill  = '#A3A3A3',
                     y = 2)
  
  intron_lwd = 1
  utr_lwd = 10
  cds_lwd = 20
  
  default_fill = '#000000'
  utr_height = 2
  cds_height = 4
  prnd_height = 1
  
  par(mar=c(0,4,0,2))
  ylims = c(-2, 6)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  segments(x0=landmarks$exon1_start, x1=landmarks$exon2_end, y0=landmarks$y, lwd=intron_lwd, lend=1)
  rect(xleft=landmarks$exon1_start, xright=landmarks$exon1_end, ybottom=landmarks$y-utr_height/2, ytop=landmarks$y+utr_height/2, col=default_fill, border=NA)
  rect(xleft=landmarks$exon2_start, xright=landmarks$exon2_end, ybottom=landmarks$y-utr_height/2, ytop=landmarks$y+utr_height/2, col=default_fill, border=NA)
  rect(xleft=landmarks$cds_start, xright=landmarks$cds_end, ybottom=landmarks$y-cds_height/2, ytop=landmarks$y+cds_height/2, col=default_fill, border=NA)
  text(x=(landmarks$exon2_end + landmarks$exon1_start)/2, y=landmarks$y-1.4, pos=1, labels=landmarks$gene, font=3, cex=0.8)
  #text(x=(landmarks$exon2_end + landmarks$exon1_start)/2, y=landmarks$y, pos=1, labels=landmarks$gene, font=3, cex=0.8)
  
}


silence_is_golden = dev.off() ### end Fig S5 humanized mouse #### 
