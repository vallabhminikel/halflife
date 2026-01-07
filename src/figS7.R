

## Figure S7 #### 
tell_user('done.\nCreating Figure S7...')
resx=300
png('display_items/figure-s7.png',width=6.5*resx,height=3*resx,res=resx)

layout_matrix = matrix(1:4, byrow=T, nrow=1)
layout(layout_matrix)

par(mar = c(6,4,2,1))
panel = 1

ffi_all %>%
  filter(protein=='PRNP' & peptide %in% use_peptides) %>%
  mutate(pepnick = substr(peptide, 1, 4)) %>%
  inner_join(leg, by='genotype') -> abun_data

for (this_age in c('young','aged')) {
  for (this_peptide in use_peptides) {
    abun_data %>%
      filter(age==this_age & peptide==this_peptide) -> subs
    xlims = c(0.5, 3.5)
    if (this_peptide == 'VVEQMCVTQYQK') {
      ylims = c(0, 2e6)
      ybigs = 0:2*1e6
      ybiglabs = gsub('\\+0','',formatC(ybigs, format='e', digits=0))
      yats = 0:20*1e5
    } else if (this_peptide == 'GENFTETDVK') {
      ylims = c(0, 4e7)
      ybigs = 0:4*1e7
      ybiglabs = gsub('\\+0','',formatC(ybigs, format='e', digits=0))
      yats = 0:40*1e6
    }
    plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
    mtext(side=3, text=paste0(this_age,' ',substr(this_peptide,1,4)), cex=0.8)
    axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
    mtext(side=1, at=leg$xgeno, text=leg$disp, col=leg$color, line=0.25, cex=0.6, las=2)
    axis(side=2, at=yats, labels=NA, tck=-0.02)
    axis(side=2, at=ybigs, labels=NA, tck=-0.05)
    axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.4)
    mtext(side=2, line=2.2, text='total intensity', cex=0.8)
    subs %>%
      group_by(xgeno, genotype, color) %>%
      summarize(.groups='keep',
                n = n(),
                mean_total = mean(total),
                l95_total = lower(total),
                u95_total = upper(total)) %>%
      ungroup() -> smry
    barwidth = 0.4
    rect(xleft=smry$xgeno-barwidth, xright=smry$xgeno+barwidth, ybottom=rep(0, nrow(smry)), ytop=smry$mean_total, col=alpha(smry$color, ci_alpha), border=NA)
    arrows(x0=smry$xgeno, y0=smry$l95_total, y1=smry$u95_total, col=smry$color, code=3, angle=90, length=0.05)
    set.seed(1)
    points(jitter(subs$xgeno,amount=.25), subs$total, pch=21, bg='#FFFFFF', col=subs$color)
    mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1
  }  
}
silence_is_golden = dev.off()
### end Figure S7 #### 

