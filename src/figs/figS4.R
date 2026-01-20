
## Figure S4 PK for ASO N ####
tell_user('done.\nCreating Figure S4...')
resx=300
png('display_items/figure-s4.png',width=resx*3.25,height=3.5*resx, res=resx)
par(mar=c(4,4,3,1))
xlims = c(0, 45)
ylims = c(0.7, 10)
yats = rep(1:9, 4) * rep(10^(-1:2), each=9)
ybigs = 10^(-1:2)
plot(NA, NA, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
axis(side = 1, at = c(0, 7, 14, 21, 28, 42, 84))
mtext(side=1, line=2.5, text='days post-dose')
points(ason$day, ason$conc_ug_g, col='#FF56AD', pch=19)
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, line=-0.25, lwd=0, labels=ybigs, las=2)
mtext(side=2, line=2.0, text='drug concentration (µg/g)')
barwidth=1
segments(x0=ason_smry$day-barwidth, x1=ason_smry$day+barwidth, y0=ason_smry$pk_mean, col=ason_smry$dose_color)
arrows(x0=ason_smry$day, y0=ason_smry$pk_l95, y1=ason_smry$pk_u95, col=ason_smry$dose_color, code=3, angle=90, length=0.02)

silence_is_golden = dev.off() ### end Figure S4 ####
