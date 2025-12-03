# Startup #### 

overall_start_time = Sys.time()
tell_user = function(...) { cat(file=stderr(), paste0(...)); flush.console() }

# Dependencies ####

tell_user('Loading required packages...')

suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(openxlsx))
suppressMessages(library(smoother))
suppressMessages(library(plotrix))
suppressMessages(library(magick))
suppressMessages(library(minpack.lm))

if (interactive()) setwd('~/d/sci/src/halflife')


# Output streams #### 

tell_user('done.\nCreating output streams...')

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T
write_stats = function(...) {
  write(paste(list(...),collapse='',sep=''),text_stats_path,append=T)
  write('\n',text_stats_path,append=T)
}

supplement_path = 'display_items/supplement.xlsx'
supplement = createWorkbook()
# options("openxlsx.numFmt" = "0.00") # this looks better for residuals but terrible for p values and weeks post-dose
supplement_directory = tibble(name=character(0), title=character(0))
write_supp_table = function(tbl, title='') {
  # write Excel sheet for supplement
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=2,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,tbl,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  
  # also write tab-sep version for GitHub repo
  write_tsv(tbl,paste0('display_items/table-',table_name,'.tsv'), na='')
  
  # and save the title in the directory tibble for later
  assign('supplement_directory',
         supplement_directory %>% add_row(name=table_name, title=title),
         envir = .GlobalEnv)
}


tell_user('done.\nDefining constants and functions...')

# Constants #### 

dt = 0.01
t = seq(0, 8, dt)

# Functions ####

## General functions ####
upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x, na.rm=T) + sds*sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
}

lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x, na.rm=T) - sds*sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
}

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x > 0 & signed, '+', ''),formatC(100*x,format='f',digits=digits),'%'))

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}
ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

parxpdt = function(expr) {
  par(xpd=T)
  expr
  par(xpd=F)
}

rbind_files = function(path, grepstring) {
  all_files = list.files(path, full.names=T)
  these_files = all_files[grepl(grepstring,all_files)]
  if (exists('rbound_table')) rm('rbound_table')
  for (this_file in these_files) {
    this_tbl = read_delim(this_file, col_types=cols()) %>% clean_names()
    this_tbl$file = gsub('.*\\/','',gsub('\\.[tc]sv','',this_file))
    if (exists('rbound_table')) {
      rbound_table = rbind(rbound_table, this_tbl)
    } else {
      rbound_table = this_tbl
    }
  }
  return (rbound_table)
}


## Functions for Fig 2 ####

# formula from Fornasiero 2018 Figure S3H ####
free_lysine = function(t) {
  pmax(0,1-0.503*(exp(-t*log(2)*0.799))-0.503*exp(-t*log(2)/39.423))
}

always_all = function(t) {
  1
}

# avails is a function telling you what lysine is available -
# either free_lysine or always_all
proportion_labeled = function(thalf, t, avails) {
  lambda = log(2) / thalf
  proportion_heavy = numeric(length(t))
  proportion_heavy[1] = 0
  for (i in 2:length(t)) {
    protein_turned_over = lambda * dt
    original = (1 - protein_turned_over) * proportion_heavy[i - 1]
    nascent = protein_turned_over * avails(t[i])
    proportion_heavy[i] = original + nascent
  }
  return (proportion_heavy)
}


proportion_labeled_t_unary = function(thalf, t, which_t, avails) {
  return (proportion_labeled(thalf, t, avails)[t==which_t])
}

proportion_labeled_t = function(thalf, t, which_t, avails) {
  result = numeric(length(thalf))
  for (i in 1:length(result)) {
    result[i] = proportion_labeled(thalf[i], t, avails)[t==which_t]
  }
  return (result)
}


function_for_uniroot = function(thalf, prop_labeled, t, which_t, avails) {
  proportion_labeled_t(thalf, t, which_t, avails) - prop_labeled
}
find_thalf_unary = function(prop_labeled, t=seq(0,8,.01), which_t=8, avails=free_lysine) {
  max_possible = proportion_labeled_t(thalf=0.01, t=t, which_t=which_t, avails=avails)
  if (prop_labeled > max_possible) {
    return (as.numeric(NA))
  }
  uniroot(function_for_uniroot, lower=0.1, upper=100, tol=0.01, 
          prop_labeled = prop_labeled, t = t, which_t = which_t, avails=avails)$root
}
# vectorized version
find_thalf = function(prop_labeled, t=seq(0,8,.01), which_t=8, avails=free_lysine) {
  result = numeric(length(prop_labeled))
  for (i in 1:length(prop_labeled)) {
    result[i] = find_thalf_unary(prop_labeled[i], t, which_t, avails)
  }
  return (result)
}

## Functions for Fig 3 ####
interpolate_rna = function(day, rna, xout) {
  tibble(day, rna) %>%
    group_by(day) %>%
    summarize(.groups='keep', rna=mean(rna)) %>%
    ungroup() %>%
    add_row(day=0, rna=1, .before=1) -> actual_data_points
  approx(x=actual_data_points$day, y=actual_data_points$rna, xout=xout)$y
}

Pt = function(R, t, lambda) {
  P = numeric(length(t))
  dP = numeric(length(t))
  P[1] = 1
  for (i in 2:length(t)) {
    dt = (t[i]- t[i-1])
    dP[i] = lambda * dt * R[i-1] - lambda * dt * P[i-1]
    P[i] = P[i-1] + dP[i]
  }
  return (P)
}

calculate_residuals = function(par, data, dt=0.01) {
  lambda = par[['lambda']]
  t = seq(0,max(data$day),dt)
  R = interpolate_rna(data$day, data$rna, t)
  P_pred = Pt(R, t, lambda)[match(data$day, t)]
  residuals = data$protein - P_pred
  return (residuals)
}


# Data ####

tell_user('done.\nReading in data...')

## IQ Proteomics ####
pivot_tab  = read_tsv("data/iqp/combined_pivot_tab.tsv", col_types=cols())
iqp_meta   = read_tsv("data/iqp/meta.tsv", col_types=cols())
lloq_stats = read_tsv('data/iqp/lloq_stats.tsv', col_types=cols())
lit_half   = read_tsv('data/iqp/forna_halflife.tsv', col_types=cols()) 
name_map   = read_tsv('data/iqp/name_map.tsv', col_types=cols())


# Figures ####

## Figure 1 #### 
tell_user('done.\nCreating Figure 1...')
resx=300
png('display_items/figure-1.png',width=6.5*resx,height=5.0*resx,res=resx)

layout_matrix = matrix(c(1,2,2,
                         1,2,2,
                         3,3,4), nrow=3, byrow=T)
layout(layout_matrix)

panel = 1


### A GTEx #### 
tpm = read.table('data/gtex/prnp_tpm_t.txt',header=F,skip=2)
colnames(tpm) = c('sampid','tpm')
tpm = as_tibble(tpm)

samp = read.table('data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',sep='\t',header=T,quote='',comment.char='')
colnames(samp) = gsub('[^a-z0-9_]','_',tolower(colnames(samp)))
samp = as_tibble(samp)

meta = read_tsv('data/gtex/tissue_metatissue.tsv', col_types=cols())

# brain color is too light, try fixing it
meta %>%
  mutate(color = case_when(color == '#EEEF4D' ~ '#CCCF00',
                           TRUE ~ color)) -> meta


tpm %>%
  inner_join(samp, by='sampid') %>%
  select(smts, smtsd, tpm) %>%
  mutate(tpm = as.numeric(tpm)) -> tpm_by_tissue

tpm_by_tissue %>%
  group_by(smts, smtsd) %>%
  summarize(.groups='keep', 
            n_samples = n(),
            median_tpm = median(tpm, na.rm=T)) %>%
  ungroup() %>%
  arrange(desc(median_tpm)) %>% 
  inner_join(meta, by=c('smtsd'='dispname')) %>%
  select(metatissue, smtsd, n_samples, median_tpm, color) -> tpm_by_tissue_smry

tpm_by_tissue_smry %>%
  group_by(metatissue) %>%
  summarize(.groups='keep', 
            n_tissues = n(),
            color = first(color),
            median_median_tpm = median(median_tpm, na.rm=T)) %>%
  ungroup() %>%
  arrange(desc(median_median_tpm)) %>%
  mutate(x = row_number()) %>%
  mutate(y = max(x) - x + 1) -> tpm_by_metatissue_smry

tpm_by_metatissue_smry %>%
  distinct(x, y, metatissue, color) -> xykey

tpm_by_tissue_smry %>%
  select(-color) %>%
  inner_join(xykey, by=c('metatissue')) -> tpm_by_tissue_smry 

par(mar=c(3,6,3,1))
xlims = c(0, max(tpm_by_tissue_smry$median_tpm)*1.1)
ylims = range(xykey$y) + c(-0.5, 0.5)
xbigs = 0:10*100
xats = 0:100*10
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, lwd=0, line=-0.5)
axis(side=1, at=xats, labels=NA, tck=-0.02)
mtext(side=1, line=1.6, text='Median TPM (GTEx v8)', cex=0.8)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, line=0.25, las=2, at=xykey$y, text=xykey$metatissue, cex=0.65)
barwidth=0.8
rect(xleft=rep(0,nrow(tpm_by_metatissue_smry)), xright=tpm_by_metatissue_smry$median_median_tpm, 
     ybottom = tpm_by_metatissue_smry$y - barwidth/2, ytop = tpm_by_metatissue_smry$y + barwidth/2, 
     col=alpha(tpm_by_metatissue_smry$color, 0.7), border=NA)
points(x=tpm_by_tissue_smry$median_tpm, y=tpm_by_tissue_smry$y, pch=21, col=tpm_by_tissue_smry$color, bg='#FFFFFF', lwd=1.5)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(tpm_by_metatissue_smry, 'GTEx median PRNP TPM by tissue.')

### B Western #### 
par(mar=c(0.5,0,3,0.5))
western_panel = image_convert(image_read('data/fig1b.png'),'png')
plot(as.raster(western_panel))
mtext(side=3, adj=-0.0, text=LETTERS[panel], line=0.5); panel = panel + 1

### C large 1:100 ELISA screen #### 

elisa_raw = read_tsv('data/058.tsv', col_types=cols())
elisa_meta = tibble(genotype = c('KO','WT'),
                    offset = c(.2, -.2),
                    color = c('#FBC74A','#00A7CD')) %>% arrange(desc(genotype))
use_dilution = 100
llq = elisa_raw %>% filter(dilution==use_dilution, flag=='LLQ') %>% slice(1) %>% pull(ngml_trunc)
elisa_raw %>%
  filter(dilution==use_dilution, !grepl('QC',detail)) %>%
  mutate(genotype = substr(detail, 1, 2),
         tissue = substr(detail, 4, 20)) %>% 
  inner_join(elisa_meta, by='genotype') %>%
  mutate(ngml_use = pmax(llq, case_when(flag=='LLQ' ~ ngml_trunc,
                                        TRUE ~ ngml))) %>%
  group_by(genotype, tissue, offset, color) %>%
  summarize(.groups='keep', 
            ngml_av = mean(ngml_use)) %>%
  ungroup() -> elisa_ready

elisa_ready %>%
  filter(genotype=="WT") %>%
  arrange(desc(ngml_av)) %>%
  mutate(x = row_number()) %>%
  mutate(y = max(x) - x + 1) %>%
  distinct(x, y, tissue) -> tissue_meta

elisa_ready %>%
  inner_join(tissue_meta, by='tissue') -> elisa_ready


par(mar=c(3,3,3,3))
ylims = c(0, max(elisa_ready$ngml_av)*1.1)
xlims = range(tissue_meta$x) + c(-0.5, 0.5)
ybigs = 0:10*10
yats = 0:100
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, lwd=0, las=2, line=-0.5)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=1.6, text='PrP (ng/mL)', cex=0.8)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
par(xpd=T)
text(x=tissue_meta$x, y=-2, adj=1, srt=45, labels=tolower(tissue_meta$tissue), cex=0.8)
par(xpd=F)
abline(h=llq, lty=3)
mtext(side=4, at=llq, las=2, text='LLQ', cex=0.8)
barwidth = 0.4
rect(xleft=elisa_ready$x + elisa_ready$offset - barwidth/2, 
     xright=elisa_ready$x + elisa_ready$offset + barwidth/2,
     ybottom = rep(0, nrow(elisa_ready)),
     ytop = elisa_ready$ngml_av, col=elisa_ready$color, border=NA)
legend('topright', elisa_meta$genotype, pch=15, col=elisa_meta$color, bty='n', cex=0.8)
mtext(side=3, adj=-0.0, text=LETTERS[panel], line=0.5); panel = panel + 1





### D refined 1:25 ELISA screen #### 
elisa_meta = tibble(animal = c('87488.1','91831.1'),
                    genotype = c('KO','WT'),
                    offset = c(.2, -.2),
                    color = c('#FBC74A','#00A7CD')) %>% arrange(desc(genotype))
alt_colors = c('#780909','#0001CD')
tissue_meta = tibble(tissue = c('Colon','Uterus','Heart','Spleen','Quad'),
                     x = 1:5) %>%
  mutate(y = max(x) - x + 1)
elisa_raw = read_tsv('data/061.tsv', col_types=cols())
use_dilution = 25
llq = elisa_raw %>% filter(dilution==use_dilution, flag=='LLQ') %>% slice(1) %>% pull(ngml_trunc)
elisa_raw %>%
  filter(dilution==use_dilution) %>%
  mutate(animal = substr(detail, 1, 7),
         tissue = substr(detail, 9, 20)) %>% 
  inner_join(elisa_meta, by='animal') %>%
  mutate(ngml_use = pmax(llq, case_when(flag=='LLQ' ~ ngml_trunc,
                                        TRUE ~ ngml))) %>%
  group_by(genotype, tissue, offset, color) %>%
  summarize(.groups='keep', 
            ngml_av = mean(ngml_use)) %>%
  inner_join(tissue_meta, by='tissue') -> elisa_ready


par(mar=c(3,3,3,3))
ylims = c(0, max(elisa_ready$ngml_av)*1.1)
xlims = range(tissue_meta$x) + c(-0.5, 0.5)
ybigs = 0:10
yats = 0:100/10
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=2, at=ybigs, labels=NA, tck=-0.05)
axis(side=2, at=ybigs, lwd=0, las=2, line=-0.5)
axis(side=2, at=yats, labels=NA, tck=-0.02)
mtext(side=2, line=1.6, text='PrP (ng/mL)', cex=0.8)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
par(xpd=T)
text(x=tissue_meta$x, y=-0.25, adj=1, srt=45, labels=tolower(tissue_meta$tissue), cex=0.8)
par(xpd=F)
abline(h=llq, lty=3)
mtext(side=4, at=llq, las=2, text='LLQ', cex=0.8)
barwidth = 0.4
rect(xleft=elisa_ready$x + elisa_ready$offset - barwidth/2, 
     xright=elisa_ready$x + elisa_ready$offset + barwidth/2,
     ybottom = rep(0, nrow(elisa_ready)),
     ytop = elisa_ready$ngml_av, col=elisa_ready$color, border=NA)
legend('topright', elisa_meta$genotype, pch=15, col=elisa_meta$color, bty='n', cex=0.8)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1


silence_is_golden = dev.off()
### end Figure 1 #### 



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


## Figure 2 #### 
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



## Figure 3 ####
tell_user('done.\nCreating Figure 3...')
resx=300
png('display_items/figure-3.png',width=resx*6.5,height=5.5*resx, res=resx)

layout_matrix = matrix(c(1,2,
                         3,4), nrow=2, byrow=T)
layout(layout_matrix)
panel = 1

xats = c(0,1,3,7,10,14,21,28,35,42)

### A ASO 6, naive WT mice #### 

rbind(cbind(read_tsv('data/103_summary.tsv', col_types=cols()),plate=103),
      cbind(read_tsv('data/104_summary.tsv', col_types=cols()),plate=104)) -> elisa

read_tsv('data/PRP220307.tsv', col_types=cols()) %>%
  clean_names() %>%
  mutate(animal = gsub('-','',subject_id)) %>%
  mutate(day = necropsy_study_day_computed) %>%
  mutate(tx = case_when(test_article_from_group=='1087171' ~ 'ASO6',
                        TRUE ~ 'none')) %>%
  mutate(rna = m_prnp_percent_pbs_from_hemisphere/100) %>%
  select(animal, day, tx, rna) -> meta

elisa %>%
  select(animal=sample, ngml_av, plate) %>%
  inner_join(meta, by='animal') %>%
  group_by(plate) %>%
  mutate(ctl_mean = mean(ngml_av[tx=='none'])) %>%
  ungroup() %>%
  mutate(rel = ngml_av/ctl_mean) -> hl


par(mar=c(4,4,3,1))
xlims = c(0, 45)
ylims = c(0, 1.25)
rna_color = '#00FB31'
prot_color = '#0001CD'
ctl_color = '#A9A9A9'
plot(NA, NA, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side = 1, at = xats, tck=-0.025, labels=NA)
axis(side = 1, at = xats[xats %in% hl$day], lwd=0, cex.axis=0.8, line=-0.5)
mtext(side=1, line=1.6, text='days post-dose', cex=0.8)
axis(side = 2, at= 0:5/4, labels=NA)
axis(side = 2, at= 0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='residual', cex=0.8)
abline(h=1, lty=3)
points(hl$day, hl$rel, col=alpha(ifelse(hl$tx=='none',ctl_color,prot_color),ci_alpha), pch=19)
points(hl$day, hl$rna, col=alpha(ifelse(hl$tx=='none',ctl_color,rna_color),ci_alpha), pch=19)
hl %>%
  filter(tx == 'ASO6') %>%
  group_by(tx, day) %>%
  summarize(.groups='keep',
            rna_mean = mean(rna),
            rna_l95 = lower(rna),
            rna_u95 = upper(rna),
            prp_mean = mean(rel),
            prp_l95 = lower(rel),
            prp_u95 = upper(rel)) %>%
  ungroup() -> smry
barwidth=0.5
segments(x0=smry$day-barwidth, x1=smry$day+barwidth, y0=smry$prp_mean, col=prot_color)
segments(x0=smry$day-barwidth, x1=smry$day+barwidth, y0=smry$rna_mean, col=rna_color)
arrows(x0=smry$day, y0=smry$prp_l95, y1=smry$prp_u95, col=prot_color, code=3, angle=90, length=0.02)
arrows(x0=smry$day, y0=smry$rna_l95, y1=smry$rna_u95, col=rna_color, code=3, angle=90, length=0.02)
dt = 0.01
t = seq(min(hl$day),max(hl$day),dt)
model_data = as.list(hl %>% filter(tx=='ASO6') %>% select(rna, protein=rel, day))
nlsfit = nls.lm(par=c(lambda=log(2)/5), fn=calculate_residuals, data=model_data, dt=0.01)
fit_lambda = as.numeric(nlsfit$par['lambda'])
thalf = log(2)/fit_lambda
interpolated_rna = interpolate_rna(model_data$day, model_data$rna, t)
predicted_protein = Pt(interpolated_rna, t, fit_lambda)
points(t, interpolated_rna, type='l', lwd=0.5, col=rna_color)
points(t, predicted_protein, type='l', lwd=1, col=prot_color)
mtext(side=3, line=0, text='naive WT mice')
legend('bottomright',c('protein','RNA'),col=c(prot_color,rna_color),text.col=c(prot_color,rna_color),bty='n')
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(smry, 'Kinetics of ASO knockdown of PrP and Prnp RNA in WT mice with ASO 6.')



### B ASO N, HuKI mice #### 

ason_prot = rbind_files('data/','_summary.tsv') %>%
  mutate(plate=as.integer(substr(file,1,3))) %>%
  filter(plate %in% c(175,177))

ason_pk_rna = read_tsv('data/PRP230911.tsv', col_types=cols()) %>%
  mutate(animal = gsub('-','',animal))

ason_prot %>%
  rename(animal = sample) %>%
  inner_join(ason_pk_rna, by='animal') %>%
  select(plate, animal, tx, dose_ug, day, ngml_av, rna, pk_ug_g) %>%
  mutate(rna = rna / mean(rna[tx=='PBS'])) %>%
  group_by(plate) %>%
  mutate(rel = ngml_av / mean(ngml_av[tx=='PBS'])) %>%
  ungroup() %>%
  mutate(day = as.numeric(gsub('D','',day))) %>%
  select(animal, tx, dose_ug, day, rel, rna, pk_ug_g) -> ason

ason %>%
  group_by(tx, dose_ug, day) %>%
  summarize(.groups='keep',
            rna_mean = mean(rna, na.rm=T), rna_l95 = lower(rna), rna_u95 = upper(rna),
            prp_mean = mean(rel, na.rm=T), prp_l95 = lower(rel), prp_u95 = upper(rel),
            pk_mean = mean(pk_ug_g, na.rm=T), pk_l95 = lower(pk_ug_g), pk_u95 = upper(pk_ug_g)) %>%
  ungroup() -> ason_smry

ason %>%
  select(animal, tx, dose_ug, day, rna, prp=rel, pk=pk_ug_g) -> ason_out
ason_smry %>%
  arrange(dose_ug, day) -> ason_smry_out
write_supp_table(ason_smry_out, 'Kinetics of ASO knockdown of PrP and PRNP RNA in ki817 mice with ASO N.')

hl = ason
par(mar=c(4,4,3,1))
xlims = c(0, 45)
ylims = c(0, 1.25)
rna_color = '#00FB31'
prot_color = '#0001CD'
ctl_color = '#A9A9A9'
plot(NA, NA, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side = 1, at = xats, tck=-0.025, labels=NA)
axis(side = 1, at = xats[xats %in% hl$day], lwd=0, cex.axis=0.8, line=-0.5)
mtext(side=1, line=1.6, text='days post-dose', cex=0.8)
axis(side = 2, at= 0:5/4, labels=NA)
axis(side = 2, at= 0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='residual', cex=0.8)
abline(h=1, lty=3)
points(hl$day, hl$rel, col=alpha(ifelse(hl$tx=='PBS',ctl_color,prot_color),ci_alpha), pch=19)
points(hl$day, hl$rna, col=alpha(ifelse(hl$tx=='PBS',ctl_color,rna_color),ci_alpha), pch=19)
hl %>%
  filter(tx != 'PBS') %>%
  group_by(tx, day) %>%
  summarize(.groups='keep',
            rna_mean = mean(rna),
            rna_l95 = lower(rna),
            rna_u95 = upper(rna),
            prp_mean = mean(rel),
            prp_l95 = lower(rel),
            prp_u95 = upper(rel)) %>%
  ungroup() -> smry
barwidth=0.5
segments(x0=smry$day-barwidth, x1=smry$day+barwidth, y0=smry$prp_mean, col=prot_color)
segments(x0=smry$day-barwidth, x1=smry$day+barwidth, y0=smry$rna_mean, col=rna_color)
arrows(x0=smry$day, y0=smry$prp_l95, y1=smry$prp_u95, col=prot_color, code=3, angle=90, length=0.02)
arrows(x0=smry$day, y0=smry$rna_l95, y1=smry$rna_u95, col=rna_color, code=3, angle=90, length=0.02)
dt = 0.01
t = seq(min(hl$day),max(hl$day),dt)
model_data = as.list(hl %>% filter(tx!='PBS') %>% select(rna, protein=rel, day))
nlsfit = nls.lm(par=c(lambda=log(2)/5), fn=calculate_residuals, data=model_data, dt=0.01)
fit_lambda = as.numeric(nlsfit$par['lambda'])
thalf = log(2)/fit_lambda
interpolated_rna = interpolate_rna(model_data$day, model_data$rna, t)
points(t, interpolated_rna, type='l', lwd=0.5, col=rna_color)
points(t, Pt(interpolated_rna, t, fit_lambda), type='l', lwd=1, col=prot_color)
points(t, Pt(interpolated_rna, t, log(2)/5), type='l', lwd=1, col=prot_color, lty=3)
mtext(side=3, line=0, text='naive ki817 mice')
# mtext(side=3, line=0, text=paste0('t1/2 = ',formatC(thalf,digits=1,format='f'), ' days'))
legend('bottomright',c('protein','RNA'),col=c(prot_color,rna_color),text.col=c(prot_color,rna_color),bty='n')
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

### C ASO 6 RML WT half-life ####

rbind(cbind(read_tsv('data/226_summary.tsv', col_types=cols()),plate=226),
      cbind(read_tsv('data/237_summary.tsv', col_types=cols()),plate=237)) -> elisa

read_tsv('data/PRP240407.tsv', col_types=cols()) %>%
  mutate(animal = as.character(animal)) -> meta

read_tsv('data/PRP240407_rna.tsv', col_types=cols()) %>%
  mutate(animal = as.character(mouse)) %>%
  select(-mouse) -> rna

elisa %>%
  select(animal=sample, ngml_av, plate) %>%
  inner_join(meta, by='animal') %>%
  group_by(plate) %>%
  mutate(ctl_mean = mean(ngml_av[tx=='none'])) %>%
  ungroup() %>%
  mutate(rel = ngml_av/ctl_mean) %>%
  inner_join(rna, by='animal') -> hl


par(mar=c(4,4,3,1))
xlims = c(0, 45)
ylims = c(0, 1.25)
rna_color = '#00FB31'
prot_color = '#0001CD'
ctl_color = '#A9A9A9'
plot(NA, NA, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side = 1, at = xats, tck=-0.025, labels=NA)
axis(side = 1, at = xats[xats %in% hl$day], lwd=0, cex.axis=0.8, line=-0.5)
mtext(side=1, line=1.6, text='days post-dose', cex=0.8)
axis(side = 2, at= 0:5/4, labels=NA)
axis(side = 2, at= 0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='residual', cex=0.8)
abline(h=1, lty=3)
par(xpd=T)
points(hl$day, hl$rel, col=alpha(ifelse(hl$tx=='none',ctl_color,prot_color),ci_alpha), pch=19)
points(hl$day, hl$rna, col=alpha(ifelse(hl$tx=='none',ctl_color,rna_color),ci_alpha), pch=19)
par(xpd=F)
hl %>%
  filter(tx == 'ASO6') %>%
  group_by(tx, day) %>%
  summarize(.groups='keep',
            rna_mean = mean(rna),
            rna_l95 = lower(rna),
            rna_u95 = upper(rna),
            prp_mean = mean(rel),
            prp_l95 = lower(rel),
            prp_u95 = upper(rel)) %>%
  ungroup() -> smry
barwidth=0.5
segments(x0=smry$day-barwidth, x1=smry$day+barwidth, y0=smry$prp_mean, col=prot_color)
segments(x0=smry$day-barwidth, x1=smry$day+barwidth, y0=smry$rna_mean, col=rna_color)
arrows(x0=smry$day, y0=smry$prp_l95, y1=smry$prp_u95, col=prot_color, code=3, angle=90, length=0.02)
arrows(x0=smry$day, y0=smry$rna_l95, y1=smry$rna_u95, col=rna_color, code=3, angle=90, length=0.02)
dt = 0.01
t = seq(min(hl$day),max(hl$day),dt)
model_data = as.list(hl %>% filter(tx=='ASO6') %>% select(rna, protein=rel, day))
nlsfit = nls.lm(par=c(lambda=log(2)/5), fn=calculate_residuals, data=model_data, dt=0.01)
fit_lambda = as.numeric(nlsfit$par['lambda'])
thalf = log(2)/fit_lambda
interpolated_rna = interpolate_rna(model_data$day, model_data$rna, t)
predicted_protein = Pt(interpolated_rna, t, fit_lambda)
points(t, interpolated_rna, type='l', lwd=0.5, col=rna_color)
points(t, predicted_protein, type='l', lwd=1, col=prot_color)
mtext(side=3, line=0, text='RML prion-infected WT mice')
# mtext(side=3, line=0, text=paste0('t1/2 = ',formatC(thalf,digits=1,format='f'), ' days'))
legend('bottomright',c('protein','RNA'),col=c(prot_color,rna_color),text.col=c(prot_color,rna_color),bty='n')
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(smry, 'Kinetics of ASO knockdown of PrP and Prnp RNA in RML prion-infected WT mice treated with ASO 6.')

### D ASO 6 rat CSF #### 

elisa = rbind_files('data/','11[789]_summary.tsv') %>%
  mutate(plate = as.integer(substr(file,1,3)))

plate_meta = tibble(plate=c(119,118,117), 
                    tissue=c('CSF','cerebellum','cerebrum'),
                    color=c('#39B7CD','#EE7600','#2E0854'))
animal_meta = read_tsv('data/PRP221222.tsv', col_types=cols())
tx_meta = tibble(tx=c('PBS','ASO 6'), pch=c(1, 19))

elisa %>%
  inner_join(plate_meta, by='plate') %>%
  inner_join(animal_meta, by=c('sample'='animal')) %>%
  inner_join(tx_meta, by=c('tx')) %>%
  rename(animal = sample) -> rhl

# do the controls look different by day?
#summary(aov(ngml_av ~ days, data=subset(rhl, tx=='PBS' & tissue=='cerebrum'))) # P = 0.041
#summary(aov(ngml_av ~ days, data=subset(rhl, tx=='PBS' & tissue=='cerebellum'))) # P = 0.0039
#summary(aov(ngml_av ~ days, data=subset(rhl, tx=='PBS' & tissue=='CSF'))) # = 0.27

rhl %>%
  filter(tx=='PBS') %>%
  group_by(plate, days, tissue) %>%
  summarize(.groups='keep', saline_mean=mean(ngml_av)) %>%
  ungroup() -> saline_means

rhl %>%
  inner_join(saline_means, by=c('plate','tissue','days')) %>%
  mutate(rel = ngml_av / saline_mean) %>%
  select(plate, animal, tissue, days, ngml_av, rel, tx, color, pch) -> rhl_indivs

rhl_indivs %>%
  group_by(tissue, days, tx, color) %>%
  summarize(.groups='keep',
            mean = mean(rel),
            l95 = lower(rel),
            u95 = upper(rel)) %>%
  ungroup() -> rhl_smry

par(mar=c(4,4,3,1))
xlims = c(14, 60)
ylims = c(0, 1.25)
yats = 0:8/4
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=1, at=unique(rhl_smry$days), tck=-0.025, labels=NA)
axis(side=1, at=unique(rhl_smry$days), lwd=0, line=-0.5, cex.axis=0.8)
mtext(side=1, line=1.5, text = 'days post-dose', cex=0.8)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
axis(side = 2, at= 0:5/4, labels=NA)
axis(side = 2, at= 0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='residual', cex=0.8)
abline(h=1, lty=3)
set.seed(1)
points(x=jitter(rhl_indivs$days,amount=2), y=rhl_indivs$rel, pch=rhl_indivs$pch, col=alpha(rhl_indivs$color, ci_alpha))
for (this_tissue in unique(rhl_smry$tissue)) {
  rhl_smry %>%
    filter(tx=='ASO 6') %>%
    filter(tissue == this_tissue) -> rhl_active
  points(x=rhl_active$days, y=rhl_active$mean, col=rhl_active$color, type='l', lwd=2)
  arrows(x0=rhl_active$days, y0=rhl_active$l95, y1=rhl_active$u95, col=rhl_active$color, code=3, angle=90, length=0.05, lwd=1.5)
}
mtext(side=3, line=0, text='naive rats')
legend('bottomleft', plate_meta$tissue, col=plate_meta$color, text.col=plate_meta$color, pch=19, bty='n', cex=0.8)
legend('bottom', tx_meta$tx, col='#000000', pch=tx_meta$pch, bty='n', cex=0.8)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(rhl_smry, 'Kinetics of PrP in rat brain and CSF following ASO 6 dosing.')

silence_is_golden = dev.off() ### End Figure 3 ####



ason_pk_rna %>%
  mutate(conc_ug_g = pk_ug_g) %>%
  mutate(day = as.integer(gsub('D','',day))) %>%
  mutate(dose_color = '#FF56AD') -> ason
ason %>%
  filter(tx != 'PBS') %>%
  mutate(day = as.integer(gsub('D','',day))) %>%
  group_by(day, dose_color) %>%
  summarize(.groups='keep',
            pk_mean = mean(conc_ug_g),
            pk_l95 = lower(conc_ug_g),
            pk_u95 = upper(conc_ug_g)) %>%
  ungroup() -> ason_smry

write_supp_table(ason_smry, 'Pharmacokinetics of ASO N in ki817 mice.')

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





## Figure 4 #### 
tell_user('done.\nCreating Figure 4...')
resx=300
png('display_items/figure-4.png',width=6.5*resx,height=4.5*resx,res=resx)

layout_matrix = matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,
                         5,5,5,6,6,6,7,7,7,8,8,8), nrow=2, byrow=T)
layout(layout_matrix)

panel = 1

ffi_all = read_tsv('data/iqp/ffi_turnover.tsv', col_types=cols()) %>%
  mutate(total = light + heavy) %>%
  mutate(prop_labeled = heavy / (total))

leg = tibble(genotype=c("129(TT-3F4-FFI)HOZ","129(TT-3F4WT)","B6/N"),
             disp = c('ki-3F4-FFI','ki-3F4-WT','C57BL/6N WT'),
             color=c('#D95F02','#22127A','#77127A'),
             xgeno = c(1,2,3),
             ttest_grouping = c('test','control','control'))
# 
# ffi_all %>% 
#   inner_join(leg, by='genotype') %>%
#   group_by(protein, peptide) %>%
#   mutate(n_ages_detected = length(unique(age))) %>%
#   filter(n_ages_detected == 2) %>%
#   ungroup() %>%
#   select(-n_ages_detected) %>%
#   group_by(protein, peptide) %>%
#   summarize(.groups='keep',
#             n = n(),
#             pval_total = t.test(total[genotype=='129(TT-3F4-FFI)HOZ'], total[genotype!='129(TT-3F4-FFI)HOZ'])$p.value,
#             ratio_total = mean(total[genotype=='129(TT-3F4-FFI)HOZ']) / mean(total[genotype!='129(TT-3F4-FFI)HOZ']),
#             label_obj = label_difference(prop_labeled, chow_days, genotype),
#             ratio_thalf = label_obj$thalf_ratio,
#             br_ffi_pval = label_obj$br_ffi_pval,
#             br_interaction_pval = label_obj$br_interaction_pval) %>%
#   ungroup() %>%
#   select(-label_obj) -> genotypic_diffs_regardless_age
# genotypic_diffs_regardless_age$bonf_total = pmin(1, genotypic_diffs_regardless_age$pval_total * nrow(genotypic_diffs_regardless_age))
# genotypic_diffs_regardless_age$bonf_label = pmin(1, genotypic_diffs_regardless_age$pval_label * nrow(genotypic_diffs_regardless_age))



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
            pval_total = t.test(total[genotype=='129(TT-3F4-FFI)HOZ'], total[genotype!='129(TT-3F4-FFI)HOZ'])$p.value,
            ratio_total = mean(total[genotype=='129(TT-3F4-FFI)HOZ']) / mean(total[genotype!='129(TT-3F4-FFI)HOZ']),
            label_obj = label_difference(prop_labeled, chow_days, genotype),
            ratio_thalf = label_obj$thalf_ratio, 
            pval_label_ffi = label_obj$br_ffi_pval, 
            pval_label_interaction = label_obj$br_interaction_pval) %>%
  ungroup() %>%
  select(-label_obj) %>%
  mutate(log2_ratio_total = log2(ratio_total),
         log2_ratio_thalf = log2(ratio_thalf)) %>%
  mutate(total_color = case_when(log2_ratio_total < 0 ~ '#FF0000',
                                 log2_ratio_total > 0 ~ '#0000FF',
                                 log2_ratio_total ==0 ~ '#000000')) %>%
  mutate(label_ffi_color = case_when(log2_ratio_thalf < 0 ~ '#FF0000',
                                 log2_ratio_thalf > 0 ~ '#0000FF',
                                 log2_ratio_thalf ==0 ~ '#000000'))%>%
  mutate(pepnick = substr(peptide, 1, 4)) -> genotypic_diffs_by_age


for (this_age in c('young','aged')) {
  ### A-D young FFI mice ####
  
  # volcano of total protein abundance
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
  # axis(side=2, at=ybigs, labels=NA, tck=-0.05)
  # axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, las=2, line=-0.5)
  axis(side=2, at=yats, labels=NA, tck=-0.02)
  #mtext(side=2, line=1.6, text='-log10(P)')
  abline(h=-log10(0.05 / nrow(ffi)))
  par(xpd=T)
  points(x=ffi$log2_ratio_thalf, y=pmin(-log10(ffi$pval_label_ffi), max(ylims)), pch=20, col=ffi$label_ffi_color)
  ffi %>% 
    filter(protein=='PRNP') -> to_label
  # amazingly, bquote is not vectorized? the only way to do this is in a loop??
  # text(x=to_label$log2_ratio_total, y=pmin(-log10(to_label$pval_total), max(ylims)), pos=4, labels=bquote(italic(.(to_label$protein)) * ' ' * .(to_label$pepnick)), cex=0.6)
  points(x=to_label$log2_ratio_thalf, y=pmin(-log10(to_label$pval_label_ffi)), pch=1)
  for (i in 1:nrow(to_label)) {
    text(x=to_label$log2_ratio_thalf[i], y=pmin(-log10(to_label$pval_label_ffi[i]), max(ylims)), 
         pos=ifelse(to_label$pepnick[i]=='VVEQ' & panel==6, 2, 4), # cheat to avoid label overlap 
         labels=bquote(italic(.(to_label$protein[i])) * ' ' * .(to_label$pepnick[i])), cex=0.6)
  }
  par(xpd=F)
  mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1
  
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
    mtext(side=1, line=1.6, text='day')
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
      
      points(x=this_smry$chow_days, y=this_smry$mean, type='l', lwd=2, col=this_smry$color)
      polygon(x=c(this_smry$chow_days, rev(this_smry$chow_days)), y=c(this_smry$l95, rev(this_smry$u95)), col=alpha(this_smry$color, ci_alpha), border=NA)
      points(x=subs$chow_days, y=subs$prop_labeled, col=subs$color, pch=1, bg='#FFFFFF')
      # 
      # thalf_8day_estimate = thalf_ests %>%
      #   filter(peptide==this_peptide & genotype==this_genotype) %>% 
      #   pull(thalf_8day_estimate)
      # theoreticals = proportion_labeled(thalf=thalf_8day_estimate, t=t, avails=free_lysine)
      # 
      # points(x=t, y=theoreticals, type='l', col='#000000', lty=1, lwd=0.5)
      # 
      #  mtext(side=4, las=2, line=0.25, at=this_smry$mean[this_smry$chow_days==8], text=formatC(thalf_8day_estimate,format='f',digits=1), cex=0.6)
      
      legend('topleft', leg$disp, col=leg$color, lwd=1, text.col=leg$color, bty='n', cex=0.7)
      
    }
    mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1
    
  }
  
  #write_supp_table(smry, paste0('Mean labeling of peptides by isotopic chow days in ',this_age,' FFI mice and comparators.'))
  #write_supp_table(thalf_ests, paste0('PrP half-life estimates for ',this_age,' FFI mice and comparators.'))

}

silence_is_golden = dev.off()
### end Fig 4 ####




ffi_turnover %>%
  filter(protein=='PRNP' & !grepl('QHT',peptide)) %>%
  group_by(protein, peptide, genotype, age, chow_days) %>%
  summarize(.groups='keep',
            mean_light=mean(light),
            mean_heavy=mean(heavy),
            mean_labeled=mean(heavy/(heavy+light))) %>%
  ungroup() %>%
  group_by(protein, peptide, genotype, age, chow_days) %>%
  summarize(.groups='keep',
            thalf_estimate = find_thalf(mean_labeled, t, which_t=chow_days, avails=free_lysine)) %>%
  ungroup() -> prnp_ffi_turnover


ffi_turnover %>%
  group_by(protein, peptide, genotype, age, chow_days) %>%
  summarize(.groups='keep',
            mean_light=mean(light),
            mean_heavy=mean(heavy),
            mean_labeled=mean(heavy/(heavy+light))) %>%
  ungroup() %>%
  group_by(protein, peptide, genotype, age, chow_days) %>%
  summarize(.groups='keep',
            thalf_estimate = find_thalf(mean_labeled, t, which_t=chow_days, avails=free_lysine)) %>%
  ungroup() -> all_ffi_turnover





# Supplement #### 

tell_user('done.\nFinalizing supplementary tables...')

# write the supplement directory / table of contents
supplement_directory %>% rename(table_number = name, description=title) -> contents
addWorksheet(supplement,'contents')
bold_style = createStyle(textDecoration = "Bold")
writeData(supplement,'contents',contents,headerStyle=bold_style,withFilter=T)
freezePane(supplement,'contents',firstRow=T)
# move directory to the front
original_order = worksheetOrder(supplement)
n_sheets = length(original_order)
new_order = c(n_sheets, 1:(n_sheets-1))
worksheetOrder(supplement) = new_order
activeSheet(supplement) = 'contents'
# now save
saveWorkbook(supplement,supplement_path,overwrite = TRUE)

elapsed_time = Sys.time() - overall_start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',round(as.numeric(elapsed_time),1),' ',units(elapsed_time),'.\n'))


