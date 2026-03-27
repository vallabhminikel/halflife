
## Figure 1 #### 
tell_user('done.\nCreating Figure 1...')
resx=300
png('display_items/figure-1.png',width=6.5*resx,height=4.0*resx,res=resx)

layout_matrix = matrix(c(1,2,3,3,
                         1,2,3,3,
                         4,4,4,5), nrow=3, byrow=T)
layout(layout_matrix, widths=c(.75, .75, 1, 1))

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

par(mar=c(2,4,3,1))
xlims = c(0, max(tpm_by_tissue_smry$median_tpm)*1.1)
ylims = range(xykey$y) + c(-0.5, 0.5)
xbigs = 0:10*100
xats = 0:100*10
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, lwd=0, line=-0.5)
axis(side=1, at=xats, labels=NA, tck=-0.02)
mtext(side=1, line=1.6, text='Median TPM', cex=0.8)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, line=0.25, las=2, at=xykey$y, text=xykey$metatissue, cex=0.4)
mtext(side=3, line=1.0, text='Human', cex=0.6)
barwidth=0.8
rect(xleft=rep(0,nrow(tpm_by_metatissue_smry)), xright=tpm_by_metatissue_smry$median_median_tpm, 
     ybottom = tpm_by_metatissue_smry$y - barwidth/2, ytop = tpm_by_metatissue_smry$y + barwidth/2, 
     col=alpha(tpm_by_metatissue_smry$color, 0.7), border=NA)
points(x=tpm_by_tissue_smry$median_tpm, y=tpm_by_tissue_smry$y, pch=21, col=tpm_by_tissue_smry$color, bg='#FFFFFF', lwd=1, cex=0.6)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(tpm_by_metatissue_smry, 'Human median PRNP TPM by tissue (GTEx).')

### B Mouse TPMs ####
mouse_color_key = read_tsv('data/misc/sollner-2017-color-key.tsv', col_types=cols())
mouse_tpm_raw = read_tsv('data/misc/sollner-2017-mouse_tpm.txt', col_types=cols(), skip = 1, col_names = c("gene", scan('data/misc/sollner-2017-mouse_tpm.txt', what="", nlines=1, quiet = T)))
mouse_tpm_raw %>%
  filter(gene=='ENSMUSG00000079037') %>%
  pivot_longer(cols=-gene) %>%
  separate(name, into = c("animal", "id", "tissue_name"), sep = "_") %>%
  rename(tpm = value) %>%
  mutate(tissue_name = tolower(tissue_name)) %>%
  left_join(mouse_color_key, by=c('tissue_name'='mouse_tissue')) -> prnp_tpm

prnp_tpm %>%
  group_by(tissue_name, color) %>%
  summarize(.groups='keep',
            n = n(),
            median_tpm = median(tpm)) %>%
  ungroup() %>%
  arrange(desc(median_tpm)) %>%
  mutate(x = row_number()) %>%
  mutate(y = max(x) - x + 1) -> prnp_tpm_smry

prnp_tpm$y = prnp_tpm_smry$y[match(prnp_tpm$tissue_name, prnp_tpm_smry$tissue_name)]


par(mar=c(2,4,3,1))
xlims = c(0, 100)
ylims = range(prnp_tpm_smry$y) + c(-0.5, 0.5)
xbigs = 0:10*50
xats = 0:100*10
plot(NA, NA, xlim=xlims, ylim=ylims, ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=xbigs, labels=NA, tck=-0.05)
axis(side=1, at=xbigs, lwd=0, line=-0.5)
axis(side=1, at=xats, labels=NA, tck=-0.02)
mtext(side=1, line=1.6, text='Median TPM', cex=0.8)
axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
mtext(side=2, line=0.25, las=2, at=prnp_tpm_smry$y, text=prnp_tpm_smry$tissue_name, cex=0.5)
mtext(side=3, line=1.0, text='Mouse', cex=0.6)
barwidth=0.8
rect(xleft=rep(0,nrow(prnp_tpm_smry)), xright=prnp_tpm_smry$median_tpm, 
     ybottom = prnp_tpm_smry$y - barwidth/2, ytop = prnp_tpm_smry$y + barwidth/2, 
     col=alpha(prnp_tpm_smry$color, 0.7), border=NA)
points(x=prnp_tpm$tpm, y=prnp_tpm$y, pch=21, col=prnp_tpm$color, bg='#FFFFFF', lwd=1, cex=0.6)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

write_supp_table(prnp_tpm_smry, 'Mouse median Prnp TPM by tissue (Sollner 2017).')


### C Western #### 
par(mar=c(0.25,0,1,0.25))
western_panel = image_convert(image_read('data/fig1b.png'),'png')
plot(as.raster(western_panel))
mtext(side=3, adj=-0.0, text=LETTERS[panel], line=-0.75); panel = panel + 1

### D large 1:100 ELISA screen #### 

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
text(x=tissue_meta$x, y=-2.2, adj=1, srt=45, labels=tolower(tissue_meta$tissue), cex=0.8)
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





### E refined 1:25 ELISA screen #### 
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
