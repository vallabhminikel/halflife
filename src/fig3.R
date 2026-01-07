
## Figure 3 ####
tell_user('done.\nCreating Figure 3...')
resx=300
png('display_items/figure-3.png',width=resx*6.5,height=5.5*resx, res=resx)

layout_matrix = matrix(c(1,2,
                         3,4), nrow=2, byrow=T)
layout(layout_matrix)
panel = 1

xats = c(0,1,3,7,10,14,21,28,35,42)

# Load mixture model parameters from Julia analysis
if (file.exists('data/aso/mixture_params.tsv')) {
  aso_mixture_params = read_tsv('data/aso/mixture_params.tsv', col_types=cols())
  use_mixture_models = TRUE
  tell_user('\nUsing mixture model parameters from Julia analysis...')
} else {
  use_mixture_models = FALSE
  tell_user('\nWarning: Mixture parameters not found (data/aso/mixture_params.tsv). Using single-rate models...')
  tell_user('Run: julia --project=. src/fig3_mixture_analysis.jl to generate mixture parameters.')
}

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

# Always plot RNA
interpolated_rna = interpolate_rna(model_data$day, model_data$rna, t)
points(t, interpolated_rna, type='l', lwd=0.5, col=rna_color)

# Always fit and plot single-rate model (thin dashed line)
nlsfit = nls.lm(par=c(lambda=log(2)/5), fn=calculate_residuals, data=model_data, dt=0.01)
fit_lambda = as.numeric(nlsfit$par['lambda'])
thalf_single = log(2)/fit_lambda
predicted_single = Pt(interpolated_rna, t, fit_lambda)
points(t, predicted_single, type='l', lwd=1, lty=2, col=prot_color)

# If mixture model available, plot it too (solid line)
if (use_mixture_models) {
  params = aso_mixture_params %>% filter(panel == 'ASO6_WT')
  if (nrow(params) > 0 && !is.na(params$thalf_fast)) {
    predicted_mixture = Pt_mixture(
      interpolated_rna, t,
      lambdas = c(log(2)/params$thalf_fast, log(2)/params$thalf_slow),
      weights = c(params$prop_fast, params$prop_slow)
    )$total
    points(t, predicted_mixture, type='l', lwd=2, col=prot_color)

    # Add delta AIC annotation
    delta_aic_text = sprintf('ΔAIC = %.1f', params$delta_aic)
    text(x=max(hl$day)*0.95, y=0.95, labels=delta_aic_text, pos=2, cex=0.7, col='#666666')
  }
}
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

# Always plot RNA
interpolated_rna = interpolate_rna(model_data$day, model_data$rna, t)
points(t, interpolated_rna, type='l', lwd=0.5, col=rna_color)

# Always fit and plot single-rate model (thin dashed line)
nlsfit = nls.lm(par=c(lambda=log(2)/5), fn=calculate_residuals, data=model_data, dt=0.01)
fit_lambda = as.numeric(nlsfit$par['lambda'])
thalf_single = log(2)/fit_lambda
predicted_single = Pt(interpolated_rna, t, fit_lambda)
points(t, predicted_single, type='l', lwd=1, lty=2, col=prot_color)

# If mixture model available, plot it too (solid line)
if (use_mixture_models) {
  params = aso_mixture_params %>% filter(panel == 'ASON_ki817')
  if (nrow(params) > 0 && !is.na(params$thalf_fast)) {
    predicted_mixture = Pt_mixture(
      interpolated_rna, t,
      lambdas = c(log(2)/params$thalf_fast, log(2)/params$thalf_slow),
      weights = c(params$prop_fast, params$prop_slow)
    )$total
    points(t, predicted_mixture, type='l', lwd=2, col=prot_color)

    # Add delta AIC annotation
    delta_aic_text = sprintf('ΔAIC = %.1f', params$delta_aic)
    text(x=max(hl$day)*0.95, y=0.95, labels=delta_aic_text, pos=2, cex=0.7, col='#666666')
  }
}
mtext(side=3, line=0, text='naive ki817 mice')
# mtext(side=3, line=0, text=paste0('t1/2 = ',formatC(thalf,digits=1,format='f'), ' days'))
legend('bottomright',c('protein','RNA'),col=c(prot_color,rna_color),text.col=c(prot_color,rna_color),bty='n')
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

### C ASO 6 RML WT half-life ####

# Use only plate 226 (plate 237 missing)
cbind(read_tsv('data/226_summary.tsv', col_types=cols()),plate=226) -> elisa

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

# Always plot RNA
interpolated_rna = interpolate_rna(model_data$day, model_data$rna, t)
points(t, interpolated_rna, type='l', lwd=0.5, col=rna_color)

# Always fit and plot single-rate model (thin dashed line)
nlsfit = nls.lm(par=c(lambda=log(2)/5), fn=calculate_residuals, data=model_data, dt=0.01)
fit_lambda = as.numeric(nlsfit$par['lambda'])
thalf_single = log(2)/fit_lambda
predicted_single = Pt(interpolated_rna, t, fit_lambda)
points(t, predicted_single, type='l', lwd=1, lty=2, col=prot_color)

# If mixture model available, plot it too (solid line)
if (use_mixture_models) {
  params = aso_mixture_params %>% filter(panel == 'ASO6_RML')
  if (nrow(params) > 0 && !is.na(params$thalf_fast)) {
    predicted_mixture = Pt_mixture(
      interpolated_rna, t,
      lambdas = c(log(2)/params$thalf_fast, log(2)/params$thalf_slow),
      weights = c(params$prop_fast, params$prop_slow)
    )$total
    points(t, predicted_mixture, type='l', lwd=2, col=prot_color)

    # Add delta AIC annotation
    delta_aic_text = sprintf('ΔAIC = %.1f', params$delta_aic)
    text(x=max(hl$day)*0.95, y=0.95, labels=delta_aic_text, pos=2, cex=0.7, col='#666666')
  }
}
mtext(side=3, line=0, text='RML prion-infected WT mice')
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
