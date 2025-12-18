suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(openxlsx))
suppressMessages(library(minpack.lm))

dt = 0.01
t = seq(0, 8, dt)

runs_test = function(x) {
  x = x[!is.na(x) & x != 0]
  n = length(x)
  n_pos = sum(x > 0)
  n_neg = sum(x < 0)
  runs = rle(sign(x))
  n_runs = length(runs$lengths)
  expected_runs = 1 + (2 * n_pos * n_neg) / n
  var_runs = (2 * n_pos * n_neg * (2 * n_pos * n_neg - n)) / (n^2 * (n - 1))
  z_stat = (n_runs - expected_runs) / sqrt(var_runs)
  p_value = 2 * pnorm(-abs(z_stat))
  list(statistic = z_stat, p.value = p_value, n_runs = n_runs, n_pos = n_pos, n_neg = n_neg)
}

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
ci_alpha = 0.35

supplement_path = 'display_items/supplement_mixture.xlsx'
supplement = createWorkbook()
write_supp_table = function(tbl, title='') {
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=2,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,tbl,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  write_tsv(tbl,paste0('display_items/table-mixture-',table_name,'.tsv'), na='')
}

free_lysine = function(t) {
  pmax(0,1-0.503*(exp(-t*log(2)*0.799))-0.503*exp(-t*log(2)/39.423))
}

always_all = function(t) {
  1
}

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

proportion_labeled_mixture = function(halflives, weights, t, avails) {
  K = length(halflives)
  lambdas = log(2) / halflives
  N = numeric(length(t))
  N_components = matrix(0, nrow=length(t), ncol=K)
  N_components[1,] = weights
  N[1] = 1
  for (i in 2:length(t)) {
    for (k in 1:K) {
      protein_turned_over = lambdas[k] * dt
      original = (1 - protein_turned_over) * N_components[i-1, k]
      nascent = protein_turned_over * weights[k] * avails(t[i])
      N_components[i, k] = original + nascent
    }
    N[i] = sum(N_components[i,])
  }
  return(list(total=N, components=N_components))
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

find_thalf = function(prop_labeled, t=seq(0,8,.01), which_t=8, avails=free_lysine) {
  result = numeric(length(prop_labeled))
  for (i in 1:length(prop_labeled)) {
    result[i] = find_thalf_unary(prop_labeled[i], t, which_t, avails)
  }
  return (result)
}

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

Pt_mixture = function(R, t, lambdas, weights) {
  K = length(lambdas)
  P_components = matrix(0, nrow=length(t), ncol=K)
  P_components[1,] = weights
  P = numeric(length(t))
  P[1] = 1
  for (i in 2:length(t)) {
    dt_val = (t[i] - t[i-1])
    for (k in 1:K) {
      dP = lambdas[k] * dt_val * R[i-1] * weights[k] - lambdas[k] * dt_val * P_components[i-1, k]
      P_components[i, k] = P_components[i-1, k] + dP
    }
    P[i] = sum(P_components[i,])
  }
  return(list(total=P, components=P_components))
}

calculate_residuals = function(par, data, dt=0.01) {
  lambda = par[['lambda']]
  t = seq(0,max(data$day),dt)
  R = interpolate_rna(data$day, data$rna, t)
  P_pred = Pt(R, t, lambda)[match(data$day, t)]
  residuals = data$protein - P_pred
  return (residuals)
}

calculate_residuals_mixture = function(par, data, dt=0.01, K=2) {
  lambdas = par[paste0('lambda',1:K)]
  weights_raw = par[paste0('w',1:(K-1))]
  weights = c(weights_raw, 1-sum(weights_raw))
  if(any(weights < 0) | any(weights > 1)) {
    return(rep(1e10, length(data$day)))
  }
  t = seq(0,max(data$day),dt)
  R = interpolate_rna(data$day, data$rna, t)
  P_pred = Pt_mixture(R, t, lambdas, weights)$total[match(data$day, t)]
  residuals = data$protein - P_pred
  return (residuals)
}

pivot_tab  = read_tsv("data/iqp/combined_pivot_tab.tsv", col_types=cols())
iqp_meta   = read_tsv("data/iqp/meta.tsv", col_types=cols())
lloq_stats = read_tsv('data/iqp/lloq_stats.tsv', col_types=cols())
lit_half   = read_tsv('data/iqp/forna_halflife.tsv', col_types=cols())
name_map   = read_tsv('data/iqp/name_map.tsv', col_types=cols())

pivot_tab %>%
  rename(animal_id = identifier) %>%
  inner_join(iqp_meta, by='animal_id') %>%
  mutate(total = heavy + light) %>%
  mutate(prop_labeled = heavy / total) %>%
  inner_join(lloq_stats, by=c('peptide','protein'='gene')) %>%
  mutate(above_lloq = pmin(light,heavy) >= lloq) -> iqp_data

tissue_meta = tibble(tissue=c('brain','colon'),
                     tissue_x = c(0, 5),
                     color=c('#6D00C5','#B7760A'))

genotype_meta = tibble(genotype = c('C57BL/6N','Tg25109 het; ZH3/ZH3','Tg25109 hom; ZH3/ZH3'),
                       explanation = c('C57BL/6N','Tga20 het','Tga20 hom'),
                       xgeno=1:3)

prnp_pep_meta = tibble(peptide=c('VVEQMCVTQYQK','YPNGGWNTGGSR'),
                       short=c('VVEQ','YPNG'),
                       pep_color=c('#0001CD','#B00000'),
                       xpep=1:2)

resx=300
png('display_items/figure-2a-mixture.png',width=6.5*resx,height=3.0*resx,res=resx)

par(mar = c(3, 3, 3, 3))

xlims = c(0, 12)
ylims = c(0, 1)

plot(NA, NA, xlim = xlims, ylim = ylims, xaxs = 'i', yaxs = 'i', axes = F, ann = F)
axis(1, at = xlims, labels = NA, tck = 0)
axis(1, at = c(0:10 * 2), labels = c(0:10 * 2), cex.axis=0.7, lwd=0, line=-0.5)
axis(1, at = c(0:10 * 2), tck = -0.05, labels=NA)
axis(1, at = c(0:12 * 1), tck = -0.025, labels = NA)
mtext(side = 1, line = 1.6, text = 'protein half life', cex = 0.7)

axis(side = 2, at = 0:5 / 4, labels = scales::percent(0:5 / 4), las = 2, lwd = 0, line=-0.5, cex.axis=0.8)
mtext(side = 2, line = 2.25, text = 'proportion labeled', cex = 0.7)
axis(2, at = ylims, labels = NA)
axis(2, at = c(0:4 * 0.25), labels = NA, cex = 0.8, tck = -0.05)
axis(2, at = c(0:8 * 0.125), labels = NA, cex = 0.8, tck = -0.025)
mtext(side = 3, line = 0, text = 'brain', cex = 0.8)

thalf_values = seq(0, 12, .1)
lysine_meta = tibble(model = c('100%','plasma free'),
                     color = c("#429898", "#A90101"))
points(thalf_values, proportion_labeled_t(thalf_values, t=t, which_t = 8, avails=free_lysine), type = 'l', col = lysine_meta$color[2], lwd = 2)

lloq_meta = tibble(above_lloq = c(1,0,0),
                   color = c('#000000','#C9C9C9','#BFB7EF'),
                   disp = c('OK','below LLOQ','controls'))

reported = lit_half %>%
  mutate(reported_halflife = case_when(is.na(cerebellum) ~ cortex,
                                       is.na(cortex) ~ cerebellum,
                                       TRUE ~ (cortex + cerebellum)/2)) %>%
  inner_join(name_map, by=c('gene'='fornasiero_name')) %>%
  filter(iqp_name %in% iqp_data$protein)

iqp_data %>%
  filter(genotype=='C57BL/6N',
         tissue == 'brain') %>%
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

points(smry$reported_halflife, smry$mean, pch=1, cex = 1.2, col =smry$color)

iqp_data %>%
  filter(labeled == 'yes') %>%
  filter(protein == 'PRNP') %>%
  filter(!(genotype=='Tg25109 het; ZH3/ZH3' & peptide=='VVEQMCVTQYQK')) %>%
  filter(above_lloq) %>%
  inner_join(prnp_pep_meta, by='peptide') %>%
  inner_join(tissue_meta, by='tissue') %>%
  inner_join(genotype_meta, by='genotype') -> all_prnp_plab

all_prnp_plab %>%
  group_by(genotype, explanation, tissue, peptide, color, pep_color, short) %>%
  summarize(.groups='keep',
            n = n(),
            mean = mean(prop_labeled),
            l95 = lower(prop_labeled),
            u95 = upper(prop_labeled)) %>%
  ungroup() %>%
  mutate(estimated_thalf = find_thalf(mean),
         estimated_thalf_if_100 = find_thalf(mean, avails=always_all)) -> prnp_iqp_smry

prnp_iqp_smry %>%
  filter(tissue=='brain') %>%
  filter(genotype=='C57BL/6N') %>%
  arrange(desc(mean)) -> subs
prnp_color = '#0001CD'

avails = free_lysine
label_x = 1.5
segments(x0=rep(0, nrow(subs)), x1=find_thalf(subs$mean,avails=avails), y0=subs$mean, col=prnp_color, lty=3)
text(x=rep(label_x,2), y=subs$mean+c(-0.02,0.02), pos=c(3,1), labels=subs$short, col="#0001CD", cex=0.5)
segments(x0=find_thalf(subs$mean,avails=avails), y0=rep(0,nrow(subs)), y1=subs$mean, col=prnp_color, lty=3)

legend("topright", legend = lloq_meta$disp, col =lloq_meta$color, pch = 1, cex = 0.5, bty='n')

mtext(side=3, adj=0, text='A', line=0.5)

dev.off()

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
  inner_join(meta, by=c('animal')) %>%
  group_by(plate) %>%
  mutate(ctl_mean = mean(ngml_av[tx=='none'])) %>%
  ungroup() %>%
  mutate(rel = ngml_av/ctl_mean) -> hl_aso6

model_data_aso6 = as.list(hl_aso6 %>% filter(tx=='ASO6') %>% select(rna, protein=rel, day))

nlsfit_single = nls.lm(par=c(lambda=log(2)/5), fn=calculate_residuals, data=model_data_aso6, dt=0.01)
fit_lambda_single = as.numeric(nlsfit_single$par['lambda'])
thalf_single = log(2)/fit_lambda_single

nlsfit_mixture = nls.lm(
  par=c(lambda1=log(2)/2, lambda2=log(2)/10, w1=0.7),
  fn=calculate_residuals_mixture,
  data=model_data_aso6,
  dt=0.01,
  K=2,
  lower=c(log(2)/20, log(2)/50, 0.01),
  upper=c(log(2)/0.5, log(2)/1, 0.99)
)

fit_lambdas = as.numeric(nlsfit_mixture$par[c('lambda1','lambda2')])
fit_w1 = as.numeric(nlsfit_mixture$par['w1'])
fit_weights = c(fit_w1, 1-fit_w1)
fit_halflives = log(2)/fit_lambdas

t_fig = seq(0, max(model_data_aso6$day), dt)
interpolated_rna = interpolate_rna(model_data_aso6$day, model_data_aso6$rna, t_fig)
predicted_protein_single = Pt(interpolated_rna, t_fig, fit_lambda_single)
predicted_protein_mixture = Pt_mixture(interpolated_rna, t_fig, fit_lambdas, fit_weights)$total

residuals_single = model_data_aso6$protein - predicted_protein_single[match(model_data_aso6$day, t_fig)]
residuals_mixture = model_data_aso6$protein - predicted_protein_mixture[match(model_data_aso6$day, t_fig)]

runs_test_data = sign(residuals_single)
runs_test_result = runs_test(runs_test_data)

runs_test_table = tibble(
  model = 'Single Exponential',
  statistic = as.numeric(runs_test_result$statistic),
  p_value = as.numeric(runs_test_result$p.value),
  n_runs = sum(rle(runs_test_data)$lengths > 0),
  n_positive = sum(runs_test_data > 0),
  n_negative = sum(runs_test_data < 0)
)

write_supp_table(runs_test_table, 'Wald-Wolfowitz Runs Test on log residual signs for single exponential fit vs observed data (ASO6 study)')

mixture_params_table = tibble(
  component = 1:2,
  half_life_days = fit_halflives,
  decay_rate = fit_lambdas,
  mixing_proportion = fit_weights
) %>%
  arrange(half_life_days)

write_supp_table(mixture_params_table, 'Mixture model parameters (K=2) for ASO6 PrP knockdown data')

resx=300
png('display_items/figure-3-mixture.png',width=resx*6.5,height=5.5*resx, res=resx)

layout_matrix = matrix(c(1,2,
                         3,4), nrow=2, byrow=T)
layout(layout_matrix)
panel = 1

xats = c(0,1,3,7,10,14,21,28,35,42)

par(mar=c(4,4,3,1))
xlims = c(0, 45)
ylims = c(0, 1.25)
rna_color = '#00FB31'
prot_color = '#0001CD'
ctl_color = '#A9A9A9'
mix_color = '#FF0066'

plot(NA, NA, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side = 1, at = xats, tck=-0.025, labels=NA)
axis(side = 1, at = xats[xats %in% hl_aso6$day], lwd=0, cex.axis=0.8, line=-0.5)
mtext(side=1, line=1.6, text='days post-dose', cex=0.8)
axis(side = 2, at= 0:5/4, labels=NA)
axis(side = 2, at= 0:5/4, labels=percent(0:5/4), las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='residual', cex=0.8)
abline(h=1, lty=3)
points(hl_aso6$day, hl_aso6$rel, col=alpha(ifelse(hl_aso6$tx=='none',ctl_color,prot_color),ci_alpha), pch=19)
points(hl_aso6$day, hl_aso6$rna, col=alpha(ifelse(hl_aso6$tx=='none',ctl_color,rna_color),ci_alpha), pch=19)

hl_aso6 %>%
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

points(t_fig, interpolated_rna, type='l', lwd=0.5, col=rna_color)
points(t_fig, predicted_protein_single, type='l', lwd=1, col=prot_color)
points(t_fig, predicted_protein_mixture, type='l', lwd=2, col=mix_color, lty=2)

mtext(side=3, line=0, text='ASO6 WT mice - Single vs Mixture')
legend('bottomright',c('protein (single)','protein (mixture)','RNA'),
       col=c(prot_color,mix_color,rna_color),
       text.col=c(prot_color,mix_color,rna_color),
       lty=c(1,2,1),lwd=c(1,2,0.5),bty='n',cex=0.7)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

par(mar=c(4,4,3,1))
xlims = c(0, 45)
ylims = c(-0.3, 0.3)
plot(NA, NA, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side = 1, at = xats, tck=-0.025, labels=NA)
axis(side = 1, at = xats[xats %in% model_data_aso6$day], lwd=0, cex.axis=0.8, line=-0.5)
mtext(side=1, line=1.6, text='days post-dose', cex=0.8)
axis(side = 2, las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='residuals', cex=0.8)
abline(h=0, lty=3)
points(model_data_aso6$day, residuals_single, col=prot_color, pch=19, cex=1.2)
points(model_data_aso6$day, residuals_mixture, col=mix_color, pch=17, cex=1.2)
legend('topright',c('single exp','mixture K=2'),col=c(prot_color,mix_color),pch=c(19,17),bty='n',cex=0.7)
mtext(side=3, line=0, text='Residuals')
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

par(mar=c(4,4,3,1))
xlims = c(0, 45)
ylims = c(0.001, 1)
plot(NA, NA, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
axis(side = 1, at = xats, tck=-0.025, labels=NA)
axis(side = 1, at = xats[xats %in% model_data_aso6$day], lwd=0, cex.axis=0.8, line=-0.5)
mtext(side=1, line=1.6, text='days post-dose', cex=0.8)
axis(side = 2, at=10^(-3:0), las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='log protein (residual)', cex=0.8)
points(t_fig, predicted_protein_single, type='l', lwd=1, col=prot_color)
points(t_fig, predicted_protein_mixture, type='l', lwd=2, col=mix_color, lty=2)
points(model_data_aso6$day, model_data_aso6$protein, col='#000000', pch=19, cex=1.2)
mtext(side=3, line=0, text='Log Scale')
legend('bottomleft',c('observed','single exp','mixture K=2'),
       col=c('#000000',prot_color,mix_color),
       lty=c(NA,1,2),pch=c(19,NA,NA),lwd=c(NA,1,2),bty='n',cex=0.7)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

predicted_components = Pt_mixture(interpolated_rna, t_fig, fit_lambdas, fit_weights)$components
frac_fast = predicted_components[,1] / predicted_protein_mixture
frac_slow = predicted_components[,2] / predicted_protein_mixture

par(mar=c(4,4,3,1))
xlims = c(0, 45)
ylims = c(0, 1)
plot(NA, NA, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side = 1, at = xats, tck=-0.025, labels=NA)
axis(side = 1, at = xats[xats %in% model_data_aso6$day], lwd=0, cex.axis=0.8, line=-0.5)
mtext(side=1, line=1.6, text='days post-dose', cex=0.8)
axis(side = 2, at= 0:4/4, labels=NA)
axis(side = 2, at= 0:4/4, labels=percent(0:4/4), las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='population composition', cex=0.8)

fast_color = '#FF6B00'
slow_color = '#9B00FF'
polygon(c(t_fig, rev(t_fig)), c(rep(0, length(t_fig)), rev(frac_fast)), col=alpha(fast_color, 0.3), border=NA)
polygon(c(t_fig, rev(t_fig)), c(frac_fast, rev(rep(1, length(t_fig)))), col=alpha(slow_color, 0.3), border=NA)
lines(t_fig, frac_fast, col=fast_color, lwd=2)
lines(t_fig, frac_slow, col=slow_color, lwd=2)
mtext(side=3, line=0, text='Component Proportions')
legend('right',
       c(paste0('Fast (t1/2=',formatC(fit_halflives[1],digits=1,format='f'),'d)'),
         paste0('Slow (t1/2=',formatC(fit_halflives[2],digits=1,format='f'),'d)')),
       col=c(fast_color,slow_color),lwd=2,bty='n',cex=0.7)
mtext(side=3, adj=-0.2, text=LETTERS[panel], line=0.5); panel = panel + 1

dev.off()

resx=300
png('display_items/figure-4-mixture.png',width=resx*6.5,height=4.5*resx, res=resx)

layout_matrix = matrix(c(1,1,2,2,
                         3,3,4,4), nrow=2, byrow=T)
layout(layout_matrix)
panel = 1

halflives_sim = c(2, 15)
weights_sim = c(0.7, 0.3)
rates_sim = log(2) / halflives_sim
t_sim = seq(0, 40, 0.05)

y_mix = weights_sim[1] * exp(-rates_sim[1]*t_sim) + weights_sim[2] * exp(-rates_sim[2]*t_sim)

sq_error = function(k) sum((y_mix - exp(-k * t_sim))^2)
opt_result = optimize(sq_error, c(0, 5))
k_best = opt_result$minimum
y_fit = exp(-k_best * t_sim)

avg_rate = sum(weights_sim * rates_sim)
y_tangent = exp(-avg_rate * t_sim)

n_unstable = weights_sim[1] * exp(-rates_sim[1]*t_sim)
n_stable = weights_sim[2] * exp(-rates_sim[2]*t_sim)
frac_unstable = n_unstable / y_mix
frac_stable = n_stable / y_mix

par(mar=c(3,4,3,1))
xlims = c(0, max(t_sim))
ylims = c(0, 1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(0,40,10), tck=-0.025, labels=seq(0,40,10), lwd=0, line=-0.5)
mtext(side=1, line=1.6, text='time (days)', cex=0.8)
axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='fraction remaining', cex=0.8)
lines(t_sim, y_mix, col='#000000', lwd=3)
lines(t_sim, y_fit, col='#FF0000', lwd=2, lty=2)
legend('topright',c('True Mixture','Best Single Fit'),col=c('#000000','#FF0000'),lwd=c(3,2),lty=c(1,2),bty='n')
mtext(side=3, line=0, text='Decay Comparison')
mtext(side=3, adj=0, text=LETTERS[panel], line=0.5); panel = panel + 1

par(mar=c(3,4,3,3))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(0,40,10), tck=-0.025, labels=seq(0,40,10), lwd=0, line=-0.5)
mtext(side=1, line=1.6, text='time (days)', cex=0.8)
axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='population composition', cex=0.8)
lines(t_sim, frac_stable, col='#9B00FF', lwd=3, lty=1)
lines(t_sim, frac_unstable, col='#FF6B00', lwd=3, lty=1)
legend('right',
       c(paste0('Stable (t1/2=',halflives_sim[2],'d)'),
         paste0('Unstable (t1/2=',halflives_sim[1],'d)')),
       col=c('#9B00FF','#FF6B00'),lwd=3,bty='n')
mtext(side=3, line=0, text='Population Dynamics')
mtext(side=3, adj=0, text=LETTERS[panel], line=0.5); panel = panel + 1

par(mar=c(3,4,3,1))
ylims = c(0.001, 1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
axis(side=1, at=seq(0,40,10), tck=-0.025, labels=seq(0,40,10), lwd=0, line=-0.5)
mtext(side=1, line=1.6, text='time (days)', cex=0.8)
axis(side=2, at=10^(-3:0), las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='log fraction remaining', cex=0.8)
lines(t_sim, y_mix, col='#000000', lwd=3)
lines(t_sim, y_fit, col='#FF0000', lwd=2, lty=2)
lines(t_sim, y_tangent, col='#00AA00', lwd=2, lty=3)
legend('bottomleft',c('Mixture','Best Fit (Secant)','Avg Rate (Tangent)'),
       col=c('#000000','#FF0000','#00AA00'),lwd=c(3,2,2),lty=c(1,2,3),bty='n',cex=0.7)
mtext(side=3, line=0, text='Log-Convexity')
mtext(side=3, adj=0, text=LETTERS[panel], line=0.5); panel = panel + 1

diff = y_fit - y_mix
par(mar=c(3,4,3,1))
ylims_diff = range(diff) * 1.2
plot(NA, NA, xlim=xlims, ylim=ylims_diff, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=seq(0,40,10), tck=-0.025, labels=seq(0,40,10), lwd=0, line=-0.5)
mtext(side=1, line=1.6, text='time (days)', cex=0.8)
axis(side=2, las=2, lwd=0, line=-0.25)
mtext(side=2, line=2.5, text='error (fit - truth)', cex=0.8)
lines(t_sim, diff, col='#000000', lwd=2)
abline(h=0, lty=3, col='#999999')
text(5, 0.03, 'Fit > Mix\n(too slow)', col='#FF0000', cex=0.7)
text(25, -0.02, 'Fit < Mix\n(too fast)', col='#0000FF', cex=0.7)
mtext(side=3, line=0, text='Systematic Error')
mtext(side=3, adj=0, text=LETTERS[panel], line=0.5); panel = panel + 1

dev.off()

cat('Analysis complete.\n')
cat('Single exponential half-life:', formatC(thalf_single, digits=2, format='f'), 'days\n')
cat('Mixture component 1 half-life:', formatC(fit_halflives[1], digits=2, format='f'), 'days (weight:', formatC(fit_weights[1], digits=3, format='f'), ')\n')
cat('Mixture component 2 half-life:', formatC(fit_halflives[2], digits=2, format='f'), 'days (weight:', formatC(fit_weights[2], digits=3, format='f'), ')\n')
cat('Runs test p-value:', formatC(runs_test_result$p.value, digits=4, format='f'), '\n')
