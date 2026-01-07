# Startup #### 

overall_start_time = Sys.time()
tell_user = function(...) { cat(file=stderr(), paste0(...)); flush.console() }

# Dependencies ####

tell_user('Loading required packages...')

suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(openxlsx))
# suppressMessages(library(smoother))  # Not available for R 4.5.1
suppressMessages(library(plotrix))
# suppressMessages(library(magick))  # Installation failed
suppressMessages(library(minpack.lm))



# Output streams #### 

tell_user('done.\nCreating output streams...')

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T


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
ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot



# Lt is proportion labeled at time t
Lt = function(avails, t, lambda) {
  L = numeric(length(t))
  dL = numeric(length(t))
  L[1] = 0
  for (i in 2:length(t)) {
    dt = (t[i]- t[i-1])
    # avails[i] * lambda * dt *  (1-L[i-1]) : accumulation of label. note because this is a differential (d), only the replacement of unlabeled with labeled is accounted for - replacement of labeled with labeled is invisible
    # (1-avails) * lambda * dt * (1-L[i-1]): decay of labeled and its replacement with (1-avails) worth of unlabeled
    dL[i] = avails[i] * lambda * dt * (1-L[i-1]) - (1-avails[i]) * lambda * dt * (L[i-1]) 
    L[i] = L[i-1] + dL[i]
  }
  return (L)
}


calculate_residuals_Lt = function(par, data, dt=0.01) {
  lambda = par[['lambda']]
  t = seq(0,max(data$day),dt)
  L_pred = Lt(avails, t, lambda)[match(data$day, t)]
  residuals = data$prop_labeled - L_pred
  return (residuals)
}


fit_isotopic_thalf = function(chow_days, prop_labeled, start_lambda=log(2)/5, avails_function=free_lysine, dt=0.01) {
  t = seq(0, max(chow_days), dt)
  avails = avails_function(t) 
  nlsfit = nls.lm(par=c(lambda=start_lambda), fn=calculate_residuals_Lt, data=tibble(day=chow_days, prop_labeled), dt=dt)
  fit_lambda = as.numeric(nlsfit$par['lambda'])
  thalf_found = log(2)/fit_lambda
  return(thalf_found)
}



# Data ####

tell_user('done.\nReading in data...')

## IQ Proteomics ####
pivot_tab  = read_tsv("data/iqp/combined_pivot_tab.tsv", col_types=cols())
iqp_meta   = read_tsv("data/iqp/meta.tsv", col_types=cols())
lloq_stats = read_tsv('data/iqp/lloq_stats.tsv', col_types=cols())
lit_half   = read_tsv('data/iqp/forna_halflife.tsv', col_types=cols()) 
name_map   = read_tsv('data/iqp/name_map.tsv', col_types=cols())


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


write_stats = function(...) {
  write(paste(list(...),collapse='',sep=''),text_stats_path,append=T)
  write('\n',text_stats_path,append=T)
}