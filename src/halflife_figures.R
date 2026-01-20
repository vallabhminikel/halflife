source("src/utils.R")
source("src/funcs.R")


# Figures ####

source("src/fig1.R")
source("src/figS1.R")
source("src/figS2.R")
source("src/figS3.R")
source("src/fig2.R")
source("src/fig3.R")
source("src/figS4.R")
source("src/figS5.R")
source("src/figS6.R")
source("src/fig4.R")
source("src/figS7.R")

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

## Supp tables related to Figure 4 ####
write_supp_table(genotypic_diffs_by_age, 'Differences in protein abundance and half-life between FFI mice and controls by peptide and age.')

ffi_all %>%
  group_by(protein, peptide, age, genotype, chow_days) %>%
  summarize(.groups='keep',
            n = n(),
            mean_prop_labeled = mean(prop_labeled),
            l95_prop_labeled = lower(prop_labeled),
            u95_prop_labeled = upper(prop_labeled)) %>%
  ungroup() -> ffi_labeled_smry

write_supp_table(ffi_labeled_smry, 'Proportion labeled by peptide, genotype, age, and days of isotopic chow in FFI and control mice.')

ffi_all %>%
  group_by(protein, peptide, age, genotype) %>%
  summarize(.groups='keep',
            n = n(),
            mean_total = mean(total),
            l95_total = lower(total),
            u95_total = upper(total)) %>%
  ungroup() -> ffi_total_smry

write_supp_table(ffi_total_smry, 'Abundance by peptide, genotype, age in FFI and control mice.')

ffi_all %>% 
  inner_join(leg, by='genotype') %>%
  group_by(protein, peptide) %>%
  mutate(n_ages_detected = length(unique(age))) %>%
  filter(n_ages_detected == 2) %>%
  ungroup() %>%
  select(-n_ages_detected) %>%
  group_by(protein, peptide, genotype, age) %>%
  summarize(.groups='keep',
            n = n(),
            thalf_single_rate = fit_isotopic_thalf(chow_days, prop_labeled)
            ) %>%
  ungroup() -> thalves_by_genotype_and_age

# Join with mixture model results from Julia
mixture_model_results = read_tsv('data/fig4/table-mixture-model-results.tsv', col_types=cols())

thalves_by_genotype_and_age = thalves_by_genotype_and_age %>%
  inner_join(mixture_model_results %>% 
               select(protein, peptide, genotype, age, 
                      thalf_fast, prop_fast, thalf_slow, prop_slow, 
                      thalf_mixture_effective = thalf_effective,
                      aic_single, aic_mixture, delta_aic),
             by = c("protein", "peptide", "genotype", "age"))

write_supp_table(thalves_by_genotype_and_age, 'Protein half-life estimates by peptide, genotype, and age in FFI and control mice.')



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


