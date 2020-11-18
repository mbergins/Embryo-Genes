library(tidyverse)
library(here)

# parnell_data = read_rds(here('Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_forShinyApp.rds'))
# parnell_data = parnell_data %>%
#   mutate(Timepoint = case_when(
#     Timepoint == "GD7" ~ "E7",
#     Timepoint == "GD7.25" ~ "E7.25",
#     Timepoint == "GD7.5" ~ "E7.5",
#     TRUE ~ "Eh?"
#   ))
# 
# parnell_data = parnell_data %>%
#   mutate(Treatment = case_when(
#     Treatment == "Control" ~ "Vehicle",
#     TRUE ~ Treatment
#   )) %>%
#   # mutate(Gene = as.factor(Gene),
#   #        Timepoint = as.factor(Timepoint),
#   #        Strain = as.factor(Strain),
#   #        Treatment = as.factor(Treatment)) %>%
#   identity()
# write_rds(parnell_data,
#           here('Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_forShinyApp.rds'),
#           compress="gz")

# parnell_data = read_rds(here('Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_forShinyApp.rds'))
# parnell_grouped = parnell_data %>% group_by(Gene) %>% nest()
# parnell_data_list = list()
# 
# for (i in 1:dim(parnell_grouped)[1]) {
#   parnell_data_list[[parnell_grouped$Gene[i]]] = parnell_grouped$data[i][[1]] %>%
#     mutate(Gene = parnell_grouped$Gene[i])
# }
# write_rds(parnell_data_list,here('Parnell_data_split.rds'),compress="gz")

parnell_data_full = read_rds(here('Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_table.rds')) %>%
  pivot_longer(!Gene,names_to= "data_set", values_to = "values") %>%
  mutate(Timepoint = case_when(
    str_detect(data_set,"GD7_") ~ "E7",
    str_detect(data_set,"GD7.25_") ~ "E7.25",
    str_detect(data_set,"GD7.5_") ~ "E7.5",
    TRUE ~ "Eh?"
  )) %>% 
  mutate(Strain = case_when(
    str_detect(data_set,"C57BL6J") ~ "C57BL6J",
    str_detect(data_set,"C57BL6N") ~ "C57BL6N",
    TRUE ~ "Eh?")) %>%
  mutate(Treatment = case_when(
    str_detect(data_set,"Control") ~ "Vehicle",
    str_detect(data_set,"EtOH") ~ "EtOH",
    TRUE ~ "Eh?")) %>%
  select(-data_set)

parnell_data_full_grouped = parnell_data_full %>% group_by(Gene) %>% nest()
parnell_data_list = list()

for (i in 1:dim(parnell_data_full_grouped)[1]) {
  parnell_data_list[[parnell_data_full_grouped$Gene[i]]] = parnell_data_full_grouped$data[i][[1]] %>%
    mutate(Gene = parnell_data_full_grouped$Gene[i])
}
write_rds(parnell_data_list,here('Parnell_data_full_split.rds'),compress="gz")