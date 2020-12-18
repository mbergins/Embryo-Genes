library(tidyverse)
library(here)

parnell_data = read_rds(here('Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_forShinyApp.rds'))
parnell_data = parnell_data %>%
  mutate(Timepoint = case_when(
    Timepoint == "GD7" ~ "E7.0",
    Timepoint == "GD7.25" ~ "E7.25",
    Timepoint == "GD7.5" ~ "E7.5",
    TRUE ~ "Eh?"
  )) %>%
  mutate(Treatment = case_when(
    Treatment == "Control" ~ "Vehicle",
    Treatment == "EtOH" ~ "PAE",
    TRUE ~ Treatment
  )) %>%
  identity() %>%
  write_rds(here('Parnell_data_reformated.rds'), compress = "gz")

parnell_grouped = parnell_data %>% group_by(Gene) %>% nest()
parnell_data_list = list()

for (i in 1:dim(parnell_grouped)[1]) {
  parnell_data_list[[parnell_grouped$Gene[i]]] = parnell_grouped$data[i][[1]] %>%
    mutate(Gene = parnell_grouped$Gene[i])
}
write_rds(parnell_data_list,here('Parnell_data_split.rds'),compress="gz")

parnell_data_full = read_rds(here('Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_table.rds')) %>%
  pivot_longer(!Gene,names_to= "data_set", values_to = "values") %>%
  mutate(Timepoint = case_when(
    str_detect(data_set,"GD7_") ~ "E7.0",
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
    str_detect(data_set,"EtOH") ~ "PAE",
    TRUE ~ "Eh?")) %>%
  select(-data_set)

parnell_data_full_grouped = parnell_data_full %>% group_by(Gene) %>% nest()
parnell_data_list = list()

for (i in 1:dim(parnell_data_full_grouped)[1]) {
  parnell_data_list[[parnell_data_full_grouped$Gene[i]]] = parnell_data_full_grouped$data[i][[1]] %>%
    mutate(Gene = parnell_data_full_grouped$Gene[i])
}
write_rds(parnell_data_list,here('Parnell_data_full_split.rds'),compress="gz")

lower_case_match = list()
for (gene in names(parnell_data_list)) {
    lower_case_match[[tolower(gene)]] = gene
}
write_rds(lower_case_match,here('lower_case_match.rds'))

# library(tictoc)
# tic()
# cached_graphs = list()
# for (this_gene in names(Parnell_data_split)) {
#   cached_graphs[[this_gene]] = ggplot(Parnell_data_split[[this_gene]], aes(x=Strain,y=Mean,fill=Treatment)) +
#     geom_bar(stat="identity",position="dodge",alpha=0.5) +
#     geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),width=0.2,position=position_dodge(width=0.9)) +
#     geom_point(data=Parnell_data_full_split[[this_gene]], 
#                mapping=aes(x=Strain,y=values,color=Treatment), 
#                position=position_jitterdodge(jitter.width=0.4)) +
#     facet_wrap(~Timepoint) +
#     ggtitle(this_gene) +
#     theme(text = element_text(size=20)) +
#     ylim(c(-3,11))
# }
# write_rds(cached_graphs, here('cached_graphs.rds'), compress="gz");
# toc()

# cols <- c("C57BL6J\nPAE" = "grey", 
#           "C57BL6J\nVehicle" = "darkgreen", 
#           "C57BL6N\nPAE" = "grey", 
#           "C57BL6N\nVehicle" = "darkblue")
# 
# parnell_data_list[["Wdfy1"]] %>%
#   mutate(str_treat = paste0(Strain,"\n",Treatment)) %>%
#   ggplot(aes(x=Strain,y=Mean,color=str_treat, fill=str_treat)) +
#       geom_bar(stat="identity",position="dodge",alpha=0.5) +
#       geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),width=0.2,position=position_dodge(width=0.9)) +
#       geom_point(data=parnell_data_full_list[["Wdfy1"]] %>% mutate(str_treat = paste0(Strain,"\n",Treatment)),
#                  mapping=aes(x=Strain,y=values,color=str_treat, fill = str_treat),
#                  position=position_jitterdodge(jitter.width=0.25)) +
#       labs(y="Mean-centered expression", 
#            color = "Strain\nTreatment",
#            fill = "Strain\nTreatment") +
#       facet_wrap(~Timepoint) +
#       ggtitle('Wdfy1') +
#       theme(text = element_text(size=20)) + 
#       theme(panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank(),
#             panel.background = element_blank(), 
#             axis.line = element_line(colour = "black")) +
#       scale_colour_manual(values = cols) +
#       scale_fill_manual(values = cols) +
#       ylim(c(-3,11))
# # ggsave(here('sample_four_color_plot.png'))