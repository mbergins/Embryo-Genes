library(crayon,lib="/proj/jmsimon/Rlibs36")
library(vctrs,lib="/proj/jmsimon/Rlibs36")
library(backports,lib="/proj/jmsimon/Rlibs36")
library(rstudioapi,lib="/proj/jmsimon/Rlibs36")
library(cli,lib="/proj/jmsimon/Rlibs36")
library(withr,lib="/proj/jmsimon/Rlibs36")
library(tidyverse,lib="/proj/jmsimon/Rlibs36")
library(magrittr,lib="/proj/jmsimon/Rlibs36")
library(patchwork,lib="/proj/jmsimon/Rlibs36")
library(labeling,lib="/proj/jmsimon/Rlibs36")
library(farver,lib="/proj/jmsimon/Rlibs36")
library(digest,lib="/proj/jmsimon/Rlibs36")
library(utf8,lib="/proj/jmsimon/Rlibs36")
library(rsample,lib="/proj/jmsimon/Rlibs36")
library(ellipsis,lib="/proj/jmsimon/Rlibs36")
library(readxl,lib="/proj/jmsimon/Rlibs36")
library(rematch,lib="/proj/jmsimon/Rlibs36")


setwd("/proj/jmsimon/Parnell/GD7_GD725_GD75_B6N_B6J_EtOH_forOnlineTool")

dodge = position_dodge(width=0.9)

GD7_Veh_J_vs_N = read_excel("Parnell_GD7_rerun_DESeq2_results_summary.xlsx",sheet=1,col_types=c("text","skip","skip","skip","skip","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

global_mean_GD7 = GD7_Veh_J_vs_N %>% 
	select(starts_with("GD7")) %>% 
	pivot_longer(cols=starts_with("GD7")) %>% 
	summarize(Mean = mean(value)) %>%
	pull(Mean)

GD7_avgs = GD7_Veh_J_vs_N %>% 
	mutate_at(vars(starts_with("GD7")),function(x){x - global_mean_GD7}) %>%
	pivot_longer(cols=starts_with("GD7")) %>%
	mutate(name=str_replace(name,"\\.","_")) %>%
	separate(name,into=c("Timepoint","Replicate","Strain","Pool"),sep="_") %>%
	select(!Pool) %>%
	mutate(Replicate = as.integer(Replicate)) %>%
	mutate(Replicate=ifelse(Replicate>6,yes=Replicate-6,no=Replicate)) %>%
	add_column(Treatment="Control") %>%
	group_by(Gene,Timepoint,Strain,Treatment) %>%
	summarize(Mean=mean(value),SE=sd(value)/sqrt(n()))



GD725_all = read_excel("Parnell_GD7_rerun_DESeq2_results_summary.xlsx",sheet=2,col_types=c("text","skip","skip","skip","skip",rep("numeric",24)))

global_mean_GD725 = GD725_all %>% 
	select(starts_with("C57")) %>% 
	pivot_longer(cols=starts_with("C57")) %>% 
	summarize(Mean = mean(value)) %>%
	pull(Mean)

GD725_avgs = GD725_all %>% 
	mutate_at(vars(starts_with("C57")),function(x){x - global_mean_GD725}) %>%
	pivot_longer(cols=starts_with("C57")) %>%
	separate(name,into=c("Strain","Treatment","Replicate"),sep="_") %>%
	mutate(Replicate = str_replace(Replicate,"Rep","")) %>%
	mutate(Replicate = as.integer(Replicate)) %>%
	add_column(Timepoint="GD7.25") %>%
	mutate(Treatment = str_replace(Treatment,"LRS","Control")) %>%
	group_by(Gene,Timepoint,Strain,Treatment) %>%
	summarize(Mean=mean(value),SE=sd(value)/sqrt(n()))




GD75_all = read_excel("Parnell_GD7_rerun_DESeq2_results_summary.xlsx",sheet=5,col_types=c("text","skip","skip","skip","skip",rep("numeric",21)))

global_mean_GD75 = GD75_all %>% 
	select(starts_with("C57")) %>% 
	pivot_longer(cols=starts_with("C57")) %>% 
	summarize(Mean = mean(value)) %>%
	pull(Mean)

GD75_avgs = GD75_all %>% 
	mutate_at(vars(starts_with("C57")),function(x){x - global_mean_GD75}) %>%
	pivot_longer(cols=starts_with("C57")) %>%
	separate(name,into=c("Strain","Treatment","Replicate","Date"),sep="_") %>%
	select(!Date) %>%
	mutate(Replicate = str_replace(Replicate,"Rep","")) %>%
	mutate(Replicate = as.integer(Replicate)) %>%
	add_column(Timepoint="GD7.5") %>%
	mutate(Treatment = str_replace(Treatment,"LRS","Control")) %>%
	group_by(Gene,Timepoint,Strain,Treatment) %>%
	summarize(Mean=mean(value),SE=sd(value)/sqrt(n()))

	

combined = bind_rows(GD7_avgs,GD725_avgs,GD75_avgs)	

saveRDS(combined,"Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_forShinyApp.rds")

genelist = c("Shh","Rara","Fgf5","Cdk4")

pdf("Parnell_GD7_GD725_GD75_EtOH_Control_B6JN_barplots_testgenes.pdf",onefile=TRUE)
for(i in 1:length(genelist)) {
	pdf = combined %>%
		filter(Gene==genelist[i]) %>%
		group_by(Timepoint,Strain,Treatment) %>%
		ggplot(aes(x=Strain,y=Mean,fill=Treatment)) +
		geom_bar(stat="identity",position="dodge") +
		geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),width=0.2,position=dodge) +
		facet_wrap(~Timepoint) +
		ggtitle(genelist[i]) +

		ylim(-3,11)
	print(pdf)
}
dev.off()





# Explore mean expression distribution
GD7_Veh_J_vs_N %>% 
	pivot_longer(cols=starts_with("GD")) %>% 
	group_by(Gene) %>% 
	summarize(Mean = mean(value)) %>% 
	ggplot(aes(x=log10(Mean))) + 
	geom_density() + 
	geom_vline(xintercept = log10(global_mean_GD7))

GD725_all %>% 
	pivot_longer(cols=starts_with("C57")) %>% 
	group_by(Gene) %>% 
	summarize(Mean = mean(value)) %>% 
	ggplot(aes(x=log10(Mean))) + 
	geom_density() + 
	geom_vline(xintercept = log10(global_mean_GD725))

GD75_all %>% 
	pivot_longer(cols=starts_with("C57")) %>% 
	group_by(Gene) %>% 
	summarize(Mean = mean(value)) %>% 
	ggplot(aes(x=log10(Mean))) + 
	geom_density() + 
	geom_vline(xintercept = log10(global_mean_GD75))




# Export one re-scaled/normalized table with all values

GD7_scaled = GD7_Veh_J_vs_N %>% 
	mutate_at(vars(starts_with("GD7")),function(x){x - global_mean_GD7}) %>%
	pivot_longer(cols=starts_with("GD7")) %>%
	mutate(name=str_replace(name,"\\.","_")) %>%
	separate(name,into=c("Timepoint","Replicate","Strain","Pool"),sep="_") %>%
	select(!Pool) %>%
	mutate(Replicate = as.integer(Replicate)) %>%
	mutate(Replicate=ifelse(Replicate>6,yes=Replicate-6,no=Replicate)) %>%
	add_column(Treatment="Control")

GD725_scaled = GD725_all %>% 
	mutate_at(vars(starts_with("C57")),function(x){x - global_mean_GD725}) %>%
	pivot_longer(cols=starts_with("C57")) %>%
	separate(name,into=c("Strain","Treatment","Replicate"),sep="_") %>%
	mutate(Replicate = str_replace(Replicate,"Rep","")) %>%
	mutate(Replicate = as.integer(Replicate)) %>%
	add_column(Timepoint="GD7.25") %>%
	mutate(Treatment = str_replace(Treatment,"LRS","Control"))
	
GD75_scaled = GD75_all %>% 
	mutate_at(vars(starts_with("C57")),function(x){x - global_mean_GD75}) %>%
	pivot_longer(cols=starts_with("C57")) %>%
	separate(name,into=c("Strain","Treatment","Replicate","Date"),sep="_") %>%
	select(!Date) %>%
	mutate(Replicate = str_replace(Replicate,"Rep","")) %>%
	mutate(Replicate = as.integer(Replicate)) %>%
	add_column(Timepoint="GD7.5") %>%
	mutate(Treatment = str_replace(Treatment,"LRS","Control"))
	
combined_scaled = bind_rows(GD7_scaled,GD725_scaled,GD75_scaled) %>%
	unite(Sample,c(Timepoint,Strain,Treatment,Replicate),sep="_") %>%
	arrange(Sample) %>%
	pivot_wider(values_from="value",names_from="Sample")
	
write_csv(combined_scaled,"Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_table.csv")
saveRDS(combined_scaled,"Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_table.rds")
