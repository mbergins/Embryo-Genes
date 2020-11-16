# library(crayon)
# library(vctrs)
# library(backports)
# library(rstudioapi)
# library(cli)
# library(withr)
library(tidyverse)
# library(magrittr)
# library(patchwork)
# library(labeling
# library(farver,lib="/proj/jmsimon/Rlibs36")
# library(digest,lib="/proj/jmsimon/Rlibs36")
# library(utf8,lib="/proj/jmsimon/Rlibs36")
# library(rsample,lib="/proj/jmsimon/Rlibs36")
# library(ellipsis,lib="/proj/jmsimon/Rlibs36")
# library(readxl,lib="/proj/jmsimon/Rlibs36")
# library(rematch,lib="/proj/jmsimon/Rlibs36")

library(ggplot2)
library(here)

dodge <- position_dodge(width=0.9)

combined <- readRDS(here("Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_forShinyApp.rds"))

combined %>%
	filter(Gene=="Shh") %>%
	group_by(Timepoint,Strain,Treatment) %>%
	ggplot(aes(x=Strain,y=Mean,fill=Treatment)) +
	geom_bar(stat="identity",position="dodge") +
	geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),width=0.2,position=dodge) +
	facet_wrap(~Timepoint) +
	ggtitle("Shh") +
	ylim(-3,11)


# For download of full gene expression table

# If full table is needed as an R object
full <- readRDS("Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_table.rds")

# Otherwise a csv.gz file exists in the same directory named
# Parnell_GD7_GD7.25_GD7.5_EtOH_scaledVST_table.csv.gz
