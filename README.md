# Embryo Transcriptomes Data Browser

This repository contains the code used to create the data browser available through the web at: http://parnell-lab.med.unc.edu/Embryo-Transcriptomics/. The browser was written in R Shiny and can be run locally in RStudio with this command:

```
shiny::runApp('.')
```

In addition to shiny, this application requires the tidyverse, here, ggplot2, reactable, dqshiny and shinylogs packages. The dqshiny package is only available through github: https://github.com/daqana/dqshiny.


## Citation

If you use the code or data assocaited with this repository, please cite:

Boshen et al. Transcriptomic analyses of two closely related substrains of gastrulation-stage mouse embryos with differential susceptibility to prenatal alcohol exposure, Currently Under Review
