# Somatic-IV analysis: Unveiling Candidate Drivers in Cancer Progression
Identifying causal drivers of cancer progression is crucial for developing effective anti-cancer therapies. However, disentangling causality from correlation remains challenging due to confounding factors within the complex genetic and transcriptomic landscapes of cancer.
To address this challenge, we introduce the Somatic Instrumental Variable analysis (Somatic-IV), which integrates genetic, transcriptomic, and clinical outcome data to identify candidate driver genes likely playing a causal role in disease progression. Somatic-IV estimates genetic-exposure and genetic-outcome associations, utilizing MR-Egger regression to estimate bias-reduced causal effects.
  
## Vignettes for running the Somatic-IV analysis
[Somatic_IV_vignettes.ipynb](https://github.com/wzhang1984/Somatic_IV_package/blob/main/Somatic_IV_vignettes.ipynb)  

## Running environment    
#### Docker container:
- quay.io/jupyter/datascience-notebook:2024-05-27  
#### python packages:  
- networkx  
- openpyxl  
- adjustText  
- lxml  
#### R packages  
- timereg
