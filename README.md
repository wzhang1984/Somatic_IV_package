# Somatic-IV analysis: Unveiling Candidate Drivers in Cancer Progression
Identifying causal drivers of cancer progression is crucial for developing effective anti-cancer therapies. However, disentangling causality from correlation remains challenging due to confounding factors within the complex genetic and transcriptomic landscapes of cancer.
To address this challenge, we introduce the Somatic Instrumental Variable analysis (Somatic-IV), which integrates genetic, transcriptomic, and clinical outcome data to identify candidate driver genes likely playing a causal role in disease progression. Somatic-IV estimates genetic-exposure and genetic-outcome associations, utilizing MR-Egger regression to estimate bias-reduced causal effects. ![image](https://github.com/user-attachments/assets/666de51a-e6fe-40fd-a74d-d98cda202573)
  
## Viegnettes for running the Somatic-IV analysis
Somatic_IV_vignettes.ipynb  

## Running environment    
__Docker container:__  
- quay.io/jupyter/datascience-notebook:2024-05-27  
__python packages:__  
- networkx  
- openpyxl  
- adjustText  
- lxml  
__R packages__  
- timereg
