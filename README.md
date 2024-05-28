# master_thesis_ua
## Modelling the cost-effectiveness of RSV interventions in England and Wales. 
RSV causes acute lower respiratory tract infections (ALRI) in people of all ages, but especially in children under 5 years and the elderly. New monoclonal antibodies, offering rapid immune protection to infants, maternal vaccines, conferring protection to infants through placental transfer of antibodies, and elderly vaccines are becoming available. Recently, a maternal vaccine (Abrysvo) and a monoclonal antibody (Nirsevimab) have demonstrated higher efficacy and duration of protection than the current prophlaxis, Palivizumab in phase 3 trials and have been approved for use against RSV in infants by the FDA and EMA since September 2023. 

In this project, I used mathematical modelling to assess the impact and cost-effectiveness of these new interventions taking into account waning immunity and herd immunity using a whole-population-wide age-startified compartmental model of selected countries in England and Wales. This work is an extension of previous study by [Hodgson et al. 2020](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8) and for fullfilment of my Master of Epidemiology at the University of Antwerp. Supervision by [Professor Lander Willem](https://www.uantwerpen.be/en/staff/lander-willem/) at the Centre for Health Economics and Modelling of Infectious Diseases (CHERMID) at the University of Antwerp, Belgium. 

![RSV_model_schema.pdf](https://github.com/arumadri/master_thesis_ua/files/15475728/RSV_model_schema.pdf)

# Implementation guide 
* Quick instructions*
  Clone the repository and look through the RMarkdown files in the /vingettes folder for implementation steps.

# Overview of files 
## /data folder 
+ posteriors.Rda, is dataframe which contains posterior distributions from Hodgson et al 2020.
+ uk_data_sum.RData, is a dataframe which contains information on the population of England and Wales used for generating the initial states.
+ states.csv, is a dataframe which contains the generated initial states used in this study. The procedure is described in the /vignettes folder.

## /figs folder 
+ contains main and supplementary plots used in this study.

## /functions folder 
+ RunInterventions.R, contains all the funtions used for the analysis, generating tables and figures.
+ rsv_model_merged.R, is the age-startified RSV model used in the analysis.

## /src folder 
+ aging_check.R, a walkthrough to assess the implementation of the aging process in the model.
+ check_health.R, a walkthrough to assess the implementation of the epidemiological process of the model.
+ figure_3.R, plotting function for figure 3 in the main outputs.
+ generate_initial_states.R, a walkthrough of generating the initial states.
+ incidence_comparison.R, plotting function for figure 1 in the main outputs.
+ scenarios_cost_effect.R, a walkthrough of estimating the impact of the intervention strategies and their cost-effectiveness
+ tables.R, a walkthrough of generating all tables in the study.

## Contact details 
Please email [varumadri@gmail.com](varumadri@gmail.com) with any queries relating to this code.
