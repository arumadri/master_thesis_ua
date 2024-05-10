library(knitr)
library(kableExtra)
library(grid)
library(gridExtra)
library(patchwork)

# table 1
table_1 <- data.frame(
   x = c("μ", "pR","ξ","ω","σ","γ0","g1","g2","α","b1","φ","ψ","δ1","δ2","δ3","p","",""),
  Parameter = c("Daily number of live births", "Proportion of neonates born with protection", "Average duration of maternally derived immunity (days)", 
                "Average of disease-induced/post-infection immunity (days)", "Average duration of exposure (days)", "Average duration of primary infection (days)", 
                "Decrease in secondary infection duration relative to primary ", "Decrease in subsequent infection duration relative to secondary", 
                "Relative reduction in infectiousness for asymptomatic infections ", "Relative amplitude of transmission during peak", "Seasonal shift in transmission",
                "Seasonality wavelength constant","Reduction in secondary infection (relative to primary)", "Reduction in tertiary infection (relative to secondary)",
                "Reduction in subsequent infections (relative to tertiary)", "Probability infection is asymptomatic [0-4 years]:","[5-64 years]: ","[65+]:"),
  Value = c("1863", "0.5", "133.5", "358.9", "4.98", "6.062", "0.878", "0.812", "0.634", "1.998", "0.614","0.236","0.89","0.81","0.33","0.255", "0.635", "0.753"),
  Reference = c("<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a> ", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a> ",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a> ", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>")
)

table_1 <- kable(table_1, "html", booktabs = TRUE,escape = FALSE, na.action = "insert", col.names = c("", "Parameter", "Value", "Reference")) %>%
  kable_styling(font_size = 11, full_width = FALSE ,bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold = TRUE, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(2,4,6,8,10,12,14,16,18), color = "black", background = "#fde4e2") %>% 
  row_spec(18, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(1,3,5,7,9,11,13,15,17), color = "black", background = "white") %>% 
  column_spec(1, width = "10%", border_left = TRUE) %>% 
  column_spec(2, width = "50%") %>%  
  column_spec(3, width = "10%") %>%  
  column_spec(4, width = "10%", border_right = TRUE) %>%
  footnote(general = "", 
           threeparttable = TRUE,
           general_title = "Table 1: Model parameters used in the transmission model of RSV",
           title_format = "bold") 

# save by copying to clipboard 
table_1


# table 2
table_2 <- data.frame(
  Parameter = c("<b>Costs per outcome</b>", "GP consultation", "Hospital admission episode", " ",
                "<b>Life years gained</b>", "0-4 years", "65+ years ",
                "<b>Probability of clinical outcomes</b>", "Per-infection probability of GP consultation","" ,"Per-infection probability of death", "",
                "Per-infection probability of hospital admissions","",""),
  Value = c(" ", "£36.00 (all ages)", " <5 years: £1100.23", " ≥15 years: £652.29 ", " ", "85", " 13", " ", "<5 years: 0.006 ", "≥5 years: 0.016",
            "<5 years: 8.197 x 10^-6", "≥5 years: 5.697 x 10^-6","<5 years: 0.004", "5-64 years: 4.688 x 10^-5","65+ years: 6.197 x 10^-5"),
  Reference = c(" ", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>", "<a href= 'https://www.thelancet.com/journals/lanepe/article/PIIS2666-7762(23)00248-X/fulltext'>22</a> ", 
                "<a href= 'https://www.thelancet.com/journals/lanepe/article/PIIS2666-7762(23)00248-X/fulltext'>22</a>","", "<a href= 'https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/articles/lifeexpectancycalculator/2019-06-07'>21</a> ","<a href= 'https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/articles/lifeexpectancycalculator/2019-06-07'>21</a> ",
                "", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>")
)
table_2$Value <- gsub("<", "&lt;", table_2$Value, fixed = TRUE) 


table_2 <- kable(table_2, "html", booktabs = TRUE,escape = FALSE, na.action = "insert") %>%
  kable_styling(font_size = 11, full_width = FALSE ,bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold = TRUE, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(15, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(2,4,6,8,10,12,14), color = "black", background = "#fde4e2") %>% 
  row_spec(c(1,3,5,7,9,11,13,15), color = "black", background = "white") %>% 
  column_spec(1, width = "33%", border_left = TRUE) %>%  
  column_spec(2, width = "33%", ) %>%  
  column_spec(3, width = "33%", border_right = TRUE) %>%
  footnote(general = "", 
           threeparttable = TRUE,
           general_title = "Table 2: Intervention model and Economic parameters used in this study",
           title_format = "bold") 
table_2

# table 3
table_2_1 <- data.frame(
  x = c("T", "<b><i>Palivizumab</i></b>", "ωpal","epal","cpal", "padmin", "pcost","<b><i>Nirservimab</i></b>","ωlmab","elmab","clmab","lmabadmin", 
        "lmabcost","<b><i>Abrysvo/maternal vaccine</i></b> ","ωmat","emat","cmat", "matadmin", "matcost","<b><i>Abrysvo/elderly vaccine</i></b>",
        "ω65+","e65+","c65+","65+admin","65+cost"),
  Parameter = c("Time horizon (years)", "", "Average period of protection (days)", "Efficacy against all infections (%)", "Coverage (%)", "Administration cost (per course)","Purchasing price (per course)","", 
                "Average period of protection (days)", "Efficacy against all infections (%)", "Coverage (%)","Administration cost (per course)","Purchasing price (per course)", "", 
                "Average period of protection (days)", "Efficacy against all infections (%)","Coverage (%)","Administration cost (per course)","Purchasing price (per course)", "",
                "Average period of protection (days)","Efficacy against all infections (%)", "Coverage (%)","Administration cost (per course)","Purchasing price (per course)"),
  Value = c("10","", "30", "33.80", "90", "£57.50","£4035.50", "", "150", "78.10", "90", "£11.00", "£238","","180","I: 51.30 M: 71.70","60","£9.00","£238","","180","66.70","70","£9.00","£238"),
  Reference = c("", "","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>13</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>13</a>", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>13</a>", 
              "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>13</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>13</a>", "",
                "<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2110275'>9</a>","<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2110275'>9</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>13</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>13</a>","<a href= 'https://hhpharmacy.co.uk/price-list/'>23</a>",
                "","<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2216480'>10</a>","<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2216480'>10</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>13</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>12</a>","<a href= 'https://hhpharmacy.co.uk/price-list/'>23</a>",
                "","<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2213836'>11</a>","<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2213836'>11</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>13</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>13</a>","<a href= 'https://hhpharmacy.co.uk/price-list/'>23</a>"),
  stringsAsFactors = FALSE
)

table_2_1$Value <- gsub("<", "&lt;", table_2_1$Value, fixed = TRUE) 

table_2_1 <- kable(table_2_1, "html", booktabs = TRUE,escape = FALSE, na.action = "insert", col.names = c("", "Parameter", "Value", "Reference")) %>%
  kable_styling(font_size = 11, full_width = FALSE ,bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold = TRUE, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(2,4,6,8,10,12,14,16,18,20,22,24), color = "black", background = "#fde4e2") %>% 
  row_spec(25, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(1,3,5,7,9,11,13,15,17,19,21,23,25), color = "black", background = "white") %>% 
  column_spec(1, width = "20%", border_left = TRUE) %>% 
  column_spec(2, width = "40%") %>%  
  column_spec(3, width = "20%") %>%  
  column_spec(4, width = "20%", border_right = TRUE) %>%
  footnote(general = "", 
           threeparttable = TRUE,
           general_title = "Table 2: Intervention model and Economic parameters used in this study",
           title_format = "bold") 
table_2_1

## table 1 supplementary
table_1_sup <- data_frame(
  Intervention = c("<b><i>Palivizumab</i></b>", "","", "", "<b><i>Niservimab</i></b>", "", "",
                   "", "<b><i>Abrysvo/maternal vaccine</i></b>", "", "", "","<b><i>Abrysvo/elderly vaccine</i></b>",
                   "", "", ""),
  'Age group' = c("","0-4 years", "5-64 years","65+", 
                  "","0-4 years", "5-64 years", "65+",
                  "","0-4 years","5-64 years", "65+",
                  "","0-4 years","5-64 years",  "65+"),
  
  'Cases averted' = c("","9256", "1794", "44", 
                      "","102841", "17951","272",
                      "","37258", "3077", "114",
                      "","84", "58","332"),
  
  'GP consultation' = c("","56",  "29",  "1", 
                        "","617", "287", "4",
                        "","224",  "48", "2",
                        "","0", "1", "5"),
  
  'Hospital admissions' = c("","37","1", "0",
                            "","411", "8", "0",
                            "","149", "1", "0",
                            "","0", "0", "0"),
  
  Death = c("","1", "0", "0",
            "","8", "1", "0",
            "","3", "0", "0",
            "","0", "0","0") 
)

table_1_sup <- kable(table_1_sup, "html", booktabs = TRUE,escape = FALSE, na.action = "insert", col.names = c("Intervention", "Age group", "Cases averted", "GP consultation", "Hospital admissions", "Death")) %>%
  kable_styling(font_size = 11, full_width = FALSE ,bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold = TRUE, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(2,4,6,8,10,12,14,16), color = "black", background = "#fde4e2") %>% 
  row_spec(16, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(1,3,5,7,9,11,13,15), color = "black", background = "white") %>% 
  column_spec(1, width = "20.67%", border_left = TRUE) %>% 
  column_spec(2, width = "20.67%") %>% 
  column_spec(3, width = "16.67%") %>% 
  column_spec(4, width = "16.67%") %>% 
  column_spec(5, width = "16.67%%") %>%  
  column_spec(6, width = "10%", border_right = TRUE) %>%
  footnote(general = "", 
           threeparttable = TRUE,
           general_title = "Table 3: Number of RSV cases, GP consulations, hospital admissions and deaths averted per intervention by age group",
           title_format = "bold") 

table_1_sup

## averted symptomatic/GP probability per-infection 
# <5 years: 0.006
# ≥5 years: 0.016

## averted hospitalization probability per-infection 
# <5 years: 0.004
# 5-64 years: 4.688 x 10^-5
# 65+ years: 6.197 x 10^-5

## averted death probability per-infection 
# <5 years: 8.197 x 10^-6
# ≥5 years: 5.697 x 10^-6

