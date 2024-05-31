# packages

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)

pacman::p_load(knitr, kableExtra, grid, gridExtra, patchwork)

# table 1
table_1 <- data.frame(
   x = c("μ", "pR","ξ","ω","σ","γ<sub>0</sub>","γ<sub>1</sub>","γ<sub>2</sub>","α","b<sub>1</sub>","φ","ψ","δ<sub>0</sub>","δ<sub>1</sub>","δ<sub>2</sub>","δ<sub>3</sub>","p","",""),
  Parameter = c("Daily number of live births", "Proportion of neonates born with protection", "Average duration of maternally derived immunity (days)", 
                "Average of disease-induced/post-infection immunity (days)", "Average duration of exposure (days)", "Average duration of primary infection (days)", 
                "Decrease in secondary infection duration relative to primary ", "Decrease in subsequent infection duration relative to secondary", 
                "Relative reduction in infectiousness for asymptomatic infections ", "Relative amplitude of transmission during peak", "Seasonal shift in transmission",
                "Seasonality wavelength constant","Reduction in primary infection (relative to primary)","Reduction in secondary infection (relative to primary)", "Reduction in tertiary infection (relative to secondary)",
                "Reduction in subsequent infections (relative to tertiary)", "Probability infection is asymptomatic [0-4 years]:","Probability infection is asymptomatic [5-64 years]: 
                ","Probability infection is asymptomatic [65+ years]:"),
  Value = c("1863", "0.5", "133.5", "358.9", "4.98", "6.062", "0.878", "0.812", "0.634", "1.998", "0.614","0.236","1","0.89","0.81","0.33","0.255", "0.635", "0.753"),
  Reference = c("<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a> ", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a> ",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a> ", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>")
)

table_1 <- kable(table_1, "html", booktabs = TRUE,escape = FALSE, na.action = "insert", col.names = c("", "Parameter", "Value", "Reference")) %>%
  kable_styling(font_size = 11, full_width = FALSE ,bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold = TRUE, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(2,4,6,8,10,12,14,16,18), color = "black", background = "#fde4e2") %>% 
  row_spec(19, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(1,3,5,7,9,11,13,15,17,19), color = "black", background = "white") %>% 
  column_spec(1, width = "10%", border_left = TRUE) %>% 
  column_spec(2, width = "50%") %>%  
  column_spec(3, width = "10%") %>%  
  column_spec(4, width = "10%", border_right = TRUE) %>%
  footnote(general = "", 
           threeparttable = TRUE,
           general_title = "Table 1: Model parameters used in the transmission model of RSV",
           title_format = "bold") 

# save by copying to clipboard or by zooming out in viewer and taking a screenshot 
table_1

# table 2
table_2 <- data.frame(
  Parameter = c("<b>Costs per outcome</b>", "GP consultation", "Hospital admission episode", "Hospital admission episode ",
                "<b>Life years gained</b>", "0-4 years", 
                "<b>Probability of clinical outcomes</b>", "Per-infection probability of GP consultation","Per-infection probability of GP consultation" ,"Per-infection probability of death", "Per-infection probability of death",
                "Per-infection probability of hospital admissions","Per-infection probability of hospital admissions","Per-infection probability of hospital admissions"),
  Value = c(" ", "£36.00 (all ages)", " <5 years: £1100.23", " ≥15 years: £652.29 ", " ", "85", " ", "<5 years: 0.006 ", "≥5 years: 0.016",
            "<5 years: 8.197 x 10^-6", "≥5 years: 5.697 x 10^-6","<5 years: 0.004", "5-64 years: 4.688 x 10^-5","65+ years: 6.197 x 10^-5"),
  Reference = c(" ", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", "<a href= 'https://www.thelancet.com/journals/lanepe/article/PIIS2666-7762(23)00248-X/fulltext'>25</a> ", 
                "<a href= 'https://www.thelancet.com/journals/lanepe/article/PIIS2666-7762(23)00248-X/fulltext'>25</a>","", "<a href= 'https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/articles/lifeexpectancycalculator/2019-06-07'>24</a> ",
                "", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>")
)
table_2$Value <- gsub("<", "&lt;", table_2$Value, fixed = TRUE) 


table_2 <- kable(table_2, "html", booktabs = TRUE,escape = FALSE, na.action = "insert") %>%
  kable_styling(font_size = 11, full_width = FALSE ,bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold = TRUE, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(14, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(2,4,6,8,10,12,14), color = "black", background = "#fde4e2") %>% 
  row_spec(c(1,3,5,7,9,11,13), color = "black", background = "white") %>% 
  column_spec(1, width = "33%", border_left = TRUE) %>%  
  column_spec(2, width = "33%", ) %>%  
  column_spec(3, width = "33%", border_right = TRUE) %>%
  footnote(general = "", 
           threeparttable = TRUE,
           general_title = "Table 2: Intervention model and Economic parameters used in this study",
           title_format = "bold") 
# save by copying to clipboard or by zooming out in viewer and taking a screenshot 
table_2

# table 3
table_2_1 <- data.frame(
  x = c("T<sub>cea</sub>", "T<sub>mort</sub>","<b><i>Palivizumab</i></b>", "ω<sub>pal</sub>","e<sub>pal</sub>","uptake<sub>pal</sub>", "admin<sub>pal</sub>", "cost<sub>pal</sub>","<b><i>Nirservimab</i></b>","ω<sub>lmab</sub>","e<sub>lmab</sub>","uptake<sub>lmab</sub>", "admin<sub>lmab</sub>", "cost<sub>lmab</sub>", 
     "<b><i>Abrysvo/maternal vaccine</i></b>", "ω<sub>mat</sub>","e<sub>mat</sub>","uptake<sub>mat</sub>", "admin<sub>mat</sub>", "cost<sub>mat</sub>","<b><i>Abrysvo/elderly vaccine</i></b>",
     "ω<sub>65+</sub>","e<sub>mat</sub>","uptake<sub>65+</sub>", "admin<sub>65+</sub>", "cost<sub>65+</sub>"),
  Parameter = c("Time horizon (years) [CEA]", "Time Horizon (years) [Mortality]","", "Average period of protection (days)", "Efficacy against all infections (%)", "Coverage (%)", "Administration cost (per course)","Purchasing price (per course)","", 
                "Average period of protection (days)", "Efficacy against all infections (%)", "Coverage (%)","Administration cost (per course)","Purchasing price (per course)", "", 
                "Average period of protection (days)", "Efficacy against all infections (%)","Coverage (%)","Administration cost (per course)","Purchasing price (per course)", "",
                "Average period of protection (days)","Efficacy against all infections (%)", "Coverage (%)","Administration cost (per course)","Purchasing price (per course)"),
  Value = c("1", "Life long","", "30", "33.80", "90", "£57.50","£4035.50", "", "150", "78.10", "90", "£11.00", "£238","","180","I: 51.30 M: 71.70","60","£9.00","£238","","180","66.70","70","£9.00","£238"),
  Reference = c("", "","","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", 
              "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", 
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>", "",
                "<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2110275'>9</a>","<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2110275'>9</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://hhpharmacy.co.uk/price-list/'>26</a>",
                "","<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2216480'>10</a>","<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2216480'>10</a>",
                "<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://hhpharmacy.co.uk/price-list/'>26</a>",
                "","<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2213836'>11</a>","<a href= 'https://www.nejm.org/doi/10.1056/NEJMoa2213836'>11</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01802-8'>15</a>","<a href= 'https://hhpharmacy.co.uk/price-list/'>26</a>"),
  stringsAsFactors = FALSE
)

table_2_1$Value <- gsub("<", "&lt;", table_2_1$Value, fixed = TRUE) 

table_2_1 <- kable(table_2_1, "html", booktabs = TRUE,escape = FALSE, na.action = "insert", col.names = c("", "Parameter", "Value", "Reference")) %>%
  kable_styling(font_size = 11, full_width = FALSE ,bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold = TRUE, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(2,4,6,8,10,12,14,16,18,20,22,24,26), color = "black", background = "#fde4e2") %>% 
  row_spec(26, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(1,3,5,7,9,11,13,15,17,19,21,23,25), color = "black", background = "white") %>% 
  column_spec(1, width = "20%", border_left = TRUE) %>% 
  column_spec(2, width = "40%") %>%  
  column_spec(3, width = "20%") %>%  
  column_spec(4, width = "20%", border_right = TRUE) %>%
  footnote(general = "", 
           threeparttable = TRUE,
           general_title = "Table 2: Intervention model and Economic parameters used in this study",
           title_format = "bold") 
# save by copying to clipboard or by zooming out in viewer and taking a screenshot 
table_2_1

## table 3

table_3_averted <- data.frame(
  Intervention = c("<b><i>Base scenario<i></b>", "","","","<b><i>Palivizumab</i></b>", "","", "", "<b><i>Niservimab</i></b>", "", "",
                   "", "<b><i>Abrysvo/maternal vaccine</i></b>", "", "", "","<b><i>Abrysvo/elderly vaccine</i></b>",
                   "", "", ""),
  'Age group' = c("","0-4 years", "5-64 years","65+",
                  "","0-4 years", "5-64 years","65+", 
                  "","0-4 years", "5-64 years", "65+",
                  "","0-4 years","5-64 years", "65+",
                  "","0-4 years","5-64 years",  "65+"),
  'Base cases' = c("","1621452", "420636", "9426", 
                     "","", "", "",
                     "","", "", "",
                     "","", "", "",
                     "","","",""),
  
  'Cases averted' = c("","0","0","0",
                      "","1004", "66", "0", 
                      "","12958", "657","-8",
                      "","4410", "-41", "-1",
                      "","-6", "-2","60"),
  
  'GP consultations averted' = c("","0","0","0",
                                 "","6",  "1",  "0", 
                                 "","78", "11", "0",
                                 "","27",  "0", "0",
                                 "","0", "0", "1"),
  
  'Hospital admissions averted' = c("","0","0","0",
                                    "","4","0", "0",
                                    "","52", "0", "0",
                                    "","18", "0", "0",
                                    "","0", "0", "0"),
  
  'Deaths averted' = c("","0","0","0",
                       "","0", "0", "0",
                       "","1", "0", "0",
                       "","0", "0", "0",
                       "","0", "0","0") 
)

table_3_averted <- kable(table_3_averted, "html", booktabs = TRUE,escape = FALSE, na.action = "insert", col.names = c("Intervention", "Age group", "Base cases", "Cases averted", "GP consultations averted", "Hospital admissions averted", "Deaths averted")) %>%
  kable_styling(font_size = 11, full_width = FALSE ,bootstrap_options = c("striped", "hover")) %>%
  row_spec(0, bold = TRUE, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(2,4,6,8,10,12,14,16,18,20), color = "black", background = "#fde4e2") %>% 
  row_spec(20, color = "black", background = "#fde4e2", extra_css = "border-bottom: 1px solid #000;") %>% 
  row_spec(c(1,3,5,7,9,11,13,15,17), color = "black", background = "white") %>% 
  column_spec(1, width = "20.67%%", border_left = TRUE) %>% 
  column_spec(2, width = "20.67%%") %>% 
  column_spec(3, width = "16.67%%") %>%  
  column_spec(4, width = "16.67%%") %>% 
  column_spec(5, width = "16.67%%") %>% 
  column_spec(6, width = "16.67%%") %>%  
  column_spec(7, width = "16.67%%%", border_right = TRUE) %>%
  footnote(general = "", 
           threeparttable = TRUE,
           general_title = "Table 3: Seasonal number of RSV cases, GP consulations, hospital admissions and deaths averted per intervention by age group",
           title_format = "bold") 

# save by copying to clipboard or by zooming out in viewer and taking a screenshot 
table_3

