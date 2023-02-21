#########################################################
#####            AMR Systematic Review             #####
#########################################################

# This script categorises and summarises the risk factors

## Setup

pacman::p_load(
  tidyverse,
  rio,
  janitor,
  gtsummary,
  writexl
)

setwd("C:/Users/ANSD/OneDrive - Folkehelseinstituttet/AMR/ResCan/Oversiktsstudie/GitHub/AMRreview")

## Infection

### Data management

Infection <- import("Table 1.xlsx", which = "Table 1a - Infection") %>%
  select(Authors, `Microbial etiology and resistance`, `Factors included in the final model`) %>%
  clean_names() %>%
  separate(factors_included_in_the_final_model, sep = ",", into = c("factor1", "factor2", "factor3", "factor4", "factor5", "factor6",
                                                                      "factor7", "factor8", "factor9", "factor10", "factor11", "factor12",
                                                                      "factor13", "factor14", "factor15", "factor16", "factor17", "factor18",
                                                                    "factor19", "factor20", "factor21", "factor22", "factor23", "factor24",
                                                                    "factor25", "factor26", "factor27", "factor28", "factor29", "factor30"),
           extra = "warn", fill = "warn") %>%
  pivot_longer(cols = c("factor1", "factor2", "factor3", "factor4", "factor5", "factor6",
                        "factor7", "factor8", "factor9", "factor10", "factor11", "factor12",
                        "factor13", "factor14", "factor15", "factor16", "factor17", "factor18",
                        "factor19", "factor20", "factor21", "factor22", "factor23", "factor24",
                        "factor25", "factor26", "factor27", "factor28", "factor29", "factor30")) %>%
  drop_na(value) %>%
  mutate(category = case_when(
    str_detect(value, regex("Antibiot|piperacillin|quinolone|(?<!during|by).carba|azole|TMP/SMX|cephalosporin|glycopeptide|
                            |(CDI treat)|Ciprofloxacin|Levofloxacin|Meropenem|penicillin|lactam|cefepime|Vancomycin|Daptomycin|
                            |quinolones|tazobactam|TMP-SMX|Cefpodoxime|Amoxicillin|antifungals|anti SA|antimicrobial|aminoglycoside|
                            |combined|antifungal", ignore_case = T)) ~ "Antibiotic use",

    str_detect(value, regex("Comorbid|Diabetes|Pleural|diseases|Hematological disease|distortion|Carlson|
                            |vascular disease|renal disease|Atrial fibrillation|Hypertension|Chronic liver|Elixhauser|
                            |Renal disease|heart|Charleston|ischemia|Charlson|CCI|underlying disease", ignore_case = T)) ~ "Comorbidities/clinical condition",

    str_detect(value, regex("BMI|race|(?<!\\w)Age|Gender|sex|Male|(?<!prior 1).Year|Obesity|household income|^Year|
                            |smoking", ignore_case = T)) ~ "Basic characteristics",

    str_detect(value, regex("Histology|albumin|hemoglobin|FBG|Trigylceride|cholesterol|HDL|Hemoglobin|\\b(OB)\\b|Hypoproteinemia|
                            |Procalcitonin|C Reactive Protein|Amyloid|leukocytes", ignore_case = T)) ~ "Laboratory findings (non-microbiological)",

    str_detect(value, regex("NSAID|PPI|Radiation|Proton pump|corticosteroid|PP inhibitor|antacid|glucocorticoid|Hemodialysis|
                            |steroids|radiotherapy", ignore_case = T)) ~ "Other treatments/medications",

    str_detect(value, regex("(?<!from).chemotherap|^Chemotherap|^ Chemo|Conditioning|Immunosuppress|clorafabine|G-CSF|VGPR|Rituximab|Early de-escalation|
                            |decitabine|CHOP|Clofarabine", ignore_case = T)) ~ "Chemotherapy/immunosuppressants",

    str_detect(value, regex("neutropenia|neutropaenia|Typhilitis|\\b(ANC)\\b|neutropenic", ignore_case = T)) ~ "Neutropenia",

    str_detect(value, regex("hospitali|\\bLOS\\b|ICU|ward|Hospital stay|Transfer|Hospital bed|Hospital teaching|Expected primary|
                            |prior Hospital|Lengh of stay|Excessive bed|inpatient|admission|region|location|of stay", ignore_case = T)) ~ "Hospital-related",

    str_detect(value, regex("malignancy|Cancer|tumour|AML|HCI-CI|KPS|Karnofsky|Neoplasia|lymphoma|Leukemia|Lymphomia|Site|Morphology|
                            |Differentiation|T stage|N stage|neoplasm|tumor|myeloma|fatal|remission", ignore_case = T)) ~ "Cancer-related",

    str_detect(value, regex("Gastrointestinal(?!\\s+cancer)|stool|biliary|Diarrhea|Inflammatory bowel|IBD", ignore_case = T)) ~ "Gastrointestinal",

    str_detect(value, regex("transplant|HSCT|ASCT|HCST|donor|GVHD", ignore_case = T)) ~ "Transplantation-related",

    str_detect(value, regex("transfusion|CVC|(?<!portal).vein|intravenous|Central CVP|venous|parenteral", ignore_case = T)) ~ "Intravascular access",

    str_detect(value, regex("Mechanical Ventilation|COPD|respiratory fail", ignore_case = T)) ~ "Respiratory",

    str_detect(value, regex("surg|Operation|blood loss|portal vein|invasive", ignore_case = T)) ~ "Surgery-related",

    str_detect(value, regex("urinary|ileal|indiana|neobladder|diversion type", ignore_case = T)) ~ "Urinary",

    str_detect(value, regex("bacteremia|difficile|enterococcal|maltophilia|CDI disease|infection|pneumonia|epidermidis|colonization|
                            |Culture|CMV|Herpesviridae|PTZ-R|Outpatient CDI|Type of Organism|Prior CRE|Severe CDI|CDI exposure|BSBL|ESBL|
                            |Fever days|VRE|BSI|shock|bacteriemia|pathogen|E.coli|Klebsiella|Pseudomonas|Acinetobacter|Candida|
                            |tropicalis", ignore_case = T)) ~ "Infection-related",

    TRUE ~ "Other"))

### What does the "other" category contain now?

Other_inf <- Infection %>%
  filter(category == "Other or unknown") %>%
  group_by(value) %>%
  summarise(Freq=n()) %>%
  arrange(desc(Freq))
print(Other_inf)

### Create categories of microbes

Infection <- Infection %>%
  mutate(microbe=case_when(
    str_detect(microbial_etiology_and_resistance, "negative|aeruginosa|baumannii|Escherichia coli|Stenotrophomonas|Enterobacteriaceae|Klebsiella|ESBL") ~ "Gram-negative bacteria",
    str_detect(microbial_etiology_and_resistance, "Several|several|Gram negative|47 species|Enterobacter spp.|MDRO") ~ "Several bacteria and/or fungi",
    str_detect(microbial_etiology_and_resistance, "Clostridioides difficile|C. difficile") ~ "Clostridioides difficile",
    str_detect(microbial_etiology_and_resistance,"aureus|MRSA|enterococc|vancomycin|Bacillus|epidermidis") ~ "Gram-positive bacteria",
    str_detect(microbial_etiology_and_resistance, "Aspergillus|molds|aspergillus|Molds|fungi|andida") ~ "Fungi",
    TRUE ~ "Other or unknown"
    ))

### Table with results

Infection_table <- Infection %>%
  group_by(microbe, category) %>%
  summarize(count=n()) %>%
  arrange(microbe, desc(count)) %>%
  group_by(microbe) %>%
  slice_head(n = 5)
write_xlsx(Infection_table, "Infection_table.xlsx")

### List the different risk factors in each category

Infection_riskfactors <- Infection %>%
  group_by(category) %>%
  summarise(riskfactors = paste(value, collapse = ", ")) %>%
  mutate(outcome = "Infection")

## Mortality

### Data management

Mortality <- import("Table 1.xlsx", which = "Table 1b - Death") %>%
  select(Authors, `Microbial etiology and resistance`, `Factors included in the final model`) %>%
  clean_names() %>%
  separate(factors_included_in_the_final_model, sep = ",", into = c("factor1", "factor2", "factor3", "factor4", "factor5", "factor6",
                                                                    "factor7", "factor8", "factor9", "factor10", "factor11", "factor12",
                                                                    "factor13", "factor14", "factor15", "factor16", "factor17", "factor18",
                                                                    "factor19", "factor20", "factor21", "factor22", "factor23", "factor24",
                                                                    "factor25", "factor26", "factor27", "factor28", "factor29", "factor30"),
           extra = "warn", fill = "warn") %>%
  pivot_longer(cols = c("factor1", "factor2", "factor3", "factor4", "factor5", "factor6",
                        "factor7", "factor8", "factor9", "factor10", "factor11", "factor12",
                        "factor13", "factor14", "factor15", "factor16", "factor17", "factor18",
                        "factor19", "factor20", "factor21", "factor22", "factor23", "factor24",
                        "factor25", "factor26", "factor27", "factor28", "factor29", "factor30")) %>%
  drop_na(value) %>%
  mutate(category = case_when(
    str_detect(value, regex("(?<!inactive|ongoing).Antibiot|piperacillin|quinolone|(?<!during|by).carba(!penemase|!penem resist)|azole|TMP/SMX|cephalosporin(!s resist)|glycopeptide|
                            |(CDI treat)|Ciprofloxacin|Levofloxacin|Meropenem|penicillin|lactam(!ase)|cefepime|Vancomycin(!e)|Daptomycin|
                            |quinolones|tazobactam|TMP-SMX|Cefpodoxime|Amoxicillin|antifungals|anti SA|antimicrobial|aminoglycoside|
                            |combined|antifungal|within 4 weeks|appropriate|echinocandin pre|4th generation|empirical|IIAT|anti-staph|
                            |initial adequate therapy", ignore_case = T)) ~ "Antibiotic use",

    str_detect(value, regex("Comorbid|Diabetes|Pleural|diseases|Hematological disease|distortion|Carlson|
                            |vascular disease|renal disease|Atrial fibrillation|Hypertension|Chronic liver|Elixhauser|
                            |Renal disease|heart|Charleston|ischemia|Charlson|\\bCCI|SAPS II|SOFA score|
                            |APACHE|functional capacity|pulmonary disease|liver disease|comorbidities|renal failure|
                            |kidney|multiorgan|end organ|cirhhosis|Child–Pugh|Signs of severity|hypotension|
                            |CRI non hemorragic|renal Replacement|ECOG|number of organ|renal insufficiency|
                            |consciousness|hepatic failure|organ failure|hypoxia|dementia|underlying disease", ignore_case = T)) ~ "Comorbidities/clinical condition",

    str_detect(value, regex("BMI|race|(?<!\\w)Age|Gender|sex|Male|(?<!prior 1).Year|Obesity|household income|^Year|adult|
                            |Ethnicity|smoking", ignore_case = T)) ~ "Basic characteristics",

    str_detect(value, regex("Histology|albumin|hemoglobin|FBG|Trigylceride|cholesterol|HDL|Hemoglobin|\\b(OB)\\b|Hypoproteinemia|
                            |Procalcitonin|C Reactive Protein|Amyloid|leukocytes|Intravascular Coagulation|blood pressure|Hemoglobin|
                            |platelet|Alpha-fetoprotein|Lymphocytes|PCT|T cell|prothromb|D-Dimer|creatinin|alkaline|bilirubin|
                            |blood lactate|leucocyte|CRP|\\b(AST)\\b|glucose|Nadir WBC", ignore_case = T)) ~ "Laboratory findings (non-microbiological)",

    str_detect(value, regex("NSAID|PPI|Radiation|Proton pump|corticosteroid|PP inhibitor|antacid|glucocorticoid|Hemodialysis|
                            |steroids|radiotherapy|steroid|aciclovir|anthracycline|vasopressor", ignore_case = T)) ~ "Other treatments/medications",

    str_detect(value, regex("(?<!from).chemotherap|^Chemotherap|^ Chemo|Conditioning|Immunosuppress|clorafabine|G-CSF|VGPR|Rituximab|Early de-escalation|
                            |decitabine|CHOP|Clofarabine|blast clearance|sorafenib", ignore_case = T)) ~ "Chemotherapy/immunosuppressants",

    str_detect(value, regex("neutropenia|neutropaenia|Typhilitis|\\b(ANC)\\b|neutropenic|neutrophil|monocyt", ignore_case = T)) ~ "Neutropenia",

    str_detect(value, regex("hospitali|\\bLOS\\b|ICU|ward|Hospital stay|Transfer|Hospital bed|Hospital teaching|Expected primary|
                            |prior Hospital|Lengh of stay|Excessive bed|inpatient|admission|region|location|of stay|
                            |care unit|teaching status", ignore_case = T)) ~ "Hospital-related",

    str_detect(value, regex("malignanc|Cancer|tumour|AML|HCI-CI|KPS|Karnofsky|Neoplasia|lymphoma|Leukemia|Lymphomia|Site|Morphology|
                            |Differentiation|T stage|N stage|neoplasm|tumor|myeloma|fatal|remission|refractory|recent diagnosis|
                            |Performance status|BCLC stage|ALBI grade|metastas|Myelodysplastic|blastoma|sarcoma|cytology|MDS|
                            |relapse|cholangiocarcinoma", ignore_case = T)) ~ "Cancer-related",

    str_detect(value, regex("Gastrointestinal(?!\\s+cancer)|stool|biliary|Diarrhea|Inflammatory bowel|IBD", ignore_case = T)) ~ "Gastrointestinal",

    str_detect(value, regex("transplant|HSCT|ASCT|HCST|donor|GVHD|Antigen compatibility|allograft", ignore_case = T)) ~ "Transplantation-related",

    str_detect(value, regex("transfusion|CVC|(?<!portal).vein|intravenous|Central CVP|venous|central catheterization|catheter removal|parenteral", ignore_case = T)) ~ "Intravascular access",

    str_detect(value, regex("Mechanical Ventilation|COPD|respiratory fail|respiratory|lung disease|intubation|
                            |bronchial", ignore_case = T)) ~ "Respiratory",

    str_detect(value, regex("surg|Operation|blood loss|resection|invasive", ignore_case = T)) ~ "Surgery-related",

    str_detect(value, regex("urin|ileal|indiana|neobladder|diversion type|foley catheter", ignore_case = T)) ~ "Urinary",

    str_detect(value, regex("bacteremia|bacteraemia|difficile|enterococc|maltophilia|CDI disease|infection|pneumonia|epidermidis|colonization|
                            |Culture|CMV|Herpesviridae|PTZ-R|Outpatient CDI|Organism|Prior CRE|Severe CDI|CDI exposure|BSBL|ESBL|
                            |Fever days|VRE|BSI|shock|bacteriemia|pathogen|E.coli|Klebsiella|Pseudomonas|Acinetobacter|Candida|
                            |tropicalis|bacteria|CRGNB|Enterobacter|carbapenem resistance|fungi|polymic|organisms|
                            |Infected|sepsi|albicans|Biofilm|resistance|phenotype|fever|caspofungin resistance|serotype|P.aeruginosa|
                            |\\bcons\\b|Gram-positive|MDR|clinical source|fungemia|MASCC|staphylococc|bacilli|
                            |beta lactamase|source|Candidemia|\\bHIV\\b|Amphotericin|aspergillosis|Pitt|endocarditis|A.baumanii|
                            |parapsilosis|\\bCRAB\\b|carba-non sensitive|CR-PA|septic", ignore_case = T)) ~ "Infection-related",

    TRUE ~ "Other"))


### Create categories of microbes

Mortality <- Mortality %>%
  mutate(microbe=case_when(
    str_detect(microbial_etiology_and_resistance, "Several|several|Gram negative|47 species|Enterobacter spp.|MDRO|organisms") ~ "Several bacteria and/or fungi",
    str_detect(microbial_etiology_and_resistance, "maltophilia|sobria|negative|aeruginosa|baumannii|Escherichia coli|Stenotrophomonas|Enterobacteriaceae|Klebsiella|ESBL") ~ "Gram-negative bacteria",
    str_detect(microbial_etiology_and_resistance, "Clostridioides difficile|C. difficile") ~ "Clostridioides difficile",
    str_detect(microbial_etiology_and_resistance,"aureus|MRSA|enterococc|vancomycin|Bacillus|epidermidis|pneumo") ~ "Gram-positive bacteria",
    str_detect(microbial_etiology_and_resistance, "Aspergillus|molds|aspergillus|Molds|fungi|andida") ~ "Fungi",
    TRUE ~ "Other or unknown"
  ))

### Table with results

Mortality_table <- Mortality %>%
  group_by(microbe, category) %>%
  summarize(count=n()) %>%
  arrange(microbe, desc(count)) %>%
  group_by(microbe) %>%
  slice_head(n = 5)
write_xlsx(Mortality_table, "Mortality_table.xlsx")

### List the different risk factors in each category

Mortality_riskfactors <- Mortality %>%
  group_by(category) %>%
  summarise(riskfactors = paste(value, collapse = ", ")) %>%
  mutate(outcome = "Mortality")

## List of risk factors in each category

RiskFactors <- Infection_riskfactors %>%
  bind_rows(Mortality_riskfactors)
write_xlsx(RiskFactors, "RiskFactor_table.xlsx")

