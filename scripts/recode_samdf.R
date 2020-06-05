samdf <- read.csv("sample_data/Sample_info.csv", header=TRUE) %>%
  dplyr::select(-X)


## Add dates
samdf <- samdf %>% 
  mutate(collection_date = case_when(
    str_detect(ExtractID, "CT1-") ~ "2017-10-5",
    str_detect(ExtractID, "CT2") ~ "2017-10-5",
    str_detect(ExtractID, "CT3") ~ "2017-10-13",
    str_detect(ExtractID, "CT4") ~ "2017-10-20",
    str_detect(ExtractID, "CT5") ~ "2017-10-20",
    str_detect(ExtractID, "CT6") ~ "2017-11-9",
    str_detect(ExtractID, "CT7") ~ "2017-11-17",
    str_detect(ExtractID, "CT8") ~ "2017-11-24",
    str_detect(ExtractID, "CT9") ~ "2017-12-15",
    str_detect(ExtractID, "CT10") ~ "2018-1-19",
    str_detect(ExtractID, "CT11") ~ "2017-12-17",
    str_detect(ExtractID, "CT12") ~ "2018-1-4",
    # Dros - mornington
    str_detect(ExtractID, "^DM0") ~ "2017-12-13",
    str_detect(ExtractID, "^M2") ~ "2017-12-27",
    str_detect(ExtractID, "^M4") ~ "2018-1-10",
    str_detect(ExtractID, "^M6") ~ "2018-1-24",
    str_detect(ExtractID, "^M8") ~ "2018-2-7",
    str_detect(ExtractID, "^M10") ~ "2018-2-28",
    # Dros tat
    str_detect(ExtractID, "^T2") ~ "2018-1-15",
    str_detect(ExtractID, "^T4") ~ "2018-1-29",
    str_detect(ExtractID, "^T6") ~ "2018-2-12",
    str_detect(ExtractID, "^T8") ~ "2018-2-26",
    str_detect(ExtractID, "^T10") ~ "2018-3-12",
    TRUE ~ "NA"
  ))


## Add orchard
samdf <- samdf %>% 
  mutate(geo_loc_name = case_when(
    str_detect(ExtractID, "CT1-") ~ "Wandown",
    str_detect(ExtractID, "CT2") ~ "Wandown",
    str_detect(ExtractID, "CT3") ~ "Wandown",
    str_detect(ExtractID, "CT4") ~ "Wandown",
    str_detect(ExtractID, "CT5") ~ "Wandown",
    str_detect(ExtractID, "CT6") ~ "Wandown",
    str_detect(ExtractID, "CT7") ~ "Wandown",
    str_detect(ExtractID, "CT8") ~ "Wandown",
    str_detect(ExtractID, "CT9") ~ "Wandown",
    str_detect(ExtractID, "CT10") ~ "Carina",
    str_detect(ExtractID, "CT11") ~ "Lake Powell",
    str_detect(ExtractID, "CT12") ~ "Lake Powell",
    TRUE ~ as.character(geo_loc_name)
  ))


#Add material column

samdf <- samdf %>%
  mutate(material = case_when(
    #Carpophilus
    str_detect(ExtractID, "CM[:digit:]") ~ "Carpophilus Adults",
    str_detect(ExtractID, "CT") ~ "Carpophilus Adults",
    str_detect(ExtractID, "CML[:digit:]") ~ "Carpophilus Larvae",
    #Drosophila
    str_detect(ExtractID, "D100M") ~ "Drosophila Adults",
    str_detect(ExtractID, "D250M") ~ "Drosophila Adults",
    str_detect(ExtractID, "D500M") ~ "Drosophila Adults",
    str_detect(ExtractID, "D1000M") ~ "Drosophila Adults",
    str_detect(ExtractID, "DM") ~ "Drosophila Adults",
    str_detect(ExtractID, "DLarv") ~ "Drosophila Larvae", 
    #Traps
    str_detect(ExtractID, "FF") ~ "Mixed Larvae", 
    str_detect(ExtractID, "SPD") ~ "Mixed Adults", 
    str_detect(ExtractID, "ACV") ~ "Mixed Adults", 
    str_detect(ExtractID, "DC") ~ "Mixed Adults", 
    str_detect(ExtractID, "Sach") ~ "Mixed Adults", 
    #Synthetics
    str_detect(ExtractID, "SynMock") ~ "Synthetic",    
    str_detect(ExtractID, "POS") ~ "Synthetic",  
    str_detect(ExtractID, "NTC") ~ "Blank",    
    str_detect(ExtractID, "BLANK") ~ "Blank",
    str_detect(ExtractID, "extblank") ~ "Blank",
    str_detect(ExtractID, "pcrblank") ~ "Blank"
  )) 

samdf %>% filter(is.na(material))


#Add treatment column
samdf <- samdf %>%
  mutate(treatment = case_when(
    #Carpophilus
    str_detect(ExtractID, "CM") ~ "Mock",
    str_detect(ExtractID, "CT") ~ "Trap",
    str_detect(ExtractID, "CML") ~ "Mock",
    #Drosophila
    str_detect(ExtractID, "D100M") ~ "Mock",
    str_detect(ExtractID, "D250M") ~ "Mock",
    str_detect(ExtractID, "D500M") ~ "Mock",
    str_detect(ExtractID, "D1000M") ~ "Mock",
    str_detect(ExtractID, "DM") ~ "Mock",
    str_detect(ExtractID, "DLarv") ~ "Mock", 
    #Traps
    str_detect(ExtractID, "FF") ~ "Fruit crush", 
    str_detect(ExtractID, "SPD") ~ "SPD", 
    str_detect(ExtractID, "ACV") ~ "ACV", 
    str_detect(ExtractID, "DC") ~ "DC", 
    str_detect(ExtractID, "Sach") ~ "Sachet", 
    #Synthetics
    str_detect(ExtractID, "SynMock") ~ "Synthetic",    
    str_detect(ExtractID, "POS") ~ "Synthetic",  
    str_detect(ExtractID, "NTC") ~ "Blank",    
    str_detect(ExtractID, "BLANK") ~ "Blank",
    str_detect(ExtractID, "extblank") ~ "Blank",
    str_detect(ExtractID, "pcrblank") ~ "Blank"
  ))

samdf %>% filter(is.na(treatment))

write.csv(samdf, "sample_data/Sample_info2.csv")
