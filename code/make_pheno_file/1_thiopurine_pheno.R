library(tidyverse)
library(readr)
library(data.table)
library(readxl)

## read shells file
thio <- read_xlsx('/mnt/share6/FOR_Takeo/WES2/thiopurine/thiopurine_sub_from_shell_with_WES.xlsx')


dim(thio)
thio <- thio %>% 
  filter(!is.na(Genetic.ID)) %>% 
  distinct(Genetic.ID,.keep_all =T)

## reaad original file
original <- read_xlsx('/mnt/share6/FOR_Takeo/phenotypdata/Thiopurine toxicityforDalin_v2.xlsx')
colnames(original) <- make.names(colnames(original))


original <- original %>% 
  drop_na(LAB.ID) %>% 
  distinct(LAB.ID,.keep_all = T) 



panc <-original %>% 
  select(LAB.ID,Pancreatitis)


## merge 2 files

panc1 <- panc %>% 
  full_join(thio,by = c("LAB.ID" = "Genetic.ID"))


panc1 <- panc1 %>% 
  mutate(pancreatitis = ifelse(!is.na(all_thiopurine),all_thiopurine,
                               ifelse(is.na(all_thiopurine),Pancreatitis.x,NA))) %>% 
  mutate(pancreatitis = ifelse(pancreatitis == "1","Yes",
                               ifelse(pancreatitis == "0","No",pancreatitis))) %>% 
  mutate(pancreatitis = ifelse(pancreatitis == "Yes",2,
                               ifelse(pancreatitis == "No",1,pancreatitis))) %>% 
  select(LAB.ID,pancreatitis)


panc1 <- panc1 %>% 
  mutate(GeneticID = str_replace_all(LAB.ID,"-","0"))
### read fam file of GSA and keyfile to add genetic ID
fam = fread('/mnt/isilon_data/For_Takeo2/From_nick/imputed_ichip/plinkfiles/ichip1t7/chr1_bfile.fam')
fam <- fam %>% 
  rename(IID1 = V2)

fam <- fam %>% 
  mutate(IID2 = IID1) %>% 
  separate(IID1,into = c("GeneticID","no_need"),sep = "_")
  





length(intersect(fam$GeneticID,panc1$GeneticID))

## merge 2 datasets

panc2 <- panc1 %>% 
  left_join(fam,by = "GeneticID") %>% 
  drop_na(no_need)


### should be removed subjects

removed_sub = fread("/mnt/isilon_data/For_Takeo2/WES2/thiopurine/shoudl_be_removed_subjects_from_shell",
                    header = F)


panc3 <- panc2 %>% 
  filter(!LAB.ID %in% removed_sub$V1)



### export as new data
panc3 %>% 
  select(V1,IID2,pancreatitis) %>% 
  rename(FID = V1,
         IID = IID2) %>% 
  write_tsv("/mnt/isilon_data/For_Takeo2/From_nick/imputed_ichip/pancreatitis_pheno_imputed_ichip")

