library(tidyverse)
library(readr)
library(data.table)
library(readxl)

# read pihat file, pheno and missing files.
# use pihat cut off 0.25
relate <- read_xlsx("/mnt/share6/FOR_Takeo/ICHIP/IBDichip1to8D_washU_HG19B_BBC_c_filtered3_genome.genome.xlsx")


pheno <- fread("/mnt/isilon_data/For_Takeo2/From_nick/imputed_ichip/pancreatitis_pheno_imputed_ichip") %>% 
  dplyr::select(IID,pancreatitis) %>% 
  mutate(IID2 = IID) %>% 
  separate(IID2,into= c("FID","non_need"),sep = "_")

freeze <-  read_xls("/mnt/share6/FOR_Takeo/phenotypdata/Phenotype _TRICS _2019_07_26_rev.xls",
                    col_types = "text")
colnames(freeze) <- make.names(colnames(freeze))
freeze <- freeze %>% 
  mutate(FID = str_replace_all(Genetic.ID,"-","0")) %>% 
  dplyr::select(FID,Age)

missing <- fread("/mnt/share6/FOR_Takeo/ICHIP/no_filter/subclinical/cedars_qced.imiss")



relate1 <- relate %>% 
  filter(PI_HAT >= 0.25) %>% 
  dplyr::select(FID1,FID2,PI_HAT)
dim(relate1)



# add columns of phenotype and missing rate to relatedness file
# if Miss_one and Miss_two are both  missing, these pairs are not European --> remove
relate2 <- relate1 %>% 
  left_join(pheno,by = c("FID1" = "FID") ) %>% 
  left_join(pheno,by = c("FID2" = "FID")) %>% 
  dplyr::rename(IBD_one = pancreatitis.x) %>% 
  dplyr::rename(IBD_two = pancreatitis.y) %>% 
  left_join(missing,by = c("FID1" = "FID")) %>% 
  left_join(missing,by = c("FID2" = "FID")) %>% 
  dplyr::rename(Miss_one = F_MISS.x ) %>% 
  dplyr::rename(Miss_two = F_MISS.y ) %>% 
  dplyr::select(FID1,IBD_one,Miss_one,FID2,IBD_two,Miss_two) %>% 
  mutate(minimum = pmin(Miss_one,Miss_two)) %>% 
  mutate(Miss_one = Miss_one + 0.1,
         Miss_two = Miss_two + 0.1) %>% 
  filter(!is.na(Miss_one) | !is.na(Miss_two))


## input "missing" to NA for filtering 

relate2 <- relate2 %>% 
  mutate_at(.vars = vars("IBD_one","Miss_one","IBD_two","Miss_two","minimum"),
            list(~ifelse(is.na(.),"missing",.)))

relate3 <- relate2 %>% 
  filter(IBD_one != "missing" & IBD_two != "missing")



relate3 <- relate3 %>% 
  mutate(jogai = ifelse(IBD_one == 2 & IBD_two == 1, FID2,
                        ifelse(IBD_one == 1 & IBD_two == 2, FID1,
                               ifelse(Miss_one > Miss_two, FID2,
                                      ifelse(Miss_one < Miss_two, FID1,NA))))) %>% 
  select(jogai) 


fam = fread('/mnt/isilon_data/For_Takeo2/From_nick/imputed_ichip/plinkfiles/ichip1t7/chr1_bfile.fam')
fam <- fam %>% 
  rename(IID1 = V2)

fam <- fam %>% 
  mutate(IID2 = IID1) %>% 
  separate(IID1,into = c("GeneticID","no_need"),sep = "_")


relate4 <- relate3 %>% 
  left_join(fam,by= c("jogai" = "GeneticID")) %>%
  select(V1,IID2) %>% 
  rename(FID = V1,
         IID = IID2)
relate4 %>% 
  distinct(IID,.keep_all = T) %>% 
  write_tsv("/mnt/isilon_data/For_Takeo2/From_nick/imputed_ichip/related_imputed_ichip")


