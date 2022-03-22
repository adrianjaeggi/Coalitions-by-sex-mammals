## set working directory
setwd("//files.iem.uzh.ch/Data/Institute/Human_Ecology/ajaegg/Private/My Documents/GitHub/Coalitions-by-sex-mammals")

## read excel file
d<- readxl::read_xlsx("Sex bias in intragroup coalitions across mammals_3-5-22.xlsx", sheet = 1)

# select variables to be used
d<- d[,c("Genus_species", "Order", "Philopatry", "SexComposition", "Sex_Dim", "Food_resource_defendable", "Within-Conflict_FormationBySex_zero_center_exclude_domestic")]

# rename and recode coalitions
colnames(d)[7]<- "coalitions"
d$coalitions[d$coalitions=="-1"]<- "Males"
d$coalitions[d$coalitions=="0"]<- "Both"
d$coalitions[d$coalitions=="1"]<- "Females"
d$coalitions<- as.factor(d$coalitions)

# remove species with missing data on coalitions
d<- d[d$coalitions!="NA",]

# fix Philopatry categories
d$Philopatry[d$Philopatry=="B/N"]<- "B" # this species has records of both male and female philopatry (see spreadsheet, hence B)




## load phylogeny
library(ape)
trees<- read.nexus("TreeBlockSexBiasCooperation.nex") # list of 100 trees


# compare species names and resolve discrepancies
setdiff(d$Genus_species, trees[[1]]$tip.label)
d$Genus_species[d$Genus_species=="Stenella_longirostris_longirostris"]<- "Stenella_longirostris"
d$Genus_species[d$Genus_species=="Panthera_leo_leo"]<- "Panthera_leo"
d$Genus_species[d$Genus_species=="Equus_ferus_przewalskii"]<- "Equus_ferus"
d$Genus_species[d$Genus_species=="Equus_burchellii_quagga"]<- "Equus_quagga"
setdiff(d$Genus_species, trees[[1]]$tip.label)

write.csv(d, "sexbiascoalitions.csv", row.names = FALSE)


