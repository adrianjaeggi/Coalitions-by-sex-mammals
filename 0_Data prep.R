## set working directory
setwd("//files.iem.uzh.ch/Data/Institute/Human_Ecology/ajaegg/Private/My Documents/GitHub/sex-bias-cooperation")

## read excel file
d<- readxl::read_xlsx("Sex bias in intragroup coalitions across mammals_2-11-22_PM.xlsx", sheet = 1)

# select variables to be used
d<- d[,c("Genus_species", "Order", "Philopatry", "SexComposition", "Sex_Dim", "Food_resource_defendable", "Within-Conflict_FormationBySex_zero_center")]

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
d$Philopatry[d$Philopatry=="F?"]<- "F" 



## load phylogeny
library(ape)
trees<- read.nexus("VertlifeTreeBlock2.nex") # list of 100 trees


# compare species names and resolve discrepancies
setdiff(d$Genus_species, trees[[1]]$tip.label)
d$Genus_species[d$Genus_species=="Stenella_longirostris_longirostris"]<- "Stenella_longirostris"
d$Genus_species[d$Genus_species=="Tursiops_spp"]<- "Tursiops_truncatus"
d$Genus_species[d$Genus_species=="Panthera_leo_leo"]<- "Panthera_leo"
d$Genus_species[d$Genus_species=="Equus_ferus_przewalskii"]<- "Equus_ferus"
d$Genus_species[d$Genus_species=="Equus_burchellii_quagga"]<- "Equus_quagga"
d$Genus_species[d$Genus_species=="Cercopithecus_kandti"]<- "Cercopithecus_mitis"
# remove dogs
d<- d[which(d$Genus_species!="Canis_lupus_familiaris"),]
setdiff(d$Genus_species, trees[[1]]$tip.label)

write.csv(d, "sexbiascoalitions.csv", row.names = FALSE)

## drop unused species
trees.mov<- trees
for(i in 1:100){
  trees.mov[[i]]<- drop.tip(trees.mov[[i]], setdiff(trees.mov[[i]]$tip.label, d$Genus_species))
}

write.tree(trees.mov, "TreeBlockSexBiasCooperation.tre")


