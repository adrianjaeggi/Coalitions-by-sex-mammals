## set working directory
setwd("//files.iem.uzh.ch/Data/Institute/Human_Ecology/ajaegg/Private/My Documents/GitHub/Coalitions-by-sex-mammals")

## read excel file
d0<- readxl::read_xlsx("Sex bias in intragroup coalitions across mammals_6-15-22.xlsx", sheet = 1)

# select variables to be used - use unique() to get only one record for species-level data
d<- unique(d0[,c("Genus_species", "Order", "Philopatry", "SexComposition", "Sex_Dim", "Food_resource_defendable", "Within-Conflict_FormationBySex_zero_center_exclude_domestic", "Male mass (kg)", "Female mass (kg)")])

# rename and recode coalitions
colnames(d)[7]<- "coalitions"
d$coalitions[d$coalitions=="-1.0"]<- "Males"
d$coalitions[d$coalitions=="0.0"]<- "Both"
d$coalitions[d$coalitions=="1.0"]<- "Females"
d$coalitions<- as.factor(d$coalitions)

# remove species with missing data on coalitions
d<- d[d$coalitions!="NA",]

# fix Philopatry categories
d$Philopatry[d$Philopatry=="B/N"]<- "B" # this species has records of both male and female philopatry (see spreadsheet, hence B)

# generate bibliometric frequency measures
d.bib<- d0[,c("Genus_species", "Adult males intervene in an ongoing fight within the group or form an intragroup coalition? (y =1, n = 0)", "Adult females intervene in an ongoing fight within the group or form an intragroup coalition? (y =1, n = 0)")]
colnames(d.bib)[2]<- "freq.males"
colnames(d.bib)[3]<- "freq.females"
d.bib$freq.males<- as.numeric(d.bib$freq.males)
d.bib$freq.females<- as.numeric(d.bib$freq.females)

# aggregate
d.bib<- aggregate(d.bib[,2:3], by=list(Genus_species=d.bib$Genus_species), sum, na.rm=T)

# merge
d<- merge(d, d.bib, all.x=TRUE)

# recode frequencies >2 to be 2
d$freq.males[d$freq.males>2]<- 2
d$freq.females[d$freq.females>2]<- 2

# recode sex dim following Smith 1999
d[,"Male mass (kg)"]<- as.numeric(d[,"Male mass (kg)"])
d[,"Female mass (kg)"]<- as.numeric(d[,"Female mass (kg)"])

d$Sex_Dim[which(d$'Female mass (kg)'>d$'Male mass (kg)')] <- 2- d$'Female mass (kg)'[which(d$'Female mass (kg)'>d$'Male mass (kg)')] / d$'Male mass (kg)'[which(d$'Female mass (kg)'>d$'Male mass (kg)')]



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


