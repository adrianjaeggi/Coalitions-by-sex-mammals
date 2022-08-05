This repository contains all the data, code and markdown to reproduce the analyses, results section and supplemental results for the paper 
'Sex differences in cooperative coalitions: A mammalian perspective'
by Jennifer Smith, Adrian Jaeggi, Rose Holmes, and Joan Silk, published in the Philosophical Transactions of the Royal Society B. 

Note that two versions of the analyses were done, one using a consensus phylogeny - code: 1_Main_models.R - and one using a sample of 100 phylogenies - code: 2_Main_models_looped.R. The former was used for the initial submission of the paper due to time constraints, the latter for the final published version. The results did not differ. If you'd like access to the saved posterior samples - which were too large to push - please contact me (Adrian Jaeggi, adrian.jaeggi@iem.uzh.ch) directly.

List of files and explanation:

Sex bias in intragroup coalitions across mammals_6-15-22.xlsx - Dataset compiled by Jennifer Smith and Rose Holmes

TreeBlockSexBiasCooperation.nex - Sample of 100 phylogenies for the included species, downloaded from vertlife.org

0_Data prep.R - Code for reading excel dataset and phylogenies, some data cleaning

sexbiascoalitions.csv - Clean dataset produced by 0_Data prep.R

1_Main_models.R - Code for running analyses using consensus phylogeny (not updated after initial submission)

2_Main_models_looped.R - Code for running analyses using sample of phylogenies

Sex-bias coalitions results.Rmd - Markdown script that produces the results section of the paper from the saved posterior samples

Sex diff coalitions supplement.Rmd - Markdown script that produces the supplemental results file from the saved posterior samples

PhyForFigure.tre - Consensus phylogeny produced in the process of making figure 2
