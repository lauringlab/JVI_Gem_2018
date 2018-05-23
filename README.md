# The complexities of viral mutation rates (Peck and Lauring 2018)

Data for generating Figure 1 for [The complexities of viral mutation rates (Peck and Lauring 2018)](http://jvi.asm.org/content/early/2018/04/26/JVI.01031-17.short).

## Figure_1_mu_and_K_data.csv
Contains data compiled from [Sanjuan 2012](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1002685) (for evolutionary rates), [Sanjuan and Domingo-Calap 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5075021/) (mutation rates), and additional references for estimates of evolutionary and mutation rates. Data file also includes group, family, and genome size for each viral species.

## Figure_1_mu_and_K_references.docx
Full references for the mutation rate and evolutionary rate data in Figure_1_mu_K_data.csv. References are for values reported in Sanjuan 2012, Domingo-Calap and Sanjuan 2016, as well as a few additional reference (post-2016 estimates).

## Lynch_2016_mu_data.csv
Data from [Lynch et al. 2016](https://www.nature.com/articles/nrg.2016.104) ([Supplemental Table 1](https://media.nature.com/original/nature-assets/nrg/journal/v17/n11/extref/nrg.2016.104-s1.xlsx)) which includes the genome size and mutation rate estimates for a number of multicellular, unicellular, and eubacteria taxa.

## Figure_1_plots.R
R code to read in the above .csv files and output: 
* Figures 1A: viral evolutionary rates, K, vs. mutation rates, mu, for each Baltimore class (data from Figure_1_mu_and_K_data.csv)
* Figure 1B: viral evolutionary rates, K, vs. mutation rates, mu, for individual viruses (ones that include an estimate for both K and mu) (data from Figure_1_mu_and_K_data.csv)
* Figure 1C: mutation rate, mu, vs. genome size, G, for viruses (data from Figure_1_mu_and_K_data.csv) as well as eubacteria, unicellular, and multicellular taxa (data from Lynch_2016_mu_data.csv)

## Shiny app (server.R and ui.R)
To run the shiny app displaying Figures 1A-1C, open R or RStudio and run:

install.packages("shiny")

library(shiny)

runGitHub("JVI_Gem_2018", "lauringlab")

Alternatively, download the files, open them in RStudio and click on "Run App".
