#----------------
#
# Code to generate Figure 1 (A-C) for Peck and Lauring 2018, JVI
#
# Developed by Kayla Peck - Last updated 18.05.23
#
#----------------

# Read in data 
#
## This is data compiled from Sanjuan and Domingo-Calap 2016 (for mutation rates)
## and Sanjuan 2012 (for evolution rates, geometric mean calculated) plus added estimates that I found in the literature 
## for both mutation rates and evolution rates since those publications. The geometric mean was used when incorporating additional estimates.
#
# Compiled data: Figure_1_mu_and_K_data.csv
# Compiled references: Figure_1_mu_and_K_references.docx

dat <- read.csv("Figure_1_mu_and_K_data.csv")
#dat <- dat[,1:6] #remove unneccessary columns with reference info

#Read in data from Lynch et al. 2016 (for Figure 1C)
datL <- read.csv("Lynch_2016_mu_data.csv")
names(datL)[1] <- "group" #rename group column

#Add in column to store genomic mutation rates
dat$U <- dat$G*1000*dat$mu 
datL$U <- datL$G*1000000*datL$mu 

#Calculate K and mu for each Baltimore class
#Note that geometric means are used for both metrics
K <- mu <- G <- 0
for(i in 1:length(unique(dat$group))){
  sub <- subset(dat, group==unique(dat$group)[i] & Included.in.2018.publication == T)
  K[i] <- 10^mean(log10(sub$K), na.rm=T)
  mu[i] <- 10^mean(log10(sub$mu), na.rm=T)
  G[i] <- 10^mean(log10(sub$G), na.rm=T)
}

#Set up color palettes
cols <- c("gold","forestgreen","orangered","firebrick","dodgerblue", "darkorchid")
colsL <- c("black", "white","grey70")
dat$colors <- "dodgerblue"
for(i in 1:length(dat$group)){
  if(dat$group[i]=="ss_neg_RNA"){
    dat$colors[i] <- "firebrick"
  }
  if(dat$group[i]=="dsRNA"){
    dat$colors[i] <- "forestgreen"
  }
  if(dat$group[i]=="retro"){
    dat$colors[i] <- "orangered"
  }
  if(dat$group[i]=="ssDNA"){
    dat$colors[i] <- "darkorchid"
  }
  if(dat$group[i]=="dsDNA"){
    dat$colors[i] <- "gold"
  }
}
datL$colors <- "black"
for(i in 1:length(datL$group)){
  if(datL$group[i]=="eubacteria"){
    datL$colors[i] <- "grey70"
  }
  if(datL$group[i]=="unicellular"){
    datL$colors[i] <- "white"
  }
}

#Figure 1A: Plot mean K vs. mean mu for Baltimore classes (save as 700 x 700)
pdf(file="K_vs_mu.pdf")
par(mar=c(5,6,1,1))
plot(log10(K)~log10(mu), ylim=c(-5,-2), xlim=c(-7.5,-3.8), pch=16, col=cols, cex=2,
     ylab=expression(paste("Evolutionary rate (s/n/y)")), xlab=expression(paste("Mutation rate (s/n/c)")), 
     xaxt='n', yaxt='n', cex.lab=2)
points(log10(K)~log10(mu), pch=1, cex=2)
axis(1, cex.axis=1.25, at=c(-7,-6,-5,-4), labels=c(expression(paste(10^{-7})), expression(paste(10^{-6})), expression(paste(10^{-5})), expression(paste(10^{-4}))))
axis(2, cex.axis=1.25, las=2, at=c(-5,-4,-3,-2), labels=c(expression(paste(10^{-5})), expression(paste(10^{-4})), expression(paste(10^{-3})), expression(paste(10^{-2}))))
legend("topleft", ncol=2, c("dsDNA","dsRNA","retro","(-)ssRNA", "(+)ssRNA", "ssDNA"), 
       pch=21, pt.bg=cols, cex=1.5, bty='n')
dev.off()

#Figure 1B: Plot K vs. mu for individual viruses
pdf(file="K_vs_mu_individual.pdf")
par(mar=c(5,6,1,1))
dat.sub <- na.omit(dat[dat$Included.in.2018.publication==T,])
plot(log10(K)~log10(mu), dat=dat.sub, ylim=c(-5,-1), xlim=c(-7.5,-3.5), pch=c(21,22,23,24,25,21,22,21,23,21),
     bg=c(rep("dodgerblue",5),rep("firebrick",2),rep("orangered",2),"gold"),
     ylab=expression(paste("Evolutionary rate (s/n/y)")), xlab=expression(paste("Mutation rate (s/n/c)")), 
     xaxt='n', yaxt='n',cex=2,cex.lab=2)
axis(1, cex.axis=1.25, at=c(-7,-6,-5,-4), labels=c(expression(paste(10^{-7})), expression(paste(10^{-6})), expression(paste(10^{-5})), expression(paste(10^{-4}))))
axis(2, cex.axis=1.25, las=2, at=c(-5,-4,-3,-2,-1), labels=c(expression(paste(10^{-5})), expression(paste(10^{-4})), expression(paste(10^{-3})), expression(paste(10^{-2})), expression(paste(10^{-1}))))
legend("topleft", ncol=2, c("Tobacco mosaic virus", "Human rhinovirus", "Poliovirus 1",
                            "Human norovirus", "Hepatitis C virus", "Influenza A virus",
                            "Measles virus", "Avian HBV", "HIV-1", "Herpes simplex virus"), 
       pt.bg=c(rep("dodgerblue",5),rep("firebrick",2),rep("orangered",2),"gold"), pch=c(21,22,23,24,25,21,22,21,23,21),
       cex=1.2, bty='n')
dev.off()

#Figure 1C: Plot of mu vs. G (viral and nonviral data)
pdf(file="mu_vs_G_all.pdf")
par(mar=c(5,6,1,1))
plot(log10(mu)~log10(G*1000), dat=dat[dat$Included.in.2018.publication==T,], ylim=c(-11,-1), xlim=c(3,10),
     ylab=expression(paste("Mutation rate (s/n/c or s/n/g)")), xlab=expression(paste("G (bp)")), 
     xaxt='n', yaxt='n', bg=dat$colors, pch=21, cex=1.5, cex.lab=2)
points(log10(mu)~log10(G*1e6),dat=datL[datL$included==T,], bg=datL$colors, pch=21, cex=1.5)
axis(1, cex.axis=1.25, at=c(3,4,5,6,7,8,9,10), labels=c(expression(paste(10^{3})),expression(paste(10^{4})),expression(paste(10^{5})),expression(paste(10^{6})),expression(paste(10^{7})),expression(paste(10^{8})),expression(paste(10^{9})),expression(paste(10^{10}))))
axis(2, cex.axis=1.25, las=2, at=c(-11, -9,-7,-5,-3,-1), labels=c(expression(paste(10^{-11})), expression(paste(10^{-9})), expression(paste(10^{-7})), expression(paste(10^{-5})), expression(paste(10^{-3})), expression(paste(10^{-1}))))
legend("topright", ncol=2, c("multicellular","unicellular", "eubacteria","", "", "", "dsDNA", "dsRNA", "retro", "(-)ssRNA", "(+)ssRNA", "ssDNA"), 
       pch=21, pt.bg=c(colsL,rep("white",3),cols), col=c(rep(1,3),rep("white",3),rep(1,6)),cex=1.5, bty='n')
dev.off()


