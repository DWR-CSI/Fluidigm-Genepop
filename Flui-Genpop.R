library(caroline)
getwd()
#Set working directory
setwd("C:/Users/sbrown/OneDrive - California Department of Water Resources/Salvage Pilot Study Documents/Fluidigm Data/20240321_Salvage_JPE/R Script Processing")
#Read in data
data1 <- read.csv("Results_20240321_Salvage_JPE.csv", skip=14, nrows=97)

#pull out a key for the snp names
snps <-  as.vector(data1[1,])
nums <- colnames(data1)
snpkey <- cbind(snps, nums)

#remove the snp row
data1 <- data1[-1,-1]


new <- data.frame(ID = data1[,1], com=rep("," , times=(length=nrow(data1))))
for(i in 2:ncol(data1)){   #change the number here (2) to the column where the genotype data starts in your dataset
  temp1 = vector(length=nrow(data1))
  temp1[data1[,i]=='XX'] = "100100"
  temp1[data1[,i]=='XY'] = "100200"
  temp1[data1[,i]=='YY'] = "200200"
  temp1[data1[,i]=='No Call'] = "000000"
  temp1[data1[,i]=='NTC'] = "000000"
  temp1[data1[,i]=='Invalid'] = "000000"
  
  new <- cbind(new, temp1)  
}
colnames(new) <- c("SampleID", "com", colnames(data1[2:ncol(data1)]))
new
newest <- new[order(new$SampleID),]
#Get loci in order
loci <- read.csv("Locus_Key_New_Order_240129.csv")

corrected_loci <- data.frame(ID = loci, com=rep("," , times=(length=nrow(newest))))

daf <- matrix(nrow=96, ncol=96)
for (i in 1:length(loci[,1])) {
  a <- loci[i,2]
  b <- new[,(as.numeric(a)+2)]
  daf[,i] <- b
}
daf
daf <- cbind(new[,1:2],daf)
daf


#export the genepop file
write.delim(daf, "Genepop_20240321_salvage_JPE.gen", quote = FALSE, row.names=F, sep="\t", col.names=F)
#export snp key
write.csv(snpkey, "SNPKey.csv")


#Run a quick PCA to see if the format works
library(adegenet)
library(RColorBrewer)

adphen <- read.genepop("Genepop_20240321_salvage_JPE.gen", ncode=3)
adphen$pop
x <- scaleGen(adphen, NA.method="mean")
#quick look at data before putting populations into genepop file
#look for things like too many groups, lots of samples in the middle etc. 
pca1 <- dudi.pca(x, cent=F, scale=F, scannf=F, nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
col <- brewer.pal(4, "Dark2")
s.class(pca1$li, pop(adphen),xax=1,yax=2, col=col, axesell=FALSE,
        cstar=0, cpoint=1.5, grid=FALSE, cellipse = 0)

#filter out samples with lots of missing data
tab_ad <- adphen$tab
missing <- matrix(nrow = 96, ncol=2)
missing[,1] <- row.names(tab_ad)
samples <- row.names(tab_ad)
for (i in 1:length(tab_ad[,1])) {
  missing[i,2] <- 1-length(na.omit(tab_ad[i,]))/192
}
colnames(missing) <- c("Sample","% Missing")
missing
write.csv(missing, "missing_data.csv")

# Order using a list (advanced) -------------------------------------------

#order using a list
sites <- read.csv("Samples_site_order.csv")

by_sites <- data.frame(ID = sites, com=rep("," , times=(length=nrow(p2))))

daf <- data.frame()
for (i in 1:length(sites[,1])) {
  a <- by_sites[i,1]
  b <- new[which(new$SampleID==a),]
  daf <- rbind(daf, b)
}
daf
colnames(daf) <- c("SampleID", "com", colnames(p2[2:ncol(p2)]))
daf
row.names(daf) <- daf$SampleID
daf
daf <- daf[,-1]
daf 

#remove individuals with a ton of missing data
#remove individuals (aka individuals with low data)
#read in a csv with just the individual names
bad_samps <- read.csv('missing_data.csv', header=T)
rem_samps <- bad_samps[which(bad_samps$X..Missing>0.1),2]
length(rem_samps)
rem_samps
#remove the samples
minbad <- adphen[!row.names(adphen@tab) %in% rem_samps]
#check the right number of samples were removed
length(row.names(minbad@tab))
write.delim(minbad, "badout.gen", quote = FALSE, row.names=T, sep="\t")
genind2genpop(minbad, "badout.gen")

#bring in file with bad sequenced samples removed
badout <- read.genepop("badout.gen", ncode=1)
x <- scaleGen(minbad, NA.method="mean")
pca1 <- dudi.pca(x, cent=F, scale=F, scannf=F, nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
s.class(pca1$li, pop(minbad),xax=1,yax=2, col=col, axesell=FALSE,
        cstar=0, cpoint=1.5, grid=FALSE, cellipse = 0, clabel=0)

