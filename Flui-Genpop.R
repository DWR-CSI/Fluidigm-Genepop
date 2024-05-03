# Setup -----------
library(caroline)
library(data.table)
library(tidyverse)
library(adegenet)
library(RColorBrewer)

# Settings -----------------------
# Set working directory
# setwd("C:/Users/sbrown/OneDrive - California Department of Water Resources/Salvage Pilot Study Documents/Fluidigm Data/20240321_Salvage_JPE/R Script Processing")
setwd("/Users/bryan/Documents/DWR/fluidigm/") # set to your working directory where relevant files exist.
study_name <- "JPE salvage Fluidigm" # choose study name accordingly
fluidigm_input_loc <- "Results_20240321_Salvage_JPE.csv"
snp_keyfile_name <- "SNPKey.csv" # name it whatever
genpop_input_filename <- "~/Documents/DWR/fluidigm/genepop_test2.gen" # name it whatever
locus_keyfile_name <- "Locus_Key_New_Order_240129.csv"

# Functions ----------------------

read_fluidigm <- function(x, skip = 15, ...) { # x = filename
  data <- read_csv(x, skip = skip, show_col_types = FALSE, ...) %>%
    select(-c("...1")) %>% # skip first column
    rename(samplename = "...2") %>% # second column is samplename
    mutate(across(-1, ~ factor(.x, levels = c("XX", "XY", "YY", "No Call", "NTC", "Invalid")))) # rest of columns are factor data
  return(data)
}

gen_snpkey <- function(x) {
  data <- tibble(
    SNP = colnames(x)[-1], # remove first col which is sample names
    number = 1:length(colnames(x)[-1])
  ) # numbered sequentially
  return(data)
}

translate_fluidigm_to_genepop <- function(data) { # translates XX to 100100, etc...
  # levels_original <- c("XX","XY","YY","No Call", "NTC", "Invalid") #not used
  new_data <- data %>%
    mutate(across(
      -1, # skip first column
      ~ fct_recode(
        .x,
        "100100" = "XX", # translate data
        "100200" = "XY",
        "200200" = "YY",
        "000000" = "No Call",
        "000000" = "NTC",
        "000000" = "Invalid"
      )
    ))
  return(new_data)
}

write_genpop_file <- function(data, filename, title = "Fluidigm", ...) { # Write a table with genepop format calls to a file
  initial_data <- c(title, 1:96, "Pop") # required header data for genepop, 1:96 represents the SNPs
  write.table(initial_data,
    file = filename,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE, ...
  )
  data_with_comma <- data %>%
    add_column(comma = ",", .after = 1) # add comma between rows 1 and 2
  write.table(data_with_comma,
    file = filename,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE, # Append to the file
    quote = FALSE, ...
  )
}

resort_fl2gp <- function(data, keyfile = "Locus_Key_New_Order_240129.csv") {
  order_info <- read_csv(keyfile) %>% # must be csv
    mutate(Number = Number + 1) # leave 1 column for sample names, no commas at this point
  new_order <- c(1, order_info$Number)
  data_reordered <- data[, new_order]
  return(data_reordered)
}
# Working code ---------------------
data2 <- read_fluidigm(fluidigm_input_loc)

# pull out a key for the snp names
snpkey <- gen_snpkey(data2)

translated_data2 <- translate_fluidigm_to_genepop(data2)

# re-order loci according to keyfile

resorted_data2 <- resort_fl2gp(translated_data2, keyfile = locus_keyfile_name)

# export the genepop file
write_genpop_file(data = resorted_data2, filename = genpop_input_filename, title = study_name)
write.csv(snpkey, snp_keyfile_name)


# Run a quick PCA to see if the format works


adphen <- read.genepop(genpop_input_filename, ncode = 3)
adphen$pop
x <- scaleGen(adphen, NA.method = "mean")
# quick look at data before putting populations into genepop file
# look for things like too many groups, lots of samples in the middle etc.
pca1 <- dudi.pca(x, cent = F, scale = F, scannf = F, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))
col <- brewer.pal(4, "Dark2")
s.class(pca1$li, pop(adphen),
  xax = 1, yax = 2, col = col, axesell = FALSE,
  cstar = 0, cpoint = 1.5, grid = FALSE, cellipse = 0
)

# filter out samples with lots of missing data
tab_ad <- adphen$tab
missing <- matrix(nrow = 96, ncol = 2)
missing[, 1] <- row.names(tab_ad)
samples <- row.names(tab_ad)
for (i in 1:length(tab_ad[, 1])) {
  missing[i, 2] <- 1 - length(na.omit(tab_ad[i, ])) / 192
}
colnames(missing) <- c("Sample", "% Missing")
missing
write.csv(missing, "missing_data.csv")

# Order using a list (advanced) -------------------------------------------

# order using a list
sites <- read.csv("Samples_site_order.csv")

by_sites <- data.frame(ID = sites, com = rep(",", times = (length <- nrow(p2))))

daf <- data.frame()
for (i in 1:length(sites[, 1])) {
  a <- by_sites[i, 1]
  b <- new[which(new$SampleID == a), ]
  daf <- rbind(daf, b)
}
daf
colnames(daf) <- c("SampleID", "com", colnames(p2[2:ncol(p2)]))
daf
row.names(daf) <- daf$SampleID
daf
daf <- daf[, -1]
daf

# remove individuals with a ton of missing data
# remove individuals (aka individuals with low data)
# read in a csv with just the individual names
bad_samps <- read.csv("missing_data.csv", header = T)
rem_samps <- bad_samps[which(bad_samps$X..Missing > 0.1), 2]
length(rem_samps)
rem_samps
# remove the samples
minbad <- adphen[!row.names(adphen@tab) %in% rem_samps]
# check the right number of samples were removed
length(row.names(minbad@tab))
write.delim(minbad, "badout.gen", quote = FALSE, row.names = T, sep = "\t")
genind2genpop(minbad, "badout.gen")

# bring in file with bad sequenced samples removed
badout <- read.genepop("badout.gen", ncode = 1)
x <- scaleGen(minbad, NA.method = "mean")
pca1 <- dudi.pca(x, cent = F, scale = F, scannf = F, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))
s.class(pca1$li, pop(minbad),
  xax = 1, yax = 2, col = col, axesell = FALSE,
  cstar = 0, cpoint = 1.5, grid = FALSE, cellipse = 0, clabel = 0
)




# Greb1L Calling ----------------------------------------------------------
library(readxl)

getwd()
# Set working directory
setwd("C:/Users/smeyer/OneDrive - California Department of Water Resources/Scott/Fluidigm-git/Fluidigm-Genepop")
# Read in data
data1 <- read.csv("Results_20240321_Salvage_JPE.csv", skip = 14, nrows = 97)

# pull out a key for the snp names
snps <- as.vector(data1[1, ])
nums <- colnames(data1)
snpkey <- cbind(snps, nums)


# remove the snp row and add SNP names
colnames(data1) <- data1[1, ]
data1 <- data1[-1, -1]

# pull out greb data
grebs <- data1[, ]

grebcols <- grep("GREB", snpkey[, 1])
grebdata <- matrix(nrow = 96, ncol = 14)
for (i in 1:length(grebcols)) {
  b <- grebcols[i]
  a <- data1[, b - 1]
  grebdata[, i] <- a
}
colnames(grebdata) <- snpkey[grep("GREB", snpkey[, 1]), 1]

greb_key <- read.csv("greb1l_calls_fluidigm_guide.csv")

for (i in 1:length(greb_key[, 1])) {
  a <- greb_key[i, 1]
  b <- grebdata[, a]
  homox <- which(b == "XX")
  homoy <- which(b == "YY")
  het <- which(b == "XY")
  grebdata[homox, a] <- greb_key[i, 4]
  grebdata[homoy, a] <- greb_key[i, 6]
  grebdata[het, a] <- "Het"
}

grebdata <- as.data.frame(grebdata)
grebdata$Samples <- data1[, 1]

greblong <- pivot_longer(grebdata, cols = 1:14)

library(ggplot2)
ggplot(greblong, aes(x = name, y = Samples, fill = value)) +
  geom_tile()
