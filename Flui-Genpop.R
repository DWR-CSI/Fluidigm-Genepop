# Setup -----------
library(caroline)
library(tidyverse)
library(adegenet)
library(RColorBrewer)
library(readxl)

# Settings -----------------------
# Set working directory
# setwd("C:/Users/sbrown/OneDrive - California Department of Water Resources/Salvage Pilot Study Documents/Fluidigm Data/20240321_Salvage_JPE/R Script Processing")
setwd("C:/Users/bryann/Documents/code/fl2gp") # set to your working directory where relevant files exist.
study_name <- "JPE salvage Fluidigm" # choose study name accordingly
fluidigm_input_loc <- "Results_20240321_Salvage_JPE.csv" #input, set to corresponding file
snp_keyfile_name <- "SNPKey.csv" # output, name it whatever
genpop_filename <- "genepop_test2.gen" # An output, name it whatever
locus_keyfile_name <- "Locus_Key_New_Order_240129.csv" # Input, set accordingly
greb_key_file <- "greb1l_calls_fluidigm_guide.csv" # Input, set accordingly

greb_output_filename <-"greb1l_data.tsv" # Output, give an informative and unique name

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
    number = 1:length(colnames(x)[-1]), # numbered sequentially
    greb1l = ifelse(grepl("GREB", colnames(x)[-1]), TRUE, FALSE) # greb1l or other
  )
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

# Greb functions ----------------------------
lookup_greb <- function(allele, locus_name, greb_key) { # Look up greb allele in greb key file
  ifelse(greb_key$X_early[which(greb_key$SNP_NAME == locus_name)] & allele == "X" | !greb_key$X_early[which(greb_key$SNP_NAME == locus_name)] & allele == "Y", "Early", "Late")
}

translate_greb_locus <- function(gcall, locus, greb_key) {
  # Given a greb1l genotype, locus info and key file, return early or late
  if (is.na(gcall)) {
    return(NA)
  } else if (gcall == "No Call") {
    return("No Call")
  } else {
    allele1 <- substr(as.character(gcall), 1, 1)
    allele2 <- substr(as.character(gcall), 2, 2)
    if (allele1 == allele2) {
      lookup_greb(allele = allele1, locus_name = locus, greb_key = greb_key)
    } else if (allele1 == "X" & allele2 == "Y") {
      return("Heterozygote")
    } else {
      return("Error")
    }
  }
}

translate_grebs <- function(x, greb_key) { #separate the two alleles per locus
  data <- x %>%
    pivot_longer(cols = -samplename, names_to = "locus", values_to = "gcall") %>%
    mutate(
      Allele1 = ifelse(gcall != "No Call", substr(gcall, 1, 1), NA),
      Allele2 = ifelse(gcall != "No Call", substr(gcall, 2, 2), NA),
      greb_el = map2_chr(gcall, locus, ~ translate_greb_locus(.x, .y, greb_key))
    )
  return(data)
}

extract_grebs <- function(data, greb_key) {
  data2 <- data %>%
    select(
      "samplename",
      starts_with("GREB")
    ) %>%
    translate_grebs(greb_key)
  return(data2)
}

# Working code ---------------------
data2 <- read_fluidigm(fluidigm_input_loc)

# Pull out Greb1l loci -------
greb_key_data <- read.csv(greb_key_file) %>% # load greb locus info key file
  mutate(X_early = (Allele.X.pheno == "early"))
grebs <- extract_grebs(data2, greb_key_data)

grebs_wider <- pivot_wider( # Format: Each locus is a column
  grebs,
  id_cols = samplename,
  names_from = locus,
  values_from = greb_el
) # You can filter out grebs here if you don't want all of them

# Total up Early, Late, Het, No Call counts
grebs_wider$early_greb_count <- rowSums(select(grebs_wider, starts_with("GREB")) == "Early", na.rm = TRUE)
grebs_wider$late_greb_count <- rowSums(select(grebs_wider, starts_with("GREB")) == "Late", na.rm = TRUE)
grebs_wider$het_greb_count <- rowSums(select(grebs_wider, starts_with("GREB")) == "Heterozygote", na.rm = TRUE)
grebs_wider$nocall_greb_count <- rowSums(select(grebs_wider, starts_with("GREB")) == "No Call", na.rm = TRUE)


write_tsv(grebs_wider, file = greb_output_filename)

# Translate fluidigm data to genepop format ----------------------------
# pull out a key for the snp names
snpkey <- gen_snpkey(data2)

translated_data2 <- translate_fluidigm_to_genepop(data2)

# re-order loci according to keyfile

resorted_data2 <- resort_fl2gp(translated_data2, keyfile = locus_keyfile_name)

# export the genepop file
write_genpop_file(data = resorted_data2, filename = genpop_filename, title = study_name)
write.csv(snpkey, snp_keyfile_name)


# Run a quick PCA to see if the format works


adphen <- read.genepop(genpop_filename, ncode = 3)
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