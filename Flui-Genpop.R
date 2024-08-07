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
study_name <- "20240605 Salvage Fluidigm 80" # choose study name accordingly
n_loci_cutoff <- 80 # Number of loci, set accordingly
max_missing_data <- 0.3 # Maximum missing data allowed
# For the Fluidigm input file, use the file created by the Fluidigm analysis software when you do: File > Export
fluidigm_input_loc <- file.choose() #"Results_20240321_Salvage_JPE.csv" #input, set to corresponding file
snp_keyfile_name <- "SNPKey.csv" # output, name it whatever
genpop_filename <- "genepop_80_bn.gen" # An output, name it whatever
locus_keyfile_name <- "Locus_Key_New_Order_240129.csv" # Input, set accordingly

greb_key_file <- "greb1l_calls_fluidigm_guide.csv" # Input, set accordingly

greb_output_filename <-"greb1l_80_data.tsv" # Output, give an informative and unique name

the_good_grebs <- c("GREB1l_pos2194538", "GREB1l_pos2198644", "GREB1l_pos2199210", "GREB1l_pos2200828",	"GREB1l_pos2202893")

log_filename <- paste0(gsub("[: ]","_",Sys.time()), "_log.txt") # Output, give an informative and unique name
sink(log_filename, append = FALSE) # Log file)

# Logging -----------------------
print(Sys.Date())
print(Sys.time())
print(paste0("Study name: ", study_name))
print(paste0("Input file: ",fluidigm_input_loc))
print(paste0("Working directory: ", getwd()))
print("Settings  -------")
print(paste0("SNP key file: ",snp_keyfile_name))
print(paste0("Locus key file: ",locus_keyfile_name))
print(paste0("Number of loci: ",n_loci_cutoff))
print(paste0("Grebs key file: ",greb_key_file))
print(paste0("Missing data cutoff: ",max_missing_data))
print("Output Files -------")
print(paste0("Grebs output file: ",greb_output_filename))
print(paste0("Genepop output file: ",genpop_filename))
# Functions ----------------------

read_fluidigm <- function(x, skip = 15, ...) { # x = filename
  data <- read_csv(x, skip = skip, show_col_types = FALSE, ...) %>%
    select(-c("...1")) %>% # skip first column
    rename(samplename = "...2") %>% # second column is samplename
    mutate(across(-1, ~ factor(.x, levels = c("XX", "XY", "YY", "No Call", "NTC", "Invalid")))) # rest of columns are factor data
  first_na_row <- which(is.na(data$samplename))[1] # find the first row with NA in the samplename column
  if(!is.na(first_na_row)) { #if not trimmed, trim
    return(data[1:(first_na_row - 1), ]) # return the data up to that row
  } else {
    return(data)
  }

  
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

write_genpop_file <- function(data, filename, title = "Fluidigm", n_loci = 96, loci_names = NULL, ...) { # Write a table with genepop format calls to a file
  if (!missing(loci_names)) {
    # Replace underscores with dashes in loci names
    loci_names <- gsub("_", "-", loci_names)
    initial_data <- c(title, loci_names[1:n_loci], "Pop") # required header data for genepop, use SNP names if present
  } else {
    initial_data <- c(title, 1:n_loci, "Pop") # required header data for genepop, numbers represents the SNPs
  }
  
  
  write.table(initial_data,
    file = filename,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE, ...
  )
  data_with_comma <- data %>%
    select(1:(n_loci+1)) %>% # select n_loci + 1 columns
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
grebs_wider$error_greb_count <- rowSums(select(grebs_wider, starts_with("GREB")) == "Error", na.rm = TRUE)

grebs_wider$early_good_greb_count <- rowSums(select(grebs_wider, any_of(the_good_grebs)) == "Early", na.rm = TRUE)
grebs_wider$late_good_greb_count <- rowSums(select(grebs_wider, any_of(the_good_grebs)) == "Late", na.rm = TRUE)
grebs_wider$het_good_greb_count <- rowSums(select(grebs_wider, any_of(the_good_grebs)) == "Heterozygote", na.rm = TRUE)
grebs_wider$nocall_good_greb_count <- rowSums(select(grebs_wider, any_of(the_good_grebs)) == "No Call", na.rm = TRUE)
grebs_wider$error_good_greb_count <- rowSums(select(grebs_wider, any_of(the_good_grebs)) == "Error", na.rm = TRUE)

write_tsv(grebs_wider, file = greb_output_filename)

# Translate fluidigm data to genepop format ----------------------------
# pull out a key for the snp names
snpkey <- gen_snpkey(data2)

translated_data2 <- translate_fluidigm_to_genepop(data2)

# re-order loci according to keyfile

resorted_data2 <- resort_fl2gp(translated_data2, keyfile = locus_keyfile_name)

qc_data <- resorted_data2 %>%
  mutate(
    missing = rowSums(. == "000000") / (ncol(resorted_data2) - 1)
  ) %>%
  select(samplename, missing)

# Drop rows where missing data is too high
qced_resorted_data2 <- resorted_data2[which(qc_data$missing < max_missing_data), ]
print(paste0("# of samples with too much missing data: ", nrow(resorted_data2) - nrow(qced_resorted_data2)))

# Print samples that were removed
print("Samples removed:")
if (nrow(resorted_data2) - nrow(qced_resorted_data2) > 0) {
  print(qc_data[which(qc_data$missing >= max_missing_data), 1])
} else {
  print("None")
}
# export the genepop file
write_genpop_file(data = qced_resorted_data2, filename = genpop_filename, title = study_name, n_loci = n_loci_cutoff, loci_names = c(colnames(resorted_data2)[-1]))
write.csv(snpkey, snp_keyfile_name)