# Fluidigm Data Processing in R

This project contains R scripts for processing and analyzing Fluidigm data. The scripts include functions for reading Fluidigm data, generating SNP keys, and translating Fluidigm data to Genepop format.

## Installation

Clone the repository to your local machine:

```bash
git clone git@github.com:Scomeyer/Fluidigm-Genepop.git .
```

## Dependencies
This project requires the following R packages:

- caroline
- tidyverse
- adegenet
- RColorBrewer
- readxl

You can install these packages in R with:
```r
install.packages(c("caroline", "tidyverse", "adegenet", "RColorBrewer", "readxl"))
```

## Usage
The `Flui-Genpop.R` file contains a Settings sections. Be sure to set the working directory, the correct input file names, and the desired output file names within the Settings section, customized to your computer setup.

Set the working directory to the folder where you have gathered the required input files. The output files will also write to the working directory folder.

The primary functions of the R code are defined within the functions section. Some working code is in the bottom half of the code. Some files are not included in the GitHub directory inclusing `Samples_site_order.csv`. However, the Genepop file generation and Greb1l calling code has been tested and is working.
