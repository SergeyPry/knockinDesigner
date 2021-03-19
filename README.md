# CRISPR Knock-in Designer 
This is a repository for the **CRISPR Knock-in Designer** tool.

You can run **CRISPR Knock-in Designer** online through a web interface here: [CRISPR Knock-in Designer](https://crisprtools.shinyapps.io/knockinDesigner/)
 
### If you already have R and/or RStudio installed, you can jump to [here](https://github.com/SergeyPry/knockinDesigner#run_locally) to immmediately start running **CRISPR Knock-in Designer** locally.

# How to Run **CRISPR Knock-in Designer** locally
You will need to have the ability to install software on the computer you are using to run **CRISPR Knock-in Designer** locally; this may require administrator privileges. 

[1. Download and Install R](https://github.com/SergeyPry/knockinDesigner#1-download-and-install-r)

[2. Download and Install RStudio](https://github.com/SergeyPry/knockinDesigner#2-download-and-install-rstudio)

[3. Run **CRISPR Knock-in Designer** locally](https://github.com/SergeyPry/knockinDesigner#run_locally)


## [1. Download and Install R](#1-download-and-install-r)
**CRISPR Knock-in Designer** has been tested directly with 3.5 and 3.6 versions of R, but can be expected to work with any version of R. 

Download R for your appropriate operating system:

[CRAN website for R Download](https://cran.r-project.org/)
Once you have downloaded the R installer, run it to install R.

The most recent R version for Windows that we directly tested is available from [R 3.6.3 download page](https://cran.r-project.org/bin/windows/base/old/3.6.3/).

If you need additional help installing R, please check the installation instructions for your operating system:

Windows:    https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-Windows

Mac OS:     https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-macOS

Linux/Unix: https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-Unix_002dalikes

## [2. Download and Install RStudio](#2-download-and-install-rstudio) (optional)
**CRISPR Knock-in Designer** does not require the use of the RStudio development environment, but if you are interested in examining or modifying the **CRISPR Knock-in Designer** code, it may be easier to do it in RStudio. 

You can download RStudio for free here: https://www.rstudio.com/products/rstudio/download/#download

After downloading the RStudio installer, follow the installation instructions. If you have both R and RStudio installed, you should only do the following steps in RStudio.

## [3. Run **CRISPR Knock-in Designer** locally](#run_locally)
You can run this RShiny web app in R (or RStudio) by opening up an R or RStudio session.

You can copy and paste the code blocks below into your R/RStudio console to run them.

###Run this code ONLY THE FIRST TIME you run this tool on a computer, or when you need to update these packages:

```
#Install packages required to run CRISPR Knock-in Designer; you can also run this code to update these packages

#Install CRAN packages
install.packages(c("BiocManager","stringr", "readr"))
install.packages(c("shiny", "shinyjs", "Rcpp", "xml2") )
install.packages(c("httr", "jsonlite", "shinyFeedback", "shinycssloaders"))
install.packages(c("seqinr", "TmCalculator"))

#Install Biostrings from Bioconductor
library(BiocManager)
BiocManager::install("Biostrings")
```

### Run this code every time you want to use the tool, including the first time:

```
#Load Shiny in the R/RStudio Environment
library(shiny)

#Retrieve, load, and run main CRISPR Knock-in Designer app from GitHub
runGitHub("knockinDesigner", "SergeyPry")

#Retrieve, load, and run tp53 R144H CRISPR Knock-in Designer app from GitHub
runGitHub("knockinDesigner", "SergeyPry", subdir = "tp53_R144H")
```

You're all set!
