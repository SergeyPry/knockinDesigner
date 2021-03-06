# CRISPR Knock-in Designer 
This is a repository for the **CRISPR Knock-in Designer** tool.

You can run **CRISPR Knock-in Designer** online through a web interface here: [CRISPR Knock-in Designer](https://crisprtools.shinyapps.io/knockinDesigner/).
The instructions and explanations of what the online version of the tool does are available on its website.
 
### If you already have R and/or RStudio installed, you can immediately start running **CRISPR Knock-in Designer** locally.

# How to Run **CRISPR Knock-in Designer** locally
You will need to have the ability to install software on the computer you are using to run **CRISPR Knock-in Designer** locally; this may require administrator privileges. 

1. [Download and Install R](https://github.com/SergeyPry/knockinDesigner#1-download-and-install-r)

2. [Download and Install RStudio](https://github.com/SergeyPry/knockinDesigner#2-download-and-install-rstudio-optional)

3. [Run CRISPR Knock-in Designer locally](https://github.com/SergeyPry/knockinDesigner#3-run-crispr-knock-in-designer-locally)


## 1. Download and Install R
**CRISPR Knock-in Designer** has been tested directly with 3.5 and 3.6 versions of R, but can be expected to work with any version of R. 

Download R for your appropriate operating system: [CRAN website for R Download](https://cran.r-project.org/)
Once you have downloaded the R installer, run it to install R.

The most recent R version for Windows that we directly tested is available from [R 3.6.3 download page](https://cran.r-project.org/bin/windows/base/old/3.6.3/).

If you need additional help installing R, please check the installation instructions for your operating system:

Windows:    https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-Windows

Mac OS:     https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-macOS

Linux/Unix: https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-Unix_002dalikes

## 2. Download and Install RStudio (optional)
**CRISPR Knock-in Designer** does not require the use of the RStudio development environment, but if you are interested in examining or modifying the **CRISPR Knock-in Designer** code, it may be easier to do it in RStudio. 

You can download RStudio for free here: https://www.rstudio.com/products/rstudio/download/#download

After downloading the RStudio installer, follow the installation instructions. If you have both R and RStudio installed, you should only do the following steps in RStudio.

## 3. **Run CRISPR Knock-in Designer locally**
You can run this RShiny web app in R (or RStudio) by opening up an R or RStudio session.

You can copy and paste the code blocks below into your R/RStudio console to run them.

### Run this code ONLY THE FIRST TIME you run this tool on a computer, or when you need to update these packages:

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

The following code runs the main version of the **CRISPR Knock-in Designer** app locally.
```
#Load Shiny in the R/RStudio Environment
library(shiny)

#Retrieve, load, and run main CRISPR Knock-in Designer app from GitHub
runGitHub("knockinDesigner", "SergeyPry")
```

If you want to run one of the demo versions of the app with data embedded in the code, please try one of the following lines.
```
#Load Shiny in the R/RStudio Environment
library(shiny)

#Retrieve, load, and run tp53 R144H CRISPR Knock-in Designer app demo from GitHub
runGitHub("knockinDesigner", "SergeyPry", subdir = "tp53_R144H")

#Retrieve, load, and run lmna R471L CRISPR Knock-in Designer app demo from GitHub
runGitHub("knockinDesigner", "SergeyPry", subdir = "lmna_R471L")

#Retrieve, load, and run cacna1c G419R CRISPR Knock-in Designer app demo from GitHub
runGitHub("knockinDesigner", "SergeyPry", subdir = "cacna1c_G419R")

#Retrieve, load, and run iqgap1 R495S CRISPR Knock-in Designer app demo from GitHub
runGitHub("knockinDesigner", "SergeyPry", subdir = "iqgap1_R495S")

#Retrieve, load, and run lef1 G179V CRISPR Knock-in Designer app demo from GitHub
runGitHub("knockinDesigner", "SergeyPry", subdir = "lef1_G179V")

#Retrieve, load, and run tp53 S277R CRISPR Knock-in Designer app demo from GitHub
runGitHub("knockinDesigner", "SergeyPry", subdir = "tp53_S277R")
```

The alternative to the code above is to download the entire repository, go to the appropriate folder on your own computer, open an app.R file in RStudio and run it.


You're all set!
