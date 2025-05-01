# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Engineering

## Dependencies ----
## Amend DESCRIPTION with dependencies read from package code parsing
## install.packages('attachment') # if needed.
attachment::att_amend_desc()
usethis::use_package("ggplot2")
usethis::use_package("dplyr")
usethis::use_package("purrr")
usethis::use_package("stringr")
usethis::use_package("readr")
usethis::use_package("ggtext")
usethis::use_package("phyloseq")
usethis::use_package("phangorn")
usethis::use_package("decontam")
usethis::use_package("dada2")
usethis::use_package("Biostrings")
usethis::use_package("vegan")
usethis::use_package("patchwork")
usethis::use_package("Polychrome")
usethis::use_package("RColorBrewer")
usethis::use_package("egg")
usethis::use_package("rmarkdown")
usethis::use_package("knitr")
usethis::use_package("kableExtra")
usethis::use_package("biomeUtils")
usethis::use_package("biomformat")
usethis::use_package("pander")
usethis::use_package("microbiome")
usethis::use_package("edgeR")
usethis::use_package("vsn")
usethis::use_package("doParallel")
usethis::use_package("GUniFrac")
usethis::use_package("Rtsne")
usethis::use_package("umap")
usethis::use_package("ComplexHeatmap")
usethis::use_package("circlize")
usethis::use_package("ape")
usethis::use_package("dendextend")
usethis::use_package("ggpubr")
usethis::use_package("rstatix")
usethis::use_package("broom")
usethis::use_package("pairwiseAdonis")
usethis::use_package("DESeq2")
usethis::use_package("cowplot")
usethis::use_package("gtable")
usethis::use_package("lemon")
usethis::use_package("ggplotify")
usethis::use_package("writexl")
usethis::use_package("ggrepel")
usethis::use_package("microbiomeMarker")
usethis::use_package("parallel")
usethis::use_package("metacoder")
usethis::use_package("UpSetR")
usethis::use_package("grid")
usethis::use_package("pals")
usethis::use_package("colorspace")
usethis::use_package("data.table")
usethis::use_package("ALDEx2")
usethis::use_package("ggraph")
usethis::use_package("ggforce")
usethis::use_package("job")
usethis::use_package("hexbin")
usethis::use_package("shinydashboard")
usethis::use_package("shinythemes")
usethis::use_package("DT")
usethis::use_package("later")
usethis::use_package("spsComps")
usethis::use_package("plyr")
usethis::use_package("tibble")
usethis::use_package("scales")
usethis::use_package("metapod")
usethis::use_package("igraph")
usethis::use_package("tidyr")
usethis::use_package("fresh")
usethis::use_package("phangorn")
usethis::use_package("vegan")
usethis::use_package("GUniFrac")
usethis::use_package("Rtsne")
usethis::use_package("umap")
usethis::use_package("circlize")
usethis::use_package("dendextend")
usethis::use_package("ggpubr")
usethis::use_package("broom")
usethis::use_package("writexl")
usethis::use_package("ggrepel")
usethis::use_package("UpSetR")
usethis::use_package("pals")
usethis::use_package("ggraph")
usethis::use_package("ggforce")
usethis::use_package("job")
usethis::use_package("hexbin")
usethis::use_package("fpc")
usethis::use_package("readr")
usethis::use_package("RColorBrewer")
usethis::use_package("knitr")
usethis::use_package("grid")
usethis::use_package("colorspace")
usethis::use_package("data.table")
usethis::use_package("decontam")
usethis::use_package("dada2")
usethis::use_package("biomformat")
usethis::use_package("microbiome")
usethis::use_package("edgeR")
usethis::use_package("vsn")
usethis::use_package("ComplexHeatmap")
usethis::use_package("DESeq2")
usethis::use_package("DirichletMultinomial")
usethis::use_package("ALDEx2")
usethis::use_package("mia")
usethis::use_package("ANCOMBC")
usethis::use_package("biomeUtils")
usethis::use_package("pairwiseAdonis")
usethis::use_package("ggord")
## Add modules ----
## Create a module infrastructure in R/
golem::add_module(name = "upload_files_mod", with_test = TRUE) # Name of the module
golem::add_module(name = "name_of_module2", with_test = TRUE) # Name of the module

## Add helper functions ----
## Creates fct_* and utils_*
golem::add_fct("helpers", with_test = TRUE)
golem::add_utils("custom_theme", with_test = TRUE)

## External resources
## Creates .js and .css files at inst/app/www
golem::add_js_file("script")
golem::add_js_handler("handlers")
golem::add_css_file("custom")
golem::add_sass_file("custom")

## Add internal datasets ----
## If you have data in your package
usethis::use_data_raw(name = "my_dataset", open = FALSE)

## Tests ----
## Add one line by test you want to create
usethis::use_test("app")

# Documentation

## Vignette ----
usethis::use_vignette("golemMB")
devtools::build_vignettes()

## Code Coverage----
## Set the code coverage service ("codecov" or "coveralls")
usethis::use_coverage()

# Create a summary readme for the testthat subdirectory
covrpage::covrpage()

## CI ----
## Use this part of the script if you need to set up a CI
## service for your application
##
## (You'll need GitHub there)
usethis::use_github()

# GitHub Actions
usethis::use_github_action()
# Chose one of the three
# See https://usethis.r-lib.org/reference/use_github_action.html
usethis::use_github_action_check_release()
usethis::use_github_action_check_standard()
usethis::use_github_action_check_full()
# Add action for PR
usethis::use_github_action_pr_commands()

# Travis CI
usethis::use_travis()
usethis::use_travis_badge()

# AppVeyor
usethis::use_appveyor()
usethis::use_appveyor_badge()

# Circle CI
usethis::use_circleci()
usethis::use_circleci_badge()

# Jenkins
usethis::use_jenkins()

# GitLab CI
usethis::use_gitlab_ci()

# You're now set! ----
# go to dev/03_deploy.R
rstudioapi::navigateToFile("dev/03_deploy.R")
