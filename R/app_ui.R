#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import phyloseq
#' @import phangorn
#' @import decontam
#' @import dada2
#' @import Biostrings
#' @import vegan
#' @import patchwork
#' @import Polychrome
#' @import RColorBrewer
#' @import egg
#' @import rmarkdown
#' @import knitr
#' @import kableExtra
#' @import biomeUtils
#' @import biomformat
#' @import pander
#' @import microbiome
#' @import edgeR
#' @import vsn
#' @import doParallel
#' @import GUniFrac
#' @import Rtsne
#' @importFrom umap umap
#' @import ComplexHeatmap
#' @import circlize
#' @import ape
#' @import dendextend
#' @import ggpubr
#' @import rstatix
#' @import broom
#' @import pairwiseAdonis
#' @import cowplot
#' @import gtable
#' @import lemon
#' @import ggplotify
#' @import writexl
#' @import ggrepel
#' @import microbiomeMarker
#' @import parallel
#' @import metacoder
#' @import UpSetR
#' @import grid
#' @import pals
#' @import colorspace
#' @import data.table
#' @import ALDEx2
#' @import ggraph
#' @import ggforce
#' @import job
#' @import hexbin
#' @import shinydashboard
#' @import shinythemes
#' @import DT
#' @import later
#' @import spsComps
#' @import plyr
#' @import tibble
#' @import scales
#' @import metapod
#' @import igraph
#' @import tidyr
#' @import fresh
#' @import fpc
#' @import DirichletMultinomial
#' @import mia
#' @import ANCOMBC
#' @import ggord
#' @import ggplot2
#' @import dplyr
#' @import purrr
#' @import stringr
#' @import readr
#' @import ggtext
#' @noRd

# # load libraries
# library(shinythemes)
# library(shinydashboard)
# library(fresh)
# library(DT)
# library(shiny)
# library(spsComps)
# library(purrr)
# library(bslib)

# source text list
#source("/R/helper_functions/text_list.R")

app_ui <- function(request) {
  tagList(
    # adding external resources
    golem_add_external_resources(),
    tags$head(
      # Inline CSS for immediate testing
      tags$style(HTML("
        body {
            color: white !important;
            background-color: black !important; /* Optional for contrast */
        }
        .box-body {
            color: #ffffff !important; /* Ensure white text inside boxes */
        }
      ", "
                      div.datatables {
                      color: white !important;}"))
    ),
    # UI logic
    shinydashboard::dashboardPage(
      skin = NULL,
      dashboardHeader(title = "golemMB v1.0"), #Amplicore #CladeAid
      dashboardSidebar(sidebar_menu()),
      dashboardBody(
        use_theme(custom_theme()),
        tab_items()
      )
    )
  )
}

# sidebar menu function
sidebar_menu <- function() {
  sidebarMenu(
    menuItem("Read me!", tabName = "readme", icon = icon("info-circle")),
    menuItem("Upload Files", tabName = "upload_files", icon = icon("folder-open"),
    menuSubItem("OTU", tabName = "uOTU"),
    menuSubItem("ZOTU", tabName = "uZOTU")),
    menuItem("Set Parameters", tabName = "set_parameters", icon = icon("sliders-h")),
    menuItem("Data Pre-Processing", tabName = "pre_processing", icon = icon("filter"),
    menuSubItem("Import Data to Phyloseq", tabName = "import_phyloseq"),
    menuSubItem("Count Transform Data", tabName = "count_transform")),
    menuItem("Alpha Diversity", tabName = "alpha_diversity", icon = icon("chart-bar")),
    menuItem("Beta Diversity", tabName = "b_diversity", icon = icon("circle-nodes"),
    menuSubItem("Calculate Distance Matrices", tabName = "distance_matrices"),
    menuSubItem("Beta Diversity", tabName = "beta_diversity")),
    menuItem("Differential Abundance", tabName = "differential_abundance", icon = icon("chart-line")),
    menuItem("Functional Analyses", tabName = "functional_analyses", icon = icon("dna"),
    menuSubItem("Functional Abundance", tabName = "functional_abundance"),
    menuSubItem("Phenotype Prediction", tabName = "bugbase")),
    menuItem("Correlation Analysis", tabName = "correlations", icon = icon("table"))
        )
}

# tab items function
tab_items <- function() {
  tabItems(
    tabItem(tabName = "readme", create_readme_tab()),
    tabItem(tabName = "upload_files"),
    tabItem(tabName = "uOTU", create_OTU_tab()),
    tabItem(tabName = "uZOTU", create_ZOTU_tab()),
    tabItem(tabName = "set_parameters", create_set_parameters_tab()),
    tabItem(tabName = "import_phyloseq", create_import_tab()),
    tabItem(tabName = "alpha_diversity", create_alpha_diversity_tab()),
    tabItem(tabName = "count_transform", create_count_transform_tab()),
    tabItem(tabName = "bugbase", create_bugbase_tab()),
    tabItem(tabName = "distance_matrices", create_distance_tab()),
    tabItem(tabName = "beta_diversity", create_beta_diversity_tab()),
    tabItem(tabName = "differential_abundance", create_differential_abundance_tab()),
    tabItem(tabName = "functional_abundance", create_functional_abundance_tab()),
    tabItem(tabName = "correlations", create_correlations_tab())
        )
}

### readme tab
create_readme_tab <- function() {
  tagList(
    h2("Welcome!"),
    p("This 16S analysis app was created within the scope of my master's thesis. The purpose of this app is
              to analyse 16S data in a reproducible and contained environment, i.e. consistent RStudio versioning, package versioning, dependencies etc."),
    br(),
    h4("Requirements:"),
    tags$ul(
      tags$li("(Z)OTU table(s) containing taxonomic assigments"),
      tags$li("Tree file(s) (see ", a('https://www.arb-silva.de/'),")"),
      tags$li("fasta file(s)"),
      tags$li("meta data file in .rds format")
    ),
    br(),
    p("This app will NOT generate these files."),
    p("Once you uploaded your main files and selected your variables you can progressively analyse your data. Your results will be saved
      as a PDF report and in .rds files. Some analyses will also seperately save the resulting figures. As you analyze your data you will
      be required to upload your .rds results to continue.")

  )
}

### upload files tab
create_upload_files_tab <- function() {
    fluidRow(
      box(width = 12, title = "Upload Files", status = "primary", solidHeader = TRUE,
          column(6, fileInput("uploadOTUtab", "Upload OTU file", accept = ".tab")),
          column(6, fileInput("uploadZOTUtab", "Upload ZOTU file", accept = ".tab")),
          column(6,fileInput("uploadOTUtre", "Upload OTU Tree file", accept = ".tre")),
          column(6,fileInput("uploadZOTUtre", "Upload ZOTU Tree file", accept = ".tre")),
          column(6,fileInput("uploadOTUfasta", "Upload OTU fasta file", accept = ".fasta")),
          column(6,fileInput("uploadZOTUfasta", "Upload ZOTU fasta file", accept = ".fasta")),
          column(6, br(), fileInput("uploadmeta", "Upload meta file", accept = ".rds")),
          column(6, br(), fileInput("uploadbugbasepredicts", "Upload BugBase file", accept = ".txt")),
          column(6, br(), fileInput("uploadbugbasevsearch", "Upload BugBase vsearch stats", accept = ".txt")),
          column(12, tags$hr()),
          column(12, p("The necessary files are aquired with each analysis. You will need to upload
                       them after every download to progress.")),
          column(6,fileInput("uploadOTUphyseq", "Upload OTU phyloseq file", accept = ".rds")),
          column(6,fileInput("uploadZOTUphyseq", "Upload ZOTU phyloseq file", accept = ".rds")),
          column(6,fileInput("uploadOTUalpha", "Upload OTU alpha diversity file", accept = ".rds")),
          column(6,fileInput("uploadZOTUalpha", "Upload ZOTU alpha diversity file", accept = ".rds")),
          column(6,fileInput("uploadOTUdists", "Upload OTU distance matrices file", accept = ".rds")),
          column(6,fileInput("uploadZOTUdists", "Upload ZOTU distance matrices file", accept = ".rds")),
          column(6,fileInput("uploadOTUdistsgenus", "Upload OTU distance matrices (Genus) file", accept = ".rds")),
          column(6,fileInput("uploadZOTUdistsgenus", "Upload ZOTU distance matrices (Genus) file", accept = ".rds")),
          column(6,fileInput("uploadOTUord", "Upload OTU ordination file", accept = ".rds")),
          column(6,fileInput("uploadZOTUord", "Upload ZOTU ordination file", accept = ".rds")),
          column(6,fileInput("uploadOTUordgenus", "Upload OTU ordination (Genus) file", accept = ".rds")),
          column(6,fileInput("uploadZOTUordgenus", "Upload ZOTU ord (Genus) file", accept = ".rds")),
          column(6,fileInput("uploadOTUtrans", "Upload transformed OTUs file", accept = ".rds")),
          column(6,fileInput("uploadZOTUtrans", "Upload transformed ZOTUs file", accept = ".rds")),
          column(6,fileInput("uploadOTUtransgenus", "Upload transformed OTUs (Genus) file", accept = ".rds")),
          column(6,fileInput("uploadZOTUtransgenus", "Upload transformed ZOTUs (Genus) file", accept = ".rds")),
          column(6,fileInput("uploadpicrustfolder", "Upload your ZIPPED PICRUSt2 output folder", accept = ".zip")),
          column(6, fileInput("uploadpicruststats", "Upload your PICRUSt2 stats file", accept = ".txt"))

      )
  )
}


### upload OTU files tab
create_OTU_tab <- function() {
  tagList(
  column(12, h3("Main OTU files")),
  column(12, p("Upload your files here.")),
  column(6, fileInput("uploadOTUtab", "Upload OTU file", accept = ".tab")),
  column(6, fileInput("uploadOTUtre", "Upload OTU Tree file", accept = ".tre")),
  column(6, fileInput("uploadOTUfasta", "Upload OTU fasta file", accept = ".fasta")),
  column(6, fileInput("uploadmeta", "Upload meta file", accept = ".rds")),
  br(),
  column(12, h3("Analysis results")),
  column(12, p("Upload your result files that you generated throughout your analysis here.")),
  column(6, fileInput("uploadOTUphyseq", "Upload OTU phyloseq file", accept = ".rds")),
  column(6, fileInput("uploadOTUalpha", "Upload OTU alpha diversity file", accept = ".rds")),
  column(6, fileInput("uploadOTUdists", "Upload OTU distance matrices file", accept = ".rds")),
  column(6, fileInput("uploadOTUdistsgenus", "Upload OTU distance matrices (Genus) file", accept = ".rds")),
  column(6, fileInput("uploadOTUord", "Upload OTU ordination file", accept = ".rds")),
  column(6, fileInput("uploadOTUordgenus", "Upload OTU ordination (Genus) file", accept = ".rds")),
  column(6, fileInput("uploadOTUtrans", "Upload transformed OTUs file", accept = ".rds")),
  column(6, fileInput("uploadOTUtransgenus", "Upload transformed OTUs (Genus) file", accept = ".rds")),
  column(6, br(), fileInput("uploadbugbasepredicts", "Upload BugBase predictions file", accept = ".txt")),
  column(6, br(), fileInput("uploadbugbasevsearch", "Upload BugBase vsearch stats", accept = ".txt")),
  )

}

### upload ZOTU files tab
create_ZOTU_tab <- function() {
  tagList(
    column(12, h3("Main ZOTU files")),
    column(12, p("Upload your ZOTU files here.")),
    column(6, fileInput("uploadZOTUtab", "Upload ZOTU file", accept = ".tab")),
    column(6, fileInput("uploadZOTUtre", "Upload ZOTU Tree file", accept = ".tre")),
    column(6, fileInput("uploadZOTUfasta", "Upload ZOTU fasta file", accept = ".fasta")),
    br(),
    column(12, h3("Analysis results")),
    column(12, p("Upload your result files that you generated throughout your analysis here.")),
    column(6, fileInput("uploadZOTUphyseq", "Upload ZOTU phyloseq file", accept = ".rds")),
    column(6, fileInput("uploadZOTUalpha", "Upload ZOTU alpha diversity file", accept = ".rds")),
    column(6, fileInput("uploadZOTUdists", "Upload ZOTU distance matrices file", accept = ".rds")),
    column(6, fileInput("uploadZOTUdistsgenus", "Upload ZOTU distance matrices (Genus) file", accept = ".rds")),
    column(6, fileInput("uploadZOTUord", "Upload ZOTU ordination file", accept = ".rds")),
    column(6, fileInput("uploadZOTUordgenus", "Upload ZOTU ord (Genus) file", accept = ".rds")),
    column(6, fileInput("uploadZOTUtrans", "Upload transformed ZOTUs file", accept = ".rds")),
    column(6, fileInput("uploadZOTUtransgenus", "Upload transformed ZOTUs (Genus) file", accept = ".rds")),
    column(6,fileInput("uploadpicrustfolder", "Upload your ZIPPED PICRUSt2 output folder", accept = ".zip")),
    column(6, fileInput("uploadpicruststats", "Upload your PICRUSt2 stats file", accept = ".txt"))

  )
}

### set parameters tab
create_set_parameters_tab <- function() {
  fluidRow(
    box(width = 12, title = "Meta File", status = "primary", solidHeader = TRUE, collapsible = TRUE,
        p("You may take a look at your meta file if you're unsure which variables to select."),
        br(), br(), br(),
        DT::dataTableOutput("meta")),
    box(width = 12, title = "Save/Restore settings", status = "primary", solidHeader = TRUE,
        p("You can bookmark your session and save your settings here. Clicking
                      the button will give you a URL of the immediate state of your current session, which you can
                      then re-connect to whenever the app is running. You will find a folder named
                      'shiny_bookmarks' storing your uploaded files and inputs in your working directory. I strongly recommend keeping a .txt file
                      tracking all your URLs and to periodically clear your subdirectory as to not bloat your disk."),
        #actionButton("save_inputs", "Save settings"),
        bookmarkButton(),
        #actionButton("load_inputs", "Load settings")
    ),

    box(width = 4, title = "Set analysis parameters", status = "primary", solidHeader = TRUE,
        fluidRow(
          column(12, selectInput("analysis_type", "Analysis type", choices = c("OTU", "ZOTU", "both"))),
          column(12, selectInput("var_test", "Analysis variable", choices = NULL)),
          column(12, textInput("var_formula", "Variable formula")),
          column(12, numericInput("min_sample_counts", "Minimum sample counts", value = 5000, min = 5000, step = 1000)),
          column(12, selectInput("count_transform", "Normalization method (import only!)", choices = c("compositional", "Z", "log10", "log10p", "hellinger", "identity", "CLR", "ALR"), selected = "log10p")),
          column(12, selectInput("samples_outlier", "Excluding following outlier(s):", choices = NULL, multiple = TRUE)))),

    box(width = 4, title = "Figure Options", status = "primary", solidHeader = TRUE,
        fluidRow(
          column(12,selectInput("point_label", "Label datapoints by this variable:", choices = NULL)),
          column(12, textInput("var_name", "Analysis variable name")),
          column(12,selectInput("shape_var", "Shape datapoints by this variable:", choices = NULL)),
          column(12,selectInput("label_set", "Choose your variable order", choices = NULL, multiple = TRUE)),
          column(12,textInput("labels", "Name your variables for figures later. Separate by comma and mind the order!")),
          column(12,selectInput("color_set", "Select variable colors and mind the order!", choices =
                                  c("grey60", "black", "tomato", "red4", "skyblue", "deepskyblue4", "olivedrab3", "olivedrab",
                                    "pink", "violetred", "orchid1", "darkmagenta", "cyan", "aquamarine4", "yellow2", "yellow4",
                                    "lightsalmon", "chocolate"), multiple = TRUE))


        )),

    box(width = 4, title = "Set comparison list", status = "primary", solidHeader = TRUE,
        fluidRow(
          column(8, selectizeInput("comp_list", "Compare these two variables:", choices = NULL, multiple = TRUE, options = list(maxItems = 2))),
          actionButton("add_comp", NULL, style = "border-radius: 0px;
                                             border-width: 0px; background-color: transparent; margin-top: 25px", icon("circle-plus")),
          # actionButton("rmv_comp", NULL, style = "border-radius: 0px;
          #                        border-width: 0px; background-color: transparent; margin-top: 25px", icon("circle-minus")),
          column(8, div(id = "placeholder")),
          column(8, div(id = "placeholder2")),
          column(12, p("Tip: you are comparing something to your REFERENCE. That means the reference is always at the second position."))
        )
    )
  )
}

### import phyloseq tab
create_import_tab <- function() {
  fluidRow(
    box(width = 12, title = "Import data to phyloseq", status = "primary", solidHeader = TRUE,
        h4("Importing your data to a phyloseq object"),
        p("The first step to processing your data is to convert it to a phyloseq object. phyloseq is an R package that allows you easily wrangle,
                      filter, store and plot your phylogenetic microbial NGS data. Here, negative controls, low-performing samples and leftover eukaryotic sequences
                      will be filtered out. The resulting PDF will give you an overview of your filtered data, including a representation of your sequencing depth,
                      sample clustering, sequence length distribution and phylum prevalence in your samples."),
        br(),
        h4("Requirements to import your data:"),
        tags$ul(
          tags$li("(Z)OTU table"),
          tags$li("(Z)OTU tree file"),
          tags$li("(Z)OTU fasta file"),
          tags$li("meta file"),
          tags$li("excluded samples"),
          tags$li("analysis variable"),
          tags$li("analysis variable name"),
          tags$li("variable shapes"),
          tags$li("color set"),
          tags$li("variable label order"),
          tags$li("variable labels"),
          tags$li("minimum sample counts")
        ),
        br(),
        actionButton("run", "Run Import", style = "width: 200px"),
        downloadButton("download_import", "download phyloseq object", style = "width: 200px")

    )
  )
}

### alpha diversity tab
create_alpha_diversity_tab <- function() {
  fluidRow(
    box(width = 12, title = "Alpha Diversity", status = "primary", solidHeader = TRUE,
        p("Calculate alpha diversity of your samples in the following metrics: Chao1, ACE, Shannon,
                    Simpson, Inverse Simpson, and observed taxa."),
        h4("Requirements to calculate alpha diversity:"),
        tags$ul(
          tags$li("meta file"),
          tags$li("analysis variable"),
          tags$li("analysis variable color (see below). This will color code the factor levels of you analysis variable.")),
        br(),
        selectInput("var_color", "Color alpha diversity by variable", choices = NULL),
        actionButton("calc_alpha", "Calculate alpha diversity", style = "width: 200px"),
        downloadButton("download_alpha1", "download alpha diversity", style = "width: 200px"),
        br(), br(),
        tags$hr(),
        h4("Requirements to plot alpha diversities:"),
        tags$ul(
          tags$li("calculated alpha diversity"),
          tags$li("analysis type"),
          tags$li("analysis variable"),
          tags$li("analysis variable name"),
          tags$li("analysis variable formula"),
          tags$li("comparisons"),
          tags$li("variable label order"),
          tags$li("comparison variable names"),
          tags$li("color set")
        ),
        br(),
        actionButton("plot_alpha", "Plot alpha diversity", style = "width: 200px"),
        downloadButton("download_alpha2", "download alpha plot", style = "width: 200px")

    )
  )
}

### count transform function
create_count_transform_tab <- function() {
  fluidRow(
    box(width = 12, title = "Count Transform Data", status = "primary", solidHeader = TRUE,
        p("Transform your data with the following methods: relative log expression normalization (RLE),
          regularized log normalization (Rlog), variance stabilizing transformation (vst), centered log ratio
          transformation (CLR) and a modified CLR transformation (CLR2)."),
        h4("Requirements to count transform your data:"),
        tags$ul(
          tags$li("imported your data to phyloseq"),
          tags$li("analysis type"),
          tags$li("analysis variable"),
          tags$li("analysis variable name")
        ),
        br(),
        actionButton("count_transform", "Transform data", style = "width: 200px"),
        downloadButton("download_ct", "download transformed data", style = "width: 200px")
    )
  )
}

### bugbase function
create_bugbase_tab <- function() {
  fluidRow(
    box(width = 12, title = "Functional prediction via BugBase", status = "primary", solidHeader = TRUE,
        h4("What is BugBase?"),
        p("BugBase is a phenotype prediction algorithm for microbial sequencing data. Predictions for
                      bacterial characteristics include:"),
        tags$ul(
          tags$li("oxygen tolerance"),
          tags$li("Gram staining"),
          tags$li("pathogenicity"),
          tags$li("biofilm formation"),
          tags$li("mobile element content"),
          tags$li("oxidative stress tolerance")),
        p("For more see the preprint", a("https://doi.org/10.1101/133462")),
        p(strong("NOTE:"), "The BugBase algorithm is not implemented in this app. You can only plot your results.
                      See ", a("https://bugbase.cs.umn.edu/"), "to run your data through the BugBase pipeline.", strong("The resulting file
                      'predictions.txt' is required for plotting.")),
        br(),
        h4("Input requirements:"),
        tags$ul(
          tags$li("analysis variable"),
          tags$li("analysis variable name"),
          tags$li("analysis variable formula"),
          tags$li("comparisons"),
          tags$li("variable label order"),
          tags$li("comparison variable names"),
          tags$li("color set")),
        br(),
        actionButton("plot_bb", "plot BugBase", style = "width: 200px"),
        downloadButton("download_bugbase", "download BugBase plots", style = "width: 200px")
    )
  )
}

### distance matrices function
create_distance_tab <- function() {
  fluidRow(
    box(width = 12, title = "Calculate Distance Matrices", status = "primary", solidHeader = TRUE,
        p("Calculate distances to plot beta diversity with dimensionality reduction techniques. Distance
          metrics include Jensen-Shannon Divergence, Gower's Distance, Bray-Curtis
          Dissimilaity, Double Principle Coordinate Analysis, unweighted, weighted,
          variance-adjusted weighted and generalized UniFrac distances."),
        h4("Requirements to calculate distance matrices:"),
        tags$ul(
          tags$li("transformed your data"),
          tags$li("analysis type"),
          tags$li("analysis variable"),
          tags$li("analysis variable name"),
          tags$li("variable order AND names"),
          tags$li("variable colors")
        ),
        br(),
        actionButton("dist_ma", "Calculate distance matrices", style = "width: 200px"),
        downloadButton("download_dm", "download distance files", style = "width: 200px")
    )
  )
}

### beta diversity function
create_beta_diversity_tab <- function() {
  fluidRow(
    box(width = 12, title = "Beta Diversity", status = "primary", solidHeader = TRUE,
        h4("Requirements to calculate and plot beta diversity:"),
        tags$ul(
          tags$li("transformed your data"),
          tags$li("calculated distance matrices"),
          tags$li("analysis type"),
          tags$li("analysis variable"),
          tags$li("analysis variable name"),
          tags$li("analysis variable formula"),
          tags$li("variable colors"),
          tags$li("variable label order AND names"),
          tags$li("variable shape")
        ),
        br(),
        actionButton("plot_beta", "Plot Beta Diversity", style = "width: 200px"),
        downloadButton("download_betadiversity", "download beta diversity results", style = "width: 200px")
    )
  )
}

### differential abundance function
create_differential_abundance_tab <- function() {
  fluidRow(
    box(width = 12, title = "Differential Abundance", status = "primary", solidHeader = TRUE,
        h4("Requirements to calculate and plot differential abundance:"),
        tags$ul(
          tags$li("transformed your data"),
          tags$li("analysis type"),
          tags$li("analysis variable"),
          tags$li("variable comparison list"),
          tags$li("variable order AND names"),
          tags$li("variable colors")
        ),
        br(),
        actionButton("plot_DA", "Plot Differential Abundance", style = "width: 200px"),
        downloadButton("download_diffabund", "download results", style = "width: 200px")
    )
  )
}

### functional abundance function
create_functional_abundance_tab <- function() {
  fluidRow(
    box(width = 12, title = "Functional Abundance", status = "primary", solidHeader = TRUE,
        h4("Requirements to calculate and plot functional abundance:"),
        tags$ul(
          tags$li("imported your data to phyloseq"),
          tags$li("picrust_2_out_stratified folder (see above)"),
          tags$li("analysis type"),
          tags$li("analysis variable"),
          tags$li("variable comparison list"),
          tags$li("variable order AND names"),
          tags$li("variable colors")
        ),
        br(),
        actionButton("plot_FA", "Plot Functional Abundance", style = "width: 200px"),
        downloadButton("download_funcabund", "download results", style = "width: 200px")
    )
  )
}

### correlations function
create_correlations_tab <- function() {
  fluidRow(
    box(width = 12, title = "Correlations", status = "primary", solidHeader = TRUE,
        h4("Calculating correlations between independent meta variables and TUs"),
        p("This script allows you to calculate Pearson correlations between your
                          independent meta variables and TU abundances. Your output will be a
                          TSV file with your TU identifiers, your independent variables,
                          correlation coefficients, effect sizes, p-values and adjusted
                          p-values (after Benjamini-Hochberg) and a TSV file containing
                          the taxonomic classification of each TU."),
        br(),
        p("The abundance filter removes TUs in a sample with a relative abundance below the cutoff.
                          TUs that have zero abundance across all samples will be excluded from further calculations.
                          The default (recommended) value is 0.5%."),
        numericInput("abundance_filter", "Set abundance cutoff (%)", value = 0.5, min = 0, max = 1, step = 0.1),
        p("The prevalence filter removes TUs with an abundance below the cutoff relative to all samples.
                          The default value is 0.3 (30%)."),
        numericInput("prevalence_filter", "Set prevalence cutoff", value = 0.3, min = 0, max = 0.7, step = 0.1),
        br(),
        p(strong("Filter Pearson correlation table")),
        p("You can filter your output by setting cutoffs for correlation coefficients,
                          p-values, and effect sizes. The effect size cutoff depends on your data."),
        numericInput("r_cutoff", "Set correlation cutoff (default |0.5|)", value = 0.5, step = 0.1),
        numericInput("pval_cutoff", "Set p-value cutoff (default 0.05)", value = 0.05, min = 0, step = 0.01),
        numericInput("support_cutoff", "Set effect-size cutoff (default 4)", value = 4, min = 0, step = 1),
        h4("Requirements for calculating correlations:"),
        tags$ul(
          tags$li("imported your data to phyloseq"),
          tags$li("your meta file")
        ),
        br(),
        actionButton("calc_corr", "Calculate correlations", style = "width: 200px"),
        downloadButton("download_corr", "download results", style = "width: 200px")
    ))
}



#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "golemMB"
    ),
    # Include custom CSS
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
