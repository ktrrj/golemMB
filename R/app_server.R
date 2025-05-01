#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
#'
app_server <- function(input, output, session) {
  # Your application server logic
  # track the number of comp_list widgets
  #set upload limit to 2GB
  options(shiny.maxRequestSize = 2048 * 1024^2)

  count <- reactiveValues(comp_count = 1, comp_ids = list())


  #access rmd file folder
  rmd_file_path <- function(rmd_file) {
    file.path("R/rmd", rmd_file)
  }

  # source helper functions
  # source("R/helper_functions/download_function.R")
  # source("R/helper_functions/check_function2.R", local = TRUE)

  ###########will add properly later##################

  filter_var <- "infection"
  author <- "Dr. Roman Gerlach"
  workdir <- "U:/00_Shiny/Test"
  samples_excluded <- c()



  ####################################################


  # pull variable input options from uploaded meta file
  observeEvent(input$uploadmeta, {       #if metafile is uploaded the variable options below will update

    meta_file <- readRDS(input$uploadmeta$datapath)

    updateSelectInput(session, "var_test", choices = names(meta_file))
    updateSelectInput(session, "point_label", choices = c((names(meta_file)), "none"))
    updateSelectInput(session, "shape_var", choices = c((names(meta_file)), "none"))
    updateSelectInput(session, "samples_outlier", choices = c(meta_file[meta_file$samp_ctrl != "neg_ctrl", "sampleID"], "none"))
    updateSelectInput(session, "var_color", choices = names(meta_file))

    #update comp_list and label_set options with var_test factor levels
    observe({

      updateSelectInput(session, "comp_list", choices = levels(meta_file[[input$var_test]]))
      updateSelectInput(session, "label_set", choices = levels(meta_file[[input$var_test]]))

    })
    output$meta <- DT::renderDataTable({
      meta <- input$uploadmeta
      meta_path <- meta$datapath
      readRDS(meta_path)
    })

  })


  ### ~~~set global variables~~~ ###



  # add comparison vector as UI element
  # add counter that tracks and creates the number/ID of newly created comparison vectors
  observeEvent(input$add_comp, {

    req(input$uploadmeta)

    count$comp_count <- count$comp_count + 1
    count_id <- paste0("comp_list", count$comp_count)
    count$comp_ids <- c(count$comp_ids, count_id)

    # insert new comparison vector

    insertUI(
      selector = "#placeholder",
      where = "beforeEnd",
      ui = selectInput(count_id, paste("Compare these two variables", count$comp_count), choices = NULL, multiple = TRUE))

    # update choices with factor levels from var_test
    later::later(function() {
      observe({
        for (comp_id in count$comp_ids) {
          meta_file <- readRDS(input$uploadmeta$datapath)
          updateSelectizeInput(session, count_id, choices = levels(meta_file[[input$var_test]]), options = list(maxItems = 2))

        }
      })
    }, 0.1)





  }, ignoreInit = T)



### ~~~run import~~~ ###

observeEvent(input$run, {

  output_dir1 <- file.path(tempdir(), "phyloseq")

  # create list of all required inputs
  missing_inputs_import <- list("Metadata file" = input$uploadmeta$datapath,
                                "analysis variable" = input$var_test,
                                "analysis variable name" = input$var_name,
                                "variable shape" = input$shape_var,
                                "point label" = input$point_label,
                                "label order" = input$label_set,
                                "label labels" = input$labels,
                                "color set" = input$color_set,
                                "analysis type" = input$analysis_type,
                                "Minimum sample counts" = input$min_sample_counts) #input$count_transform) #input$samples_excluded, input$samples_outlier

  if (input$analysis_type == "OTU") {
    missing_inputs_import <- c(missing_inputs_import, list(
      "OTU table file" = input$uploadOTUtab$datapath,
      "OTU fasta file" = input$uploadOTUfasta$datapath,
      "OTU tree file" = input$uploadOTUtre$datapath))
  } else

    if (input$analysis_type == "ZOTU") {
      missing_inputs_import <- c(missing_inputs_import, list(
        "ZOTU table file" = input$uploadZOTUtab$datapath,
        "ZOTU fasta file" = input$uploadZOTUfasta$datapath,
        "ZOTU tree file" = input$uploadZOTUtre$datapath))
    }


  if (input$analysis_type == "both") {
    missing_inputs_import <- c(missing_inputs_import, list(
      "OTU table file" = input$uploadOTUtab$datapath,
      "OTU fasta file" = input$uploadOTUfasta$datapath,
      "OTU tree file" = input$uploadOTUtre$datapath,
      "ZOTU table file" = input$uploadZOTUtab$datapath,
      "ZOTU fasta file" = input$uploadZOTUfasta$datapath,
      "ZOTU tree file" = input$uploadZOTUtre$datapath))
  }

  # filter list for missing inputs
  missing_inputs_import <- lapply(missing_inputs_import, function(x) {
    if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
      NULL
    } else {
      x
    }
  })

  # pull names of missing inputs to paste into error message
  missing <- names(missing_inputs_import)[sapply(missing_inputs_import, is.null)]

  # catch and print missing inputs
  if (length(missing) > 0) {
    spsComps::shinyCatch(

      stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")

    return()
  }


  if (input$samples_outlier == "none") {
    samples_outlier <- c()
  } else {
    samples_outlier <- input$samples_outlier
  }

  # start Import script

  params_import <- list(
    output_file = list("01_import_OTU.pdf", "01_import_ZOTU.pdf"),
    TU_type = list("OTU", "ZOTU"),
    TU_file = list(input$uploadOTUtab$datapath, input$uploadZOTUtab$datapath),
    TU_tree_file = list(input$uploadOTUtre$datapath, input$uploadZOTUtre$datapath),
    TU_fasta_file = list(input$uploadOTUfasta$datapath, input$uploadZOTUfasta$datapath),
    species_assign = list(FALSE, FALSE)
  )


  params_import <- case_when(
    input$analysis_type == "OTU" ~ lapply(params_import, `[`, 1),
    input$analysis_type == "ZOTU" ~ lapply(params_import, `[`, 2),
    input$analysis_type == "both" ~ params_import
  )


  # check if color_set, labels and label_set have the same length
  color_set <- input$color_set
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- input$label_set

  check_color_set(color_set, labels, label_set)

  color_set <- setNames(input$color_set, input$label_set)
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- setNames(labels, input$label_set)

  withProgress(message = "Rendering R Markdown...", value = 0, {


    purrr::pwalk(params_import,
                 ~rmarkdown::render(
                   input = system.file("rmd", "01_phyloseq_import_v1.10.Rmd", package = "golemMB"),
                   output_format = "pdf_document",
                   output_dir = output_dir1,
                   output_file = ..1,
                   params = list(
                     author = author,
                     TU_type = ..2,
                     TU_file = ..3,
                     TU_tree_file = ..4,
                     TU_fasta_file = ..5,
                     meta_file = input$uploadmeta$datapath,
                     species_assign = ..6,
                     filter_var = filter_var,
                     samples_excluded = samples_excluded,
                     samples_outlier = samples_outlier,
                     min_sample_counts = input$min_sample_counts,
                     min_tu_counts = 1,
                     count_transform = "log10p", #input$count_transform,
                     var_test = input$var_test,
                     var_name = input$var_name,
                     point_label = input$point_label,
                     shape_var = input$shape_var,
                     color_set = color_set,
                     label_set = label_set)
                 ))

    setProgress(1, message = "Rendering complete!")
    Sys.sleep(3)})

  # download phyloseq object
  output$download_import <- download_function("outputs_import", output_dir1)

})

### ~~~end run import~~~ ###

### ~~~calculate alpha diversity~~~ ###
observeEvent(input$calc_alpha, {

  ### test directory ###
  output_dir2 <- file.path(tempdir(), "alpha_diversity")


  # error handling for missing inputs/files # FIX THIS FOR ZOTUS TOO

  # create list of all required inputs
  missing_inputs_ac <- list("Metadata file" = input$uploadmeta$datapath,
                            "alpha diversity color" = input$var_color,
                            "analysis variable" = input$var_test
  )




  if (input$analysis_type == "OTU") {
    missing_inputs_ac <- c(missing_inputs_ac, list(
      "OTU Import file" = input$uploadOTUphyseq$datapath
    ))

    # filter list for missing inputs
    missing_inputs_ac <- lapply(missing_inputs_ac, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })

  }

  if (input$analysis_type == "ZOTU") {
    missing_inputs_ac <- c(missing_inputs_ac, list(
      "ZOTU Import file" = input$uploadZOTUphyseq$datapath
    ))
    # filter list for missing inputs
    missing_inputs_ac <- lapply(missing_inputs_ac, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })

  }


  if (input$analysis_type == "both") {
    missing_inputs_ac <- c(missing_inputs_ac, list(
      "OTU Import file" = input$uploadOTUphyseq$datapath,
      "ZOTU Import file" = input$uploadZOTUphyseq$datapath))



    # filter list for missing inputs
    missing_inputs_ac <- lapply(missing_inputs_ac, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })

  }

  # pull names of missing inputs to paste into error message
  missing <- names(missing_inputs_ac)[sapply(missing_inputs_ac, is.null)]

  #catch and print missing inputs
  if (length(missing) > 0) {
    spsComps::shinyCatch(

      stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")

    return()
  }

  # set analysis variable for pwalk
    analysis_type <- case_when(
    input$analysis_type == "OTU" ~ "OTU",
    input$analysis_type == "ZOTU" ~ "ZOTU",
    input$analysis_type == "both" ~ c("OTU", "ZOTU")
  )

#######################################################################################################
#!!!!!!!    # analysis_type <- ifelse(input$analysis_type == "both", c("OTU, "ZOTU"), .) !!!!!!!!!!!!!#
#######################################################################################################

  import_Rhea  <- FALSE
  alpha_file   <- "/Rhea/2.Alpha-Diversity/alpha.tab"
  meta_file <- gsub("\\\\", "/", input$uploadmeta$datapath)
  # render Rmarkdown document

  # progress bar



  withProgress(message = "Rendering R Markdown...", value = 0, {

    purrr::walk(analysis_type,
                ~rmarkdown::render(
                  input = system.file("rmd", "02_calc_alpha_diversity_param_v1.20.Rmd", package = "golemMB"),
                  output_format = "pdf_document",
                  output_dir = output_dir2,
                  output_file = paste0("02_calc_alpha_diversity_", .x, ".pdf"),
                  intermediates_dir = output_dir2,
                  params = list(
                    author = author,
                    workdir = workdir,
                    TU_type = .x,
                    var_test = input$var_test,
                    var_color = input$var_color
                  )))
    setProgress(1, message = "Rendering complete!")
    Sys.sleep(3)

    # print(list.files(output_dir))

    ### test ###
    ### download function ###
    output$download_alpha1 <- download_function("outputs_alphadiversity", output_dir2)


  })
})

### ~~~end calculate alpha diversity~~~ ###


### ~~~count transform data~~~ ###
observeEvent(input$count_transform, {

  # create temporary directory to store the output
  output_dir3 <- file.path(tempdir(), "count_transforming")

  missing_inputs_ct <- list("Metadata file" = input$uploadmeta$datapath,
                            "analysis type" = input$analysis_type,
                            "analysis variable name" = input$var_name,
                            "analysis variable" = input$var_test)


  # for analysis type = OTU: check all missing inputs (and files)
  if (input$analysis_type == "OTU") {
    missing_inputs_ct <- c(missing_inputs_ct, list(
      "OTU Import file" = input$uploadOTUphyseq$datapath
    ))
    # filter list for missing inputs
    missing_inputs_ct <- lapply(missing_inputs_ct, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })

  }

  # for analysis type = ZOTU: check all missing inputs (and files)
  if (input$analysis_type == "ZOTU") {
    missing_inputs_ct <- c(missing_inputs_ct, list(
      "ZOTU Import file" = input$uploadZOTUphyseq$datapath
    ))

    # filter list for missing inputs
    missing_inputs_ct <- lapply(missing_inputs_ct, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })
  }

  # for analysis type = both: check all missing inputs (and files)
  if (input$analysis_type == "both") {
    missing_inputs_ct <- c(missing_inputs_ct, list(
      "OTU Import file" = input$uploadOTUphyseq$datapath,
      "ZOTU Import file" = input$uploadZOTUphyseq$datapath))

    # filter list for missing inputs
    missing_inputs_ct <- lapply(missing_inputs_ct, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })
  }


  # pull names of missing inputs to paste into error message
  missing <- names(missing_inputs_ct)[sapply(missing_inputs_ct, is.null)]

  #catch and print missing inputs
  if (length(missing) > 0) {
    spsComps::shinyCatch(

      stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")
    return()
  }



  analysis_type <- case_when(
    input$analysis_type == "OTU" ~ "OTU",
    input$analysis_type == "ZOTU" ~ "ZOTU",
    input$analysis_type == "both" ~ c("OTU", "ZOTU")
  )


  withProgress(message = "Rendering R Markdown...", value = 0, {

    purrr::walk(analysis_type,
                ~rmarkdown::render(
                  input = system.file("rmd", "03_count_transform_param_v1.21.Rmd", package = "golemMB"),
                  output_format = "pdf_document",
                  output_dir = output_dir3,
                  output_file = paste0("03_count_transform_", .x,".pdf"),
                  intermediates_dir = output_dir3,
                  params = list(
                    author = author,
                    workdir = workdir,
                    TU_type = .x,
                    var_test = input$var_test,
                    var_name = input$var_name
                  )))
    setProgress(1, message = "Rendering complete!")
    Sys.sleep(3)
  })

  output$download_ct <- download_function("outputs_count_transformation", output_dir3)

})

### ~~~end count transform data~~~ ###

### ~~~calculate distance matrices~~~ ###
observeEvent(input$dist_ma, {

  # create temporary directory to store output files
  output_dir4 <- file.path(tempdir(), "dist_matrices") %>%
    gsub("\\\\", "/", .)



  missing_inputs_dm <- list("Metadata file" = input$uploadmeta$datapath,
                            "analysis type" = input$analysis_type,
                            "analysis variable name" = input$var_name,
                            "analysis variable" = input$var_test,
                            "color set" = input$color_set,
                            "variable order" = input$label_set,
                            "variable labels" = input$labels)




  # for analysis type = OTU: check all missing inputs (and files)
  if (input$analysis_type == "OTU") {
    missing_inputs_dm <- c(missing_inputs_dm, list(
      "OTU count transformation file" = input$uploadOTUtrans$datapath,
      "OTU count transformation file (Genus)" = input$uploadOTUtransgenus$datapath
    ))
    # filter list for missing inputs
    missing_inputs_dm <- lapply(missing_inputs_dm, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })

  }

  # for analysis type = ZOTU: check all missing inputs (and files)
  if (input$analysis_type == "ZOTU") {
    missing_inputs_dm <- c(missing_inputs_dm, list(
      "ZOTU count transformation file" = input$uploadZOTUtrans$datapath,
      "ZOTU count transformation file (Genus)" = input$uploadZOTUtransgenus$datapath
    ))

    # filter list for missing inputs
    missing_inputs_dm <- lapply(missing_inputs_dm, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })
  }

  # for analysis type = both: check all missing inputs (and files)
  if (input$analysis_type == "both") {
    missing_inputs_dm <- c(missing_inputs_dm, list(
      "OTU count transformation file" = input$uploadOTUtrans$datapath,
      "OTU count transformation file (Genus)" = input$uploadOTUtransgenus$datapath,
      "ZOTU count transformation file" = input$uploadZOTUtrans$datapath,
      "ZOTU count transformation file (Genus)" = input$uploadZOTUtransgenus$datapath
      ))

    # filter list for missing inputs
    missing_inputs_dm <- lapply(missing_inputs_dm, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })
  }


  # pull names of missing inputs to paste into error message
  missing <- names(missing_inputs_dm)[sapply(missing_inputs_dm, is.null)]

  #catch and print missing inputs
  if (length(missing) > 0) {
    spsComps::shinyCatch(

      stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")
    return()
  }


  # set analysis type variable for purrr::walk

  analysis_type <- case_when(
    input$analysis_type == "OTU" ~ "OTU",
    input$analysis_type == "ZOTU" ~ "ZOTU",
    input$analysis_type == "both" ~ c("OTU", "ZOTU")
  )

  # check if color_set, labels and label_set have the same length
  color_set <- input$color_set
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- input$label_set

  check_color_set(color_set, labels, label_set)

  color_set <- setNames(input$color_set, input$label_set)
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- setNames(labels, input$label_set)

  # progress bar
  withProgress(message = "Rendering R Markdown...", value = 0, {

    purrr::walk(analysis_type,
                ~rmarkdown::render(
                  input = system.file("rmd", "04_dist_ord_calculation_param_v1.16.Rmd", package = "golemMB"),
                  output_format = "pdf_document",
                  output_dir = output_dir4,
                  output_file =  paste0("04_dist_ord_calculation_", .x, ".pdf"),
                  intermediates_dir = output_dir4,
                  params = list(
                    author = author,
                    workdir = workdir,
                    TU_type = .x,
                    var_test = input$var_test,
                    var_name = input$var_name,
                    color_set = color_set,
                    label_set = label_set
                  )))
    setProgress(1, message = "Rendering complete!")
    Sys.sleep(3)
  })

  # download output
  output$download_dm <- download_function("outputs_distance_matrices", output_dir4)

})

### ~~~end calculate distance matrices~~~ ###


### ~~~plot alpha diversity~~~ ###
observeEvent(input$plot_alpha, {


  ### test directory ###
  output_dir5 <- file.path(tempdir(), "alpha_diversity_plots") %>%
    gsub("\\\\", "/", .)


  missing_inputs_pa <- list("Metadata file" = input$uploadmeta$datapath,
                            "analysis type" = input$analysis_type,
                            "analysis variable name" = input$var_name,
                            "analysis variable formula" = input$var_formula,
                            "analysis variable" = input$var_test,
                            "color set" = input$color_set,
                            "variable order" = input$label_set,
                            "variable labels" = input$labels
                            #"OTU alpha calculation file" = alpha_otu <- paste0(workdir, "/data/alpha_list_otu.rds")
  )

  # for analysis type = OTU: check all missing inputs (and files)
  if (input$analysis_type == "OTU") {
    missing_inputs_pa <- c(missing_inputs_pa, list(
      "OTU alpha diversity file" = input$uploadOTUalpha$datapath)
    )
    # filter list for missing inputs
    missing_inputs_pa <- lapply(missing_inputs_pa, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })

  }

  # for analysis type = ZOTU: check all missing inputs (and files)
  if (input$analysis_type == "ZOTU") {
    missing_inputs_pa <- c(missing_inputs_pa, list(
      "ZOTU alpha diversity file" = input$uploadZOTUalpha$datapath
    ))

    # filter list for missing inputs
    missing_inputs_pa <- lapply(missing_inputs_pa, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })
  }

  # for analysis type = both: check all missing inputs (and files)
  if (input$analysis_type == "both") {
    missing_inputs_pa <- c(missing_inputs_pa, list(
      "OTU count transformation file" = input$uploadOTUalpha$datapath,
      "ZOTU count transformation file" = input$uploadZOTUalpha$datapath))

    # filter list for missing inputs
    missing_inputs_pa <- lapply(missing_inputs_pa, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })
  }


  # pull names of missing inputs to paste into error message
  missing <- names(missing_inputs_pa)[sapply(missing_inputs_pa, is.null)]

  #catch and print missing inputs
  if (length(missing) > 0) {
    spsComps::shinyCatch(

      stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")
    return()
  }



  comp_lists <- lapply(count$comp_ids, function(id) input[[id]])

  #test
  color_set <- input$color_set
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- input$label_set

check_color_set(color_set, labels, label_set)


  color_set <- setNames(input$color_set, input$label_set)
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- setNames(labels, input$label_set)
  var_name <- input$var_name
  var_formula <- input$var_formula
  var_test <- input$var_test
  comp_list <- list(append(list(input$comp_list), comp_lists))


  grouped_pt = FALSE
  color_param = NULL
  data_param = NULL

  var_list <- list(var_test)
  var_formula_list <- list(var_formula)
  comp_list_list <- comp_list
  color_set_list <- list(color_set)
  label_set_list <- list(label_set)
  var_name_list <- list(var_name)

  # this is for grouped plots, will add at a later point#-------------------------
  label_set_b <- list(NULL, NULL)
  grouped_list <- list(grouped_pt)
  color_param_list <- list(color_param)
  data_param_list <- list(data_param)
  #-------------------------------------------------------------------------------
  # format parameter list for pwalk
  params_alpha_bugbase <- list(
    subfolder = lapply(var_list,
                       function(x) {paste(format(Sys.time(), '%y%m%d'), x, sep = "_")}),
    var_test = var_test,
    var_formula = var_formula_list,
    var_name = var_name_list,
    comp_list = comp_list_list,
    color_set = color_set_list,
    label_set = label_set_list,
    label_set_b = list(NULL),
    grouped = grouped_list,
    color_param = color_param_list,
    data_param = data_param_list
  )


  params_alpha <- Map(c, params_alpha_bugbase, params_alpha_bugbase)
  params_alpha <- c(
    TU_type = list(rep(list("OTU", "ZOTU"),
                       each = length(var_list))),
    params_alpha
  )

  params_alpha <- Map(c, params_alpha_bugbase, params_alpha_bugbase)
  params_alpha <- c(
    TU_type = list(rep(list("OTU", "ZOTU"),
                       each = length(var_list))),
    params_alpha
  )

  params_alpha <- case_when(
    input$analysis_type == "OTU" ~ lapply(params_alpha, `[`, 1),
    input$analysis_type == "ZOTU" ~ lapply(params_alpha, `[`, 2),
    input$analysis_type == "both" ~ params_alpha
  )


  withProgress(
    message = "Rendering R Markdown...", value = 0, {

      # render Rmarkdown document
      purrr::pwalk(params_alpha,
                   ~rmarkdown::render(
                     input = system.file("rmd", "06_plot_alpha_param_v1.52.Rmd", package = "golemMB"),
                     output_format = "pdf_document",
                     output_dir = output_dir5,
                     output_file = paste("06_plot_alpha", ..3, ..1, "pdf", sep = "."),
                     intermediates_dir = output_dir5,
                     params = list(
                       author = author,
                       workdir = workdir,
                       subfolder = ..2,
                       TU_type = ..1,
                       var_test = ..3,
                       var_formula = ..4,
                       var_name = ..5,
                       comp_list = comp_list_list, #NULL for no stats, "all" for all comparisons, ..6 for specific comparisons
                       color_set = ..7, #NULL for bw plots ..7 for colors
                       label_set = ..8,
                       label_set_b = ..9,
                       grouped = ..10,
                       dodge_value = 0.9,
                       color_param = ..11,
                       data_param = ..12,
                       alpha_vars = c("Observed", "Chao1", "ACE", "Shannon",
                                      "Simpson", "InvSimpson"),
                       min_repl = 3,
                       rm_ns = TRUE,
                       padj_method = "fdr",
                       cutoff_pval = 0.05,
                       rel_step = 6
                     )))

      setProgress(1, message = "Rendering complete!")
      Sys.sleep(3)
    })

  # download rmd output
  output$download_alpha2 <- download_function("outputs_alphadiversity_plots", output_dir5)

})

### ~~~end plot alpha diversity~~~ ###


### ~~~plot BugBase results~~~ ###
observeEvent(input$plot_bb, {

  # create temporary directory to store output files
  output_dir6 <- file.path(tempdir(), "bugbase") %>%
    gsub("\\\\", "/", .)

  missing_inputs_bb <- list("Metadata file" = input$uploadmeta$datapath,
                            "OTU import file" = input$uploadOTUphyseq$datapath,
                            "analysis type" = input$analysis_type,
                            "analysis variable name" = input$var_name,
                            "analysis variable formula" = input$var_formula,
                            "analysis variable" = input$var_test,
                            "color set" = input$color_set,
                            "variable order" = input$label_set,
                            "variable labels" = input$labels,
                            "Sample outliers" = input$samples_outlier,
                            "BugBase prediction file" = input$uploadbugbasepredicts$datapath,
                            "BugBase vsearch stats" = input$uploadbugbasevsearch$datapath)


  # filter list for missing inputs
  missing_inputs_bb <- lapply(missing_inputs_bb, function(x) {
    if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
      NULL
    } else {
      x
    }
  })



  # pull names of missing inputs to paste into error message
  missing <- names(missing_inputs_bb)[sapply(missing_inputs_bb, is.null)]

  #catch and print missing inputs
  if (length(missing) > 0) {
    spsComps::shinyCatch(

      stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")
    return()
  }


  if (input$samples_outlier == "none") {
    samples_outlier <- c()
  } else {
    samples_outlier <- input$samples_outlier
  }


  meta_file <- input$uploadmeta$datapath
  comp_lists <- lapply(count$comp_ids, function(id) input[[id]])

  # check if color_set, labels and label_set have the same length
  color_set <- input$color_set
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- input$label_set

  check_color_set(color_set, labels, label_set)

  color_set <- setNames(input$color_set, input$label_set)
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- setNames(labels, input$label_set)
  var_name <- input$var_name
  var_formula <- input$var_formula
  var_test <- input$var_test
  comp_list <- list(append(list(input$comp_list), comp_lists))


  filter_var = "infection"
  grouped_pt = FALSE
  color_param = NULL
  data_param = NULL

  var_list <- list(var_test)
  var_formula_list <- list(var_formula)
  comp_list_list <- comp_list
  color_set_list <- list(color_set)
  label_set_list <- list(label_set)
  var_name_list <- list(var_name)

  # this is for grouped plots, will add at a later point#-------------------------
  label_set_b <- list(NULL, NULL)
  grouped_list <- list(grouped_pt)
  color_param_list <- list(color_param)
  data_param_list <- list(data_param)
  #-------------------------------------------------------------------------------
  # format parameter list for pwalk
  params_alpha_bugbase <- list(
    subfolder = lapply(var_list,
                       function(x) {paste(format(Sys.time(), '%y%m%d'), x, sep = "_")}),
    var_test = var_test,
    var_formula = var_formula_list,
    var_name = var_name_list,
    comp_list = comp_list_list,
    color_set = color_set_list,
    label_set = label_set_list,
    label_set_b = list(NULL),
    grouped = grouped_list,
    color_param = color_param_list,
    data_param = data_param_list
  )


  params_alpha_bugbase <- case_when(
    input$analysis_type == "OTU" ~ lapply(params_alpha_bugbase, `[`, 1),
    input$analysis_type == "ZOTU" ~ lapply(params_alpha_bugbase, `[`, 2),
    input$analysis_type == "both" ~ params_alpha_bugbase
  )


  #progress bar
  withProgress(message = "Rendering R Markdown...", value = 0, {

    purrr::pwalk(params_alpha_bugbase,
                 ~rmarkdown::render(
                   input = system.file("rmd", "06a_plot_BugBase_param_v1.53.Rmd", package = "golemMB"),
                   output_format = "pdf_document",
                   output_dir = output_dir6,
                   output_file = paste("06a_plot_BugBase", ..2, "pdf", sep = "."),
                   intermediates_dir = output_dir6,
                   params = list(
                     workdir = workdir,
                     subfolder = ..1,
                     var_test = ..2,
                     var_formula = ..3,
                     var_name = ..4,
                     filter_var = filter_var, #sample
                     samples_excluded = samples_excluded,
                     samples_outlier = samples_outlier,
                     comp_list = ..5, #NULL for no stats, "all" for all comparisons, ..5 for specific comparisons
                     color_set = ..6, #NULL for bw plots ..6 for colors
                     label_set = ..7,
                     label_set_b = ..8,
                     grouped = ..9,
                     dodge_value = 0.9,
                     color_param = ..10,
                     data_param = ..11,
                     min_repl = 3,
                     rm_ns = TRUE,
                     padj_method = "fdr",
                     cutoff_pval = 0.05,
                     rel_step = 6
                   )))

    incProgress(1, message = "Rendering complete!")
    Sys.sleep(3)
  })

  print(output_dir6)

  # download output files
  output$download_bugbase <- download_function("outputs_BugBase", output_dir6)

})

### ~~~end plot BugBase results~~~ ###


### ~~~plot beta diversity~~~ ###
observeEvent(input$plot_beta, {

  # create temporary directory to store output files
  output_dir7 <- file.path(tempdir(), "beta_diversity") %>%
    gsub("\\\\", "/", .)


  missing_inputs_bd <- list("Metadata file" = input$uploadmeta$datapath,
                            "analysis type" = input$analysis_type,
                            "analysis variable name" = input$var_name,
                            "analysis variable formula" = input$var_formula,
                            "analysis variable" = input$var_test,
                            "color set" = input$color_set,
                            "variable order" = input$label_set,
                            "variable labels" = input$labels,
                            "variable shape" = input$shape_var)


  if (input$analysis_type == "OTU") {
    missing_inputs_bd <- c(missing_inputs_bd, list(
      "OTU distance matrices files" = input$uploadOTUdists$datapath,
      "OTU Genus distance matrices files" = input$uploadOTUdistsgenus$datapath,
      "OTU count transformation files" = input$uploadOTUtrans$datapath,
      "OTU Genus count transformation files" = input$uploadOTUtransgenus$datapath,
      "OTU ordination files" = input$uploadOTUord$datapath,
      "OTU Genus ordination files" = input$uploadOTUordgenus$datapath
    ))
    # filter list for missing inputs
    missing_inputs_bd <- lapply(missing_inputs_bd, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })

  }

  # for analysis type = ZOTU: check all missing inputs (and files)
  if (input$analysis_type == "ZOTU") {
    missing_inputs_bd <- c(missing_inputs_bd, list(
      "ZOTU distance matrices files" = input$uploadZOTUdists$datapath,
      "ZOTU Genus distance matrices files" = input$uploadZOTUdistsgenus$datapath,
      "ZOTU count transformation files" = input$uploadZOTUtrans$datapath,
      "ZOTU Genus count transformation files" = input$uploadZOTUtransgenus$datapath,
      "ZOTU ordination files" = input$uploadZOTUord$datapath,
      "ZOTU Genus ordination files" = input$uploadZOTUordgenus$datapath
    ))

    # filter list for missing inputs
    missing_inputs_bd <- lapply(missing_inputs_bd, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })
  }

  # for analysis type = both: check all missing inputs (and files)
  if (input$analysis_type == "both") {
    missing_inputs_bd <- c(missing_inputs_bd, list(
      "OTU distance matrices files" = input$uploadOTUdists$datapath,
      "OTU Genus distance matrices files" = input$uploadOTUdistsgenus$datapath,
      "OTU count transformation files" = input$uploadOTUtrans$datapath,
      "OTU Genus count transformation files" = input$uploadOTUtransgenus$datapath,
      "OTU ordination files" = input$uploadOTUord$datapath,
      "OTU Genus ordination files" = input$uploadOTUordgenus$datapath,
      "ZOTU distance matrices files" = input$uploadZOTUdists$datapath,
      "ZOTU Genus distance matrices files" = input$uploadZOTUdistsgenus$datapath,
      "ZOTU count transformation files" = input$uploadZOTUtrans$datapath,
      "ZOTU Genus count transformation files" = input$uploadZOTUtransgenus$datapath,
      "ZOTU ordination files" = input$uploadZOTUord$datapath,
      "ZOTU Genus ordination files" = input$uploadZOTUordgenus$datapath
    ))

    # filter list for missing inputs
    missing_inputs_bd <- lapply(missing_inputs_bd, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })
  }

  # pull names of missing inputs to paste into error message
  missing <- names(missing_inputs_bd)[sapply(missing_inputs_bd, is.null)]

  #catch and print missing inputs
  if (length(missing) > 0) {
    spsComps::shinyCatch(

      stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")
    return()
  }


  # check if color_set, labels and label_set have the same length
  color_set <- input$color_set
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- input$label_set

  check_color_set(color_set, labels, label_set)

  var_test <- input$var_test
  var_name <- input$var_name
  var_formula <- input$var_formula
  shape_var <- input$shape_var
  color_set <- setNames(input$color_set, input$label_set)
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- setNames(labels, input$label_set)



  if (input$shape_var == "none") {
    shape_var <- NULL
  }

  var_list <- list(var_test)
  var_formula_list <- list(var_formula)
  var_name_list <- list(var_name)
  shape_var_list <- list(shape_var)
  color_set_list <- list(color_set)
  label_set_list <- list(label_set)

  trans_methods <- c("CLR", "ra")
  dist_methods <- c("guni", "bray")
  ord_methods <- c("MDS", "NMDS", "tSNE", "UMAP")

  # format parameter list for pwalk: beta diversity
  params_beta <- list(
    subfolder = lapply(
      var_list, function(x) {paste(format(Sys.time(), '%y%m%d'), x, sep = "_")}),
    var_test = var_list,
    var_comp_test = var_formula_list,
    var_name = var_name_list,
    shape_var = shape_var_list,
    color_set = color_set_list,
    label_set = label_set_list
  )



  # multiply parameters for OTU/ZOTU/Genus by four
  params_beta <- Map(c, params_beta, params_beta, params_beta, params_beta)
  params_beta <- c(
    TU_type = list(
      rep(list("ZOTU", "OTU", "ZOTU", "OTU"), each = length(var_list))),
    analysis_level = list(
      rep(list("TU", "Genus"), each = length(var_list)*2)),
    params_beta)

  params_beta <- case_when(
    input$analysis_type == "OTU" ~ lapply(params_beta, `[`, c(2,4)),
    input$analysis_type == "ZOTU" ~ lapply(params_beta, `[`, c(1,3)),
    input$analysis_type == "both" ~ params_beta
  )


  # progress bar
  withProgress(message = "Rendering R Markdown...", value = 0, {

    purrr::pwalk(params_beta,
                 ~rmarkdown::render(
                   input = system.file("rmd", "07_beta_diversity_param_v1.15.Rmd", package = "golemMB"),
                   output_format = "pdf_document",
                   output_dir = output_dir7,
                   output_file = paste("07_beta_diversity", ..4, ..1, ..2,
                                       #ifelse(..2 == "Genus", "genus", NULL),
                                       "pdf", sep = "."),
                   intermediates_dir = output_dir7,
                   params = list(
                     author = author,
                     workdir = workdir,
                     subfolder = ..3,
                     TU_type = ..1,
                     analysis_level = ..2,
                     var_test = ..4,
                     var_comp_test = ..5,
                     var_name = ..6,
                     trans_methods = trans_methods,
                     dist_methods = dist_methods,
                     ord_methods = ord_methods,
                     shape_var = ..7,
                     color_set = ..8,
                     label_set = ..9
                   )));
    setProgress(1, message = "Rendering complete!")
    Sys.sleep(3)
  })

  # download output files
  output$download_betadiversity <- download_function("outputs_betadiversity", output_dir7)

})

### ~~~end plot beta diversity~~~ ###


### ~~~plot differential abundance~~~ ###
observeEvent(input$plot_DA, {

  # create temporary directory to store output files
  output_dir8 <- file.path(tempdir(), "diff_abund") %>%
    gsub("\\\\", "/", .)

  missing_inputs_da <- list("Metadata file" = input$uploadmeta$datapath,
                            "analysis type" = input$analysis_type,
                            "analysis variable" = input$var_test,
                            "color set" = input$color_set,
                            "variable order" = input$label_set,
                            "variable labels" = input$labels,
                            "comparison list" = input$comp_list)


  if (input$analysis_type == "OTU") {
    missing_inputs_da <- c(missing_inputs_da, list(
      "OTU count transformation files" = input$uploadOTUtrans$datapath
    ))

    missing_inputs_da <- lapply(missing_inputs_da, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })
  }


  if (input$analysis_type == "ZOTU") {
    missing_inputs_bd <- c(missing_inputs_da, list(
      "ZOTU count transformation files" = input$uploadZOTUtrans$datapath
    ))

    missing_inputs_da <- lapply(missing_inputs_da, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })

  }


  if (input$analysis_type == "both") {
    missing_inputs_da <- c(missing_inputs_da, list(
      "OTU count transformation files" = input$uploadOTUtrans$datapath,
      "ZOTU count transformation files" = input$uploadZOTUtrans$datapath
    ))

    missing_inputs_da <- lapply(missing_inputs_da, function(x) {
      if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
        NULL
      } else {
        x
      }
    })

  }


  # pull names of missing inputs to paste into error message
  missing <- names(missing_inputs_da)[sapply(missing_inputs_da, is.null)]

  #catch and print missing inputs
  if (length(missing) > 0) {
    spsComps::shinyCatch(

      stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")
    return()
  }


  trans_method <- "CLR"

  # which taxonomic level to use for metacoder plots?
  metacoder_tax <- "Genus"

  comp_lists <- lapply(count$comp_ids, function(id) input[[id]])

  # check if color_set, labels and label_set have the same length
  color_set <- input$color_set
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- input$label_set

  check_color_set(color_set, labels, label_set)

  var_test <- input$var_test
  comp_list <- c(list(input$comp_list), comp_lists)
  color_set <- setNames(input$color_set, input$label_set)
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- setNames(labels, input$label_set)
  var_name <- input$var_name
  var_list <- list(var_test)
  comp_list_list <- list(comp_list)


  # format parameter list for pwalk
  params_da <- list(
    analysis_level = rep(list("all"), each = length(var_list)),
    subfolder = lapply(var_list,
                       function(x) {paste(format(Sys.time(), '%y%m%d'), x, sep = "_")}),
    var_test = var_list,
    comp_list = comp_list_list
  )
  # duplicate parameters for OTU/ZOTU
  params_da <- Map(c, params_da, params_da)
  params_da <- c(
    TU_type = list(rep(list("OTU", "ZOTU"),
                       each = length(var_list))),
    params_da
  )


  params_da <- case_when(
    input$analysis_type == "OTU" ~ lapply(params_da, `[`, 1),
    input$analysis_type == "ZOTU" ~ lapply(params_da, `[`, 2),
    input$analysis_type == "both" ~ params_da
  )


  # progress bar
  withProgress(message = "Rendering R Markdown...", value = 0, {

    # loop through all variables using pwalk
    purrr::pwalk(params_da,
                 ~rmarkdown::render(
                   input = system-file("rmd", "08_differential_abundance_param_v3.12.Rmd", package = "golemMB"),
                   output_format = "pdf_document",
                   output_dir = output_dir8,
                   output_file = paste("08_DA_analysis", ..4, ..1, ..2, metacoder_tax, "pdf", sep = "."),
                   intermediates_dir = output_dir8,
                   params = list(
                     author = author,
                     workdir = workdir,
                     subfolder = ..3,
                     TU_type = ..1,
                     useMC = 4,
                     analysis_level = ..2,
                     var_test = ..4,
                     trans_method = trans_method,
                     comp_list = ..5,
                     seed = 42,
                     min_count = 10, # counts in at least 'min_count_prev'% samples
                     min_count_prev = 10,
                     min_rel_abund = 0.1, # relative abundance in at least 'min_rel_abund_prev'% samples
                     min_rel_abund_prev = 10,
                     cutoff_pval = 0.05,
                     padj_method = "fdr",
                     cutoff_fcda = 2,
                     cutoff_lda = 2,
                     cutoff_W = 0.75,
                     cutoff_aldex = 1,
                     metacoder_tax = metacoder_tax,
                     volcano_label = 20
                   )));
    incProgress(1, message = "Rendering complete!")
    Sys.sleep(3)
  })

  output$download_diffabund <- download_function("outputs_differential_abundance", output_dir8)

})


### ~~~end plot differential abundance~~~ ###



### ~~~plot functional abundance~~~ ###
observeEvent(input$plot_FA, {



  # create temporary directory to store output files
  output_dir9 <- file.path(tempdir(), "func_abund") %>%
    gsub("\\\\", "/", .)

  missing_inputs_da <- list("Metadata file" = input$uploadmeta$datapath,
                            "analysis variable" = input$var_test,
                            "analysis variable name" = input$var_name,
                            "color set" = input$color_set,
                            "variable order" = input$label_set,
                            "variable labels" = input$labels,
                            "comparison list" = input$comp_list,
                            "ZOTU Import files" = input$uploadZOTUphyseq$datapath,
                            "PICRUSt2 out_stratified folder" = input$uploadpicrustfolder$datapath,
                            "PICRUSt2 stats file" = input$uploadpicruststats$datapath

  )


  # filter list for missing inputs
  missing_inputs_da <- lapply(missing_inputs_da, function(x) {
    if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
      NULL
    } else {
      x
    }
  })

  # pull names of missing inputs to paste into error message
  missing <- names(missing_inputs_da)[sapply(missing_inputs_da, is.null)]

  #catch and print missing inputs
  if (length(missing) > 0) {
    spsComps::shinyCatch(

      stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")
    return()
  }


  # unzip picrust folder
  zip_dir <- gsub("\\\\", "/", tempfile())
  dir.create(zip_dir)
  unzip(input$uploadpicrustfolder$datapath, exdir = zip_dir)

  # print directory
  check_files <- list.files(output_dir9, recursive = TRUE, full.names = TRUE)

  browser()

  comp_lists <- lapply(count$comp_ids, function(id) input[[id]])

  # check if color_set, labels and label_set have the same length
  color_set <- input$color_set
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- input$label_set

  check_color_set(color_set, labels, label_set)

  var_test <- input$var_test
  comp_list <- list(c(input$comp_list, comp_lists))
  color_set <- setNames(input$color_set, input$label_set)
  labels <- unlist(strsplit(input$labels, ","))
  label_set <- setNames(labels, input$label_set)
  var_name <- input$var_name




  var_list <- list(var_test)
  comp_list_list <- list(comp_list)

  params_fa <- list(
    subfolder = lapply(var_list,
                       function(x) {paste(format(Sys.time(), '%y%m%d'), x, sep = "_")}),
    var_test = var_list,
    comp_list = comp_list_list
  )

  # progress bar
  withProgress(message = "Rendering R Markdown...", value = 0, {

    # loop through all variables using pwalk
    purrr::pwalk(params_fa,
                 ~rmarkdown::render(
                   input = system.file("rmd", "09_functional_abundance_param_v1.50.Rmd", package = "golemMB"),
                   output_format = "pdf_document",
                   output_dir = output_dir9,
                   output_file = paste("09_FA_analysis", ..2, "pdf", sep = "."),
                   intermediates_dir = output_dir9,
                   params = list(
                     author = author,
                     workdir = workdir,
                     subfolder = ..1,
                     useMC = 4,
                     var_test = ..2,
                     comp_list = ..3,
                     seed = 42,
                     cutoff_pval = 0.05,
                     cutoff_aldex = 1
                   )));
    incProgress(1, message = "Rendering complete!")
    Sys.sleep(3)
  })

  # download output
  output$download_funcabund <- download_function("outputs_functional_abundance", output_dir9)


})

### ~~~end plot functional abundance~~~ ###


# save logs

# helper function for NULL (no input given) inputs
# safe_file_name <- function(input_element) {
#   if (is.null(input_element)) {
#     return("")
#   } else {
#     return(input_element$name)
#   }
# }
#


# observeEvent(input$save_inputs, {
#
#   input_values <- c(input$analysis_type, input$var_test, input$var_formula,
#                          input$min_sample_counts,
#                          input$var_name, input$shape_var,
#                          label_set <- as.vector(paste0(input$label_set, collapse = ",")),
#                          labels <- as.vector(paste0(input$labels, collapse = ",")),
#                          color_set <- as.vector(paste0(input$color_set, collapse = ",")),
#                          input$point_label, input$var_color,
#                          comp_list <- as.vector(paste0(input$comp_list, collapse = ",")),
#                          safe_file_name(input$uploadmeta),
#                          safe_file_name(input$uploadOTUtab),
#                          safe_file_name(input$uploadOTUtre),
#                          safe_file_name(input$uploadOTUfasta),
#                          safe_file_name(input$uploadZOTUtab),
#                          safe_file_name(input$uploadZOTUtre),
#                          safe_file_name(input$uploadZOTUfasta))
#
#   input_labels <- c("Analysis Type", "Analysis variable", "Analysis variable formula",
#                       "Min. sample counts", "Analysis variable name", "Shape by x", "Label order",
#                       "labels", "color set", "Datapoint label", "AD colors", "Comparison list", "meta file",
#                     "OTU table", "OTU tree file", "OTU fasta", "ZOTU table", "ZOTU tree file",
#                     "ZOTU fasta")
#
#   # input_values <- list()
#   logfile <- data.frame(input_labels, input_values)
#
#   # Save the inputs to a CSV file
#   # date
#   current_date <- format(Sys.Date(), '%y%m%d')
#   write.csv(logfile, file = paste0(current_date, "_user_inputs.csv"))
# })
# if (!file.exists("user_inputs.csv")) {
#   print("doesnt work!!")
# }



### ~~~calculate TU correlations~~~ ###
observeEvent(input$calc_corr, {

  # create temporary directory to store output files
  output_dir10 <- file.path(tempdir(), "correlation")

  # set a list of all necessary files and inputs
  missing_inputs_cor <- list("Metadata file" = input$uploadmeta$datapath
  )


if (input$analysis_type == "OTU") {
  missing_inputs_cor <- c(missing_inputs_cor, list(
    "OTU alpha diversity file" = input$uploadOTUphyseq$datapath
  ))
  # filter list for missing inputs
  missing_inputs_cor <- lapply(missing_inputs_cor, function(x) {
    if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
      NULL
    } else {
      x
    }
  })

}

#for analysis type = ZOTU: check all missing inputs (and files)
if (input$analysis_type == "ZOTU") {
  missing_inputs_cor <- c(missing_inputs_cor, list(
    "ZOTU alpha diversity file" = input$uploadZOTUphyseq$datapath
  ))

  # filter list for missing inputs
  missing_inputs_cor <- lapply(missing_inputs_cor, function(x) {
    if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
      NULL
    } else {
      x
    }
  })
}

# for analysis type = both: check all missing inputs (and files) ###do ifelse here later
if (input$analysis_type == "both") {
  missing_inputs_cor <- c(missing_inputs_cor, list(
    "OTU count transformation file" = input$uploadOTUphyseq$datapath,
    "ZOTU count transformation file" = input$uploadZOTUphyseq$datapath))

  # filter list for missing inputs
  missing_inputs_cor <- lapply(missing_inputs_cor, function(x) {
    if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
      NULL
    } else {
      x
    }
  })
}


  # add error message here
  missing <- names(missing_inputs_cor)[sapply(missing_inputs_cor, is.null)]

  #catch and print missing inputs
  if (length(missing) > 0) {
    spsComps::shinyCatch(

      stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")
    return()
  }

  #set analysis variable for pwalk
  analysis_type <- case_when(
    input$analysis_type == "OTU" ~ "OTU",
    input$analysis_type == "ZOTU" ~ "ZOTU",
    input$analysis_type == "both" ~ c("OTU", "ZOTU")
  )



  meta_file <- gsub("\\\\", "/", input$uploadmeta$datapath)

  # loop through all variables using pwalk
  purrr::walk(analysis_type,
               ~rmarkdown::render(
                 input = system.file("rmd", "10_TU_correlations_v1.5.Rmd", package = "golemMB"),
                 output_format = "pdf_document",
                 output_dir = output_dir10,
                 output_file = paste("10_correlation_analysis", .x, "pdf", sep = "."),
                 intermediates_dir = output_dir10,
                 params = list(
                   workdir = workdir,
                   TU_type = .x,
                   meta_file = meta_file,
                   abundance_cutoff = input$abundance_filter,
                   prevalence_cutoff = input$prevalence_filter,
                   support_cutoff = input$support_cutoff,
                   pval_cutoff = input$pval_cutoff,
                   r_cutoff = input$r_cutoff
                 )))

  # download output files
  output$download_corr <- download_function("outputs_correlation_analysis", output_dir10)

})


### ~~~end calculate TU correlations~~~ ###


### ~~~BOOKMARKING~~~ ###
onRestored(function(state) {
  updateSelectInput(session, "var_test", choices = names(meta_file))
  updateSelectInput(session, "point_label", choices = c((names(meta_file)), "none"))
  updateSelectInput(session, "shape_var", choices = c((names(meta_file)), "none"))
  updateSelectInput(session, "samples_outlier", choices = c(meta_file[["sampleID"]], "none"))
  updateSelectInput(session, "samples_excluded", choices = c(meta_file[["sampleID"]], "none"))
  updateSelectInput(session, "var_color", choices = names(meta_file))
  # Add more update functions for other types of inputs if needed
})




setBookmarkExclude(
  c("plot_FA", "plot_DA", "plot_beta", "dist_ma", "count_transform",
    "plot_alpha", "calc_alpha", "add_comp", "load_inputs",
    "run", "sidebarCollapsed", "sidebarItemExpanded", "meta_rows_current",
    "meta_cells_selected", "meta_columns_selected", "meta_rows_selected",
    "meta_cell_clicked", "meta_search", "meta_state", "meta_rows_all") #rmv_comp
)

### ~~~end BOOKMARKING~~~ ###

} #<- end of server function
