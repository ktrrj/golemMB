# ### create named lists with all possible inputs
#
# # input master list
# all_inputs <- list(
#   "Metadata file" = input$uploadmeta$datapath,
#   "analysis variable" = input$var_test,
#   "analysis variable name" = input$var_name,
#   "variable shape" = input$shape_var,
#   "point label" = input$point_label,
#   "label order" = input$label_set,
#   "label labels" = input$labels,
#   "color set" = input$color_set,
#   "analysis type" = input$analysis_type,
#   "Minimum sample counts" = input$min_sample_counts,
#   "alpha diversity color" = input$var_color,
#   "analysis variable formula" = input$var_formula,
#   "Sample outliers" = input$samples_outlier,
#   "BugBase prediction file" = input$uploadbugbasepredicts$datapath,
#   "BugBase vsearch stats" = input$uploadbugbasevsearch$datapath,
#   "comparison list" = input$comp_list,
#   "PICRUSt2 out_stratified folder" = input$uploadpicrustfolder$datapath,
#   "PICRUSt2 stats file" = input$uploadpicruststats$datapath
#   )
#
# # OTU specific input list
# otu_specific_inputs <- list(
#   "OTU table file" = input$uploadOTUtab$datapath,
#   "OTU fasta file" = input$uploadOTUfasta$datapath,
#   "OTU tree file" = input$uploadOTUtre$datapath,
#   "OTU Import file" = input$uploadOTUphyseq$datapath,
#   "OTU distance matrices files" = input$uploadOTUdists$datapath,
#   "OTU Genus distance matrices files" = input$uploadOTUdistsgenus$datapath,
#   "OTU count transformation files" = input$uploadOTUtrans$datapath,
#   "OTU Genus count transformation files" = input$uploadOTUtransgenus$datapath,
#   "OTU ordination files" = input$uploadOTUord$datapath,
#   "OTU Genus ordination files" = input$uploadOTUordgenus$datapath
# )
#
# # ZOTU specific input list
# zotu_specific_inputs <- list(
#   "ZOTU table file" = input$uploadZOTUtab$datapath,
#   "ZOTU fasta file" = input$uploadZOTUfasta$datapath,
#   "ZOTU tree file" = input$uploadZOTUtre$datapath,
#   "ZOTU Import files" = input$uploadZOTUphyseq$datapath,
#   "ZOTU distance matrices files" = input$uploadZOTUdists$datapath,
#   "ZOTU Genus distance matrices files" = input$uploadZOTUdistsgenus$datapath,
#   "ZOTU count transformation files" = input$uploadZOTUtrans$datapath,
#   "ZOTU Genus count transformation files" = input$uploadZOTUtransgenus$datapath,
#   "ZOTU ordination files" = input$uploadZOTUord$datapath,
#   "ZOTU Genus ordination files" = input$uploadZOTUordgenus$datapath
# )
#
# #
# missing_inputs <- case_when(
#   input$run == TRUE ~ list(all_inputs[1:10], ifelse(input$analysis_type == "OTU", otu_specific_inputs[1:3], zotu_specific_inputs[1:3]))
# )
#
# # create a list with all inputs that are either empty, 0 or NULL
# missing_inputs_import <- lapply(missing_inputs_import, function(x) {
#   if (any(length(x) == 0) || any(x == "") || any(is.null(x))) {
#     NULL
#   } else {
#     x
#   }
# })
#
# # pull the names of the missing inputs to later paste them into an error message
# missing <- names(missing_inputs_import)[sapply(missing_inputs_import, is.null)]
#
# # print error message for missing inputs and stop calculations
# if (length(missing) > 0) {
#   spsComps::shinyCatch(
#
#     stop(paste("The following inputs are missing: "), paste(missing, collapse = ", ")), blocking_level = "error")
#
#   return()
# }
