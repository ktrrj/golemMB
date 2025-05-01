#### script with accessory functions ####

#### general accessory functions ####

# function to display and format color set vectors
show_color_set <- function(x) {
  if(is.null(x)) {
    return("No colors defined - Plot all in black.")
  } else {
    suppressWarnings(require(kableExtra))
    # display table with kbl
    kbl(data.frame(
      "Original" = if( !is.null(names(x)) ) names(x) else seq(1:length(x)),
      "Color set" = unname(x)),
      align = "l", booktabs = TRUE, linesep = "") %>%
      column_spec(column = 1, width = "2in") %>%
      column_spec(column = 2, width = "1in", background = gplots::col2hex(x)) %>%
      kable_styling(font_size = 10, latex_options = c("HOLD_position")) #, "scale_down"
  }
}

# function to display and format label set vectors
show_label_set <- function(x) {
  if(length(x) == 0) {
    stop("No label set provided")
  } else {
    require(pander)
    # translate html to markdown
    label_set_md <- str_replace_all(x, "</{0,1}sup>", "^") %>%
      str_replace_all("</{0,1}sub>", "~")
    # display table with pander
    pander(data.frame("Original" = names(x), "Label set" = label_set_md),
                   justify = "left")
  }
}

# function to check validity of comparison list
# check_comp_list <- function(x, psl = physeq.list) {
#   require(phyloseq)
#   require(stringr)
#   psid <- na.omit(stringr::str_extract(names(psl), "^ori"))
#   res.list <- list()
#   unique_var <- as.character(unlist(unique(
#     as.data.frame(phyloseq::sample_data(psl[[psid]])[, var_test])
#   )))
#   check_level <- function(x) {
#     z <- x %in% unique_var
#   }
#   if(all(lengths(comp) == 2))
#     res.list[[1]] <- "Check passed: each list element is of length 2" else
#       stop("STOP: not all list elements of length 2!")
#   if(all(base::sapply(comp, check_level) == TRUE))
#     res.list[[2]] <- paste("Check passed: each comparison element is a valid level of variable:", var_test) else
#       stop(paste("STOP: not all elements of the comparisons are valid levels of variable:"), var_test)
#   return(res.list)
# }
# function to validate contrasts list
check_comp_list <- function(x, ref) {
  res.list <- list()
  if(is.null(x)) {
    res.list[[1]] <- "Doing no statistics!"
  } else if(typeof(x) == "list") {
    if(any(class(ref) == "data.frame")) {
    unique_var <- as.character(na.omit(unique(pull(ref[, var_test]))))
    } else if (class(ref) == "list") {
      require(phyloseq)
      require(stringr)
      psid <- na.omit(stringr::str_extract(names(psl), "^ori"))
      unique_var <- as.character(unlist(unique(
        as.data.frame(phyloseq::sample_data(psl[[psid]])[, var_test]))))
    }
    check_level <- function(x) {
      unique(unlist(x)[!unlist(x) %in% unique_var])
    }
    if(all(lengths(x) == 2))
      res.list[[1]] <- "Check passed: each list element is of length 2" else
        res.list[[1]] <- "STOP: NOT all list elements of length 2!"
    if(length(check_level(x)) == 0)
      res.list[[2]] <- paste("Check passed: each comparison element is a valid level of variable:", var_test) else {
        res.list[[2]] <- paste("STOP: The elements", check_level(x), "of the comparisons are NOT valid levels of variable:", var_test)}
  } else if(typeof(x) == "character" & x == "all" & any(class(ref) == "data.frame")) {
    res.list[[1]] <- "Doing all possible pairwise comparisons"
  } else if(typeof(x) == "character" & x != "all" & any(class(ref) == "data.frame")) {
    res.list[[1]] <- "please use the argument 'all' for all possible pairwise comparisons"
  }
  return(res.list)
}

# check presence and consistencies of table, tree and sequence files
# table_file <- OTU_file
# tree_file <- OTU_tree_file
# seq_file <- OTU_fasta_file
# meta_file = meta_file
check_microbiome_files <- function(table_file, tree_file, seq_file,
                                   meta_file = meta_file) {
  # list to store results
  rl <- list()

  # helper function to check file existence
  check_exist <- function(x) {
    if( !file.exists(x) ) {
      cat(red$bold("\u2717") %+%
            red$bold(x) %+% red$bold(" does NOT exist!"))
    } else {
      cat(green$bold("\n\u2713") %+%
            green(x) %+% green(" exists!"))
    }
  }

  # helper function to check file readability
  check_read <- function(x, y) {
    if( !exists(x) ) {
      red$bold("\u2717") %+%
        red$bold(y) %+% red$bold(" cannot be read! Check file format!")
    } else {
      green$bold("\n\u2713") %+%
        green(y) %+% green(" sucessfully read!")
    }
  }

  # when all files are present check content
  if( all(sapply(c(table_file, tree_file, seq_file, meta_file), file.exists)) ){

    # extract TUs from table
    tu <- read_delim(table_file, show_col_types = FALSE)
    rl[[length(rl)+1]] <- check_read("tu", table_file)
    if( exists("tu") ) {tab.tu <- dplyr::pull(tu[, 1])}

    # extract tree tip labels
    tree <- ape::read.tree(tree_file)
    rl[[length(rl)+1]] <- check_read("tree", tree_file)
    if( exists("tree") ) {tree.tu <- tree$tip.label}

    # extract fasta headers
    ref_seqs <- Biostrings::readDNAStringSet(
      file = seq_file, format = "fasta", nrec = -1L, skip = 0L,
      seek.first.rec = FALSE, use.names = TRUE)
    rl[[length(rl)+1]] <- check_read("ref_seqs", seq_file)
    if( exists("ref_seqs") ) {fasta.tu <- names(ref_seqs)}

    # import metadata
    meta <- readRDS(meta_file)
    rl[[length(rl)+1]] <- check_read("meta", meta_file)

    # compare sample names from TU table and meta data
    if( all(exists("meta"), exists("tu")) ) {
      not_in_tu <- row.names(meta)[!(row.names(meta) %in% colnames(tu))]
      if(length(not_in_tu) > 0) {
        rl[[length(rl)+1]] <- red$bold("\n\u2717") %+% red$bold(" The following sample names:\n  ") %+%
          red$bold(paste(not_in_tu, collapse = ", ")) %+%
          red$bold(" \n  are NOT present in the TU table!")
      } else {
        rl[[length(rl)+1]] <- green$bold("\n\u2713") %+%
          green(" All sample names from metadata are present in the TU table.")
      }
    }
    if( all(exists("meta"), exists("tu")) ) {
      # summarize extracted TUs
      common <- intersect(c(tree.tu, fasta.tu), tab.tu)
      if ( all(all(tree.tu %in% common) & all(tree.tu %in% common) & all(tree.tu %in% common)) ) {
        rl[[length(rl)+1]] <- green$bold("\n\u2713") %+%
          green(" Similar taxonomic units in table, tree and sequence file present!\n")
      } else {
        rl[[length(rl)+1]] <- red$bold("\n\u2717") %+%
          red$bold(" Differences in taxonomic unit descriptions in table, tree and sequence file detected!\n")
      }
    }
  } else {

    # error message if one or more files are missing
    rl[[length(rl)+1]] <- red$bold("\n\u2717") %+%
      red$bold(" At least one file is missing!\n")
  }

  # output
  purrr::walk(c(table_file, tree_file, seq_file, meta_file), check_exist)
  cat(unlist(rl))
}
#
# check_microbiome_files(OTU_file, OTU_tree_file, OTU_fasta_file, meta_file)
# check_microbiome_files(ZOTU_file, ZOTU_tree_file, ZOTU_fasta_file, meta_file)

# check consistency of taxonomy file
# e.g. whether a particular genus has similar 'parental' taxonomy
# check is done over all taxonomy levels to check entries with NA tax levels
# cave:
# No check of correct taxonomy assignment!
# No check for singular tax entries!
check_tax_consistency <- function(tabfile) {
  # settings
  require(crayon)
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)
  if ( file.exists(tabfile) ) {


    # extract file name
    fname <- str_match(tabfile, ".*\\/(Z?OTU\\/[:graph:]+\\.tab$)")[,2]

    # define tax levels from df
    taxlvl <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    if( str_detect(tabfile, "\\.tab$|\\.txt$") ) {
      x <- read.csv(tabfile, sep = "\t") %>% dplyr::select(taxonomy) %>%
        separate_wider_delim(taxonomy, delim = ";", names = taxlvl) %>%
        mutate(across(everything(), \(x) na_if(x, "")))
    } else {
      cat(red$bold("\u2717") %+% red(" OTU/ZOTU/ASV table not in text format!\n"))
      stop()}
    # x[14, "Class"] <- "Lachnospiraceae"

    # initiate list to store results
    res  <- list()

    # screen for duplicates in all tax level
    # i <- 6
    for ( i in 6:1 ) {
      current_lvl <- taxlvl[i]
      y <- x %>%
        group_by(across(all_of(taxlvl[1:i]))) %>%
        summarize(n = n()) %>% ungroup()
      # store duplicated in vector
      z <- unlist(as.vector(y[duplicated(y[, current_lvl]), current_lvl]))
      prob <- y %>% dplyr::filter(.data[[current_lvl]] %in% z[!is.na(z)])
      res[[current_lvl]] <- if ( nrow(prob) == 0 ) NULL else prob
    }
    # evaluate results
    if( length(res) == 0 ) {
      cat(green$bold("\n\u2713") %+%
            green(" Taxonomy assignment of file:\n") %+%
            green("  ../") %+% green(fname) %+%
            green("\n  seems structurally OK\n") %+%
            yellow("  Attention! Genus/Species assignment is not checked!\n"))
    } else {
      # extract unique results
      res <- bind_rows(res)
      res <- res[complete.cases(res), ]
      # helper function
      print_and_capture <- function(x) {
        cyan(paste(capture.output(print(format(x)[-3L]))[-1], collapse = "\n"))}
      # generate message text
      cat(red$bold("\n\u2717") %+% red$bold(" Warning!") %+% red(" In file") %+%
            yellow(" ../") %+% yellow(fname) %+%
            red(" the following taxonom(y/ies) contain inconsistencies:\n"),
          print_and_capture(res))
    }
    # detach("package:crayon", unload = TRUE)
  } else {
    # error message if file does not exist
    cat(red$bold("\u2717") %+%
          red$bold(" The following file does not exist:\n") %+%
          red("  ") %+% red(tabfile))
  }
}

# check all input parameters
# var_test <- var_list[[1]]
# var_formula_test <- var_formula_list[[1]]

check_analysis_parameters <- function(
    var_test, var_formula_test, shape_var, point_label, comp_list, label_set,
    color_set, grouped_pt, color_param, label_set1c, color_set1c, data_param,
    label_set1d, meta_file, filter_var, var_color) {

  # initiate empty results list
  rl <- list()

  ## check presence of variables in metadata ##

  # load metadata and extract variables
  meta <- readRDS(meta_file)
  vars <- names(meta)

  # helper function to check presence of variables
  check_vars <- function(x, y = deparse(substitute(x))) {
    if ( x %in% vars ) {
      cat(green$bold("\n\u2713") %+%
            green(" Variable '") %+%
            green(y) %+%
            green("': [") %+%
            yellow(x) %+%
            green("] is present in metadata!"))
      TRUE
    } else {
      cat(red$bold("\u2717") %+%
            red$bold(" Variable [") %+%
            red$bold(y) %+%
            red$bold("]: '") %+%
            yellow$bold(x) %+%
            red$bold("' is NOT present in metadata!"))
      FALSE
    }
  }
  # indicate new analysis set
  cat(blue$bold("\n\nNew analysis set"))

  # check 'filter_var'
  if( !is.null(samples_excluded) ) {
    rl[["filter_var"]] <- check_vars(filter_var)
  }

  # check 'var_test'
  rl[["var_test"]] <- check_vars(var_test)

  # check 'var_formula_test'
  cat("\nChecking variable 'var_formula_test':")
  formula_check <- purrr::map_lgl(unlist(str_split(var_formula_test, "[[:punct:]&&[^_]]")), check_vars)
  rl <- append(rl, formula_check)

  # check of 'samp_ctrl' variable
  rl[["samp_ctrl"]] <- check_vars("samp_ctrl", "samp_ctrl")

  # check whether factor levels 'sample' and 'neg_ctrl' are present in 'samp_ctrl'
  if( rl[["samp_ctrl"]] ) {
    if( all(c("sample", "neg_ctrl") %in% levels(meta[, "samp_ctrl"])) ) {
      cat(green$bold("\n\u2713") %+%
            green(" Factor levels 'sample' and 'neg_ctrl' are present in variable [") %+%
            yellow("samp_ctrl") %+% green("]!"))
      rl[["samp_ctrl_lvls"]] <- TRUE
    } else {
      cat(red$bold("\u2717") %+%
            red$bold(" Factor levels 'sample' and/or 'neg_ctrl' NOT present in variable [") %+%
            yellow$bold("samp_ctrl") %+% red$bold("]!"))
      rl[["samp_ctrl_lvls"]] <- FALSE
    }
  }

  # check of 'quant_reading' variable
  rl[["quant_reading"]] <- check_vars("quant_reading", "quant_reading")

  # check whether 'quant_reading' is numeric
  if( rl[["quant_reading"]] ) {
    if( is.numeric(meta[, "quant_reading"]) ) {
      cat(green$bold("\n\u2713") %+%
            green(" Variable [") %+%
            yellow("quant_reading") %+% green("] is numeric!"))
      rl[["quant_reading_num"]] <- TRUE
    } else {
      cat(red$bold("\u2717") %+%
            red$bold(" Variable [") %+%
            yellow$bold("quant_reading") %+% red$bold("] is NOT numeric!"))
      rl[["quant_reading_num"]] <- FALSE
    }
  }

  # check 'shape_var'
  if( !is.null(shape_var) ) {rl[[length(rl)+1]] <- check_vars(shape_var)
  } else cat("\nNo 'shape_var' given.")

  # check 'point_label'
  if( !is.null(point_label) ) { if( point_label == "sample_ID" ) {
    rl[["point_label"]] <- check_vars("sampleID", "point_label")} else {
      rl[["point_label"]] <- check_vars(point_label)
    } } else cat("No 'point_label' given.")

  # check 'var_color' (alpha diversity quick plot)
  rl[["var_color"]] <- check_vars(var_color)

  # check parameters for grouped plots
  if( grouped_pt ) {
    rl[["color_param"]] <- check_vars(color_param)
    rl[["data_param"]] <- check_vars(data_param)
  }

  ## check validity of comparison list ##

  # extract 'var_test' levels
  unique_var <- as.character(unique(meta[, var_test]))
  unique_var <- unique_var[!is.na(unique_var)]

  # identify vector elements which are not valid levels of variable
  comp_pres <- unique(unlist(comp_list)[!unlist(comp_list) %in% unique_var])

  # perform checks
  if( typeof(comp_list) == "list" ) {
    if( all(lengths(comp_list) == 2) ) {
      cat(green$bold("\n\u2713") %+%
            green(" Each element of the comparison list is of length 2!"))
      rl[["comp_list_2elements"]] <- TRUE
      if( length(comp_pres) == 0 ) {
        cat(green$bold("\n\u2713") %+%
              green(" Each comparison element is a valid level of variable [") %+%
              yellow(var_test) %+% green("]!"))
        rl[["comp_list_levels"]] <- TRUE
      } else {
        cat(red$bold("\u2717") %+%
              red$bold(" The element(s): ") %+% red$bold(paste(comp_pres, collapse = ", ")) %+%
              red$bold(" of the comparisons are NOT valid levels of variable: [") %+%
              yellow$bold(var_test) %+% red$bold("]"))
        rl[["comp_list_levels"]] <- FALSE
      }
    } else {
      cat(red$bold("\n\u2717") %+%
            red$bold(" NOT all list elements are of length 2!"))
      rl[["comp_list_2elements"]] <- FALSE
    }
  } else {
    cat(red$bold("\n\u2717") %+% red$bold(" 'comp_list' is not a list!"))
    rl[["comp_list_nolist"]] <- FALSE
  }

  ## check validity of label and color vectors ##
  # x = vector, y = variable
  check_color_label <- function(x, y, type = c("color", "label")) {
    # identify unique levels of variable
    unique_var <- as.character(unique(meta[, y]))
    unique_var <- unique_var[!is.na(unique_var)]

    # extract vector name
    z <- deparse(substitute(x))

    # check content
    if( type == "color" &  is.null(names(x)) ) {
      if( length(unique_var) <= length(x) ) {
        cat(green$bold("\n\u2713") %+%
              green(" Length of unnamed color vector '") %+%
              green(deparse(substitute(x))) %+%
              green("' is sufficient to cover all levels of variable [") %+%
              yellow(y) %+% green("]!"))
        rl[[paste0(z, "_length")]] <- TRUE
      } else {
        cat(red$bold("\n\u2717") %+%
              red$bold(" Length of unnamed color vector '") %+%
              red$bold(deparse(substitute(x))) %+%
              red$bold("' is NOT sufficient to cover all levels of variable [") %+%
              yellow$bold(y) %+% red$bold("]!"))
        rl[[paste0(z, "_length")]] <- FALSE
      }
    } else {

      if( all(unique_var %in% names(x)) & length(unique_var) <= length(x) ) {
        cat(green$bold("\n\u2713") %+%
              green(" All levels of variable [") %+%
              yellow(y) %+% green("] are present in vector '") %+%
              green(deparse(substitute(x))) %+% green("'!"))
        rl[[paste0(z, "_levels")]] <- TRUE
      } else {
        not_present <- unique_var[!(unique_var %in% names(x))]
        cat(red$bold("\n\u2717") %+%
              red$bold(" The following levels of variable [") %+%
              yellow$bold(y) %+% red$bold("] are NOT present in vector '") %+%
              red$bold(deparse(substitute(x))) %+% red$bold("':") %+%
              red$bold(paste(not_present, collapse = ", ")))
        rl[[paste0(z, "_levels")]] <- FALSE
      }
    }
  }

  # perform checks
  check_color_label(label_set, var_test, type = "label")
  check_color_label(color_set, var_test, type = "color")
  if( grouped_pt ) {
    check_color_label(label_set1c, color_param, type = "label")
    check_color_label(color_set1c, color_param, type = "color")
    check_color_label(label_set1d, data_param, type = "label")
  }

  ### check if 'Rhea' folder is present
  return(unlist(rl))
}

# function to generate names for subfolders
subfolder_name <- function(x = var_list, incl_date = TRUE) {
  sub_names <- lapply(x, function(x) {paste(c(if( incl_date ) {format(Sys.time(), '%y%m%d')}, x), collapse = "_")})
  sub_names[duplicated(sub_names)] <- imap_chr(sub_names[duplicated(sub_names)], \(x, idx) {paste(x, idx, sep = "_")} )
  return(sub_names)
}

# function to generate proper inputs for ungrouped and grouped data analysis
# 'ulist' = parameter list for ungrouped plots
# 'glist' = parameter list for grouped plots
generate_input_list <- function(ulist, glist, grouped = grouped_list) {
  l <- list()
  for( i in 1:length(grouped) ) {
    if( grouped[[i]] ) {
      l[[i]] <- glist[[i]]
    } else {
      l[[i]] <- ulist[[i]]
    }
  }
  return(l)
}


# function to check number of replicates and identify filtered samples
filtered_levels <- function(before, group_param = var_test,
                            min.repl = min_repl)  {

  # check and filter for replicates
  after <- before %>%
    dplyr::group_by(!!!syms(group_param), measure) %>%
    dplyr::mutate(freq = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(freq >= min.repl) %>%
    dplyr::select(-freq)

  # message which factor levels have less than required number of replicates
  # these samples are removed before pairwise comparisons
  if(nrow(before) > nrow(after)) {
    rm_levels <- as.character(unique(dplyr::pull(before[, group_param]))[!(unique(dplyr::pull(before[, group_param])) %in% unique(dplyr::pull(after[, group_param])))])
    if(length(rm_levels) > 1) {
      rm_print <- paste(rm_levels, collapse = ", ")
    } else rm_print <- rm_levels
    info_rm <- stringr::str_replace_all(paste0("Because there were less than ", min.repl ," replicates, the following levels of variable '",
                                      group_param, "' were removed from statistical comparisons: ", rm_print, "."), "_", "\\\\_")
    # make df for display table
    filter_df <- before %>%
      dplyr::filter(get(group_param) %in% rm_levels) %>%
      dplyr::group_by(!!!syms(group_param), measure) %>%
      dplyr::mutate(Frequency = n()) %>%
      dplyr::ungroup() %>%
      dplyr::select(all_of(c(group_param, "Frequency"))) %>%
      dplyr::distinct() %>%
      tibble::column_to_rownames(var = group_param)

  } else {
    rm_levels <- NULL
    filter_df <- NULL
    info_rm <- NULL
  }
  return(list(
    rm_levels,
    filter_df,
    info_rm
  ))
}

# function to check for normal distribution of data
norm_test <- function(df) {
  x <- df %>%
    dplyr::group_by(measure) %>%
    dplyr::summarize(Shapiro.test = shapiro.test(value)$p.value) %>%
    # return "Wilcoxon" for non-normal and "t-test" for normal distribution
    dplyr::mutate(
      Normal.dist = ifelse(Shapiro.test < 0.05, "no", "yes"),
      Stat.test = ifelse(Shapiro.test < 0.05, "Wilcoxon", "t-test")) %>%
    dplyr::select(Measure = measure, Shapiro.test, Normal.dist, Stat.test)
}

# function to format time difference
# https://stackoverflow.com/questions/27312292/convert-seconds-to-days-hoursminutesseconds
format_difftime <- function(first, second) {
  x <- as.numeric(as.difftime(second - first, units = "secs"), units = "secs")
  days = round(x %/% (60 * 60 * 24))
  hours = round((x - days*60*60*24) %/% (60 * 60))
  minutes = round((x - days*60*60*24 - hours*60*60) %/% 60)
  seconds = round(x - days*60*60*24 - hours*60*60 - minutes*60, digits = 1)
  days_str = ifelse(days == 0, "", paste0(days, " days "))
  hours_str = ifelse((hours == 0 & days == 0), "", paste0(hours, " hrs "))
  minutes_str = ifelse((minutes == 0 & days == 0 & hours == 0), "", paste0(minutes, " min "))
  seconds_str = paste0(seconds, " sec")
  final_str = paste0(days_str, hours_str, minutes_str, seconds_str)
  return(final_str)
}

#### manipulate / check phyloseq objects ####

# function to filter phyloseq object based on prevalence
phyloseq_filter_prev <- function(ps,
                                 min.count = min_count,
                                 min.count.prev = min_count_prev,
                                 min.rel.abund = min_rel_abund,
                                 min.rel.abund.prev = min_rel_abund_prev) {
  require(phyloseq)

  phyloseq::filter_taxa(
    ps,
    function(x){(sum(x >= min.count) > length(x)*(min.count.prev/100)) | ((sum(x >= min.count) > (length(x)*(min.rel.abund.prev/100))) & (mean(x/phyloseq::sample_sums(ps)) > (min.rel.abund/100)))},
    prune = TRUE)
}

# function to compare filtered phyloseq with original object
# return df with diff stats
id_taxa_removed <- function(before, after, min_ct = NULL) {

  require(phyloseq)
  require(dplyr)

  y <- as.data.frame(phyloseq::otu_table(before))
  z <- as.data.frame(phyloseq::otu_table(after))

  # identify removed OTUs
  prev_removed <- row.names(y)[!(row.names(y) %in% row.names(z))]

  # calculate prevalence in original phyloseq
  prev <- apply(X = y, MARGIN = 1, FUN = function(x){sum(x > 0)})

  # calculate relative abundance
  rel_abund <- apply(X = y, MARGIN = 1, FUN = function(x){mean(x/phyloseq::sample_sums(before))})

  # construct df
  df <- data.frame(Reads = rowSums(y[prev_removed, ]),
                   Prevalence = prev[prev_removed],
                   Mean_abundance = paste0(round(rel_abund[prev_removed]*100, 3), "%")) %>%
    mutate(Rel_prev = paste0(round((Prevalence * 100) / nsamples(before), 1), "%"),
           Mean_reads = round(Reads/nsamples(before), 1)) %>%
    dplyr::select(Reads, Mean_reads, Mean_abundance, Prevalence, Rel_prev) %>%
    # relocate(Mean_reads, .after = Reads) %>%
    dplyr::arrange(desc(Reads))

  # calculate prevalence above threshold counts in orig. phyloseq
  if(!is.null(min_ct)){
    prev_thresh <- apply(X = y, MARGIN = 1, FUN = function(x){sum(x > min_ct)})
    df$Prev_thresh <- prev_thresh[prev_removed]
    df$Rel_prev_thresh <- paste0(round((df$Prev_thresh * 100) / nsamples(before), 1), "%")
  }
  return(df)
}

# filter transformed physeq data based on OTUs identified in the untransformed object
filter_trans_ps <- function(df, tr_method) {
  require(phyloseq)
  phyloseq::prune_taxa(!phyloseq::taxa_names(physeq.list[[tr_method]]) %in% row.names(df), physeq.list[[tr_method]])
}
