#### script to provide functions for data import and export #####

# export statistics from comparative analyses (Alpha diversity, BugBase)
export_data_stats <- function(df, group_param = var_test, group_formula = var_formula,
                              measure_vars = measure.vars, p_threshold = cutoff_pval) {
  # prepare empty list to store data
  stat_methods <- c("Shapiro", "Kruskal", "Wilcox_pair", "ANOVA", "t-test_pair", "AOV", "Tuckey")
  stat_tests <- setNames(vector("list", length(measure_vars)*7), paste(rep(measure_vars, each = 7), stat_methods, sep = "_"))

  # transform df to wide format
  df1 <- tidyr::pivot_wider(df, id_cols = !any_of(c("measure", "value")), names_from = measure) %>%
    base::droplevels()

  # remove samples with low number of replicates
  if(!is.null(low_replicates[[1]])) {
    df2 <- dplyr::filter(df1, !get(group_param) %in% low_replicates[[1]])
  } else {df2 <- df1}

  # calculate stats based on data distribution
  for (i in measure_vars) {
    if (shapiro.test(pull(df2[, i]))$p.value < 0.05) {
      stat_tests[[paste0(i, "_Kruskal")]] <- broom::tidy(kruskal.test(as.formula(paste(i, group_param, sep = " ~ ")), data = df2))
      stat_tests[[paste0(i, "_Wilcox_pair")]] <- broom::tidy(pairwise.wilcox.test(dplyr::pull(df2[, i]), dplyr::pull(df2[, group_param]),
                                                                           p.adjust.method = p.adj.method)) %>%
        dplyr::filter(p.value < p_threshold)
    } else {
      stat_tests[[paste0(i, "_ANOVA")]] <- broom::tidy(anova(aov(as.formula(paste(i, group_param, sep = " ~ ")), data = df2)))
      stat_tests[[paste0(i, "_t-test_pair")]] <- broom::tidy(pairwise.t.test(pull(df2[, i]), dplyr::pull(df2[, group_param]),
                                                                      p.adjust.method = p.adj.method)) %>%
        dplyr::filter(p.value < p_threshold)
    }
    stat_tests[[paste0(i, "_Shapiro")]] <- broom::tidy(shapiro.test(pull(df2[, i])))
    stat_tests[[paste0(i, "_AOV")]] <- broom::tidy(aov(as.formula(paste(i, group_formula, sep = " ~ ")), data = df2))  %>%
      dplyr::filter(p.value < p_threshold)
    stat_tests[[paste0(i, "_Tuckey")]] <- broom::tidy(TukeyHSD(aov(as.formula(paste(i, group_formula, sep = " ~ ")), data = df2)))  %>%
      dplyr::filter(adj.p.value < p_threshold)
  }

  # remove empty list slots
  stat_tests <- stat_tests[lengths(stat_tests) > 0]

  # summarize similar tests to one df
  stat_tests_sum <- setNames(vector("list", length(stat_methods)), stat_methods)
  for (j in stat_methods) {
    stat_tests_sum[[j]] <- dplyr::bind_rows(stat_tests[c(str_detect(string = names(stat_tests), pattern = j))], .id = "measure")
  }

  # edit Shapiro df
  stat_tests_sum[["Shapiro"]] <- dplyr::mutate(stat_tests_sum[["Shapiro"]],
                                        measure = str_remove_all(measure, "_Shapiro"),
                                        normally_distributed = ifelse(p.value < 0.05, "no", "yes"))

  # remove empty list slots
  stat_tests_sum <- stat_tests_sum[lengths(stat_tests_sum) > 0]

  # export group-sorted raw data
  export_Prism <- function(df, df_wide, measure_vars, group_param) {
    data_set <- setNames(vector("list", length(measure_vars)), measure_vars)
    for (k in measure_vars) {
      data_set[[k]] <- df %>% dplyr::filter(measure == k) %>%
        dplyr::select(sample_ID, all_of(group_param), value) %>%
        tidyr::pivot_wider(names_from = all_of(group_param), values_from = value) %>%
        dplyr::select(sample_ID, all_of(levels(pull(df_wide[,group_param])))) #order columns according to factor levels
    }
    names(data_set) <- paste0("Raw_", names(data_set))
    return(data_set)
  }
  stat_data <- export_Prism(df = df, df_wide = df2, measure_vars = measure_vars, group_param = group_param)

  # function to calculate SEM
  sem <- function(x) sd(x)/sqrt(length(x))

  # summarize raw data
  stat_data_summary <- setNames(
    list(
      df %>%
        dplyr::select(all_of(group_param), measure, value) %>%
        dplyr::group_by(!!!syms(group_param), measure) %>%
        dplyr::summarize(mean = mean(value),
                  SD = sd(value),
                  median = median(value),
                  SEM = sem(value)) %>%
        arrange(measure)),
    "summarized_data"
  )

  # combine values with statistics and return
  stat_tests_sum <- c(stat_data_summary, stat_tests_sum, stat_data)
  return(stat_tests_sum)
}

# import data from differential or functional abundance analyses

# function to import differential abundance data to metacoder environment
# import list of microbiomeMarker analysis outputs (modified phyloseq objects)
# as a data frame compatible for insertion as a 'diff_table' in a 'metacoder'
# environment and return the final 'metacoder' environment
# ps = ps.ori.filt
# da_list = da_extract[["da_data"]]
# da_sum_list <- da_summary
# cutoff.pval = cutoff_pval
# cutoff.fcda = cutoff_fcda
# cutoff.lda = cutoff_lda
# cutoff.W = cutoff_W
# cutoff.aldex = cutoff_aldex

import_DA_metacoder <- function(ps, da_list, da_sum_list,
                                cutoff.pval = 0.05,
                                cutoff.fcda = 2,
                                cutoff.lda = 2,
                                cutoff.W = 0.75,
                                cutoff.aldex = 1
                                ) {
  require(tidyverse)
  require(metacoder)
  require(microbiomeMarker)

  # import phyloseq in metacoder
  x <- metacoder::parse_phyloseq(ps)

  # extract arbitrary taxon ids (including supertaxa) and taxonomy ('feature') from metacoder environment
  ref <- tibble(
    taxon_id = names(metacoder::classifications(x)),
    feature = metacoder::classifications(x)
  )
  # extract number of taxa
  tax_units <- tibble(
    taxon_id =  names(n_obs(x, data = "tax_data")),
    tu_number = as.numeric(n_obs(x, data = "tax_data"))
  )

  # calculate taxon abundances
  x$data$tax_abund <- calc_taxon_abund(x, "otu_table",
                                       cols = x$data$sample_data$sample_id)

  # calculate 'group' abundances
  x$data$type_abund <- calc_group_median(x, "tax_abund",
                                         cols = x$data$sample_data$sample_id,
                                         groups = x$data$sample_data[[var_test]])

  # initiate empty list to store results
  store_list <- list()
  superstore <- list()
  methods <- c("edger", "lv", "deseq2", "ancombc", "mgzig", "lefse", "aldex", "ancom", "mgziln")

  # loop through all comparisons
  for (i in names(da_list)) {
    # i <- "Cre_DSS_vs_Cre_reference"
    # loop through all methods
    for (j in methods) {
      # j <- "mgziln"
      comp_pair <- unlist(str_split(i, pattern = "_vs_"))

      # import diff abundance data using microbiomeMarker function
      if (!is.null(da_list[[i]][["mm_data"]][[j]])) {
        if (!is.null(microbiomeMarker::marker_table(da_list[[i]][["mm_data"]][[j]]))) {
          df <- microbiomeMarker::marker_table(da_list[[i]][["mm_data"]][[j]]) %>%
            as_tibble() %>%
            mutate(feature = str_replace_all(string = feature, pattern = "\\|", replacement = "\\;")) %>%
            dplyr::rename(diffabund = starts_with("ef"))

          # adjust "directionality" for lefse, mgziln and ancom
          if (j %in% c("lefse", "ancom", "mgziln")) {
            df <- df %>%
              mutate(diffabund = ifelse(enrich_group == comp_pair[2], -1*diffabund, diffabund))
          }
          if (j %in% c("mgziln")) {
            df <- df %>%
              mutate(diffabund = ifelse(enrich_group == comp_pair[2], -1*abs(diffabund), abs(diffabund)))
          }

          df <- ref %>% dplyr::left_join(df, by = "feature") %>%
            dplyr::left_join(tax_units, by = "taxon_id") %>%
            dplyr::mutate(
              treatment_1 = comp_pair[1],
              treatment_2 = comp_pair[2],
              comparison = paste(comp_pair, collapse = "_vs_"),
              method = j
            )

          if (j == "ancom") {
            df <- df %>%
              mutate(
                diffabund = ifelse(is.na(diffabund), 0, diffabund),
                W = ifelse(is.na(W) | is.null(da_list[[i]][["mm_data"]][["ancom"]]@'marker_table'), 1, W)
              ) %>%
              dplyr::select(taxon_id, treatment_1, treatment_2, comparison, method, diffabund, W, tu_number)
          } else {
            df <- df %>%
              mutate(
                diffabund = ifelse(is.na(diffabund), 0, diffabund),
                pvalue = ifelse(is.na(pvalue), 1, pvalue),
                padj = ifelse(is.na(padj), 1, padj)
              ) %>%
              dplyr::select(taxon_id, treatment_1, treatment_2, comparison, method, diffabund, pvalue, padj, tu_number)
          }
          store_list[[i]][[j]] <- df
        }} else next
    }
    superstore[[i]] <- bind_rows(store_list[[i]])
  }

  # combine all differential abundance data and add missing combinations
  all_data <- bind_rows(superstore) %>%
    tidyr::complete(taxon_id, comparison, method, fill = list(
      diffabund = 0,
      pvalue = 1,
      padj = 1,
      W = 1
    )) %>%
    fill(treatment_1, treatment_2, tu_number)
  ## remove non-significant hits
  # calculate W threshold
  ancom_counts <- all_data %>%
    group_by(method, comparison) %>%
    dplyr::summarize(n = n()) %>%
    filter(method == "ancom")
  threshold.W <- round(max(ancom_counts$n)*cutoff.W, 0)
  # filter df, set for all non-significant hits diffabund to 0
  all_data <- all_data %>%
    mutate(
      diffabund = case_when(
      (method %in% c("edger", "lv", "deseq2", "mgzig", "mgziln", "ancombc")) & ((padj >= cutoff.pval) | (abs(diffabund) < log2(cutoff.fcda))) ~ 0,
      (method == "lefse") & ((padj >= cutoff.pval) | (abs(diffabund) < cutoff.lda)) ~ 0,
      (method == "aldex") & ((padj >= cutoff.pval) | (abs(diffabund) < cutoff.aldex)) ~ 0,
      (method == "ancom") & ((W < threshold.W)  | (abs(diffabund) < log2(cutoff.fcda))) ~ 0,
      TRUE ~ diffabund),
      comparison = factor(comparison, levels = names(da_list)),
      method = factor(method, levels = methods)
    )

  # identify completed data sets
  compl_data <- all_data %>%
    anti_join(bind_rows(superstore), by = c("comparison", "method")) %>%
    group_by(comparison, method) %>%
    dplyr::summarize(n = n()) %>%
    arrange(method, comparison) %>%
    dplyr::select(method, comparison)

  if(nrow(compl_data) > 0) {
    warning(paste(nrow(compl_data), "data sets were completed due to missing data!"))
    attr(all_data, "completed") <- compl_data
  }

  x$data[["diff_table"]] <- all_data

  # summarize data of all methods
  sum_data <- all_data %>%
    filter(diffabund != 0) %>%
    dplyr::select(taxon_id, comparison, method, diffabund, tu_number) %>%
    group_by(taxon_id, comparison) %>%
    dplyr::summarize(mean_diff = mean(diffabund),
              n = n()) %>%
    dplyr::mutate(n_sig = ifelse(mean_diff < 0, -1*n, n)) %>%
    ungroup()

  sum_data <- distinct(all_data, taxon_id, comparison, tu_number) %>%
    dplyr::left_join(sum_data, by = c("taxon_id", "comparison")) %>%
    tidyr::complete(taxon_id, comparison, fill = list(
      mean_diff = 0,
      n = 0,
      n_sig = 0)) %>%
    dplyr::rename(tu_number1 = tu_number, comparison1 = comparison)

  x$data[["diff_summary"]] <- sum_data

  return(x)
  # return(superstore)
}

# function to import functional abundance data to metacoder environment
# import list of PICRUSt2/ALDEx2 analysis outputs (from pathway/PW analyses)
# as a data frame compatible for insertion as a 'diff_table' in a 'metacoder'
# environment and return the final 'metacoder' environment
import_aldex2_metacoder <- function(ax2_list, pw_hierarchy, min_eff_size, max_adj_p) {
  require(metacoder)
  # format hierarchical pathway data
  hier <- pw_hierarchy %>%
    tibble::rownames_to_column("PW") %>%
    dplyr::mutate(root = "r__Root",
           SC1 = paste0("sc1__", Superclass1),
           SC2 = ifelse(Superclass1 == Superclass2, NA, Superclass2),
           SC2 = ifelse(is.na(SC2), paste0("pw__", PW), paste0("sc2__", SC2)),
           PW2 = paste0("pw__", PW),
           PW2 = ifelse(PW2 == SC2, NA, PW2)) %>%
    tidyr::unite(taxonomy, root, SC1, SC2, PW2, sep = ";",na.rm = TRUE, remove = FALSE) %>%
    dplyr::select(PW, taxonomy, pathway)

  # import pw data from aldex2 differential analysis
  PW.list <- list()
  for (i in names(ax2_list)) {
    comp_pair <- unlist(stringr::str_split(i, pattern = "_vs_"))
    PW.list[[i]] <- ax2_list[[i]][["PW"]] %>%
      # dplyr::filter(abs(effect) > min_eff_size & we.eBH < max_adj_p) %>%
      dplyr::mutate(
        treatment_1 = comp_pair[1],
        treatment_2 = comp_pair[2],
        comparison = i
      )
  }
  # complete missing data sets with non-significant entries
  PW.df <- dplyr::bind_rows(PW.list) %>%
    # use factors to define order similar to original
    dplyr::mutate(comparison = factor(comparison,
                               levels = names(ax2_list))) %>%
    tidyr::complete(PW, comparison, fill = list(
      diff.btw = 0,
      diff.win = 0,
      effect = 0,
      we.ep = 1,
      we.eBH = 1,
      wi.ep = 1,
      wi.eBH = 1
    ))

  # identify completed data sets
  compl_data <- PW.df %>%
    dplyr::anti_join(dplyr::bind_rows(PW.list), by = c("comparison", "PW")) %>%
    dplyr::mutate(Comparison = stringr::str_replace(comparison, "_vs_", " vs. ")) %>%
    dplyr::select(Comparison, Pathway = PW)

  # show only data with an effect size > min_eff_size | < -min_eff_size AND
  # having at least in one comparison an adj. p-value of less than 'max_adj_p'
  pw.present <- pull(unique(PW.df[abs(PW.df$effect) > min_eff_size & PW.df$we.eBH < max_adj_p, "PW"]))

  if(length(pw.present) != 0) {
  # filter hierarchical PW based on data
  hier <- hier %>%
    dplyr::filter(PW %in% pw.present)
  # hier <- hier %>%
  #   filter(PW %in% unique(bind_rows(PW.list)[,"PW"]))


  # import hierarchy in metacoder
  obj <- metacoder::parse_tax_data(hier,
                        class_cols = "taxonomy", # The column in the input table
                        class_sep = ";", # What each taxon is separated by
                        class_regex = "([a-z0-9]{1,3}__)?([^;]+)", #"^([a-z0-9]{1,3})__(.*)$",
                        class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name")
  )

  # extract arbitrary taxon ids (including supertaxa) and
  # taxonomy ('feature') from metacoder environment
  ref <- tibble(
    taxon_id = names(metacoder::classifications(obj)),
    feature_path = metacoder::classifications(obj),
    PW = stringr::str_match(feature_path, ".*;(.*)$")[,2]
  ) %>%
    dplyr::select(taxon_id, PW)

  # convert to list
  PW.list <- base::split(PW.df, f = PW.df$comparison)
  names(PW.list) <- names(ax2_list)

  # combine taxon and supertaxon ids with aldex2 data
  PW.list <- lapply(PW.list, function(x) {
    y <- x %>%
      dplyr::mutate(comparison = as.character(comparison)) %>%
      dplyr::right_join(ref, by = "PW") %>%
      tidyr::fill(treatment_1, treatment_2, comparison) %>%
      dplyr::arrange(taxon_id) %>%
      dplyr::select(taxon_id, treatment_1, treatment_2, comparison, rab.all:description) #%>%
    #droplevels()
  })

  # add aldex2 results as 'diff_table' together with number of 'taxa' / pathways
  # extract number of pathways
  hier_units <- tibble(
    taxon_id =  names(n_obs(obj, data = "tax_data")),
    hier_number = as.numeric(n_obs(obj, data = "tax_data"))
  )

  # join with number of pathways and set non-significant hits to 0 (effect) and 1 (we.eBH)
  obj$data[["diff_table"]] <- bind_rows(PW.list) %>%
    dplyr::left_join(hier_units, by = "taxon_id") %>%
    dplyr::mutate(effect = ifelse(abs(effect) < min_eff_size | we.eBH >= max_adj_p, 0, effect),
           we.eBH = ifelse(we.eBH < max_adj_p, we.eBH, 1))

  # go through multiple rounds until all taxon_ids have a value
  # add missing data to supertaxa
  for (i in names(ax2_list)) {
    # sub_df <- dplyr::filter(metacoder:::get_taxmap_table(obj, "diff_table"), comparison == i)
    while(any(is.na(dplyr::filter(metacoder:::get_taxmap_table(obj, "diff_table"), comparison == i)[["effect"]]))) {
      sub1 <- obj %>%
        metacoder::filter_obs(data = "diff_table",
                   comparison == i) %>%
        metacoder::subtaxa(value = "effect", recursive = TRUE) %>%
        base::sapply(mean)
      sub2 <- tibble(
        taxon_id = names(sub1),
        value = sub1,
        comparison = i
      )
      sub2[sub2 == "NaN"] <- NA
      obj$data[["diff_table"]] <- obj$data[["diff_table"]] %>%
        dplyr::left_join(sub2, by = c("taxon_id", "comparison")) %>%
        dplyr::mutate(effect = coalesce(effect, value)) %>%
        dplyr::select(-value)
    }
  }

  # save parameters as attributes
  attr(obj$data[["diff_table"]], "sig_pw") <- length(pw.present)
  if(nrow(compl_data) > 0) {
    warning(paste(nrow(compl_data), "data sets were completed due to missing data!"))
    attr(obj$data[["diff_table"]], "completed") <- compl_data}

  # return object
  return(obj)
  } else {
    warning("No significant hits!")
    return(NULL)
  }
}

# function to import functional abundance data as an igraph object
# import list of PICRUSt2/ALDEx2 analysis outputs (from pathway/PW analyses)
# as a node/edge data frames compatible to generate an igraph graph object
# and return the final 'igraph' graph object
import_aldex2_igraph <- function(ax2_list, pw_hierarchy, min_eff_size, max_adj_p) {
  require(metapod)
  require(igraph)
  # format hierarchical pathway data
  hier <- pw_hierarchy %>%
    tibble::rownames_to_column(var = "pathway2") %>%
    dplyr::mutate(
           origin = "origin",
           SC2 = ifelse(Superclass1 == Superclass2, NA, Superclass2),
           SC2 = ifelse(is.na(SC2), pathway2, SC2),
           SC2 = dplyr::case_when(
             Superclass1 == "Generation of Precursor Metabolites and Energy" & SC2 == "Superpathways" ~ "Superpathways Metabolism and Energy",
             Superclass1 == "Biosynthesis" & SC2 == "Superpathways" ~ "Superpathways Biosynthesis",
             TRUE ~ SC2),
           PW = ifelse(pathway2 == SC2, NA, pathway2)
    ) %>%
    dplyr::select(origin, Superclass1, Superclass2 = SC2, Pathway = PW, path_ori = pathway2)

  # import pw data from aldex2 differential analysis
  PW.list <- stats::setNames(
    vector("list", length(ax2_list)),
    names(ax2_list))
  for (i in names(PW.list)) {
    PW.list[[i]] <- ax2_list[[i]][["PW"]]}

  # filter metacyc pathways
  # show only data with an effect size > min_eff_size | < -min_eff_size AND
  # having at least in one comparison an adj. p-value of less than 'max_adj_p'
  PW.df <- dplyr::bind_rows(PW.list) %>%
    dplyr::filter(abs(effect) > min_eff_size & we.eBH < max_adj_p)
  if(nrow(PW.df) != 0) {
    pw.present <- unique(PW.df[, "PW"])
    # filter hierarchy
    hier <- hier %>%
      dplyr::filter(path_ori %in% pw.present)

    # prepare edge df from hierarchy
    edge_data <- dplyr::bind_rows(
      # test2 <- bind_rows(
      dplyr::distinct(hier, origin, Superclass1) %>%
        dplyr::select(from = origin, to = Superclass1) %>%
        arrange(to),
      dplyr::distinct(hier, Superclass1, Superclass2) %>%
        dplyr::select(from = Superclass1, to = Superclass2) %>%
        dplyr::arrange(from, to),
      dplyr::distinct(hier, Superclass2, Pathway) %>%
        dplyr::select(from = Superclass2, to = Pathway) %>%
        dplyr::arrange(from, to)
    ) %>%
      dplyr::filter(!is.na(to))

    # determine maximum effect in all significant hits
    max_eff <- max(abs(PW.df$effect))
    # determine minimum adj. p-value in all data sets
    min_adj_p <- min(PW.df$we.eBH)

    # list for storing results
    ggraph.list <- stats::setNames(
      vector("list", length(ax2.results)),
      names(ax2.results))

    # loop through all conditions
    for (j in names(PW.list)) {
      # j <- "Cre_DSS_vs_WT_DSS"
      # join with hierarchy and set non-significant hits to 0 (effect) and 1 (we.eBH)
      hier2 <- dplyr::left_join(hier, PW.list[[1]], by = c("path_ori" = "PW")) %>%
        dplyr::mutate(effect = ifelse(abs(effect) < min_eff_size | we.eBH >= max_adj_p, 0, effect),
               we.eBH = ifelse(we.eBH < max_adj_p, we.eBH, 1))

      # prepare vertices df including all metadata:
      # summarized effects (mean) and p-values using 'metapod' package
      summarize_p <- function(x, grouping) {
        y <- metapod::groupedHolmMin(p.values = x, grouping = grouping)
        return(y$p.value)
      }
      vertices_data <- dplyr::bind_rows(
        hier2 %>% dplyr::group_by(origin) %>% dplyr::summarise(effect = mean(effect, na.rm = TRUE),
                                                 p.adj = summarize_p(we.eBH, origin)) %>%
          dplyr::rename(node = origin),
        hier2 %>% dplyr::group_by(Superclass1) %>% dplyr::summarise(effect = mean(effect, na.rm = TRUE),
                                                      p.adj = summarize_p(we.eBH, Superclass1)) %>%
          dplyr::rename(node = Superclass1),
        hier2 %>% dplyr::group_by(Superclass2) %>% dplyr::summarise(effect = mean(effect, na.rm = TRUE),
                                                      p.adj = summarize_p(we.eBH, Superclass2)) %>%
          dplyr::rename(node = Superclass2),
        hier2 %>% dplyr::filter(!is.na(Pathway)) %>% dplyr::select(node = Pathway, effect, p.adj = we.eBH)
      ) %>%
        dplyr::left_join(dplyr::select(hier2, path_ori, description), by = c("node" = "path_ori")) %>%
        dplyr::mutate(
          label = ifelse(node == "origin", NA, node),
          # adjust p-value to n.s. (1) for abs(effect) < min_eff_size
          p.adj = ifelse(abs(effect) < min_eff_size, 1, p.adj),
          logBH = -log10(p.adj),
          description = ifelse(is.na(description), label, description)
        ) %>%
        dplyr::select(node, label, description, effect, p.adj, logBH)

      # add grouping variable for 'leafs'
      # identify leaf entries
      myleaves = base::which(is.na(base::match(vertices_data$label, edge_data$from) ))[-1]

      # add group variable to metadata
      vertices_data[myleaves, "group"] <- hier[
        base::match(dplyr::pull(vertices_data[myleaves, "label"]), hier$Pathway), "Superclass2"]

      # define group colors
      sgroup2 <- stats::na.omit(unique(vertices_data$group))
      require(Polychrome)
      if (length(sgroup2) <= 36) {
        group_colors <- Polychrome::palette36.colors(length(sgroup2))
      } else {
        set.seed(567629)
        group_colors <- base::rev(Polychrome::createPalette(length(sgroup2), c("#5A5156", "#E4E1E3", "#F6222E")))
      }
      names(group_colors) <- sgroup2

      # generate igraph graph
      ggraph.list[[j]] <- igraph::graph_from_data_frame(edge_data, vertices = vertices_data)
      attr(ggraph.list[[j]], "group_colors") <- group_colors
    }

    # save parameters as attributes
    attr(ggraph.list,  "sig_pw") <- length(pw.present)
    attr(ggraph.list,  "max_eff") <- max_eff
    attr(ggraph.list,  "min_adj_p") <- min_adj_p

    # return list
    return(ggraph.list)
  } else {
    warning("No significant hits!")
    return(NULL)
  }
}
