#### script to provide functions for plotting #####

#### resize / fix size of ggplots ####

# function to set ggplot panel size to fixed values (in cm)
# modified from here: https://gist.github.com/infotroph/a0a7fcf034e70de8ed75 and
# https://stackoverflow.com/questions/32580946/setting-absolute-size-of-facets-in-ggplot2
# define panel size with 'p.width' and 'p.height',
# 'p.margin' is margin arround plot (on canvas)

# p = p1
# p.margin = 20
# p.width = 12
# p.height = NULL

fix_ggplot_panel <- function(
    p,
    p.margin = 5,
    p.width = NULL,
    p.height = NULL) {

  require(ggplot2)
  require(grid)
  require(ggplotify)

  # convert ggplot to grob
  if("ggplot" %in% class(p)){
    g <- ggplot2::ggplotGrob(p)
  }else if ("gtable" %in% class(p)){
    g <- p
  }else{
    stop(paste(
      "Don't know how to get sizes for object of class",
      class(p)))
  }

  # check and extract fixed aspect ratio
  if(g$respect) {
    if (!is.null(p$coordinates$ratio)) {
      asp <- p$coordinates$ratio
    } else {
      suppressWarnings(
        asp <-
          as.numeric(str_remove_all(as.character(g$heights)[str_detect(as.character(g$heights), "null")], "null")) /
          as.numeric(str_remove_all(as.character(g$widths)[str_detect(as.character(g$widths), "null")], "null")))
    }
  } else {asp <- NULL}

  # calculate plot width and height depending on aspect ratio
  if (is.null(p.width) & is.null(p.height)) {
    stop("please provide width and / or height") } else
      if (is.null(asp) & ((is.null(p.width) & !is.null(p.height)) |
                          (!is.null(p.width) & is.null(p.height)))) {
        message("no fixed aspect ratio set and only height or width given: set aspect ratio to 1")
        p.height <- ifelse(is.null(p.height), p.width, p.height)
        p.width <- ifelse(is.null(p.width), p.height, p.width) } else
          if (!is.null(asp) & !is.null(p.width) & is.null(p.height)) {
            p.height <- p.width * asp
            message(paste("plot height with fixed aspect ratio of", asp, "set to", p.height))
          } else
            if (!is.null(asp) & is.null(p.width) & !is.null(p.height)) {
              p.width <- p.height / asp
              message(paste("plot width with fixed aspect ratio of", asp, "set to", p.width))
            } else
              if (!is.null(asp) & !is.null(p.width) & !is.null(p.height)) {
                message(paste("fixed plot aspect ratio of", asp,
                              "overridden with width set to", p.width,
                              "and height set to", p.height,
                              "resulting in new aspect ratio of", p.height/p.width))
              }

  # transform panel sizes in units
  pmargin <- unit(p.margin, "mm")
  pwidth <- unit(p.width, "cm")
  pheight <- unit(p.height, "cm")

  # extract number of vertical and horizontal panels
  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)

  # fix panel sizes
  g$widths[panel_index_w] <-  rep(pwidth,  nw)
  g$heights[panel_index_h] <- rep(pheight, nh)

  # calculate plot dimensions
  p_width <- grid::convertWidth(sum(g$widths) + pmargin,
                                unitTo = "in", valueOnly = TRUE)

  p_height <- grid::convertHeight(sum(g$heights) + pmargin,
                                  unitTo = "in", valueOnly = TRUE)
  # convert to ggplot
  g <- ggplotify::as.ggplot(g)

  # save dimensions as attributes to plot:
  attr(g, "plot_width") <- p_width
  attr(g, "plot_height") <- p_height
  return(g)
}

# p3 <- fix_ggplot_panel(p1, p.height = 8)


# function to shift legend in empty grop of faceted plot
# https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
# and fix panel size
# https://stackoverflow.com/questions/32580946/setting-absolute-size-of-facets-in-ggplot2
# p <- p2
facet_shift_legend_fix_panel <- function(p,
                                         p.margin = 5,
                                         p.width = 10,
                                         p.height = 7) {
  require(ggplot2)
  require(cowplot)
  require(gtable)
  require(purrr)
  require(lemon)

  # transform panel sizes in units
  pmargin <- unit(p.margin, "mm")
  pwidth <- unit(p.width, "cm")
  pheight <- unit(p.height, "cm")

  # identify name of empty panel
  pnls <- cowplot::plot_to_gtable(p) %>% gtable::gtable_filter("panel") %>%
    with(setNames(grobs, layout$name)) %>% purrr::keep(~identical(.x,zeroGrob()))

  if( length(pnls) == 0 ) {g <- ggplot2::ggplotGrob(p)} else {

    # if present shift legend to empty panel
    g <- lemon::reposition_legend( p, "center", panel = names(pnls) )
  }

  # extract number of vertical and horizontal panels
  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)

  # fix panel sizes
  g$widths[panel_index_w] <-  rep(pwidth,  nw)
  g$heights[panel_index_h] <- rep(pheight, nh)

  # calculate plot dimensions
  p_width <- grid::convertWidth(sum(g$widths) + pmargin,
                                unitTo = "in", valueOnly = TRUE)

  p_height <- grid::convertHeight(sum(g$heights) + pmargin,
                                  unitTo = "in", valueOnly = TRUE)
  # convert to ggplot
  g <- ggplotify::as.ggplot(g)

  # save dimensions as attributes to plot:
  attr(g, "plot_width") <- p_width
  attr(g, "plot_height") <- p_height

  return(g)
}

#### function for box/scatter plots (Alpha diversity, BugBase) ####
# alpha_df <- alpha_list_otu[["no_out"]]
# low_replicates <- filtered_levels(before = alpha_df, group_param = "geno2_diet2_week", #"geno2_diet2_sex"
#                                   min.repl = 3)

# df.. = alpha_df
# group_param.. = "group" #"geno2_diet2_sex"
# group_name.. = "ASD vs. Control"
# measure_param = "ACE"
# which_stat = "Wilcoxon"
# comp_list.. = "all" #comp_list2
# p_adj_method.. = "fdr"
# rm_ns.. = TRUE
# rel_step.. = 6
# label_set.. = c(
#     "Control" = "control", "ASD" = "ASD"
#   )
# label_set_b.. = c("Control" = "control", "ASD" = "ASD")
# color_set.. = c("Control" = "chartreuse3", "ASD" = "darkorchid4")
# do_grouped_pt.. = FALSE
# dodge_value.. = 0.9
# color_param.. = NULL # "sex" #
# data_param.. = NULL #"geno2_diet2" #
# legend_rows = 3
# cutoff_pval = 0.05


plot_box_scatter <- function(df.., group_param.., group_name.., measure_param,
                             which_stat, comp_list.., p_adj_method..,
                             rm_ns.., rel_step.., label_set.., label_set_b..,
                             color_set.., do_grouped_pt.., dodge_value..,
                             color_param.., data_param.., legend_rows)
{

  # filter data set according to BugBase variable
  df <- dplyr::filter(df.., measure == measure_param)

  # remove samples with low number of replicates
  if(!is.null(low_replicates[[1]])) {
    df1 <- dplyr::filter(df, !get(group_param..) %in% low_replicates[[1]])
  } else {df1 <- df}

  # determine normal distribution and calculate pairwise comparisons
  if (!is.null(comp_list..)) {
    if (which_stat == "Wilcoxon") {
      stat.test <-
        broom::tidy(pairwise.wilcox.test(df1$value, dplyr::pull(df1[, group_param..]))) %>%
        adjust_pvalue(method = p_adj_method..) %>%
        add_significance(p.col = "p.value.adj",
                         cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "ns"))
    }
    if (which_stat == "t-test") {
      stat.test <-
        broom::tidy(pairwise.t.test(df1$value, dplyr::pull(df1[, group_param..]))) %>%
        adjust_pvalue(method = p_adj_method..) %>%
        add_significance(p.col = "p.value.adj",
                         cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "ns"))
    }
    # remove non-desired comparisons
    if(typeof(comp_list..) == "list") {
      # take care of switched groups for comparisons
      stat.list <- vector("list", nrow(stat.test)*2)
      for (j in 1:nrow(stat.test)) {
        stat.list[[j]] <- c(dplyr::pull(stat.test[j, 1]), dplyr::pull(stat.test[j, 2]))
        stat.list[[j+nrow(stat.test)]] <- c(dplyr::pull(stat.test[j, 2]), dplyr::pull(stat.test[j, 1]))
      }
      stat.test.x2 <- dplyr::bind_rows(stat.test, stat.test)
      # remove comparisons not in 'comp_list..' list
      stat.test <- stat.test.x2[stat.list %in% comp_list.., ]
      } else if(typeof(comp_list..) == "character" & comp_list.. != "all") {
          stop("please use the argument 'all' for all possible pairwise comparisons")
        }

    # optionally remove non-significant comparisons
    if (rm_ns..) {stat.test <- dplyr::filter(stat.test, p.value.adj < cutoff_pval)}

    # extract factor levels of variable to analyze
    param_levels <- base::levels(dplyr::pull(df[, group_param..]))

    # calcuate x positions for normal plots
    if (!do_grouped_pt.. & nrow(stat.test) > 0) {
      # helper to properly assign x position with custom factor levels
      # rstatix::add_xy_position does not consider custom factor levels and
      # positions are just in alphabetical order
      assign_x <- structure(seq(1,length(param_levels)),
                            names = param_levels)

      # add x coordinates to stat.test df
      stat.test <- rstatix::add_xy_position(
        data = df, test = stat.test,
        formula = as.formula(paste("value", group_param.., sep = " ~ ")),
        x = group_param..) %>%
        dplyr::mutate(xmin = assign_x[group1],
               xmax = assign_x[group2])
    }

    # calculate x positions for grouped plots
    if (do_grouped_pt.. & nrow(stat.test) > 0) {
      stat.test2 <- data.frame(
        gp = param_levels,
        x = rstatix:::get_grouped_x_position(df, x = data_param..,
                                             group = group_param..,
                                             dodge = dodge_value..))

      stat.test <- stat.test %>%
        dplyr::left_join(stat.test2, by = c("group1" = "gp")) %>%
        dplyr::rename(xmin = x) %>%
        dplyr::left_join(stat.test2, by = c("group2" = "gp")) %>%
        dplyr::rename(xmax = x) %>%
        dplyr::mutate(delta = abs(round(xmax - xmin, 2)))
    }

    # set min y position and calculate further y positions
    if (nrow(stat.test) > 0) {
      step_size <- (max(df$value)/100)*rel_step..
      min.y <- max(df$value) + step_size/2
      y.pos.vec <- seq(min.y, min.y + 15*step_size, by = step_size)

      if (!do_grouped_pt..) {
        stat.test$y.position <- y.pos.vec[1:nrow(stat.test)]
      }
      if (do_grouped_pt..) {
        stat.test <- stat.test %>%
          dplyr::mutate(y.position = ifelse(delta == dodge_value../length(pull(unique(df[, color_param..]))), min.y, NA),
                        row_count = cumsum(is.na(y.position)) + 1,
                        y.position = ifelse(is.na(y.position), y.pos.vec[row_count], y.position))
      }
    }
  }
  # make dummy 'stat.test' in case of no stats were done
  if(is.null(comp_list..)) stat.test <- tibble()
  # plotting function
  pt <- ggplot(
    df,
    aes(x = if(do_grouped_pt..) get(data_param..) else get(group_param..), y = value,
        color = if(!is.null(color_set..) & !do_grouped_pt..) get(group_param..) else NULL)) + {

          if(do_grouped_pt..) {
            geom_boxplot(aes(fill = get(color_param..)),
                         position = position_dodge(width = dodge_value..),
                         outlier.shape = NA, coef = 0, fatten = 3, width = 0.6,
                         alpha = if(is.null(color_set..)) 0 else 0.5)
          }} + {
            if(do_grouped_pt..) {
              geom_point(aes(color = get(color_param..),
                             shape = if(is.null(color_set..)) get(color_param..)),
                         size = 1.5, stroke = 1, alpha = 0.5,
                         position = position_jitterdodge(dodge.width = dodge_value.., jitter.width = 0.2))
            }} + {
              if(!do_grouped_pt..) {
                geom_boxplot(outlier.shape = NA, coef = 0, fatten = 3, width = 0.6,
                             alpha = if(is.null(color_set..)) 0 else 0.5)
              }} + {
                if(!do_grouped_pt..) {
                  geom_jitter(aes(shape = if(is.null(color_set..)) get(color_param..) else factor(1)),
                              size = 1.5, shape = 1, stroke = 1, width = 0.2, alpha = 0.5)
                }} + {

                  if(!is.null(comp_list..) & nrow(stat.test) > 0) {stat_pvalue_manual(
                    data = stat.test, label = "p.value.adj.signif",
                    tip.length = 0.01, size = 4, bracket.size = 0.3)
                  }} +
    scale_fill_manual(
      values = if(!is.null(color_set..)) color_set.. else if(do_grouped_pt..) {
        rep("black", length(unique(pull(df[, color_param..]))))
      } else if (!do_grouped_pt..) "black",
      name = group_name.., labels = if(do_grouped_pt..) label_set_b.. else label_set..,
      aesthetics = c("color", "fill")
    ) + {
      if(!is.null(color_set..)) {guides(
        color = guide_legend(title.position = "top", nrow = legend_rows, byrow = TRUE,
                             override.aes = aes(label = "")),
        shape = "none")
      }} + {
        if(is.null(color_set..)) {guides(
          colour = "none", fill = "none",
          shape = guide_legend(title = group_name..,
                               title.position = "top", nrow = legend_rows, byrow = TRUE,
                               override.aes = aes(label = "")))
        }} +
    scale_y_continuous(name = measure_param, limits = c(0, NA),
                       expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(name = "", labels = label_set..) +
    ggtitle(measure.names[measure_param]) +
    theme_article(base_size = 16) +
    theme(
      legend.position = "none",
      legend.text = element_markdown(),
      plot.title = element_markdown(hjust = 0.5, size = 12),
      plot.subtitle = element_markdown(hjust = 0.5),
      axis.title.y = element_blank(),
      axis.text.x = if(do_grouped_pt.. | (is.null(color_set..) & !do_grouped_pt..)) {
        element_markdown(angle = 45, vjust = 1, hjust = 1)
      } else element_blank())
  return(pt)
}

#### function for ordination plots ####
# plot ordination with phyloseq function
# phylo_obj <- transform_list[[trans_meth]]
# ord_obj <- ord.all[[trans_meth]][[paste(ord_meth, dist_meth, sep = ".")]]
# group_param <- var_test
# title = plot_title
# type <- "other"
# fix_asp_1 <- FALSE

# general plot settings (also used by other plot functions..)
# ord_plot_settings <- list(
#   theme_article(base_size = 12),
#   scale_color_manual(values = params$color_set, name = params$var_name,
#                      labels = params$label_set, aesthetics = c("colour", "fill")),
#   theme(plot.title = element_text(hjust = 0.5, size = 12),
#         legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.text = element_markdown())
# )


# function to plot ordinations
plot_ord <- function(phylo_obj, ord_obj, group_param, title,
                     type = c("other", "tSNE", "UMAP"), fix_asp_1 = TRUE) {

  # extract unique categorical groups
  unique_sample_data <- dplyr::pull(unique(phyloseq::sample_data(phylo_obj)[, group_param]))

  # set aspect ratio of plots
  if (type == "other") {
    ord.evals <- if (!is.null(ord_obj[["CCA"]][["eig"]]))
      ord_obj[["CCA"]][["eig"]][c(1,2)] else
        if (is.null(ord_obj$values$Eigenvalues)) c(1,1)
    else ord_obj$values$Eigenvalues[c(1,2)]
  } else ord.evals <- c(1,1)

  if (length(color_set) >= length(unique_sample_data)) {

    # use phyloseq plot_ordination function
    if (type == "other") {
      pt <- phyloseq::plot_ordination(physeq = phylo_obj, ordination = ord_obj,
                                      type = "samples", color = group_param)
    } else if(type == "tSNE") {
      df <- base::as.data.frame(ord_obj[["Y"]])
    } else if(type == "UMAP") {
      df <- base::as.data.frame(ord_obj[["layout"]])
    }

    # plot tSNE/UMAPs
    if (type %in% c("tSNE", "UMAP")) {
      pt <- df %>% tibble::rownames_to_column(var = "sample_ID") %>%
        dplyr::left_join(data.frame(phyloseq::sample_data(phylo_obj)), by = "sample_ID") %>%
        ggplot(aes(x = V1, y = V2, color = !! ggplot2::sym(group_param))) +
        {if (type == "tSNE") labs(x = "tSNE1", y = "tSNE2", title = title)} +
        {if (type == "UMAP") labs(x = "UMAP1", y = "UMAP2", title = title)}
    }

    # apply general plot settings
    pt <- pt +
      stat_ellipse(aes(fill = !! ggplot2::sym(group_param)), geom = "polygon", alpha = 0.1) +
      {if (!fix_asp_1) coord_fixed(ratio = sqrt(ord.evals[2] / ord.evals[1]))} +
      {if (is.null(shape_var)) geom_point(size = 3)} +
      {if (!is.null(shape_var)) geom_point(aes(shape = !! ggplot2::sym(shape_var)),
                                           size = 3)} +
      ggtitle(title) +
      ord_plot_settings <- list(
        egg::theme_article(base_size = 12),
        ggplot2::scale_color_manual(values = color_set, name = var_name,
                                    labels = label_set, aesthetics = c("colour", "fill")),
        ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 12),
                       legend.position = "bottom",
                       legend.title = element_blank(),
                       legend.text = element_markdown())
      ) +
      {if (!is.null(shape_var)) theme(legend.box = "vertical")} +
      {if (fix_asp_1) theme(aspect.ratio = 1)}
  } else {
    stop(paste("only", length(color_set),"colours defined in 'color_set'\nbut",
               length(unique_sample_data), "colours required"))
  }
}

# define named vectors for figure titles and legends
trans_names <- structure(c("original", "relative-transformed",
                           "log-transformed", "Hellinger-transformed",
                           "relative log expression (RLE)-transformed",
                           "variance stabilizing-transformed (vst)",
                           "mod. centered log-ratio (CLR)-transformed",
                           "alt. centered log-ratio (CLR2)-transformed"
),
names = c("ori", "ra", "log", "hell", "RLE", "vst",
          "CLR", "CLR2"))

dist_names <- structure(c("UniFrac", "weighted UniFrac",
                          "variance-adjusted weighted UniFrac",
                          "generalized UniFrac", "DPCoA", "Bray-Curtis",
                          "Gower’s", "Jensen-Shannon divergence"),
                        names = c("uni", "wuni", "vawuni", "guni", "dpcoa",
                                  "bray", "gower", "jsd"))


#### modified ALDEx2 stats plotting functions (based on ggplot2) ####
# all the plots accept a combined df and use faceting

# ALDEx2 effect histogram plots
hist.aldex.plot <- function(x, cutoff.effect = cutoff_aldex, theme.base = 8) {
  require(egg)
  require(ggtext)
  colnum <- if(length(unique(x[["comparison"]])) < 3) 2 else if(length(x) < 16) 3 else 4

  ggplot(x, aes(x = effect)) +
    geom_histogram(bins = 30, fill = "yellow", color = "black", linewidth = 0.4) +
    geom_vline(xintercept = cutoff.effect, linetype = 2, linewidth = 0.4,
               color = "darkgrey") +
    geom_vline(xintercept = -cutoff.effect, linetype = 2, linewidth = 0.4,
               color = "darkgrey") +
    labs(x = "Effect size", y = "Frequency") +
    egg::theme_article(base_size = theme.base) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.title = element_markdown(hjust = 0.5)) +
    facet_wrap(vars(comparison), ncol = colnum)
}

# Create MW and MA plots
# modified ALDEx2 plot function to use ggplot2 instead of base plot
aldex.plot2 <- function(x, ..., type = c("MW", "MA"), test = "welch",
                        rare = 0, cutoff.pval = cutoff_pval,
                        cutoff.effect = cutoff_aldex, theme.base = 8,
                        xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
                        all.col = rgb(0, 0, 0, 0.2), all.pch = 19, all.cex = 0.4,
                        called.col = "red", called.pch = 20, called.cex = 0.6,
                        rare.col = "blue", rare.pch = 20, rare.cex = 0.2,
                        thres.line.col = "darkgrey", thres.lwd = 0.4
)
{
  require(egg)
  require(dplyr)
  require(ggtext)
  colnum <- if(length(unique(x[["comparison"]])) < 3) 2 else if(length(x) < 16) 3 else 4

  type <- match.arg(type)
  if (length(x$effect) == 0)
    stop("Please run aldex.effect before plotting")
  if (test == "welch") {
    if (length(x$we.eBH) == 0)
      stop("Welch's t test results not in dataset")
    called <- x$we.eBH <= cutoff.pval
  } else if (test == "wilcox") {
    if (length(x$wi.eBH) == 0)
      stop("Wilcoxon test results not in dataset")
    called <- x$wi.eBH <= cutoff.pval
  } else if (test == "glm") {
    if (length(x$glm.eBH) == 0)
      stop("glm test results not in dataset")
    called <- x$glm.eBH <= cutoff.pval
  } else if (test == "kruskal") {
    if (length(x$kw.eBH) == 0)
      stop("Kruskall-Wallace test results not in dataset")
    called <- x$kw.eBH <= cutoff.pval
  } else if (test == "effect") {
    if (cutoff.effect <= 0.49)
      stop("Please set cutoff to at least 0.5")
    called <- abs(x$effect) >= cutoff.effect
  } else if (test == "both") {
    if (cutoff.effect <= 0.49)
      stop("Please set cutoff to at least 0.5")
    called <- abs(x$effect) >= cutoff.effect & x$wi.eBH <= cutoff.pval
  }

  # df for geom_text labels in facets
  label_data <- dplyr::distinct(x, comparison, ref, test)

  # general plot settings
  plot_settings <- list(
    geom_point(size = 0.6),
    scale_shape_manual(values = c("all" = all.pch, "rare" = rare.pch, "called" = called.pch)),
    scale_color_manual(values = c("all" = all.col, "rare" = rare.col, "called" = called.col)),
    scale_alpha_manual(values = c("all" = 0.3, "rare" = 0.3, "called" = 0.3)),
    scale_size_manual(values = c("all" = all.cex, "rare" = rare.cex, "called" = called.cex)),
    egg::theme_article(base_size = theme.base),
    theme(legend.position = "none",
          legend.key = element_rect(fill = NA),
          panel.grid = element_blank(),
          plot.title = element_markdown(hjust = 0.5)),
    # coord_fixed(),
    facet_wrap(vars(comparison), ncol = colnum)
  )

  # add data classification
  y <- x %>%
    mutate(x, class = case_when(rab.all < rare ~ "rare",
                                TRUE ~ "all"),
           class = ifelse(called, "called", class))

  if (type == "MW") {
    if (is.null(xlab))
      xlab <- expression("Median" ~ ~Log[2] ~ ~"dispersion")
    if (is.null(ylab))
      ylab <- expression("Median" ~ ~Log[2] ~ ~"difference")

    pt <- ggplot(y, aes(x = diff.win, y = diff.btw, color = class, shape = class,
                        alpha = class, size = class)) +
      plot_settings +
      labs(x = xlab, y = ylab) +
      geom_abline(intercept = 0, slope = 1, color = thres.line.col,
                  linetype = 2, linewidth = thres.lwd) +
      geom_abline(intercept = 0, slope = -1, color = thres.line.col,
                  linetype = 2, linewidth = thres.lwd) +
      geom_text(data = label_data, aes(label = ref, x = 0, y = min(x$diff.btw)),
                color = "grey", size = 2, hjust = 0, inherit.aes = FALSE) +
      geom_text(data = label_data, aes(label = test, x = 0, y = max(x$diff.btw)),
                color = "grey", size = 2, hjust = 0, inherit.aes = FALSE)
    return(pt)
  }
  if (type == "MA") {
    if (is.null(xlab))
      xlab <- expression("Median" ~ ~Log[2] ~ ~"relative abundance")
    if (is.null(ylab))
      ylab <- expression("Median" ~ ~Log[2] ~ ~"difference")

    pt <- ggplot(y, aes(x = rab.all, y = diff.btw, color = class, shape = class,
                        alpha = class, size = class)) +
      plot_settings +
      labs(x = xlab, y = ylab) +
      geom_text(data = label_data, aes(label = ref, x = max(x$rab.all), y = min(x$diff.btw)),
                color = "grey", size = 2, hjust = 1, inherit.aes = FALSE) +
      geom_text(data = label_data, aes(label = test, x = max(x$rab.all), y = max(x$diff.btw)),
                color = "grey", size = 2, hjust = 1, inherit.aes = FALSE)
    return(pt)
  }
}

# function to plot relationship between effect, difference, and P values
# Effect size and Volcano plot
effect.volcano.plot <- function(x, type = c("EFF", "VOLC"),
                                cutoff.effect = cutoff_aldex,
                                cutoff.pval = cutoff_pval,
                                theme.base = 8) {
  require(egg)
  require(tidyr)
  require(ggtext)
  colnum <- if(length(unique(x[["comparison"]])) < 3) 2 else if(length(x) < 16) 3 else 4

  type <- match.arg(type)
  m <- tidyr::pivot_longer(data = x, cols = we.ep:we.eBH, names_to = "pval")
  plot_settings <- list(
    geom_point(size = 0.6, alpha = 0.3),
    geom_hline(yintercept = cutoff.pval, linetype = 2, linewidth = 0.4,
               color = "orange"),
    scale_color_manual(name = NULL,
                       values = c("we.ep" = "blue", "we.eBH"  = "red"),
                       labels = c("we.ep" = "p-value", "we.eBH"  = "BH-adjusted")),
    scale_y_log10(),
    egg::theme_article(base_size = theme.base),
    theme(legend.position = c(1, 0.01), # position of the legend anchoring point (relative plot width and height)
          legend.justification = c(1, 0), # anchoring point of the legend (relative legend width and height)
          panel.grid = element_blank(),
          plot.title = element_markdown(hjust = 0.5),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0,0,0,0),
          legend.background = element_blank(),
          #legend.key = element_rect(fill = NA),
          legend.key.size = unit(c(0, 0), 'cm'), #if not set to 0: still separates keys
          legend.text = element_markdown(
            margin = margin(l = 0, unit = "pt"))), # no need for negative adjustment, if legend.key.size set to 0
    facet_wrap(vars(comparison), ncol = colnum)
  )
  if (type == "EFF") {
    pt <- ggplot(m, aes(x = effect, y = value, color = pval)) +
      geom_vline(xintercept = cutoff.effect, linetype = 2, linewidth = 0.4,
                 color = "darkgrey") +
      geom_vline(xintercept = -cutoff.effect, linetype = 2, linewidth = 0.4,
                 color = "darkgrey") +
      labs(x = "Effect size", y = "p-value") +
      plot_settings
    return(pt)
  }
  if (type == "VOLC") {
    pt <- ggplot(m, aes(x = diff.btw, y = value, color = pval)) +
      labs(x = "Difference", y = "p-value") +
      plot_settings
    return(pt)
  }
}

# function to compose stat plots from list of plots
compose_plots <- function(x) {
  require(patchwork)
  comb.plots <- wrap_plots(x) +
    plot_layout(
      ncol = if(length(x) < 3) 2 else if(length(x) < 16) 3 else 4,
      # heights = unit(c(8), c('cm')))
      heights = 1)
  return(comb.plots)
}

#### matrix plot functions ####
# plot the data as matrix of metacoder 'heat trees' or ggraph dendrograms

### helper functions for all matrix plots ###

# function to generate layout matrix from vector of 'treatments' (comparisons..)
# https://stackoverflow.com/questions/76772780/create-and-populate-the-upper-triangular-matrix-with-a-vector-in-r#new-answer
create_layout_matrix <- function(treat) {
  n <- length(treat)
  if(n == 0) {
    stop("There is no data for sub-plots!")
  } else {
    # calculate dimensions of the square matrix
    m <- ceiling(sqrt(0.25 + 2*n) - 0.5)
    # create empty matrix
    x <- array(NA, rep(m, 2))
    s <- base::pmin(m + 1L - row(x), col(x)) + base::pmax(m + 1L - row(x), col(x))/n + (m + 1L)*base::upper.tri(x)
    x[base::sort(base::order(s)[1:n])] <- 1:n
    t(x) }
}

# calculate coordinates of all subplots
calc_subplot_coords2 <- function(a_matrix, treat, x1 = 0, y1 = 0, x2 = 1, y2 = 1,
                                 ks, kx, ky) {
  # lowerleft = c(x1, y1), upperright = c(x2, y2)
  x_coords <- seq(from = x1, to = x2, length.out = ncol(a_matrix) + 1)[- (ncol(a_matrix) + 1)]
  y_coords <- seq(from = y1, to = y2, length.out = nrow(a_matrix) + 1)[- (nrow(a_matrix) + 1)]
  matr <- base::do.call(base::rbind, lapply(1:ncol(a_matrix), function(x) data.frame(plot_index = a_matrix[, x],
                                                                         x = x_coords[x],
                                                                         y = rev(y_coords))))
  matr <- matr[!is.na(matr$plot_index), ]
  matr <- matr[order(matr$plot_index), ]
  if(length(treat) > 1){
    return(matr)} else if(length(treat) == 1) {
      matr[1, "x"] <- ks + kx
      matr[1, "y"] <- ks + ky
      return(matr)
    }
}

# function to calculate subplot dimensions based on layout matrix
# a_matrix = layout_matrix
# b_matrix = matrix_data
# which_dim = "height"
# ks = key_size
# kx = key_x
# ky = key_y

calc_subplot_dims <- function(a_matrix, b_matrix, which_dim = c("width", "height"),
                              ks, kx, ky) {
  if(nrow(b_matrix) > 2) {
    named_row <- base::which(base::apply(a_matrix, MARGIN = 1, function(x) all(!is.na(x))))
    named_col <- base::which(base::apply(a_matrix, MARGIN = 2, function(x) all(!is.na(x))))

    if(length(named_row) == 1 & length(named_col) == 1) {
      horz_label_data <- b_matrix[match(a_matrix[named_row, ], b_matrix$plot_index), ]
      vert_label_data <- b_matrix[match(a_matrix[, named_col], b_matrix$plot_index), ]
      if(which_dim == "width") {
        return(abs(horz_label_data$x[1] - horz_label_data$x[2]))
      } else if(which_dim == "height") {
        return(abs(vert_label_data$y[1] - vert_label_data$y[2]))
      }
    } else {
      if(which_dim == "width") {return(max(abs(diff(b_matrix$x))))
      } else if(which_dim == "height") {return(max(abs(diff(b_matrix$y))))
      }
    }
  } else {
    if(which_dim == "width") {return(1 - (ks + kx))
    } else if(which_dim == "height") {return(1 - (ks + ky))
    }
  }
}


# function to plot a matrix of heat trees from differential abundance (DA) data
# modified 'heat_tree_matrix()' function: https://github.com/grunwaldlab/metacoder/blob/master/R/heat_tree_matrix.R
# - to accept a subset of all possible comparisons defined by the comparison list of the main analysis
# this is mainly a wrapper to the 'heat_tree' function to generate sub-plots based on filtering of the main data and
# a key plot labeled with the actual taxonomy
# https://github.com/grunwaldlab/metacoder/commit/6fe54ca722baab49c92b2644c0a65feb3564f979


# obj <- mc_env_agg %>%
#   mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]|^[a-z]{1}__|__$|^.*s__", replacement = ""))
# data = "diff_summary"
# label_small_trees =  FALSE
# small_trees = "summary"
# node_color_trans = "linear"
# node_size_axis_label = "Number of OTUs"
# node_color_axis_label = "Sig. methods"
# node_color_range = colorspace::diverging_hcl(7, palette = "Green-Orange")
# tax_level = "Genus"
# key_size = 0.6
# key_x = 0.05
# key_y = 0
# matrix_x1 = 0.05
# matrix_y1 = 0.05
# matrix_x2 = 0.95
# matrix_y2 = 0.95
# seed = 1
# output_files = TRUE
# export_dir = paste0(current_dir, "/DA_analysis/figures")
# TU_type = TU_type
# useMC = FALSE
# dataset = NULL

heat_tree_matrix_DA <- function(obj, data, label_small_trees =  FALSE,
                                small_trees = c("comparisons", "summary", "methods"),
                                tax_level = "Genus",
                                key_size = 0.6,
                                key_x = 0.05,
                                key_y = 0,
                                matrix_x1 = 0.05,
                                matrix_y1 = 0.05,
                                matrix_x2 = 0.95,
                                matrix_y2 = 0.95,
                                seed = 1,
                                output_files = TRUE,
                                export_dir,
                                TU_type,
                                useMC,
                                ..., dataset = NULL) {

  # Check for use of "dataset"
  if (! is.null(dataset)) {
    warning(call. = FALSE,
            'Use of "dataset" is depreciated. Use "data" instead.')
    data <- dataset
  }

  # Get plot data table
  diff_table <- metacoder:::get_taxmap_table(obj, data)

  # plot layout for 'methods' = all methods for one comparison
  # export as list with 'comparisons' as elements
  if (small_trees == "methods") {
    categories <- unique(diff_table$comparison)
    if(is.factor(categories)) {
      categories <- levels(categories)
    }
  }

  # plot layout for 'comparisons' = all comparisons for one method
  # export as list with 'methods' as elements
  if (small_trees == "comparisons") {
    categories <- unique(diff_table$method)
    if(is.factor(categories)) {
      categories <- levels(categories)
    }
  }

  # plot layout for 'summary' = summary of all methods for all comparisons
  # export as list with single entry: 'summary'
  if (small_trees == "summary") {
    categories <- "summary"
  }

  # x <- categories
  # function to loop through all comparisons
  mc_matrix <- function(x) {

    # filter diff_table
    if (small_trees == "methods") {
      df <- dplyr::filter(diff_table, comparison == x)
      # identify individual methods
      treatments <- dplyr::pull(dplyr::distinct(df, method))
      # maximum abundances
      max_vals <- dplyr::group_by(df, method) %>%
        dplyr::summarize(max_val = max(abs(diffabund))) %>%
        dplyr::mutate(max_val = plyr::round_any(max_val, 0.1, f = ceiling))
    }

    if (small_trees == "comparisons") {
      df <- dplyr::filter(diff_table, method == x)
      # identify individual comparisons
      treatments <- dplyr::pull(dplyr::distinct(df, comparison))
      # maximum abundance
      max_val <- plyr::round_any(max(abs(df$diffabund)), 0.1, f = ceiling)
    }

    if (small_trees == "summary") {
      df <- diff_table
      # identify individual comparisons
      treatments <- dplyr::pull(dplyr::distinct(df, comparison1))
      # maximum number of methods with significant results
      max_val <- max(df$n)
    }

    # calculate layout matrix based on 'treatments'
    layout_matrix <- create_layout_matrix(treatments)

    # Get subplot layout data
    # The x1, y1, x2, y2 refer to the relative size of the cowplot canvas size
    # control position and size of sub-plots by providing the lower left (x1, y1) and upper right (x2, y2) coordinates,
    # relative to the canvas size defined above (usually 0,0 (lower left) -> 1,1 (upper right))
    matrix_data <- calc_subplot_coords2(layout_matrix, treat = treatments,
                                        x1 = matrix_x1, y1 = matrix_y1,
                                        x2 = matrix_x2, y2 = matrix_y2,
                                        ks = key_size, kx = key_x, ky = key_y)

    # right join with 'treatments' data, remove rows without entry
    if(small_trees == "comparisons") {
      matrix_data <- dplyr::right_join(matrix_data,
                                       dplyr::distinct(df, comparison) %>%
                                         tibble::rownames_to_column("plot_index") %>%
                                         dplyr::mutate(plot_index = as.numeric(plot_index)),
                                       by = "plot_index")
    }

    if(small_trees == "summary") {
      matrix_data <- dplyr::right_join(matrix_data,
                                       dplyr::distinct(df, comparison1) %>%
                                         tibble::rownames_to_column("plot_index") %>%
                                         dplyr::mutate(plot_index = as.numeric(plot_index)),
                                       by = "plot_index")
    }

    if(small_trees == "methods") {
      matrix_data <- dplyr::right_join(matrix_data,
                                       dplyr::distinct(df, method) %>%
                                         tibble::rownames_to_column("plot_index") %>%
                                         dplyr::mutate(plot_index = as.numeric(plot_index)),
                                       by = "plot_index")
    }
    rownames(matrix_data) <- matrix_data$plot_index

    # calculate subgraph dimensions
    subgraph_width <- calc_subplot_dims(layout_matrix, matrix_data, "width",
                                        ks = key_size, kx = key_x, ky = key_y)

    subgraph_height <- calc_subplot_dims(layout_matrix, matrix_data, "height",
                                         ks = key_size, kx = key_x, ky = key_y)

    # correct for small number of comparisons: smaller sub-plots!
    if (length(treatments) <= 3) {
      # make plots smaller
      subgraph_width <- subgraph_width*0.7
      subgraph_height <- subgraph_height*0.7
      # push smaller plots to outer limits
      matrix_data <- matrix_data %>%
        mutate(across(c(x, y), ~ .x*1.2))
    }

    # This odd thing is used to overwrite options without evaluation
    plot_sub_plot <- ifelse(
      label_small_trees,
      ifelse(small_trees == "methods",
             function(..., make_edge_legend = FALSE, output_file = NULL) {
               metacoder::heat_tree(..., make_edge_legend = FALSE, output_file = NULL)
             },
             function(..., make_node_legend = FALSE, make_edge_legend = FALSE,
                      output_file = NULL) {
               metacoder::heat_tree(..., make_node_legend = FALSE,
                                    make_edge_legend = FALSE, output_file = NULL)
             }),
      ifelse(small_trees == "methods",
             function(..., node_label = NULL, make_edge_legend = FALSE, output_file = NULL) {
               metacoder::heat_tree(..., make_edge_legend = FALSE, output_file = NULL)
             },
             function(..., node_label = NULL, make_node_legend = FALSE,
                      make_edge_legend = FALSE, output_file = NULL) {
               metacoder::heat_tree(..., make_node_legend = FALSE,
                                    make_edge_legend = FALSE, output_file = NULL)
             }
      ))

    # initiate empty list to store sub_plots
    sub_plots <- setNames(vector("list", length(treatments)), treatments)

    # method names for subplot titles
    method_names <- c("edger" = "EdgeR", "lv"= "limma (voom)",
                      "deseq2" = "DESeq2", "ancom" = "ANCOM",
                      "ancombc" = "ANCOM-BC", "lefse" = "LEfSe",
                      "aldex" = "ALDEx2", "mgzig" = "metag (ZIG)",
                      "mgziln" = "metag (ZILN)")

    # filter sub-plots individually and fill results list
    if (length(treatments) > 0) {
      if (small_trees == "methods") {
        for (j in treatments) {
          # j <- "ancom"
          max_val <- pull(max_vals[max_vals$method == j, "max_val"])
          set.seed(seed)
          sub_plots[[j]] <- obj %>%
            filter_obs(data,
                       comparison == x,
                       method == j) %>%
            plot_sub_plot(...,
                          node_legend_title = NULL,
                          node_color_interval = c(-max_val, max_val), # symmetric interval,
                          edge_color_interval = c(-max_val, max_val), # symmetric interval,
                          title = method_names[j],
                          title_size = 0.06)
        }
      }
      if (small_trees == "comparisons") {
        for (j in treatments) {
          set.seed(seed)
          sub_plots[[j]] <- obj %>%
            filter_obs(data,
                       method == x,
                       comparison == j) %>%
            plot_sub_plot(...,
                          node_legend_title = NULL,
                          node_color_interval = c(-max_val, max_val), # symmetric interval
                          edge_color_interval = c(-max_val, max_val), # symmetric interval
                          title = stringr::str_replace(j, "_vs_", " vs. "),
                          title_size = 0.06)
        }
      }
      if (small_trees == "summary") {
        for (j in treatments) {
          set.seed(seed)
          sub_plots[[j]] <- obj %>%
            filter_obs(data,
                       comparison1 == j) %>%
            plot_sub_plot(...,
                          node_legend_title = NULL,
                          node_color_interval = c(-max_val, max_val), # symmetric interval
                          edge_color_interval = c(-max_val, max_val), # symmetric interval
                          title = stringr::str_replace(j, "_vs_", " vs. "),
                          title_size = 0.06)
        }
      }
    } else {return(NULL)}

    # Make key plot
    plot_key_plot <- ifelse(small_trees == "methods",
                            function(..., node_color = NULL, node_color_interval = NULL) {
                              heat_tree(...,
                                        node_legend_title = NULL,
                                        node_color = "grey",
                                        node_color_interval = NULL)},
                            ifelse(small_trees == "comparisons",
                                   function(..., node_color = NULL, node_color_interval = NULL,
                                            node_color_axis_label = NULL) {
                                     heat_tree(...,
                                               node_legend_title = NULL,
                                               node_color_axis_label = method_names[i],
                                               node_color = "grey",
                                               node_color_interval = c(-max_val, max_val))}, # make gradient visible!
                                   function(..., node_color = NULL, node_color_interval = NULL,
                                            node_color_axis_label = NULL) {
                                     heat_tree(...,
                                               node_legend_title = NULL,
                                               node_color_axis_label = "sig. methods",
                                               node_color = "grey",
                                               node_color_interval = c(-max_val, max_val))} # make gradient visible!
                            ))

    set.seed(seed)
    key_plot <- plot_key_plot(obj, ...)

    # Make matrix plot
    matrix_plot <- cowplot::ggdraw() +
      cowplot::draw_plot(key_plot, x = key_x, y = key_y, width = key_size, height = key_size) +
      ggplot2::theme(aspect.ratio = 1)
    # 'Comparison' title for 'methods' small plots
    if (small_trees == "methods") {
      matrix_plot <- matrix_plot +
        cowplot::draw_text(str_replace(x, "_vs_", " vs. "), x = 0.5, y = 0.97, size = 14)
    }
    # 'Methods' title for 'comparisons' small plots
    if (small_trees == "comparisons") {
      matrix_plot <- matrix_plot +
        cowplot::draw_text(method_names[x], x = 0.5, y = 0.97, size = 14)
    }
    # 'Summary' title for 'summary' small plots
    if (small_trees == "summary") {
      matrix_plot <- matrix_plot +
        cowplot::draw_text("Cumulative significant hits across all methods", x = 0.5, y = 0.97, size = 14)
    }
    # add all small plots
    for (k in seq_along(sub_plots)) {
      if (! is.null(sub_plots[[k]])) {
        matrix_plot <- matrix_plot + cowplot::draw_plot(sub_plots[[k]],
                                                        x = matrix_data[k, "x"],
                                                        y = matrix_data[k, "y"],
                                                        width = subgraph_width,
                                                        height = subgraph_height)
      }
    }

    # export combined matrix plot
    file_name <- paste0(export_dir, "/", paste(
      ifelse(TU_type, "ZOTUs", "OTUs"), "all", var_test, x, tax_level, "metacoder.png",
      sep = "."))

    if (small_trees == "methods") {
      fig_legend <- paste0(
        "Overview of differential abundance (DA) analyses comparing factors '",
        str_replace_all(str_replace(x, "_vs_" , "' vs. '"), "_", "\\\\_"),
        "' of variable '", str_replace_all(var_test, "_", "\\\\_"),
        "'. In the lower left a key plot shows the taxonomy of the data set. Smaller 'heat trees' depict the differential abundance identified by the methods indicated.",
        if(!is.null(tax_level)) {paste0(" Data was agglomerated at the ", tax_level," level.")})
    }
    if (small_trees == "comparisons") {
      fig_legend <- paste0(
        "Overview of differential abundance (DA) analyses using the method'",
        method_names[x], "'. In the lower left a key plot shows the taxonomy of the data set. Smaller 'heat trees' depict the differential abundance of the comparisons indicated.",
        if(!is.null(tax_level)) {paste0(" Data was agglomerated at the ", tax_level," level.")})
    }
    if (small_trees == "summary") {
      fig_legend <- paste0(
        "Summary of differential abundance (DA) analyses using in total ", max_val,
        " different methods. In the lower left a key plot shows the taxonomy of the data set. ",
        "Smaller 'heat trees' depict the differential abundance of the comparisons indicated. ",
        "The color scale represents the number of methods reporting significant differential abundance for the particular taxonomy.",
        if(!is.null(tax_level)) {paste0(" Data was agglomerated at the ", tax_level," level.")})
    }
    attr(matrix_plot, "file_name") <- file_name
    attr(matrix_plot, "fig_legend") <- fig_legend

    return(list(
      key = key_plot,
      subs = sub_plots,
      matrix_pt = matrix_plot,
      matrix_data = matrix_data,
      sub_width = subgraph_width,
      sub_heigth = subgraph_height
    ))
  }

  # function to format total time difference
  total_mc_difft <- function(calc, start.t = st) {
    x <- format_difftime(first = start.t, second = Sys.time())
    return(paste0("Metacoder heat tree plotting done in **", x, "** using ", calc))
  }

  if(useMC) {
    # identify number of available cores
    max_cores <- parallel::detectCores(logical = TRUE)
    # run the analysis in parallel using different functions depending on the OS
    if(.Platform$OS.type == "windows") {
      # adjust number of clusters depending on number of comparisons
      # this is mainly to save memory due to unused, but initiated clusters
      num_cores <- if(length(categories) <= max_cores) {
        length(categories)} else if(length(categories)/2 <= max_cores) {
          ceiling(length(categories)/2)} else max_cores-1

      st <- Sys.time()
      # initiate clusters
      cl <- parallel::makeCluster(num_cores, type = "PSOCK")
      # everything created before starting clusters has to be exported
      # export required libraries to the env of each cluster
      parallel::clusterEvalQ(cl, c(library(metacoder), library(dplyr),
                                   library(plyr), library(tibble), library(stringr),
                                   library(ggplot2), library(cowplot)))
      # export required objects to the env of each cluster
      parallel::clusterExport(cl, c("create_layout_matrix", "calc_subplot_coords2",
                                    "calc_subplot_dims", "seed", "var_test",
                                    "metacoder_tax", "export_dir", "TU_type",
                                    "key_size", "key_x", "key_y"),
                              envir = environment())
      export_list <- parallel::parLapply(cl = cl, X = categories, fun = mc_matrix)
      # stop clusters
      parallel::stopCluster(cl)
      attr(export_list, "total_runtime") <- total_mc_difft("parLapply (parallel).")
    } else {
      # if on Linux / MacOS: use the more efficient 'mclapply'
      st <- Sys.time()
      export_list <- parallel::mclapply(categories, mc_matrix, mc.cores = max_cores-1)
      attr(export_list, "total_runtime") <- total_mc_difft("mclapply (parallel).")
    }
  } else {
    # use serial mode (no multicore)
    st <- Sys.time()
    export_list <- lapply(categories, mc_matrix)
    attr(export_list, "total_runtime") <- total_mc_difft("lapply (serial).")
  }

  # add names to list
  names(export_list) <- categories

  # export graphs
  if (output_files) {
    lapply(categories, function(x) {ggsave(
      filename = attr(export_list[[x]][["matrix_pt"]], "file_name"),
      plot = export_list[[x]][["matrix_pt"]],
      width = 30, height = 30, units = "cm", bg = "white"
    )})
  }

  return(export_list)
}

# function to plot a matrix of heat trees from functional abundance (FA) data
# modified 'heat_tree_matrix()' function: https://github.com/grunwaldlab/metacoder/blob/master/R/heat_tree_matrix.R
# - to accept a subset of all possible comparisons defined by the comparison list of the main analysis
# this is mainly a wrapper to the 'heat_tree' function to generate sub-plots based on filtering of the main data and
# a key plot labeled with the pathway hierarchy
heat_tree_matrix_FA <- function(obj, data,
                                label_small_trees =  FALSE,
                                key_size = 0.7,
                                key_x = 0,
                                key_y = 0,
                                matrix_x1 = 0.1,
                                matrix_y1 = 0.1,
                                matrix_x2 = 1,
                                matrix_y2 = 1,
                                seed = 1,
                                output_files = TRUE,
                                workdir = workdir,
                                subdir = subdir,
                                ..., dataset = NULL) {

  require(cowplot)

  # Check for use of "dataset"
  if (! is.null(dataset)) {
    warning(call. = FALSE,
            'Use of "dataset" is depreciated. Use "data" instead.')
    data <- dataset
  }

  # Get plot data table
  df <- metacoder:::get_taxmap_table(obj, data)

  # identify individual comparisons
  treatments <- dplyr::pull(dplyr::distinct(df, comparison))
  # maximum effect
  max_val <- plyr::round_any(max(abs(df$effect)), 0.1, f = ceiling)

  # calculate layout matrix based on 'treatments'
  layout_matrix <- create_layout_matrix(treatments)

  # Get subplot layout data
  # The x1, y1, x2, y2 refer to the relative size of the cowplot canvas size
  # control position and size of sub-plots by providing the lower left (x1, y1) and upper right (x2, y2) coordinates,
  # relative to the canvas size defined above (usually 0,0 (lower left) -> 1,1 (upper right))
  matrix_data <- calc_subplot_coords2(layout_matrix, treat = treatments,
                                      x1 = matrix_x1, y1 = matrix_y1,
                                      x2 = matrix_x2, y2 = matrix_y2,
                                      ks = key_size, kx = key_x, ky = key_y)

  # right join with 'treatments' data, remove rows without entry
  matrix_data <- dplyr::right_join(matrix_data,
                            tibble(plot_index = seq(1:length(treatments)),
                                   comparison = treatments),
                            by = "plot_index")
  rownames(matrix_data) <- matrix_data$plot_index

  # calculate subgraph dimensions
  subgraph_width <- calc_subplot_dims(layout_matrix, matrix_data, "width",
                                      ks = key_size, kx = key_x, ky = key_y)
  subgraph_height <- calc_subplot_dims(layout_matrix, matrix_data, "height",
                                       ks = key_size, kx = key_x, ky = key_y)

  # correct for small number of comparisons: smaller sub-plots!
  if (length(treatments) <= 3) {
    # make plots smaller
    subgraph_width <- subgraph_width*0.7
    subgraph_height <- subgraph_height*0.7
    # push smaller plots to outer limits
    matrix_data <- matrix_data %>%
      dplyr::mutate(across(c(x, y), ~ .x*1.2))
  }

  # This odd thing is used to overwrite options without evaluation
  plot_sub_plot <- ifelse(
    label_small_trees,
    function(..., make_node_legend = FALSE, make_edge_legend = FALSE,
             output_file = NULL) {
      metacoder::heat_tree(..., make_node_legend = FALSE,
                           make_edge_legend = FALSE, output_file = NULL)
    },
    function(..., node_label = NULL, make_node_legend = FALSE,
             make_edge_legend = FALSE, output_file = NULL) {
      metacoder::heat_tree(..., make_node_legend = FALSE,
                           make_edge_legend = FALSE, output_file = NULL)
    })

  # initiate empty list to store sub_plots
  sub_plots <- setNames(vector("list", length(treatments)), treatments)

  # filter sub-plots individually and fill results list
  if (length(treatments) > 0) {
    for (j in treatments) {
      set.seed(seed)
      sub_plots[[j]] <- obj %>%
        metacoder::filter_obs(data,
                   comparison == j) %>%
        plot_sub_plot(...,
                      node_legend_title = NULL,
                      node_color_interval = c(-max_val, max_val), # symmetric interval
                      edge_color_interval = c(-max_val, max_val), # symmetric interval
                      title = str_replace(j, "_vs_", " vs. "),
                      title_size = 0.07)
    }
  } else {return(NULL)}

  # Make key plot
  plot_key_plot <- function(..., node_color = NULL, node_color_interval = NULL) {
    metacoder::heat_tree(...,
              node_legend_title = NULL,
              node_color = "grey",
              node_color_interval = c(-max_val, max_val))} # make gradient visible!

  set.seed(seed)
  key_plot <- plot_key_plot(obj, ...)

  # Make matrix plot
  matrix_plot <- cowplot::ggdraw() +
    cowplot::draw_plot(key_plot, x = key_x, y = key_y, width = key_size, height = key_size) +
    ggplot2::theme(aspect.ratio = 1)

  # add all small plots
  for (k in seq_along(sub_plots)) {
    if (! is.null(sub_plots[[k]])) {
      matrix_plot <- matrix_plot + cowplot::draw_plot(sub_plots[[k]],
                                                      x = matrix_data[k, "x"],
                                                      y = matrix_data[k, "y"],
                                                      width = subgraph_width,
                                                      height = subgraph_height)
    }
  }

  # export combined matrix plot
  file_name <- paste(
    output_dir9, paste(
      "Pathway", var_test, "metacoder.png", sep = "."), sep = "/") %>%
    gsub("\\\\", "/", .)

  fig_legend <- paste0(
    "Overview of functional abundance (FA) analyses of MetaCyc pathways using PICRUSt2 followed by ALDEx2-based differential analyses.",
    " In the lower left a key plot shows the hierarchy of all ", attr(df, "sig_pw"),
    " significantly differential abundant pathways in the data set as a tree.",
    " Smaller 'heat trees' depict the differential functional abundance of the comparisons between factors of variable '",
    str_replace_all(var_test, "_", "\\\\_"), "' as indicated.",
    " The color scale represents the effect size of the differential pathway abundance.")

  attr(matrix_plot, "file_name") <- file_name
  attr(matrix_plot, "fig_legend") <- fig_legend

  if (output_files) {
    ggplot2::ggsave(filename = file_name,
           plot = matrix_plot,
           width = 30, height = 30, units = "cm",
           bg = "white")
  }

  # fill list with results
  return(list(key = key_plot,
              subs = sub_plots,
              matrix_pt = matrix_plot,
              matrix_data = matrix_data,
              sub_width = subgraph_width,
              sub_heigth = subgraph_height))
}

### helper functions for 'ggraph_matrix' plot ###
# function to plot a single ggraph dendrogram with functional abundance data as colored dots
metacyc_ggraph <- function(x, graph_title, min.eff.size = 1, max.eff, min.adj.p, max.adj.p = 0.05) {
  pt <-
    ggraph::ggraph(x, layout = 'dendrogram', circular = TRUE) +
    ggraph::geom_edge_diagonal(colour = "grey75", width = 0.25) +

    # color scale
    # https://groups.google.com/g/ggplot2/c/7avBhoac9bI
    # https://stackoverflow.com/questions/21758175/is-it-possible-to-define-the-mid-range-in-scale-fill-gradient2
    ggplot2::scale_fill_gradientn(
      colours = c("blue", "white", "white", "white", "red"), na.value = "white",
      values = scales::rescale(c(-max.eff, -min.eff.size, 0, min.eff.size, max.eff)),
      limits = c(-max.eff, max.eff)) +

    # set size of points, do not show points below significance level p = 0.05 (max.adj.p)
    ggplot2::scale_size_continuous(range = c(0.5, 6), limits = c(-log10(max.adj.p), -log10(min.adj.p))) +

    # all data with points
    ggraph::geom_node_point(aes(fill = effect, size = logBH),
                    color = "grey75", shape = 21, stroke = 0.25) +
    ggplot2::theme_void() +
    ggplot2::guides(fill = guide_colourbar(title = "Effect", title.position = "top", order = 1),
           size = guide_legend(title = "-log10(adj.p)", title.position = "top", nrow = 1, byrow = TRUE), order = 2) +
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "top",
      legend.box = "verical",
      legend.margin = margin(),
      #legend.key.size = unit(0.5, 'cm'),
      plot.margin = unit(c(0,0,0,0),"cm"),
      panel.spacing = unit(c(0,0,0,0), "cm")) +
    ggtitle(label = str_replace(graph_title, "_vs_", " vs. ")) +
    coord_fixed()
}

# function to plot the ggraph key dendrogram in grey with or without pathways designated for the leafs
metacyc_key_ggraph <- function(x, label_leafs, color_leafs = TRUE) {
  pt <-
    ggraph::ggraph(x, layout = 'dendrogram', circular = TRUE) +
    geom_edge_diagonal(colour = "grey75", width = 0.25) +
    geom_node_text(aes(filter = !leaf,
                       label = stringr::str_wrap(label, width = 20, whitespace_only = FALSE)),
                   size = 3, repel = TRUE, max.overlaps = 30) +
    {if(label_leafs)
      geom_node_text(aes(x = x*1.01, y = y*1.01, filter = leaf, label = label,
                         hjust = "outward", color = if(color_leafs) group else NULL,
                         angle = -((-node_angle(x, y) + 90 ) %% 180) + 90),
                     size = 2.3, alpha = 1)} +
    geom_node_point(color = "grey60", shape = 16) +
    {if(color_leafs)
      scale_color_manual(values = attr(x, "group_colors"))} +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = unit(c(0,0,0,0),"cm"),
      panel.spacing = unit(c(0,0,0,0), "cm")) +
    {if(!label_leafs)
      expand_limits(x = c(-1.1, 1.1), y = c(-1.1, 1.1))} +
    {if(label_leafs)
      expand_limits(x = c(-1.45, 1.45), y = c(-1.45, 1.45))} +
    coord_fixed()
}

# function to plot a matrix of ggraph dendrograms from functional abundance (FA) data
# modified 'heat_tree_matrix()' function: https://github.com/grunwaldlab/metacoder/blob/master/R/heat_tree_matrix.R
# - to accept a subset of all possible comparisons defined by the comparison list of the main analysis
# this function is a wrapper for the 'metacyc_ggraph' and 'metacyc_key_ggraph' functions to generate sub-plots and
# a key plot labeled with the pathway hierarchy
ggraph_matrix <- function(igraph_list, label_small_trees =  FALSE,
                          min_eff_size = 1,
                          max_adj_p = 0.05,
                          label_leafs = TRUE,
                          color_leafs = TRUE,
                          key_size = 0.7,
                          key_x = 0,
                          key_y = 0,
                          matrix_x1 = 0.1,
                          matrix_y1 = 0.1,
                          matrix_x2 = 1,
                          matrix_y2 = 1,
                          seed = 1, output_files = TRUE,
                          workdir = workdir,
                          subdir = subdir) {

  require(cowplot)
  # identify individual comparisons
  treatments <- names(igraph_list)

  # calculate layout matrix based on 'treatments'
  layout_matrix <- create_layout_matrix(treatments)

  # Get subplot layout data
  # The x1, y1, x2, y2 refer to the relative size of the cowplot canvas size
  # control position and size of sub-plots by providing the lower left (x1, y1) and upper right (x2, y2) coordinates,
  # relative to the canvas size defined above (usually 0,0 (lower left) -> 1,1 (upper right))
  matrix_data <- calc_subplot_coords2(layout_matrix, treat = treatments,
                                      x1 = matrix_x1, y1 = matrix_y1,
                                      x2 = matrix_x2, y2 = matrix_y2,
                                      ks = key_size, kx = key_x, ky = key_y)

  # right join with 'treatments' data, remove rows without entry
  matrix_data <- dplyr::right_join(matrix_data,
                            tibble(plot_index = seq(1:length(treatments)),
                                   comparison = treatments),
                            by = "plot_index")
  base::rownames(matrix_data) <- matrix_data$plot_index

  # calculate subgraph dimensions
  subgraph_width <- calc_subplot_dims(layout_matrix, matrix_data, "width",
                                      ks = key_size, kx = key_x, ky = key_y)
  subgraph_height <- calc_subplot_dims(layout_matrix, matrix_data, "height",
                                       ks = key_size, kx = key_x, ky = key_y)

  # correct for small number of comparisons: smaller sub-plots!
  if (length(treatments) <= 3) {
    # make plots smaller
    subgraph_width <- subgraph_width*0.7
    subgraph_height <- subgraph_height*0.7
    # push smaller plots to outer limits
    matrix_data <- matrix_data %>%
      dplyr::mutate(across(c(x, y), ~ .x*1.2))
  }

  ### plot all the sub-plots
  set.seed(seed)
  sub_plots <- purrr::imap(igraph_list, metacyc_ggraph,
                           min.eff.size = min_eff_size,
                           max.eff = attr(igraph_list, "max_eff"),
                           min.adj.p = attr(igraph_list, "min_adj_p"),
                           max.adj.p = max_adj_p)

  # extract legend from first sub plot
  legend <- cowplot::get_legend(
    # create some space to the bottom of the legend (trbl)
    sub_plots[[1]] +
      ggplot2::theme(legend.box.margin = margin(300, 0, 0, if(label_leafs) -30 else 20))
  )

  # remove legends from the plots
  sub_plots_no_legend <- base::lapply(sub_plots, function(x) {x + ggplot2::theme(legend.position = "none")})

  ### plot the key-plot
  set.seed(seed)
  key_plot <- metacyc_key_ggraph(x = igraph_list[[1]], label_leafs = label_leafs)

  # add legend to key plot
  set.seed(seed)
  # key_plot <- cowplot::plot_grid(key_plot, legend, ncol = 1, rel_heights = c(1, .1))
  key_plot <- cowplot::plot_grid(key_plot, legend, nrow = 1, rel_widths = c(1, .3))

  ### Make matrix plot
  matrix_plot <- cowplot::ggdraw() +
    cowplot::draw_plot(key_plot, x = key_x, y = key_y, width = key_size, height = key_size) +
    ggplot2::theme(aspect.ratio = 1)

  for (k in seq_along(sub_plots_no_legend)) {
    if (! is.null(sub_plots_no_legend[[k]])) {
      matrix_plot <- matrix_plot + cowplot::draw_plot(sub_plots_no_legend[[k]],
                                                      x = matrix_data[k, "x"],
                                                      y = matrix_data[k, "y"],
                                                      width = subgraph_width,
                                                      height = subgraph_height)
    }
  }

  # export combined matrix plot
  file_name <- paste(output_dir9, paste(
    "Pathway", var_test, "ggraph.png", sep = "."), sep = "/") %>%
    gsub("\\\\", "/", .)

  fig_legend <- paste0(
    "Overview of functional abundance (FA) of MetaCyc pathways using PICRUSt2 followed by ALDEx2-based differential analyses.",
    " In the lower left a key plot shows the hierarchy of all ", attr(igraph_list, "sig_pw"),
    " significantly differential abundant pathways in the data set as a dendrogram.",
    " Smaller dendrograms depict the differential functional abundance of the comparisons between factors of variable '",
    str_replace_all(var_test, "_", "\\\\_"), "' as indicated.",
    " The color of the dots corresponds to the effect size while the diameter ",
    " represents the -log10-transformed, Benjamini-Hochberg corrected expected P-value of Welch’s t test.",
    " The P-values of higher hierarchies are combined p-values of the included pathways calculated using the minimum Holm approach.")

  attr(matrix_plot, "file_name") <- file_name
  attr(matrix_plot, "fig_legend") <- fig_legend

  # export image
  if (output_files) {
    ggplot2::ggsave(filename = file_name,
           plot = matrix_plot,
           width = 40, height = 40, units = "cm",
           bg = "white")
  }

  # return list
  return(list(key = key_plot,
              subs = sub_plots,
              matrix_pt = matrix_plot,
              matrix_data = matrix_data,
              sub_width = subgraph_width,
              sub_heigth = subgraph_height))
}
