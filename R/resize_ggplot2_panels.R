# helper script to fix plot panel size of one or multiple ggplots
library(tidyverse)
library(cowplot)
library(gtable)
library(lemon)
library(ggplotify)

# data(mtcars, package = "datasets")
# p1 <- ggplot(mtcars, aes(x=hp, y=mpg, color=factor(cyl), shape=factor(cyl))) +
#   geom_point(size=3) +
#   geom_smooth(method="lm", aes(fill=cyl)) +
#   #coord_fixed(ratio = 2)
#   theme (aspect.ratio = 0.5)
# 
# p2 <- p1 + facet_wrap(~factor(cyl), ncol = 2)
# p2

# p <- p1
# p.margin = 5
# p.width = 8
# p.height = NULL

# function to set ggplot panel size to fixed values (in cm)
# modified from here: https://gist.github.com/infotroph/a0a7fcf034e70de8ed75 and
# https://stackoverflow.com/questions/32580946/setting-absolute-size-of-facets-in-ggplot2
# define panel size with 'p.width' and 'p.height',
# 'p.margin' is margin arround plot (on canvas)
# p <- corr_plot
# p.width = 12 
# p.height = 12

fix_ggplot_panel <- function(
    p,
    p.margin = 5,
    p.width = NULL, 
    p.height = NULL) {

  # convert ggplot to grob
  if("ggplot" %in% class(p)){
    g <- ggplotGrob(p)
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
        asp <- as.numeric(str_remove_all(as.character(g$heights)[str_detect(as.character(g$heights), "null")], "null")) /
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
  
  # transform panel sizes in units
  pmargin <- unit(p.margin, "mm")
  pwidth <- unit(p.width, "cm")
  pheight <- unit(p.height, "cm")
  
  # identify name of empty panel
  pnls <- cowplot::plot_to_gtable(p) %>% gtable::gtable_filter("panel") %>%
    with(setNames(grobs, layout$name)) %>% purrr::keep(~identical(.x,zeroGrob()))
  
  if( length(pnls) == 0 ) {g <- ggplotGrob(p)} else { 
    
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

# p4 <- facet_shift_legend_fix_panel(p2)
# p4