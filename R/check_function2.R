### check if color_set, label_set and labels are the same length

check_color_set <- function(color_set, label_set, labels) {
  lengths <- c(colors = length(color_set),
               label_order = length(label_set),
               labels = length(labels))

  mismatches <- lengths[lengths != lengths[1]]

  if (length(mismatches) > 0) {
    spsComps::shinyCatch(
      stop(glue::glue("Length mismatch: colors {lengths[1]}, label order {lengths[2]}, labels {lengths[3]}.")),

      blocking_level = "error")
    return()

  }
}

### list all in


