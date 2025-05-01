#' custom_theme
#'
#' @description Sets a custom theme for the app.
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
custom_theme <- function() {
  create_theme(
    adminlte_color(
      light_blue = "#427d7b"
    ),
    adminlte_sidebar(
      dark_bg = "#343a40",   # Dark sidebar background
      dark_hover_bg = "#495057",  # Hover color for sidebar
      dark_color = "#ffffff"  # Sidebar text color
    ),
    adminlte_global(
      content_bg = "#14161e",  # Dark background for the content
      box_bg = "#28333a",      # Box background in dark mode
      info_box_bg = "#1f2c33"
    )
  )
}
