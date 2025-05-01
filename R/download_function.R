### download rmd outputs function ###
download_function <- function(folder_name, output_directory) {
  downloadHandler(
    filename = function() {
      paste(folder_name, Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      output_directory
      if (length(list.files(output_directory)) == 0) {
    spsComps::shinyCatch("Please generate (a) file(s) to download.", blocking_level = "error")
  }

      zip::zip(zipfile = file,
               files = list.files(output_directory),
               root = output_directory)
    },
    contentType = "application/zip"
  )
}
