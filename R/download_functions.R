
#' Download brazilian microdata.
#'
#'
#' @param dataset Standardized name of brazilian public microdadata. See available datasets with get_available_datasets()
#' @param i       Period(year/quarter) to download, use get_available_periods(dataset) to see available periods
#' @param unzip  (optional) logical. Should files be unzipped after download?
#' @param replace (optional) logical. Should an existing version of the data be replaced?
#' @param root_path (optional) a path to the directory where dataset should be downloaded
#'
#' @examples
#' \dontrun{
#'
#' download_sourceData("PNAD", 2014, unzip = T, root_path = "F:/Datasets/PNAD", replace = T)
#'
#'}
#'
#' @importFrom  utils download.file installed.packages  object.size read.csv2  unzip
#' @importFrom  RCurl getURL
#' @export
download_sourceData <- function(dataset, i, unzip=T , root_path = NULL, replace = FALSE){
  filename = ""
  success = FALSE
  size = NA

  dataset_list<- get_available_datasets()

  #Test if parameters are valid

  if( !(dataset %in% dataset_list ) ) {
    stop(paste0("Invalid dataset. Must be one of the following: ",paste(dataset_list, collapse=", ")) ) }

  metadata <-  read_metadata(dataset)

  i_range<- get_available_periods(metadata)
  if (!(i %in% i_range)) { stop(paste0("period must be in ", paste(i_range, collapse = ", "))) }


  ft_list  <- names(metadata)[grep("ft_", names(metadata))]
  data_path <- metadata[metadata$period == i, "path"]  %>% unlist(use.names = F) %>% paste0(".*?")
  data_file_names<- metadata[metadata$period ==i, ft_list] %>% unlist(use.names = FALSE)
  data_file_names<- paste0(data_path,
                           gsub(pattern = ".+?&", replacement = "", data_file_names[!is.na(data_file_names)]),
                           "$")

  if (!replace) {
    if(any(grepl(pattern = paste0(data_file_names,collapse = "|"), x = list.files(recursive = TRUE, path = ifelse(is.null(root_path), ".", root_path))))) {
      stop(paste0("This data was already downloaded. Check:\n - ",
                  paste(grep(pattern = paste0(data_file_names,collapse = "|"), list.files(recursive = TRUE,path = ifelse(is.null(root_path), ".", root_path), full.names = TRUE), value = T), collapse = "\n - "),
                  ")\n\nIf you want to overwrite the previous files add replace=T to the function call."))

    }
  }

  md <- metadata[metadata$period== i,]

  link <- md$download_path
  data_file_names<- md
  if(is.na(link)){stop("Can't download dataset, there are no information about the source")}

  if(!is.null(root_path)){
    if(!dir.exists(root_path)){ # file.exists fails when path has a trailing slash
      stop(paste0("Can't find ",root_path))
    }}

  if(md$download_mode == "ftp"){
    filenames <- RCurl::getURL(link, ftp.use.epsv = FALSE, ftplistonly = TRUE,
                               crlf = TRUE)
    filename <- file_dir <- gsub(link, pattern = "/+$", replacement = "", perl = TRUE) %>% gsub(pattern = ".+/", replacement = "")
    new_dir <- paste(c(root_path,file_dir), collapse = "/")
    if(!dir.exists(new_dir)) dir.create(new_dir)

    if (dataset == "CENSOPernambuco") {
      filenames <- "PE.zip"
    } else {
      filenames <- RCurl::getURL(link, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE)
      filenames<- strsplit(filenames, "\r*\n")[[1]]
    }

    file_links <- paste(link, filenames, sep = "")
    download_success <- rep(FALSE, length(filenames))

    max_loops = 20
    loop_counter = 1

    while(!all(download_success) & loop_counter < max_loops) {
      dest.files.all = sapply(filenames, function(x) {paste(c(root_path,file_dir, x),collapse = "/")})

      for(y in seq_along(filenames)[!download_success]){

        dest.files = paste(c(root_path,file_dir, filenames[y]),collapse = "/")
        download_success[y] = FALSE
        try({download.file(file_links[y],destfile = dest.files)
        download_success[y] = TRUE
      })

      }
      if(sum(file.info(dest.files.all)$size) < 100000) {

        success = F
        if(loop_counter == max_loops - 1){
          message(paste0("Downloaded files for period ", i," on the ", loop_counter, "th try were too small. Possible corruption." ))


        }else{
          message(paste0("Downloaded files for period ", i," on the ", loop_counter, "th try were too small, possible corruption, retrying download..." ))

        }
      }else{
        success = T
      }
}


    loop_counter = loop_counter + 1


    if(!all(download_success)){ message(paste0("The download of the following files failed:\n"),
                                       paste(filenames[!download_success], collapse = "\n"))}

    }else{

    filename <- link %>% gsub(pattern = ".+/", replacement = "")
    file_dir <- filename %>% gsub( pattern = "(\\.zip)|(\\.7z)|(\\.rar)", replacement = "")
    dest.files <- paste(c(root_path,filename),collapse = "/")

    print(paste("download mode", md$download_mode))
    print(link)
    print(filename)
    print(paste("file dir", file_dir))



    max_loops  = 4
    loop_counter = 1

    while(!(success) & loop_counter< max_loops){

    try({ download.file(link,destfile = dest.files)

      success = TRUE

    })
#    }
#
#
     if(success == T){
     if(sum(file.info(dest.files)$size) < 100000){
#
#       success = F
      if(loop_counter == max_loops - 1){
      message(paste0("Downloaded files for period ", i," on the ", loop_counter, "th try were too small. Possible corruption." ))


      }else{
      message(paste0("Downloaded files for period ", i," on the ", loop_counter, "th try were too small, possible corruption, retrying download..." ))

      }
    }

      loop_counter = loop_counter + 1
    }}

    if (unzip==T & success == T){
      #Won't use 'archive' in the main download function until its on CRAN
      #Unzipping main source file:
      #if(grepl(filename, pattern = "\\.7z")){

       # archive::archive_extract(paste(c(root_path,filename),collapse = "/") , paste(c(root_path,file_dir),collapse = "/"))
      #}else{
      unzip(paste(c(root_path,filename),collapse = "/") ,exdir = paste(c(root_path,file_dir),collapse = "/"))
      #}
    }

}
  if (unzip==T & success == T){


    ##unzipping the data files (in case not unziped above)
    intern_files<- list.files(paste(c(root_path,file_dir),collapse = "/"), recursive = TRUE,all.files = TRUE, full.names = TRUE)
    zip_files<- intern_files[grepl(pattern = "\\.zip$",x = intern_files)]
    rar_files<- intern_files[grepl(pattern = "\\.rar$",x = intern_files)]
    r7z_files<- intern_files[grepl(pattern = "\\.7z$",x = intern_files)]

    for(zip_file in zip_files){
      exdir<- zip_file %>% gsub(pattern = "\\.zip", replacement = "")
      unzip(zipfile = zip_file,exdir = exdir )
    }

#   for(zip_file in r7z_files){
#     exdir <- zip_file %>% gsub(pattern = "\\.7z",replacement = "")
#     archive::archive_extract(zip_file,exdir)
#   }
    # check if package "archive" is installed before trying to extract the .rar files.
    if(("archive" %in% installed.packages()[,1])){
      for(zip_file in c(rar_files,r7z_files)){
        exdir<- zip_file %>% gsub("\\.7z","",gsub(pattern = "\\.rar", replacement = ""))
        archive::archive_extract(zip_file, exdir)
        cat(paste0("Extracted ", zip_file,"\n"))
      }
    }
  }
 size=""
  try({if(all(file.info(paste(c(root_path,filename),collapse = "/"))$isdir)){

    try({size<- as.object_size(sum(file.info(list.files(paste(c(root_path,filename), collapse = "/"), recursive = TRUE, full.names = T))$size)) %>%
      format(units = "Mb")})
  }else{
    try({size = as.object_size(file.size(paste(c(root_path,filename),collapse = "/"))) %>%
      format(units = "Mb")})
  }
    })
  info.output<- data.frame(name = filename, link = link, success = success, size =  size, stringsAsFactors = F)

  return(info.output)

}



#' Wrapper for unzipping lots of .rar and .7z files with archive::archive() .
#'
#'
#' @param root_path  path to the directory where files are stored(will look for files recursively inside that folder)
#'
#' @examples
#' \dontrun{
#'
#' unzip_all_7z_rar("Datasets/micro_censo_escolar_2010")
#'
#'}
#'
#' @export
unzip_all_7z_rar <- function(root_path){


  if(!("archive" %in% installed.packages()[,1])){
    stop("The package 'archive' is needed to unzip 7z and rar files. Install it with devtools::install_github('jimhester/archive') \n More info at: https://github.com/jimhester/archive")
  }


  # #unzipping the data files (in case not unziped above)
  intern_files<- list.files(root_path, recursive = TRUE,all.files = TRUE, full.names = TRUE)
  zip_files<- intern_files[grepl(pattern = "\\.zip$",x = intern_files)]
  rar_files<- intern_files[grepl(pattern = "\\.rar$",x = intern_files)]
  r7z_files<- intern_files[grepl(pattern = "\\.7z$",x = intern_files)]

  for(zip_file in zip_files){
    exdir<- zip_file %>% gsub(pattern = "\\.zip", replacement = "")
    unzip(zipfile = zip_file,exdir = exdir )
    cat(paste0("Unzipped ", zip_file,"\n"))
  }
  for(zip_file in r7z_files){
    exdir<- zip_file %>% gsub(pattern = "\\.7z", replacement = "")
    archive::archive_extract(zip_file, exdir)
    cat(paste0("Unzipped ", zip_file,"\n"))

  }
  for(zip_file in rar_files){
    exdir<- zip_file %>% gsub(pattern = "\\.rar", replacement = "")
    archive::archive_extract(zip_file, exdir)
    cat(paste0("Unzipped ", zip_file,"\n"))

  }

}


