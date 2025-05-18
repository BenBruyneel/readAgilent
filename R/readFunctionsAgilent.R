#' @title chromatogramInfo.AgilentExport
#'
#' @description translates the chromatogram description line for Agilent
#'  Masshunter chromatograms. Supported chromatogram descriptions are mass
#'  spectrometry descriptions (TIC, BPC & EIC), UV (DAD) & traces for the pump
#'  (Pressure, Flow, etc)
#'
#' @param commentString The description character vector to be converted. Note
#'  that the characters '\' and '#', which are present after the initial file
#'  read (readLines), need to be removed before using this function
#' @param defaultCollapse only used in case of unknown traces. For these, the
#'  resulting elements are joined together into a single character vector which
#'  can be split by this character
#'
#' @returns a data.frame
#'
#' @examples
#' chromatogramInfo.AgilentExport("+ESI TIC Scan Frag=125.0V Data0001.d ")
#' chromatogramInfo.AgilentExport("+ESI BPC Scan Frag=125.0V Data0001.d ")
#' chromatogramInfo.AgilentExport("+ESI EIC(372.8974) Scan Frag=125.0V Data0001.d ")
#' chromatogramInfo.AgilentExport("DAD1 - A:Sig=280.0,4.0  Ref=550.0,100.0 Data0001.d")
#' chromatogramInfo.AgilentExport("BinPump1 - A: Pressure Data0001.d")
#'
#' @export
chromatogramInfo.AgilentExport <- function(commentString, defaultCollapse = ";"){
  tempResult <- unlist(stringr::str_split(commentString, pattern = "\\s"))
  tempResult <- tempResult[nchar(tempResult) !=0 & tempResult != "-"]
  if (grepl(tempResult[1], pattern = "ESI") | grepl(tempResult[1], pattern ="APCI")){
    result <- data.frame(
      polarity = substr(tempResult[1],1,1),
      signal = substr(tempResult[1],2,nchar(tempResult[1])),
      scan = substr(tempResult[2],1,3),
      mz = as.numeric(stringr::str_extract(tempResult[2], pattern = "[:digit:]*\\.[:digit:]*")),
      MS = tempResult[3],
      fragmentor = as.numeric(stringr::str_extract(tempResult[4], pattern = "[:digit:]*\\.[:digit:]*")),
      filename = tempResult[length(tempResult)]
    )
  } else {
    if (grepl(tempResult[1], pattern = "DAD")){
      result <- data.frame(
        signal = tempResult[1],
        scan = tempResult[2],
        wavelength = as.numeric(stringr::str_extract(tempResult[2], pattern = "[:digit:]*\\.[:digit:]")),
        filename = tempResult[length(tempResult)]
      )
    } else {
      if (grepl(tempResult[1], pattern = "BinPump")){
        result <- data.frame(
          signal = tempResult[1],
          scan = paste(c(tempResult[3:(length(tempResult)-1)]), collapse = " "),
          filename = tempResult[length(tempResult)]
        )
      } else {
        result <- data.frame(
          unknown = paste(tempResult, collapse = defaultCollapse)
        )
      }
    }
  }
  return(result)
}

#' @title readChromatogram.AgilentExport.memory
#'
#' @description function factory that returns a function that takes the data
#'  from a character vector which is in the format of an Agilent chromatogram
#'  export (single chromatogram, not multiple) and returns an object (list) with
#'  two elements: data (data.frame) and info (list)
#'
#' @note this function gets called by the \code{link[readAgilent]{readChromatogram.AgilentExport}}
#'  function
#'
#' @param textLines character vector of the data in Agilent chromatogram export
#'  format: first line = description, second line = Agilent names for x & y (not
#'  used) and the rest of the lines are c(rownumbers, x, y). The rownumbers are
#'  ignored.
#' @param sep defines the separator for the rownumber-x-y data. Default is ','
#' @param translateComment defines the function to be used to translate the
#'  description line. Default is NA. \code{link[readAgilent]{readChromatogram.AgilentExport}}
#'  is this package's function that can be used.
#' @param removePatterns defines which characters to remove from the description
#'  line (defore translation)
#' @param dataColumns defines which columns to extract from the data. Default is
#'  c(2,3). The rownumbers are ignored
#' @param columnNames defines which names to give to the extracted columns.
#'  Default is 'rt' (retention time) for x and 'intensity' for y
#'
#' @returns a function that generates a list with two elements: data (data.frame)
#'  and info (list)
#'
#' @examples
#' demoFile <- fs::path_package("extdata", "Data0001.CSV", package = "readAgilent")
#' result <- readLines(demoFile, n = 9092)
#' result |> head()
#' result <- readChromatogram.AgilentExport.memory(result,
#'  translateComment = chromatogramInfo.AgilentExport)()
#' result$data |> head()
#' result$info
#'
#' @export
readChromatogram.AgilentExport.memory <- function(textLines,
                                                  sep = ",",
                                                  translateComment = NA,
                                                  removePatterns = c("\\\"", "#"),
                                                  dataColumns = c(2,3),
                                                  columnNames = c("rt", "intensity")){
  force(textLines)
  force(sep)
  force(translateComment)
  force(removePatterns)
  force(dataColumns)
  force(columnNames)
  function(...){
    tempComment <- dataInfo::strReplaceAll(textLines[1],
                                 pattern = removePatterns,
                                 replacement = "")
    tempdf <-purrr:: map_df(textLines[3:length(textLines)],
                            ~as.data.frame(t(stringr::str_split(.x,
                                                                pattern = sep)[[1]]))) %>%
      dplyr::select(tidyr::all_of(dataColumns))
    colnames(tempdf) <- columnNames
    tempdf$rt <- as.numeric(tempdf$rt)
    tempdf$intensity <- as.numeric(tempdf$intensity)
    dataInfo::readData(dataframe = tempdf,
             info = append(list(source = "Agilent",
                                comment = tempComment),
                           dataInfo::ifelseProper(
                             !identical(translateComment, NA),
                             as.list(translateComment(tempComment)),
                             NULL)
             )
    )()
  }
}

#' @title readChromatogram.AgilentExport
#'
#' @description function factory that returns a function that takes the data
#'  from a filename which is in the format of an Agilent chromatogram export
#'  (single & multiple chromatograms) and returns a list of objects (list) with
#'  each two elements: data (data.frame) and info (list)
#'
#' @note this function gets calls the \code{link[readAgilent]{readChromatogram.AgilentExport.memory}}
#'  function
#'
#' @param filename name of the Agilent Chromatogram export file. One export may
#'  contain more than one chromatogram
#' @param sep defines the separator for the rownumber-x-y data. Default is ','
#' @param seekStart default is "#". This defines the character the function needs
#'  to use to find the separate chromatograms in the export file
#' @param translateComment defines the function to be used to translate the
#'  description line. Default is NA. \code{link[readAgilent]{readChromatogram.AgilentExport}}
#'  is this package's function that can be used.
#' @param removePatterns defines which characters to remove from the description
#'  line (defore translation)
#' @param dataColumns defines which columns to extract from the data. Default is
#'  c(2,3). The rownumbers are ignored
#' @param columnNames defines which names to give to the extracted columns.
#'  Default is 'rt' (retention time) for x and 'intensity' for y
#'
#' @returns a function that generates a list of lists with each two elements:
#'  data (data.frame) and info (list)
#'
#' @examples
#' demoFile <- fs::path_package("extdata", "Data0001.CSV", package = "readAgilent")
#' result <- readChromatogram.AgilentExport(demoFile)()
#' length(result)
#' purrr::map_df(result, ~as.data.frame(.x$info))
#' result <- readChromatogram.AgilentExport(demoFile,
#'  translateComment = chromatogramInfo.AgilentExport)()
#' length(result)
#' purrr::map_df(result, ~as.data.frame(.x$info))
#' result[[3]]$data |> head()
#' plot(result[[3]]$data, type = "l")
#' result[[5]]$data |> head()
#' plot(result[[5]]$data, type = "l")
#'
#' @export
readChromatogram.AgilentExport <- function(filename, sep = ",", seekStart = "#",
                                           translateComment = NA,
                                           removePatterns = c("\\\"", "#"),
                                           dataColumns = c(2,3),
                                           columnNames = c("rt", "intensity")){
  force(filename)
  force(sep)
  force(seekStart)
  force(translateComment)
  force(removePatterns)
  force(dataColumns)
  force(columnNames)
  function(...){
    result <- list()
    tempLines <- readLines(filename)
    # note: first line is start & description, second is header
    starts <- which(grepl(tempLines, pattern = seekStart))[c(TRUE,FALSE)]
    starts <- append(starts, length(tempLines)+1)
    for (counter in 1:(length(starts)-1)){
      result[[counter]] <- readChromatogram.AgilentExport.memory(tempLines[starts[counter]:(starts[counter+1]-1)],
                                                                 translateComment = translateComment,
                                                                 removePatterns = removePatterns,
                                                                 dataColumns = dataColumns,
                                                                 columnNames = columnNames)()
    }
    return(result)
  }
}

#' @title peakListInfo.AgilentExport
#'
#' @description translates the information string which comes with export of
#'  integration info files (peaklists)
#'
#' @param string character vector to be translated
#' @param defaultCollapse only used in case of unknown traces. For these, the
#'  resulting elements are joined together into a single character vector which
#'  can be split by this character
#'
#' @returns a data.frame
#'
#' @examples
#' info <- "D:/MassHunter/Data/Ben/Data0001.d, +ESI TIC Scan Frag=125.0V Data0001.d "
#' peakListInfo.AgilentExport(info)
#' info <- "D:/MassHunter/Data/Ben/Data0001.d, +ESI EIC(591.2807) Scan Frag=125.0V Data0001.d "
#' peakListInfo.AgilentExport(info)
#'
#' @export
peakListInfo.AgilentExport <- function(string, defaultCollapse = ";"){
  string <- unlist(stringr::str_split(string, pattern = ","))
  result <- chromatogramInfo.AgilentExport(string[2], defaultCollapse = defaultCollapse)
  result$location <-string[1]
  return(result)
}

#' @title readPeaklist.AgilentExport.memory
#'
#' @description function factory that returns a function that takes the data
#'  from a character vector which is in the format of an Agilent peaklist
#'  export (single or multiple peaklist(s)) and returns list of lists with
#'  two elements: data (data.frame) and info (list)
#'
#' @note a peaklist = the result from the integration of traces/chromatograms by
#'  the Agilent Masshunter software
#'
#' @param textLines character vector of the data in Agilent peaklist export
#'  format: first line = description and the rest of the lines is the integration
#'  data. Rownumbers are ignored.
#' @param sep defines the separator for the peaklist data. Default is ','
#' @param startsString character vector that defines which lines are the start
#'  of a peaklist lines in textLines
#'
#' @returns function that generates a list of two objects: first element (info)
#'  is data.frame (with info) and second (data) is a list of data.frame's
#'
#' @examples
#' demoFile <- fs::path_package("extdata", "Data0001_MS.csv", package = "readAgilent")
#' result <- readLines(demoFile)
#' result <- readPeaklist.AgilentExport.memory(result, startsString = "Data0001")()
#' result$info
#' result$data[[1]]
#' result$data[[3]]
#' @export
readPeaklist.AgilentExport.memory <- function(textLines, sep = ",", startsString){
  force(textLines)
  force(sep)
  force(startsString)
  function(...){
    textLines <- purrr::map_chr(textLines, ~stringr::str_replace_all(.x, pattern = "\\\"", replacement = ""))
    textLines <- purrr::map_chr(textLines, ~stringr::str_replace_all(.x, pattern = "\\\\", replacement = "/"))
    starts <- which(purrr::map_lgl(textLines, ~nrow(stringr::str_locate_all(.x, pattern = startsString)[[1]]) == 2))
    results <- purrr::map_df(textLines[starts], ~peakListInfo.AgilentExport(.x))
    results$start <- starts+1
    results$end <- c(starts[-1]-1, length(textLines))
    pkList <- list()
    for (counter in 1:nrow(results)){
      pkList[[counter]] <- stringr::str_replace_all(textLines[results$start[counter]:results$end[counter]],
                                           pattern = paste(c(".*",startsString,"\\s,"), collapse = ""),
                                           replacement = "")
      if (length(pkList[[counter]]) == 1){
        pkList[[counter]] <- NA
      } else {
        pkList[[counter]][1] <- stringr::str_remove_all(
          stringr::str_replace_all(pkList[[counter]][1],
                                   pattern = "%",
                                   replacement = "Perc"),
          pattern = "\\s")
        pkList[[counter]] <- utils::read.csv(text = pkList[[counter]],
                                             sep = ",",
                                             header = T)
      }
    }
    return(list(info = results, data = pkList))
  }
}

#' @title readPeaklist.AgilentExport
#'
#' @description function factory that returns a function that takes the data
#'  from a file which is in the format of an Agilent peaklist
#'  export (single or multiple peaklist(s)) and returns list of dataElements
#'
#' @note a peaklist = the result from the integration of traces/chromatograms by
#'  the Agilent Masshunter software
#'
#' @param filename name of the peaklist file
#' @param sep defines the separator for the peaklist data. Default is ','
#' @param startsString character vector that defines which lines are the start
#'  of a peaklist lines in the file
#'
#' @returns a list of dataElements
#'
#' @examples
#' demoFile <- fs::path_package("extdata", "Data0001_MS.csv", package = "readAgilent")
#' result <- readPeaklist.AgilentExport(demoFile, startsString = "Data0001")()
#' length(result)
#' result[[1]]
#' result[[1]]$info
#' result[[1]]$data
#' result[[3]]$info
#' result[[3]]$data
#'
#' @export
readPeaklist.AgilentExport <- function(filename, sep = ",",
                                        startsString = paste0(stringr::str_remove(tools::file_path_sans_ext(basename(filename)),
                                                                         pattern = "_.*"), ".d")){
  force(filename)
  force(sep)
  force(startsString)
  function(...){
    result <- list()
    textLines <- readLines(filename)
    tempResult <- readPeaklist.AgilentExport.memory(textLines = textLines,
                                                    sep = sep,
                                                    startsString = startsString)()
    for (counter in 1:nrow(tempResult[[1]])){
      result[[counter]] <- dataInfo::readData(dataframe = tempResult[[2]][[counter]],
                                    info = as.list(tempResult[[1]][counter,]))()
    }
    return(result)
  }
}

