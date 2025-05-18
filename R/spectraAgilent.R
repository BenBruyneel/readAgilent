#' @title mzData
#'
#' @description
#'  R6 Class dealing with mass spectrum data as output by Agilent Masshunter.
#'   Main drawback of this class is that currently all data needs to be loaded
#'   into memory. The .mzData.xml files can be quite big, especially when using
#'   high quality profile data. Profile data may also give issues with peak
#'   detection as peaks on some equipment have a tendency to be assymmetric.
#'   Consider using centroided data and/or higher background cutoff when
#'   acquiring
#'   
#' @examples
#' demoFile <- fs::path_package(
#'   "extdata",
#'   "WKL-0001.mzdata.xml",
#'   package = "readAgilent"
#' )
#' demoFile
#' spectrData <- mzData$new(filename = demoFile)
#' spectrData
#' spectrData$info(centroided = TRUE)
#' 
#' MS.Analysis::plotSpectrum(
#'   spectrData$spectrum(number = 1),
#'   centroidPlot = TRUE,
#'   labelCutOff = 0.05,
#'   labelColor = "red",
#'   labelAngle = 15
#' )
#' MS.Analysis::plotSpectrum(
#'   spectrData$spectrum(number = 2),
#'   centroidPlot = TRUE,
#'   labelCutOff = 0.05,
#'   labelColor = "red",
#'   labelAngle = 15
#' )
#' 
#' @export
mzData <- R6::R6Class(
  "mzData",
  private = list(
    # this specifies whether the base64 encoded data in the file is 'little' or 'big' endian.
    endian_ = "little",
    # filename from which the data comes
    filename_ = as.character(NA),
    # private variable intended to hold the .mzData.xml file (which is essentially a text file)
    data_ = NA,
    # loads the data from an .mzData.xml file. This is achieved
    # via the XML library function xmlToList.
    # 
    # filename character vector specifying the file name of the
    #  .mzData.xml file from which the data is to be loaded. If no file with
    #  the specified path/filename exists, the function returns FALSE, otherwise
    #  the file is loaded into the private field data_ and returns TRUE
    #  Default value is NA
    # 
    # returns logical vector
    load_ = function(filename = as.character(NA)) {
      if (!is.na(filename)) {
        if (file.exists(filename)) {
          private$filename_ <- filename
          private$data_ <- XML::xmlToList(filename)
          return(TRUE) # succesfull
        }
      }
      return(FALSE) # unsuccesfull
    }
  ),
  public = list(
    #' @description initializes the object
    #'
    #' @param filename name of the .mzData.xml file (exported by Agilent Software)
    #' @param endian specifies whether the base64 encoded data in the .mzData.xml
    #'  file is 'little' or 'big' endian
    #' @param message logical vector. If TRUE (default) a message will be printed
    #'  when the data has been loaded. Files are sometimes big and may take a while
    #'  to load
    #'  
    #' @returns a new 'mzData' object
    #' 
    #' @export
    initialize = function(
      filename = as.character(NA),
      endian = "little",
      message = TRUE
    ) {
      private$endian_ <- endian
      msg <- " No data loaded\n"
      if (!is.na(filename)) {
        if (file.exists(filename)) {
          if (private$load_(filename = filename)) {
            msg <- " Data loaded succesfully\n"
          }
        }
      }
      if (message) {
        cat(msg)
      }
      invisible(self)
    },
    #' @description
    #' For printing purposes: prints info on the data loaded into the object
    #'
    print = function() {
      if (self$empty) {
        cat(NA)
        cat("\n")
      } else {
        cat(paste(c(" Filename      : ", self$filename, "\n"), collapse = ""))
        cat(paste(c(" # of Spectra  : ", self$length, "\n"), collapse = ""))
      }
    },
    #' @description
    #'  retrieves a spectrum from the (internal) .mzData.xml data
    #' @note This may take a while since data needs to be decoded and organized
    #'  into a proper data.frame
    #'  
    #' @param number number of spectrum to retrieve
    #' 
    #' @returns data.frame, with columns 'mz' and 'intensity', or NA
    spectrum = function(number = 1) {
      result <- NA
      if (!self$empty) {
        if (number <= self$length) {
          mzs <- base64enc::base64decode(
            self$data$spectrumList[[number]]$mzArrayBinary$data$text
          )
          ints <- base64enc::base64decode(
            self$data$spectrumList[[number]]$intenArrayBinary$data$text
          )
          mzs <- purrr::map_dbl(
            1:(length(mzs) / 8),
            ~ readBin(
              mzs[(1 + ((.x - 1) * 8)):(.x * 8)],
              what = "numeric",
              size = 8,
              endian = self$endian
            )
          )
          ints <- purrr::map_dbl(
            1:(length(ints) / 4),
            ~ readBin(
              ints[(1 + ((.x - 1) * 4)):(.x * 4)],
              what = "numeric",
              size = 4,
              endian = self$endian
            )
          )
          result <- data.frame(mz = mzs, intensity = ints)
        }
      }
      return(result)
    },
    #' @description
    #'  returns a data.frame with info on spectra in the object
    #'  
    #' @param number integer vector. One or more numbers of the spectra for
    #'  which the info is to be retrieved
    #' @param centroided logical vector. Adds the centroided column to the
    #'  info data.frame. Sometimes needed for plotting and analysis purposes.
    #'
    #' @returns data.frame or NA
    info = function(number = NA, centroided = FALSE) {
      if ((!(identical(number, NA)) & (length(number) == 1))) {
        result <- NA
        if (!self$empty) {
          if (number <= self$length) {
            result <- as.data.frame(t(
              self$data$spectrumList[[
                number
              ]]$spectrumDesc$spectrumSettings$acqSpecification$.attrs
            ))
            for (counter in 1:(length(
              self$data$spectrumList[[
                number
              ]]$spectrumDesc$spectrumSettings$spectrumInstrument
            ) -
              1)) {
              result <- dplyr::bind_cols(
                result,
                data.frame(
                  x = self$data$spectrumList[[
                    number
                  ]]$spectrumDesc$spectrumSettings$spectrumInstrument[[
                    counter
                  ]][[4]]
                )
              )
              colnames(result)[ncol(result)] <- self$data$spectrumList[[
                number
              ]]$spectrumDesc$spectrumSettings$spectrumInstrument[[counter]][[
                3
              ]]
            }
            result <- dplyr::bind_cols(
              result,
              as.data.frame(t(
                self$data$spectrumList[[
                  number
                ]]$spectrumDesc$spectrumSettings$spectrumInstrument$.attrs
              ))
            )
          }
        }
        if ("TimeInMinutes" %in% colnames(result)) {
          result$Start <- purrr::map_dbl(
            result$TimeInMinutes,
            ~ as.numeric(stringr::str_split(.x, pattern = "-")[[1]][1])
          )
          result$End <- purrr::map_dbl(
            result$TimeInMinutes,
            ~ as.numeric(stringr::str_split(.x, pattern = "-")[[1]][2])
          )
        }
        result$centroided <- centroided
        return(result)
      } else {
        if (identical(number, NA)) {
          number <- 1:self$length
        }
        return(purrr::map_df(number, ~ self$info(number = .x, centroided = centroided)))
      }
    },
    #' @description function factory that generates a function to get a spectrum
    #'  in the object and to return it as a \link[dataInfo]{dataElement} (for use with
    #'  \link[MS.Analysis]{msInfo})
    #'  
    #' @param number integer value (single value!) specifying which spectrum to
    #'  retrieve from the object
    #'  
    #' @returns dataElement object
    read = function(number = 1) {
      force(number)
      function(...) {
        if (length(number) > 1){
          number <- number[1]
        }
        result <- list(dataInfo::readData())
        if ((number > 0) & (number <= self$length)) {
          result <- list(dataInfo::readData(
            dataframe = self$spectrum(number),
            info = append(
              append(list(source = "mzData"), as.list(self$info(number))),
              list(filename = self$filename)
            )
          )())
        }
        return(result)
      }
    }
  ),
  active = list(
    #' @field empty if TRUE, then data has been read. If FALSE, then not.
    #'  Read only
    empty = function(value) {
      if (missing(value)) {
        return(is.na(private$filename_))
      } else {
        # nothing, read-only
      }
    },
    #' @field filename shows filename after reading data into the object
    #'  Read only
    filename = function(value) {
      if (missing(value)) {
        return(private$filename_)
      } else {
        # nothing, read-only
      }
    },
    #' @field endian gets/sets the type of endian ('little' or 'big') to be used
    #'  when dealing with the .mzdata.xml data
    endian = function(value) {
      if (missing(value)) {
        return(private$endian_)
      } else {
        if (value %in% c("little", "big")) {
          private$endian_ <- value
        } else {
          warning("Endian not set, value must be 'little' or ' big' ")
        }
      }
    },
    #' @field data provides access to the 'raw' data after it has been read from
    #'  the file.
    data = function(value) {
      if (missing(value)) {
        return(private$data_)
      } else {
        private$data_ <- value
      }
    },
    #' @field length displays the number of spectra present in the object.
    #'  Read only
    length = function(value) {
      if (missing(value)) {
        if (!self$empty) {
          return(as.numeric(private$data_$spectrumList[length(
            private$data_$spectrumList
          )]$.attrs[["count"]]))
        } else {
          return(0)
        }
      } else {
        # nothing, read-only
      }
    }
  )
)
