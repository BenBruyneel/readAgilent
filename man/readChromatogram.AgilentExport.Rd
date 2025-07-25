% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readFunctionsAgilent.R
\name{readChromatogram.AgilentExport}
\alias{readChromatogram.AgilentExport}
\title{readChromatogram.AgilentExport}
\usage{
readChromatogram.AgilentExport(
  filename,
  sep = ",",
  seekStart = "#",
  translateComment = NA,
  removePatterns = c("\\\\\\"", "#"),
  dataColumns = c(2, 3),
  columnNames = c("rt", "intensity")
)
}
\arguments{
\item{filename}{name of the Agilent Chromatogram export file. One export may
contain more than one chromatogram}

\item{sep}{defines the separator for the rownumber-x-y data. Default is ','}

\item{seekStart}{default is "#". This defines the character the function needs
to use to find the separate chromatograms in the export file}

\item{translateComment}{defines the function to be used to translate the
description line. Default is NA. \code{link[readAgilent]{readChromatogram.AgilentExport}}
is this package's function that can be used.}

\item{removePatterns}{defines which characters to remove from the description
line (defore translation)}

\item{dataColumns}{defines which columns to extract from the data. Default is
c(2,3). The rownumbers are ignored}

\item{columnNames}{defines which names to give to the extracted columns.
Default is 'rt' (retention time) for x and 'intensity' for y}
}
\value{
a function that generates a list of lists with each two elements:
data (data.frame) and info (list)
}
\description{
function factory that returns a function that takes the data
from a filename which is in the format of an Agilent chromatogram export
(single & multiple chromatograms) and returns a list of objects (list) with
each two elements: data (data.frame) and info (list)
}
\note{
this function gets calls the \code{link[readAgilent]{readChromatogram.AgilentExport.memory}}
function
}
\examples{
demoFile <- fs::path_package("extdata", "Data0001.CSV", package = "readAgilent")
result <- readChromatogram.AgilentExport(demoFile)()
length(result)
purrr::map_df(result, ~as.data.frame(.x$info))
result <- readChromatogram.AgilentExport(demoFile,
 translateComment = chromatogramInfo.AgilentExport)()
length(result)
purrr::map_df(result, ~as.data.frame(.x$info))
result[[3]]$data |> head()
plot(result[[3]]$data, type = "l")
result[[5]]$data |> head()
plot(result[[5]]$data, type = "l")

}
