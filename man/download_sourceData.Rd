% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_functions.R
\name{download_sourceData}
\alias{download_sourceData}
\title{Download brazilian microdata.}
\usage{
download_sourceData(dataset, i, unzip = T, root_path = NULL, replace = FALSE)
}
\arguments{
\item{dataset}{Standardized name of brazilian public microdadata. See available datasets with get_available_datasets()}

\item{i}{Period(year/quarter) to download, use get_available_periods(dataset) to see available periods}

\item{unzip}{(optional) logical. Should files be unzipped after download?}

\item{root_path}{(optional) a path to the directory where dataset should be downloaded}

\item{replace}{(optional) logical. Should an existing version of the data be replaced?}
}
\description{
Download brazilian microdata.
}
\examples{
\dontrun{

download_sourceData("PNAD", 2014, unzip = T, root_path = "F:/Datasets/PNAD", replace = T)

}

}
