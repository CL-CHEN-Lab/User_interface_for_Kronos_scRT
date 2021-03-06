% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Load_multiple_df.R
\name{load_multiple_df}
\alias{load_multiple_df}
\title{Loads and binds a dataframe divided in multiple files}
\usage{
load_multiple_df(dirs, delim = "\\t", ...)
}
\arguments{
\item{dirs, }{list of directories}

\item{delim, }{delimiter for read_delim}

\item{..., }{other parameters to pass to read_delim}
}
\value{
a string informing the user file format is the required one
}
\description{
Loads and binds a dataframe divided in multiple files
}
