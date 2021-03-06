% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BinRepProbPlot.R
\name{Prepare_S_phase_cells_forBinRepProb}
\alias{Prepare_S_phase_cells_forBinRepProb}
\title{Prepare S phase cells for BinRepProb plots}
\usage{
Prepare_S_phase_cells_forBinRepProb(
  S,
  RT,
  Early.cells = c(0, 30),
  Mid.cells = c(40, 60),
  Late.cells = c(70, 100)
)
}
\arguments{
\item{S, }{S single cell CNV dataframe}

\item{RT, }{pseudobulk replication timing dataframe}

\item{Early.cells, }{percentage of replication range defying Early replicating cells (0-100)}

\item{Mid.cells, }{percentage of replication range defying Mid replicating cells (0-100)}

\item{Late.cells, }{percentage of replication range defying Late replicating cells (0-100)}
}
\description{
Prepare S phase cells for BinRepProb plots
}
