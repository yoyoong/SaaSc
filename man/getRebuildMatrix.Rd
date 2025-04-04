% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SaaSc.R
\name{getRebuildMatrix}
\alias{getRebuildMatrix}
\title{Get rebuilt matrix}
\usage{
getRebuildMatrix(object, method = NULL, dim = NULL, avc = NULL)
}
\arguments{
\item{object}{A Seurat object processed by rebuildMatrix().}

\item{method}{Dimensionality reduction and reconstruction method, in a character. The choices are NULL(do not rebuild), "mca" , "nmf" and "pca".}

\item{dim}{The number of selected dimensions, 1 <= dim <= 50. Enter only one of avc and dim.}

\item{avc}{Accumulative variance contribution threshold, 0 <= avc <= 1. Enter only one of avc and dim.}
}
\value{
A matrix of rebuilt matrix by reconstruction method.
}
\description{
Get rebuilt matrix
}
