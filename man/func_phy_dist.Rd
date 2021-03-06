% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_FPD.R
\name{func_phy_dist}
\alias{func_phy_dist}
\title{Functional-Phylogenetic Distances (FPD)}
\usage{
func_phy_dist(FDist, PDist, phyloWeight, p, ...)
}
\arguments{
\item{FDist}{A matrix of functional trait distances}

\item{PDist}{A matrix of phylogenetic distances}

\item{phyloWeight}{The phylogenetic scaling parameter \emph{a} as
described in Cadotte et al 2013. A value of 1 denotes only phylogenetic
information while a value of 0 indicates only functional trait 
information. The value must be between 0 and 1.}

\item{p}{the power to raise each matrix to when combining the two 
distance matrices.}

\item{...}{Ignored here}
}
\value{
A matrix with weighted interspecific distances
}
\description{
Calculate hybrid functional phylogenetic distance matrices
based on the method of Cadotte et al 2013
}
\note{
The function does not check to make sure names in 
\code{FDist} and \code{PDist match}. Additionally, if
NAs or NaNs in the distance matrix are generated by missing data,
they will be replaced with matrix level mean distance values and
a warning message will be generated.
}
\seealso{
\code{\link{rarefy_FPD}}
}
\author{
Sam Levin
}
