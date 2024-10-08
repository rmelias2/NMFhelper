% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmf_pull_components.R
\name{nmf_pull_components}
\alias{nmf_pull_components}
\title{nmf_pull_components}
\usage{
nmf_pull_components(nmf.object)
}
\arguments{
\item{nmf.object}{NMF object outputed from either the NMF package or the CoGAPs package}
}
\value{
Returns a list of NMF outputs, including the A and P matrix, Scaled values, Maximum P values per column, colors assigned to each pattern.
}
\description{
This function takes the output of an NMF run which are the A and P matrix, and performs common transformations that are used for downstram visualization.
}
\examples{
nmf_pull_components(NMF_object)
}
