% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_hits.R
\name{plot_hits}
\alias{plot_hits}
\title{Genespace plotting routines}
\usage{
plot_hits(
  hits,
  onlyOG = TRUE,
  onlyAnchors = TRUE,
  useBlks = TRUE,
  onlyArrayReps = TRUE,
  onlyBuffer = FALSE,
  reorderChrs = TRUE,
  minGenes2plot = 10,
  gapProp = 0.01,
  cols = NULL,
  alpha = 1,
  axisTitleCex = 0.6,
  darkChrFill = "grey60",
  lightChrFill = "grey85",
  emptyChrFill = "grey97",
  missingHitCol = "grey40",
  chrLabCex = 0.4,
  returnSourceData = F,
  chrLabFun = function(x) gsub("^0", "", gsub("^chr|^scaffold|^lg|_", "", tolower(x)))
)
}
\arguments{
\item{hits}{data.table of hits}

\item{onlyOG}{logical, should the plot only show orthogroup-constrained hits?}

\item{onlyAnchors}{logical, should the plot only show syntenic anchors?}

\item{useBlks}{logical, should the plot only show hits in blocks?}

\item{onlyArrayReps}{logical, should the plot only show only array reps?}

\item{onlyBuffer}{logical, should the plot only show hits in buffers?}

\item{reorderChrs}{logical, should the chromosomes be re-ordered based on
synteny?}

\item{minGenes2plot}{integer specifying the minimum number of genes on a
chr to plot}

\item{gapProp}{numeric (0-1) specifying the proportional size of gaps
relative to the length of the largest genome}

\item{cols}{vector of colors to use for points}

\item{alpha}{numeric (0-1) specifying transparency of the points}

\item{axisTitleCex}{character expansion for the axes and genome labels}

\item{darkChrFill}{color of the most dense chr backgrounds}

\item{lightChrFill}{color of the least populated chr backgrounds}

\item{emptyChrFill}{color of the empty chr backgrounds}

\item{missingHitCol}{character string or number coercible to an R color. The
point color of hits that are not in blocks}

\item{chrLabCex}{character expansion for the chromosome labels}

\item{returnSourceData}{logical, should the source data to build the plot
be returned?}

\item{chrLabFun}{function to parse chr IDs to make them more readible}
}
\description{
\code{plot_hits} Genespace plotting routines
}
\details{
...
}
\examples{
\dontrun{
# coming soon
}

}
