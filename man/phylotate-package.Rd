\name{phylotate-package}
\alias{phylotate-package}
\alias{phylotate}
\docType{package}
\title{Phylogenies with Annotations}
\description{
  \pkg{phylotate} provides functions that allow you to read and write
  NEXUS and Newick trees containing annotations, including those
  produced by MrBayes.

  It does this by extending \pkg{ape}'s \code{phylo} object type with
  extra data members containing per-node annotation information. This
  information is stored in such a way that it can be manipulated easily
  and will survive most manipulations using standard \pkg{ape} functions
  (e.g. \code{\link{reorder}}, \code{chronopl}).

  See the documentation for the \code{\link{parse_annotated}} function
  for more information on how annotations are stored internally.

  The functions you probably want to use for most things are
  \code{\link{read_annotated}} and \code{\link{write_annotated}}.
}
\author{Daniel Beer \email{dlbeer@gmail.com}}
\keyword{phylo}
