\name{read_annotated}
\alias{read_annotated}
\title{Read an annotated phylogenetic tree}
\usage{
read_annotated(filename, format="nexus")
}
\arguments{
  \item{filename}{a file to read tree data from}
  \item{format}{a format specifier; either \code{"nexus"} or \code{"newick"}}
}
\description{
  This function takes the given file, containing data in either NEXUS or
  Newick format, and returns annotated phylogenetic trees.
}
\value{
  an object of type \code{"phylo"} or \code{"multiPhylo"}, augmented with
  node annotations.
}
\details{
  The given file text is parsed and a tree object is constructed which
  can be used with the functions in the \pkg{ape} package. Annotations
  of the kind produced by, for example, MrBayes, are parsed and
  preserved in the returned object.

  See \code{\link{parse_annotated}} for more information about the
  structure of the returned value.
}
\author{Daniel Beer \email{dlbeer@gmail.com}}
\seealso{
  \code{\link{print_annotated}},
  \code{\link{parse_annotated}},
  \code{\link{write_annotated}}
}
\keyword{phylo}
