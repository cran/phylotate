\name{write_annotated}
\alias{write_annotated}
\title{Write an annotated phylogenetic tree to a file}
\usage{
write_annotated(tree, filename, format="nexus")
}
\arguments{
  \item{tree}{a phylogentic tree, with optional annotations}
  \item{filename}{a file to write no}
  \item{format}{a format specifier; either \code{"nexus"} or \code{"newick"}}
}
\description{
  This function takes the given tree object and returns a string
  representing the tree in the requested format.
}
\details{
  The tree object should be either a \code{"phylo"} or
  \code{"multiPhylo"} object. It may optionally be augmented with
  annotations, as described in the documentation for the
  \code{\link{parse_annotated}} function.
}
\author{Daniel Beer \email{dlbeer@gmail.com}}
\seealso{
  \code{\link{parse_annotated}},
  \code{\link{print_annotated}},
  \code{\link{read_annotated}}
}
\keyword{phylo}
