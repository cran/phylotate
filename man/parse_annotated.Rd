\name{parse_annotated}
\alias{parse_annotated}
\title{Parse an annotated phylogenetic tree}
\usage{
parse_annotated(str, format="nexus")
}
\arguments{
  \item{str}{a text string, containing tree data}
  \item{format}{a format specifier; either \code{"nexus"} or \code{"newick"}}
}
\description{
  This function takes the given text string, containing data in either
  NEXUS or Newick format, and returns annotated phylogenetic trees.
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

  In addition to \code{edge}, \code{edge.length} and \code{tip.label},
  two additional vectors are added. These are \code{node.comment} and
  \code{node.distance.comment}. These contain annotations associated
  with nodes and their distance values. These arrays are indexed by node
  number, not by edge. The reason for this is that this ensures that the
  object will remain in a valid state after a call to \code{reorder}
  which might change the ordering of the edge arrays without being aware
  of annotations. If you need to obtain annotations in edge-order,
  subset by the second column of the \code{edge} array.
}
\references{
  Paradis, E. Definition of Formats for Coding Phylogenetic Trees in R.
  \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}
}
\author{Daniel Beer \email{dlbeer@gmail.com}}
\examples{
  # Parse the example data included with this package
  data(finches)
  t <- parse_annotated(finches, format="nexus")

  # Obtain annotations in edge-order, rather than node-order
  edge.comment <- t$node.comment[t$edge[,2]]
}
\seealso{
  \code{\link{print_annotated}},
  \code{\link{read_annotated}},
  \code{\link{write_annotated}},
  \code{\link{finches}}
}
\keyword{phylo}
