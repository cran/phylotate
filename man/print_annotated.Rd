\name{print_annotated}
\alias{print_annotated}
\title{Serialize an annotated phylogenetic tree}
\usage{
print_annotated(tree, format="nexus")
}
\arguments{
  \item{tree}{a phylogentic tree, with optional annotations}
  \item{format}{a format specifier; either \code{"nexus"},
    \code{"newick"}, or \code{"newick.named"}}
}
\description{
  This function takes the given tree object and returns a string
  representing the tree in the requested format. The difference between
  the \code{"newick"} and \code{"newick.named"} formats is that the
  former uses only node numbers in its output, whereas the latter uses
  the tip labels (sanitized and deduplicated if necessary).
}
\value{
  a string containing a serialized tree.
}
\details{
  The tree object should be either a \code{"phylo"} or
  \code{"multiPhylo"} object. It may optionally be augmented with
  annotations, as described in the documentation for the
  \code{\link{parse_annotated}} function.

  The output is a string suitable for writing to a file.
}
\author{Daniel Beer \email{dlbeer@gmail.com}}
\seealso{
  \code{\link{parse_annotated}},
  \code{\link{read_annotated}},
  \code{\link{write_annotated}}
}
\keyword{phylo}
