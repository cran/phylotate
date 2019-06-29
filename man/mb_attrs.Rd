\name{mb_attrs}
\alias{mb_attrs}
\title{Parse MrBayes-supplied attributes from a NEXUS file}
\usage{
mb_attrs(tree)
}
\arguments{
  \item{tree}{an object of type \code{"phylo"}}
}
\description{
  This function takes a tree object and produces a dataframe containing
  attributes attached to each node by MrBayes.
}
\value{
  A dataframe of attributes.
}
\details{
  The returned dataframe contains one row per node, and one column per
  attribute. The attributes parsed are \code{prob}, \code{prob_stddev},
  \code{length_mean}, \code{length_median}, \code{length_95_HPD_low},
  and \code{length_95_HPD_high}.

  Attributes which are derivable from the others are not parsed (for example,
  the \code{prob_percent} attribute is not parsed, since it's \code{prob}
  times 100.
}
\author{Anusha Beer \email{anbeer29@gmail.com}}
\seealso{
  \code{\link{parse_annotated}},
  \code{\link{read_annotated}}
}
\examples{
  # Parse the example data included with this package
  data(finches)
  t <- parse_annotated(finches, format="nexus")

  # Obtain a table of MrBayes attributes for each node
  attrs <- mb_attrs(t)
}
\keyword{phylo}
