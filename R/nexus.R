# phylotate
# Copyright (c) 2017 Daniel Beer
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

########################################################################
# NEXUS parser
########################################################################

parse.nexus <- function (tok) {
    if ((length(tok) < 1) || (tok$token[1] != "W") ||
	(toupper(tok$text[1]) != "#NEXUS")) {
	stop("missing NEXUS header")
    }

    end.of.line <- index.of.next(tok$token == ";")
    output <- list()
    trees <- list()
    id.map <- NA

    i <- 2
    while (i <= nrow(tok)) {
	if (tok$token[i] != "W") { token.error(tok, i, "missing directive") }
	dir <- toupper(tok$text[i])

	j <- end.of.line[i]
	if (is.na(j)) { token.error(tok, i, "missing semicolon") }

	i <- i + 1

	if (dir == "BEGIN") {
	    if (tok$token[i] != "W") {
		token.error(tok, i, "missing block name")
	    }

	    id.map <- NULL
	    trees <- list()
	} else if (dir == "TRANSLATE") {
	    tab <- tok[i:(j-1),]
	    tab <- tab[tab$tok != ",",]

	    if (nrow(tab) %% 2 != 0) {
		token.error(tok, i, "expected an even number of items")
	    }

	    if (sum(tab$tok != "W") > 0) {
		token.error(tok, i, "invalid tokens in list")
	    }

	    tab <- tab$text
	    n <- (1:(length(tab) / 2)) * 2
	    id.map <- tab[n]
	    names(id.map) <- tab[n - 1]
	} else if (dir == "TREE") {
	    if (tok$token[i] == "*") {
		i <- i + 1
	    }

	    if ((tok$token[i] != "W") || (tok$token[i+1] != "=")) {
		token.error(tok, i, "expected name and equals sign")
	    }
	    tree.name <- tok$text[i]
	    i <- i + 2

	    trees[[length(trees) + 1]] <- parse.newick(tok[i:(j-1),])
	    names(trees)[[length(trees)]] <- tree.name
	} else if (dir == "END") {
	    if (!is.null(id.map)) {
		trees <- lapply(trees,
		    function (t) {
			t$tip.label <- id.map[t$tip.label]
			names(t$tip.label) <- NULL
			t
		    })
	    }

	    a <- length(output)
	    output[1:length(trees) + a] <- trees
	    names(output)[1:length(trees) + a] <- names(trees)

	    id.map <- NULL
	    trees <- list()
	}

	i <- j + 1
    }

    # Create multiPhylo object or simplify
    if (length(output) == 1) {
	output <- output[[1]]
    } else {
	tl <- NULL

	output <- lapply(output, function (t) {
	    if (is.null(tl)) {
		tl <<- t$tip.label
	    } else {
		if ((length(t$tip.label) != length(tl)) ||
		    (sum(t$tip.label != tl) > 0)) {
		    stop("inconsistent tip labels")
		}
	    }

	    t$tip.label <- NULL
	    t
	})

	class(output) <- "multiPhylo"
	attr(output, "TipLabel") <- tl
    }

    output
}

########################################################################
# NEXUS printer
########################################################################

print.nexus <- function (phy, printer) {
    if (class(phy) == "phylo") {
	tip.label <- phy$tip.label
	phy <- list(phy)
    } else if (class(phy) == "multiPhylo") {
	tip.label <- attr(phy, "TipLabel")
    } else {
	stop("To make a NEXUS file, we require a phylo or multiPhylo object")
    }

    tree.names <- names(phy)
    if (is.null(tree.names)) { tree.names <- rep("untitled", length(phy)) }
    tree.names[is.na(tree.names)] <- "untitled"

    printer("#NEXUS\n")
    printer("[R-package phylotate]\n")

    # Taxa section
    printer("begin taxa;\n")
    printer(sprintf("\tdimensions ntax = %d;\n", length(tip.label)))
    printer("\ttaxlabels\n")
    for (t in tip.label) {
	printer(sprintf("\t\t%s\n", t))
    }
    printer("\t;\n")
    printer("end;\n")

    # Trees section
    printer("begin trees;\n")
    printer("\ttranslate\n")
    for (i in 1:length(tip.label)) {
	printer(sprintf("\t\t%d\t%s", i, tip.label[i]))
	if (i < length(tip.label)) { printer(",") }
	printer("\n")
    }
    printer("\t;\n")
    for (i in 1:length(phy)) {
	printer(sprintf("\ttree %s = [&U] ", tree.names[i]))
	print.newick(phy[[i]], printer, Ntip=length(tip.label))
	printer(";\n")
    }
    printer("end;\n")
}
