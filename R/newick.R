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
# Newick parser
########################################################################

parse.newick <- function (tok) {
    # Strip commas and semicolons (we don't need them any more)
    tok <- tok[(tok$token != ",") & (tok$token != ";"),]

    # Attach distances to their preceeding tokens
    is.dist <- (tok$token == ":") & c((tok$token == "W")[2:nrow(tok)], FALSE)
    dist <- which(is.dist)
    tok$distance <- rep(NA, nrow(tok))
    tok$distance.comment <- rep(NA, nrow(tok))
    tok$distance[dist - 1] <- as.numeric(tok$text[dist + 1])
    tok$distance.comment[dist - 1] <- tok$comment[dist + 1]
    tok <- tok[!(is.dist | c(FALSE, is.dist)[1:nrow(tok)]),]

    # Count nodes
    Ntip <- sum(tok$token == "W")
    Nnode <- sum(tok$token == "(")
    Ntotal <- Nnode + Ntip

    # Assign node numbers. Preserve original numbering on tips if
    # possible. Node numbers for leaves are attached to W tokens. Node
    # numbers for non-leaves are attached to ( tokens. The remaining
    # type of token, ), has no node number attached.
    m <- tok$token == "W"
    m.tn <- suppressWarnings(as.integer(tok$text[m]))
    if (sum(duplicated(m.tn)) || (min(m.tn) < 1) || (max(m.tn) > Ntip)) {
	warning("renumbering non-contiguous tips")
	m.tn <- 1:Ntip
    }
    tok$node <- rep(NA, nrow(tok))
    tok$node[m] <- m.tn
    tok$node[tok$token == "("] <- 1:Nnode + Ntip

    # Make a map of node positions
    m <- !is.na(tok$node)
    node.index <- rep(NA, Ntotal)
    node.index[tok$node[m]] <- (1:nrow(tok))[m]

    # Assign parentage and copy attributes from ) tokens (distance,
    # etc.) to the matching ( tokens.
    node.parent <- rep(NA, Ntotal)
    p <- NA
    for (i in 1:nrow(tok)) {
	if (tok$token[i] == "(") {
	    node.parent[tok$node[i]] <- p
	    p <- tok$node[i]
	} else if (tok$token[i] == ")") {
	    if (is.na(p)) { token.error(tok, i, "unbalanced parentheses") }
	    cols <- c("distance", "comment", "distance.comment")
	    tok[node.index[p],cols] <- tok[i,cols]
	    p <- node.parent[p]
	} else if (tok$token[i] == "W") {
	    node.parent[tok$node[i]] <- p
	} else {
	    token.error(tok, i, "unexpected token")
	}
    }

    if (!is.na(p)) { token.error(tok, nrow(tok), "unbalanced parentheses") }

    # List of children in pre-order
    po <- tok$node[!is.na(tok$node)]
    po.children <- po[!is.na(node.parent[po])]

    # Make a phylo object
    p <- list(
	edge			= matrix(as.integer(c(node.parent[po.children],
						      po.children)),
					 nrow = length(po.children), ncol = 2),
	edge.length		= tok$distance[node.index[po.children]],
	Nnode			= Nnode,
	tip.label		= tok$text[node.index[1:Ntip]],
	node.comment		= tok$comment[node.index[1:Ntotal]],
	node.distance.comment	= tok$distance.comment[node.index[1:Ntotal]]
    )

    class(p) <- "phylo"
    attr(p, "order") <- "cladewise"
    p
}

########################################################################
# Newick printer
########################################################################

print.newick <- function (phy, printer, Ntip=NA) {
    if (class(phy) != "phylo") {
	stop("print.newick requires a phylo object")
    }

    if (is.na(Ntip)) {
	if (is.null(phy$tip.label)) {
	    stop("no tip labels found: you need to specify Ntip")
	}

	Ntip <- length(phy$tip.label)
    }

    Ntotal <- Ntip + phy$Nnode

    # Map: node -> distance from parent
    node.distance <- rep(NA, Ntotal)
    if (!is.null(phy$edge.length)) {
	node.distance[phy$edge[,2]] <- phy$edge.length
    }

    # Map: node -> parent
    node.parent <- rep(NA, Ntotal)
    node.parent[phy$edge[,2]] <- phy$edge[,1]

    # Nodes ordered by parentage
    m <- order(node.parent)
    m.parent <- node.parent[m]

    # Map: parent -> first child
    f <- m.parent != c(NA, m.parent[1:(length(m) - 1)])
    f[is.na(f)] <- TRUE
    f[is.na(m.parent)] <- FALSE
    node.first.child <- rep(NA, Ntotal)
    node.first.child[m.parent[f]] <- m[f]

    # Map: node -> next sibling
    f <- m.parent == c(m.parent[2:length(m)], NA)
    f[is.na(f)] <- FALSE
    node.sibling <- rep(NA, Ntotal)
    node.sibling[m[f]] <- m[c(FALSE, f[1:(length(f) - 1)])]

    # Select root
    n <- m[is.na(m.parent)]
    if (length(n) != 1) {
	stop("Tree must contain exactly one root node")
    }

    # Emit attributes for a given node
    finish.node <- function (n) {
	if (!is.null(phy$node.comment)) {
	    c <- phy$node.comment[n]
	    if (!is.na(c)) { printer(sprintf("[%s]", c)) }
	}

	c <- node.distance[n]
	if (!is.na(c)) {
	    printer(sprintf(":%g", c))
	    if (!is.null(phy$node.distance.comment)) {
		c <- phy$node.distance.comment[n]
		if (!is.na(c)) { printer(sprintf("[%s]", c)) }
	    }
	}
    }

    # Traverse tree
    depth <- 0
    p <- NA
    while (!(is.na(n) && is.na(p))) {
	if (is.na(n)) {
	    printer(")")
	    finish.node(p)
	    n <- node.sibling[p]
	    p <- node.parent[p]
	    if (!is.na(n)) { printer(",") }
	    depth <- depth - 1
	} else if (n <= Ntip) {
	    printer(sprintf("%d", n))
	    finish.node(n)
	    n <- node.sibling[n]
	    if (!is.na(n)) { printer(",") }
	} else {
	    printer("(")
	    p <- n
	    n <- node.first.child[n]
	    depth <- depth + 1
	    if (depth > Ntotal) { stop("loop in tree") }
	}
    }
}
