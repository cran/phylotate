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
# Lexical analyzer for NEXUS file data
#
# Input: a text string
# Output: a token table. Each row is a token. Columns are:
#   * pos: index in the original string
#   * row: row of the token in the original file
#   * col: column of the token in the original file
#   * token: token symbol
#   * text: token text (only for words and comments)
#   * comment: comment following token, if any
#
# Token symbols are single-character strings. Either a symbol like ";"
# or "W" (for a keyword).
########################################################################

tokenize.nexus <- function (stext.in) {
    # Trying to disassemble and reassemble strings using substr() is
    # extremely slow. It's much faster to split the string into
    # individual characters and calculate character classes in a bulk
    # operation.
    stext <- strsplit(stext.in, "")[[1]]

    # Character classifications
    is.whitespace <- grepl("[\t\r\n ]", stext)
    is.word <- grepl("[-+A-Za-z0-9_\\.#?]", stext)
    is.sol <- c(TRUE, (stext == "\n")[1:(length(stext)-1)])

    # Row/column decode
    row <- cumsum(is.sol)
    col <- c(1:length(stext)) - which(is.sol)[row] + 1

    # Skip-ahead maps
    whitespace.end <- index.of.next(!is.whitespace)
    comment.end <- index.of.next(stext == "]")
    word.end <- index.of.next(!is.word)

    # Tokenizer loop
    token <- rep(NA, length(stext))
    text <- rep(NA, length(stext))
    pos <- rep(NA, length(stext))
    comment <- rep(NA, length(stext))
    i <- 1
    c <- 0
    while (i <= length(stext)) {
	i <- whitespace.end[i]
	if (is.na(i)) { break }

	if (stext[i] == "[") {
	    j <- comment.end[i]

	    if (is.na(j)) {
		stop(sprintf("Line %d, column %d: unterminated comment",
			     row[i], col[i]))
	    }

	    if ((c > 0) && is.na(comment[c])) {
		comment[c] <- paste0(stext[(i+1) : (j-1)], collapse="")
	    }

	    i <- j + 1
	} else if (is.word[i]) {
	    j <- word.end[i]
	    if (is.na(j)) { j <- length(stext) + 1 }
	    c <- c + 1
	    pos[c] <- i
	    token[c] <- "W"
	    text[c] <- paste0(stext[i : (j-1)], collapse="")
	    i <- j
	} else {
	    c <- c + 1
	    pos[c] <- i
	    token[c] <- stext[i]
	    i <- i + 1
	}
    }

    # Trim vectors
    token <- token[1:c]
    text <- text[1:c]
    pos <- pos[1:c]
    comment <- comment[1:c]

    data.frame(
	pos = pos,
	row = row[pos],
	col = col[pos],
	token = token,
	text = text,
	comment = comment,
	stringsAsFactors = FALSE
    )
}

token.error <- function (tok, i, msg) {
    stop(sprintf("Line %d, column %d: %s", tok$row[i], tok$col[i], msg))
}
