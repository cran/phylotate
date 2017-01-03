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

format.parsers <- list(
    "newick" = parse.newick,
    "nexus"  = parse.nexus
)

format.printers <- list(
    "newick" = print.newick,
    "nexus"  = print.nexus
)

parse_annotated <- function (str, format="nexus") {
    p <- format.parsers[[format]]
    if (is.null(p)) {
	stop(sprintf("unknown format: %s", format))
    }

    p(tokenize.nexus(str))
}

print_annotated <- function (tree, format="nexus") {
    p <- format.printers[[format]]

    if (is.null(p)) {
	stop(sprintf("unknown format: %s", format))
    }

    out <- NULL
    p(tree, function (t) { out[length(out) + 1] <<- t })
    paste0(out, collapse="")
}

read_annotated <- function (filename, format="nexus") {
    parse_annotated(readChar(filename, file.info(filename)$size), format)
}

write_annotated <- function (tree, filename, format="nexus") {
    writeChar(print_annotated(tree, format), filename, eos=NULL)
}
