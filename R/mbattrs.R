# phylotate
# Copyright (c) 2017 Daniel Beer
# Copyright (c) 2017, 2018 Anusha Beer
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

# Given an annotated tree produced by MrBayes, return a dataframe
mb_attrs <- function (tree) {
    lrange <- gsub("^.*[&,]length_95%HPD=\\{([^}]*)\\}.*$", "\\1",
                   tree$node.distance.comment)
    data.frame(
        prob = as.numeric(gsub(
                "^.*[&,]prob=([^,]*).*$", "\\1",
                tree$node.comment)),
        prob_stddev = as.numeric(gsub(
                "^.*[&,]prob_stddev=([^,]*).*$", "\\1",
                tree$node.comment)),
        length_mean = as.numeric(gsub(
                "^.*[&,]length_mean=([^,]*).*$", "\\1",
                tree$node.distance.comment)),
        length_median = as.numeric(gsub(
                "^.*[&,]length_median=([^,]*).*$", "\\1",
                tree$node.distance.comment)),
        length_95_HPD_low = gsub(",.*$", "", lrange),
        length_95_HPD_high = gsub("^[^,]*,", "", lrange)
    )
}
