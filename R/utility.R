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
# Given a boolean vector, return an integer vector of the same size
# where each element maps its index to the location of the previous/next
# TRUE item in the original vector. For example:
#
#    Index:      1     2     3     4      5     6     7     8     9
#    Input:  FALSE FALSE FALSE  TRUE  FALSE FALSE  TRUE FALSE FALSE
#    Prev:      NA    NA    NA     4      4     4     7     7     7
#    Next:       4     4     4     4      7     7     7    NA    NA
########################################################################

index.of.prev <- function (m) { c(NA, which(m))[cumsum(m) + 1] }
index.of.next <- function (m) { length(m) - rev(index.of.prev(rev(m))) + 1 }
