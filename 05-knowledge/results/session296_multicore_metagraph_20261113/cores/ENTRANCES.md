# Independently checked ready entrances

A word (possibly empty) describes the transitive tail of a
one-start funnel. Start from a rooted singleton. To apply the
word to a rooted core U, let U dominate the new transitive tail.
A fresh root q beats exactly the old root inside U, and beats
exactly the zero bits of the tail. An empty word is used only on
a nontrivial core; a nonempty word ends in 1. These guards make
the new tournament strong. The first SCC of the deletion is U
and the tail vertices are later singleton SCCs; thus a full
path starting at q must enter at the old root, take its unique
start path, and traverse the tail. The unique-start guard
propagates at every nested funnel.

After the listed words append one vertex z beaten only by the
last root q. The new root is z, with sole predecessor q. The
identity delta_z=s_q gives readiness one. Independently count
all paths, rather than the unique marked sector, by standard
subset/last-vertex DP.

| k | 405·2^k+71 | words from singleton (ε empty) |
|---:|---:|---|
|1|881|11,1,11|
|2|1691|1011,101|
|3|3311|1,ε,01,ε,1|
|4|6551|101,ε,01101|
|5|13031|010101111011|
|6|25991|011,ε,0011001|
|7|51911|1,ε,11,ε,ε,1,ε,ε|
|8|103751|1,001111,1,ε|
|9|207431|11,ε,ε,00011,1|

The k=1 chain has more entries:

| H | words from singleton |
|---:|---|
|2645|111,11,1|
|7937|01,01111101|
|23813|1011,ε,ε,ε,ε,1|
|71441|000011010110001|

A second ready 881 is 1,ε,00101. A second ready 71441 is
1,1100001110001. These comparisons exhibit distinct endpoint
responses even when H and delta coincide: the two ready 71441
tournaments have deletion counts 35249 and 16869 respectively.
The last ready source reaches 214325 by the standard two births;
the return at 214325 has not been supplied.

Each new target is attained twice: first by the raw two births
on the preceding ready source, then by the stated same-count
ready replacement. A precise named replacement is to identify
the two new roots and their predecessors and then align the
remaining old/new names in increasing construction order until
one side is exhausted. Delete the unmatched old names, insert
the unmatched new names, and set each nonroot arc to the target
template. Every retained root-star arc agrees, because the root
is beaten only by its named predecessor on both sides. This is
a legal finite graph surgery, not an induced extension, a short
native complement certificate, or an assertion of count
preservation at intermediate individual edits.

The addresses are finite positive certificates, not a promise
that all larger k or all deeper targets admit analogous
funnel words. These exits already cross the mod-45 barrier of
the old one-template recursion.
