# Reopening a second square rail: m=19·2^v

The previous reset permits widths r=v+4 just as well as r=v+3,
with a=3·2^v. Put s=2v+2 and m=19·2^v. Then the same transitive
carry gives counts 4m²−3, 4m²−1, 4m²+1. Two new rooted seeds give
the following exact continuation at every v≥2:

| Offset B | A | h,g,R,S | Outgoing rows (root 8) |
|--:|--:|:--|:--|
| +3 | 1444 | 1343,311,774,981 | (266,68,25,306,291,135,317,351,38) |
| +5 | 722 | 663,177,393,506 | (266,276,33,22,69,283,431,319,28) |

The [rooted formula](BANK.md) gives H=B+A·2^ell. Use
ell=s−2 and ell=s−1, respectively, to reach 361·2^s+B.
This yields five consecutive odds on a second all-height rail.

The reopen can again retain all old vertices and preserve every
arc incident to a long common spacer of length s−2. In the reset,
the maximal unmarked middle interval has length s, with five old
tape vertices and the head beating it, and four tape vertices
and the sentinel losing to it. Keep the first s−2 interval
vertices as active near-sinks. Assign the six old predecessors
to nonroot core vertices and add the other two. Of the old
successors pick one as the root and use the rest as neutral
suffixes, together with the last two old spacer vertices.
All remaining old changes can be confined to the 13 old
vertices outside the retained interval. This is a bounded
local rewrite; the number of flips within that boundary is
not optimized. The +5 reader promotes one neutral successor
to be one more active near-sink after the shared interval.

No selection of these two coefficients (121 or 361) as a
function of an arbitrary requested count is supplied. It does
show the reader change is not limited to the one square of 11.
