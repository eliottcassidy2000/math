# A concrete growth grammar and its tested branching boundary

**FINITE-EXACT through order11; all-order definition and response laws PROVED.**
This is a restricted constructor census, not an extension of the repository's
existing ordinary-H coverage through99,999,999.

Start at one named vertex0. At stage z, choose one old port p_z in0,...,z-1,
and append z beating precisely p_z and losing to every other old vertex.
Thus the backward arcs in creation order form the recursive rooted tree
z->p_z. A word(p_1,...,p_(n-1)) reconstructs every arc of the actual
tournament. Different words can yield isomorphic tournaments; the universe
has exactly(n-1)! presentations, not that many isomorphism classes.

Every presentation through order11 was counted by forward endpoint subset
DP. The cumulative distinct H counts and first missing allowed odd count are:

| Maximum order | Presentations at that order | Cumulative H values | First missing allowed odd | Largest H |
| ---: | ---: | ---: | ---: | ---: |
| 5 | 24 | 7 | 17 | 15 |
| 6 | 120 | 16 | 35 | 37 |
| 7 | 720 | 55 | 63 | 123 |
| 8 | 5040 | 192 | 311 | 413 |
| 9 | 40320 | 706 | 991 | 1509 |
| 10 | 362880 | 2816 | 4585 | 6083 |
| 11 | 3628800 | 11888 | 20407 | 25827 |

Consequently every positive odd h<=20405 other than7 and21 has a saved
literal constructor word in this narrow grammar. All11888 saved witnesses
were independently checked with a top-down starting-endpoint recurrence;
Python separately checks a deterministic sample and the critical endpoints.
The normal and optimized Python receipts agree. The two independent full
scans through10 and11 give identical earlier layers.

The smallest word retained is of minimum order **within this grammar and
this complete cutoff**: the scan proceeds by increasing order. It is not
necessarily a minimum-order tournament among all constructions. The result
does not prove that the growing coverage windows continue to overlap.

## Why choosing a port is more than choosing a binary digit

The complete labelled source census through order6, reduced to ordinary
isomorphism representatives, tests whether one-port births can realize both
2H-1 and2H+1 from the SAME source. The numbers of such classes at orders
2,...,6 are1,1,2,3,4. At order6 their H values are1,5,9,37.

The two ordinary six-vertex H17 classes illustrate the missing joint
condition. In lexicographic edge-mask encoding:

    mask16: child counts33,27,41,59,85,47;
    mask83: child counts49,29,31,35,77,71.

One supplies2*17-1, the other2*17+1; neither supplies both. A scalar
H17 label loses exactly which extension is available. This is not a
nonexistence theorem for larger H17 representatives or other operations.

The inherited [binary-growth obstruction](../session142_tournament_odd_spectrum_20261102/constructions/BINARY_NO_GO.md)
concerns a more restricted choice of the current or freshly added mark.
Here all old ports may be selected. The new finite counts do not contradict
that all-depth obstruction, and they do not establish binary closure here.

## What an inductive proof would still need

An actual family of source states must be proved closed under enough
target-producing moves to cover both allowed prime lanes. Existence of a
minus branch somewhere and a plus branch somewhere at one finite order
does not give a branch-closed family or an induction. The growth-faithfulness
theorem shows why a universal quotient for all named ports cannot simply
forget the source orientation. A restricted sufficient construction language
could still have a much smaller response state, as the current singleton/D
grammar does. Its all-target entry theorem remains the decisive obligation.
