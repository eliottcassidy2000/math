# Independent binary-clock span audit

2026-11-05. **PASS after two local wording corrections.**

Final audited proof hash:
`8e7a7fca634179d970795b08fc6a79a21390a2c037d50fa05dbb1c82fd7c8a01`.
Verifier hash:
`63f6a8032a6359b2ab35122858a88037bd09e7175dd4eea47056111eae95d578`.
Validation hash:
`e0f01233a44384948cd3c1754a9972dbfaecf9f8c0e9def8168ef4b87622de80`.

## Corrections made before signoff

In the original root-before-a case, the statements e<=2 and l>=n-1 were
false: ordering r,a,all B gives first-B position 3, while all B,r,a gives
last-B position n-2. The producer replaced them with the exact counting
fact actually needed: n-2 B vertices occupy distinct positions, so
`l-e>=n-3`. The theorem and subsequent inequality are unchanged.

The fixed-width maximum-count statement now explicitly assumes a positive
integer width w>=1. Width zero has no member of this n>=3 cyclic family,
so it has no largest clock count.

## All-order proof checks

When r is not last, its unique outneighbor a must follow immediately.
Choosing an arbitrary subset of B for the initial increasing run determines
one path, and the complementary increasing run follows a. These are
exactly 2^(n-2) paths. When r is last, a must be first and the remaining
path is unique. Thus H=2^(n-2)+1 and the returning scalar step is 2H-1.

For an arbitrary vertex order with a before r, the three forced backward
gaps sum to n-1, so 3w>=n-1. Zero boundary gaps are harmless. With r before
a, the inequalities `w>=p-e` and `w>=l-q` remain valid even when either
right side is negative. Their sum is at least
`(p-q)+(l-e)>=1+(n-3)=n-2`. The directed triangle ensures w>=1, which gives
`n-1<=2w+1<=3w`. This argument allows every permutation of B.

The balanced construction has two nonnegative exterior gaps and one
positive middle gap, all at most ceil((n-1)/3), summing to n-1. Its B
vertices are put in transitive order, so no internal B edge is backward.
Every possible backward edge is one of the three listed kinds. At least
one gap reaches the ceiling, since three smaller integer gaps could not
sum to n-1. The minimum span is therefore exactly ceil((n-1)/3).

For positive width w this is equivalent to n<=3w+1, yielding the stated
largest family count. It is a graph-family obstruction, not a count
obstruction or an obstruction to replacing the graph.

## Independent computations

`audit.cpp` independently constructs the tournament and inspects every
ordering through order 10, a total of 4,037,910 orderings. It checks each
lower-bound inequality in the appropriate position case, not only the
minimum. Every minimum agrees with the claimed ceiling.

It also independently minimizes over every pair of a/r positions through
order 500. This second exact method is valid because sorting B preserves
all exterior incidences and removes, rather than adds, internal backward
edges. The balanced gap construction is checked through order 10,000.

Literal permutation counting confirms the clock counts through order 10.
The alternative word 022333332 is converted to adjacency by the explicit
distance-two bit rule, independently of the source graph helper. Literal
enumeration of all 9! vertex orders gives exactly 513 Hamiltonian paths,
and its displayed order has backward span exactly two. Meanwhile C_11
has 513 paths and minimum span four by the proved formula. These graphs
are distinct even in order (nine versus eleven).

The copied source verifier also passes, including its subset-DP path
counts through order 12. No remaining mathematical or implementation
issue was found. Results are in `AUDIT.json` and `source_verification.log`.
Reproduce with `g++ -O3 -std=c++17 audit.cpp -o audit && ./audit`.
No author source files were edited by this auditor.
