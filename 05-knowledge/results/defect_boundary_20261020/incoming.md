# Incoming boundary work and the 105 matching dictionary

**PROVED:** the explicit matching dictionary and fixed-budget negligibility
below. **FINITE-EXACT:** the 105 labelled dictionary checks.
The mid-session incoming commit was `fe848abe68`, *exact strip audit and
metagraph degree asymptotics*. It supplied a useful connection rather than
a replacement for this packet's counting mechanism.

## 1. What was actually imported

The incoming [degree-asymptotics note](../boundary_support_20261020_asymptotics.md)
defines the complement-line multigraph with loops counted twice and retains
parallel lines. It derives ordinary class-count asymptotics by bounding
nonidentity permutation fixed sets, and distinguishes maximum fibre degree
from maximum Hamiltonian-path count.

Its external Hamiltonian-path upper bound is a **CITED** input to that note;
this packet does not use it as a new dependency. The elementary loss
distinction is relevant here: a count of nodes, a total number of endpoints,
and the actual matching of endpoints are three different objects. Our
one-chord family supplies specific leaf incidences; our class compiler supplies
node counts; neither silently supplies all metagraph wiring.

The incoming [research board](../boundary_support_20261020_board.md) also
recovers ordered-port colouring supports and the 105-matching atlas. Its
proposed full metagraph companion was not present in the fetched commit,
so no unprovided wiring or bridge result is assumed.

The checked strip example is particularly useful as a hostile: on its
three-state blocks the second and fourth powers have diagonal/off-diagonal
weights `(2,1)` and `(6,5)`. Their nonzero supports agree, but their weights
are not proportional. Accordingly a support grammar can preserve existence
without preserving the counts needed here.

## 2. The same matching object, with an explicit map

At feedback budget four, the leading sector of our compiler consists of
perfect matchings on eight **ordered defect endpoints**. It has 105 members.
For each matching, put a heavy vertex of capped weight five in every one
of the nine gaps, and give each endpoint weight one. This gives a pair-free
core of order 17 and total capped weight 53. Distinct matchings give distinct
capped weighted-core classes.

The [October 18 atlas classification](../hamming_matching_atlases_20261018e.md)
identifies optimal four-chart partitions of the code `H8+I2` with perfect
matchings on eight slopes in `F2^3`. Fix the map from our endpoint positions
`0,...,7` to the binary vectors of `F2^3`. Send a defect pair to the chart
whose two slopes are those endpoint vectors. Four matching pairs give the
four disjoint charts.

Here the slope `l` has four codewords

`f(X)=s+l dot X` on the eight physical positions, followed by `(t,t)`,
with `s,t in F2`.

The decoder of a chart reads its two slopes, giving the inverse map.
Thus this is a literal bijection between structured objects, not a pairing
of two lists merely because both have size 105. Pair differences split the
templates into `7+42+56`, according as there are one, two, or four distinct
nonzero XOR differences. The inherited atlas proof establishes the same
classification; [incoming_check.py](incoming_check.py) independently forms all
charts and partitions from their 32 codewords.

Preserved: the endpoint pairing, disjointness, pair differences under the
chosen binary address, and the full matching. Lost by forgetting the address:
the interpretation of the XOR direction classes. Required sidecars:
the total endpoint order, the binary identification and the code's physical
coordinate basis.

This is **not** an equivariant quotient of tournament classes by the full
affine group on the eight slopes. The nine heavy vertices impose an intrinsic
total order, so the 105 tournament templates remain distinct. Even though
permutations act transitively on abstract perfect matchings, quotienting by
that action would destroy the ordered tournament data.

This dictionary suggests a controlled next question: express low-degree
correction operations—omitting or shortening a gap, or sharing an endpoint—
on the marked matching/atlas side. The A/B/C classification in
[next_coefficient.md](next_coefficient.md) supplies exact test cases.
No inherited geometric realization or unit-distance interpretation is
transported by this map.

## 3. The sparse grammar is necessarily a thin part of A000568

The full number `a_n` of tournament classes is at least
`2^binom(n,2)/n!`, since each class contains at most `n!` labelled tournaments.
For fixed `k`, our polynomial theorem therefore gives

\[
0\leq\frac{a_{n,\leq k}}{a_n}
\leq O_k(n^{2k})\,n!\,2^{-\binom n2}\longrightarrow0.
\]

No external asymptotic theorem is needed for this conclusion. It identifies
the remaining scale problem precisely: an efficient full count must let the
defect budget grow or use a complementary dense grammar. The live divisor
transfer is one such complementary exact representation; it compresses
permutation-cycle histories rather than sparse ordered templates.

This is why the two approaches are useful together. A fixed sparse grammar
gives structural formulas, automorphism support bounds and exact local
families. Future-profile elimination handles the full counting sum. A
uniform transition law between defect budgets is still missing.
