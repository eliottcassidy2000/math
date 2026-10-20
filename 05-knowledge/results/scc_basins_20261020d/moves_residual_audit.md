# Independent audit of the alternating-score exclusion

**PASS: proof audit and independent complete certificate replay.**
[moves.md, section 11](moves.md) proves the stated nonregular self-converse
Hamiltonian-cycle theorem relative to its explicitly CITED Thomassen input.
Together with its previously established regular and nonstrong cases, the
claimed consequence is valid: every blue line at a mixed singleton-padded
source `1 -> C -> 1` has a black detour of length at most five. This does not
prove global black-core connectivity or cover every larger outer SCC.

## Connectivity reductions

If every Hamiltonian cycle has only complementary-score edges, following
one cycle makes the scores alternate. Nonregularity forces even order
`n=2m`, two distinct complementary scores, and equal low/high classes.
An edge within the low class is unbalanced. The published 3-strong arc
theorem would put it on a Hamiltonian cycle, so connectivity is at most two.
This audit uses exactly that cited statement and does not claim to audit
Thomassen's original proof.

The cut-vertex exclusion is also correct. For a minimum-score deleted
vertex v, the predecessor u of v in any Hamiltonian cycle lies in the
last SCC B of `T-v`. If `t=|B|>1`, u has at most `t-2` internal outneighbors
(B is strong), plus its arc to v, hence `d(u)<=t-1`. A vertex in an earlier
component beats all t vertices of B, so the maximum score is at least t.
Thus the closing edge uv is unbalanced. If t=1, u has degree one, below
the maximum score at every relevant order. The maximum-score deletion
case is its converse. Every vertex is extremal, so no cut vertex exists.
The only putative obstruction therefore has connectivity exactly two.

## The two-vertex cut

For a cut `T-S=A -> B`, removing the two vertices of S from an alternating
Hamiltonian cycle leaves at most two directed paths. Each path's A vertices
are an initial segment. Each alternating segment has low/high imbalance at
most one, giving `|H(A)-L(A)|<=2`. Empty segments cause no exception.

For score gap one, both the outdegrees on A and the indegrees on B are at
most m. Summation gives `2m>=a-1+2b` and `2m>=b-1+2a`; since `a+b=2m-2`,
these imply `4m>=6m-8`, excluding m>=5.

For gap at least three, a low vertex in A and a high vertex in B would
give `delta>=b` and `Delta<=b+1`, an immediate contradiction. Consequently
A is entirely high or B entirely low; taking converse makes B entirely
low in either case. If A is entirely high, its at least `m-2>=3` vertices
already violate the imbalance bound. Otherwise A contains a low vertex,
so `b<=delta`. Averaging B's outdegrees gives `b>=2delta-3`, hence
`delta<=3`. Strongness requires an arc from the full low class to the full
high class, so `m delta>m(m-1)/2`. For m>=5 this forces delta=3,
m in {5,6}, and b=3. A low vertex of A spends its entire outdegree on B;
there can be at most one such vertex. With `t=|L intersect S|`, therefore,
`t>=m-4>=1`, whereas the imbalance in A equals `1+2t>=3`. This is the
required contradiction. All inequalities and the converse exchange agree.

## Completeness and certificate replay

The blue-orbit argument is sound. Redei parity makes H(T) odd. The
automorphism group has odd order (an involutory tournament automorphism
cannot transpose the endpoints of an arc), and acts freely on ordered
Hamiltonian paths; hence the path-orbit count is odd. For an
anti-automorphism alpha, `P -> reverse(alpha(P))` induces an involution on
these orbits: alpha normalizes the automorphism group and its square is an
automorphism. A fixed orbit supplies an automorphism g such that
`reverse(alpha(P))=g(P)`. Then `g^{-1} alpha` is exactly the positional
anti-reversal of P. Thus every self-converse class has a fixed-blue
presentation. This is an existence argument, not a claim that all paths
are blue.

The independent checker [moves_residual_audit.py](moves_residual_audit.py)
imports no producer modules. It reconstructs the free reflection-pair
coordinates from scratch. At order `2m`, there are `m(m-1)` free bits
after forcing the standard path forward. It obtains respectively 4, 64,
and 4096 blue matrices at orders 4, 6, and 8. Independent forward/backward
reachability and the two-score filter leave exactly 3, 21, and 537
matrices. Their sets agree exactly with the producer receipt, with no
duplicates or missing retained matrices.

All **561** certificates replay: each cycle visits every vertex exactly
once, every consecutive and closing arc exists, the designated edge
actually belongs to the cycle, and its independently calculated degree
imbalance is nonzero and equals the recorded value. The producer receipt
SHA256 is `ad8e2b15671eddec7b6a5c6d0c41236af96ee5f7ee22dfd6555d499c2720d108`.
The [independent receipt](moves_residual_audit.json) records its own source
digest, exact universes, and imbalance histograms.

```text
python 05-knowledge/results/scc_basins_20261020d/moves_residual_audit.py
python -O 05-knowledge/results/scc_basins_20261020d/moves_residual_audit.py
```

Both modes pass with identical JSON. The prior regular-center finite
checks are inherited inputs to the final singleton-padded corollary;
this new independent replay covers exactly the 561 nonregular residual
certificates above.
