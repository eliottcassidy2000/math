# Independent attachment-lane audit of the boundary congruence and strong sharpness

**AUDIT PASS**, 2026-11-02. The author of `attachment.md` did not author
`BOUNDARY_BUDGET.md` or its computation. No mathematical objection found
to B1, B4–B7, including their strong-target and fixed-old-source quantifiers.

For B1 the essential invariant is the exact set of old blocks touched by
new-old **path adjacencies**, not the set of blocks adjacent to a new vertex
in the ambient graph. On that path stratum every acted-on old block is
disjoint from the endpoints of every new-old step. Thus independently
supported old automorphisms preserve all steps of the word even when they
are not automorphisms of the enlarged tournament. Faithfulness and a full
vertex order make this action free. The gcd and prime-wise valuation
statements follow with the displayed quantifiers.

For B4, a path affected by a changed new-new edge uses two incidences of
the new-vertex path degrees on that edge. At most `2r-2` incidences remain
available for old-new adjacencies; touched old blocks can only be fewer.
The same stratum action fixes the prescribed new-new edge and acts freely.
Apply the divisor separately to paths lost and paths gained. This proves
the signed difference congruence; it does not assert that the two affected
sets have a common bijection. Repeating a sequence of new-edge reversals
preserves the congruence because every intermediate extension has the same
named old source and old-new ports.

The strong-sharpness construction has exactly r new vertices and exactly
k old directed triangles. With `t=k-2r`, its count before rewiring is
`41^(r-1) H(W_t)`, with exact valuation t because
`45*4^t-1=-1 mod3` for t>=1 and `H(W_0)=41`. The exceptional t=0 formula
is necessary. For r>=2 reversing the first-to-last apex edge changes the
count by a multiple of `3^(t+2)`, hence cannot change its valuation t.

Every block remains internally strong. All forward connections between
blocks still exist through old vertices; the one new backward edge links
the final apex to the first. These paths make the whole tournament strong.
No old-old or old-new edge is changed by the rewiring. The r=1 case is
already strong, and no r=0 assertion is silently added. The construction
proves sharpness of the necessary divisibility budget inside the strong
locus, not production of a prime value.

An independent implementation placed **all apices after all old vertices**,
built the graph from explicit arc sets, and counted Hamiltonian paths by
a dictionary subset DP. It imported no producer code or result table.
Separate reachability searches checked strongness from every vertex, and
every old edge was compared before and after the new-edge reversal.

| r | t | Order | Before | After |
|---:|---:|---:|---:|---:|
|1|0|7|41|41|
|1|1|10|537|537|
|1|2|13|6471|6471|
|2|0|14|1681|7009|
|2|1|17|22017|107985|

Each before value equals the claimed formula, each after value has
valuation t, and each difference is divisible by `3^(t+2)`. These checks
support the written all-order argument; the table is not its substitute.

The relation to this lane's triangle chirality detector is consistent.
B1 budgets how many independently rotating blocks a path can touch, while
Theorem B in `attachment.md` determines how two actual new-vertex port
signatures detect one old triangle's reversed orientation. Neither may
replace named attachment signatures with their counts or discard the
actual source weights.
