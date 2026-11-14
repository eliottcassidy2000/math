# A target-first selector from a single broken cyclic seam

Take the cyclic tournament with three nonempty blocks
`p → {b₀,b₁} → {x₁,…,x_n} → p`.
The two nonsingleton blocks are transitive in the listed order.
For a selected (j), reverse exactly the arc (b_i→x_j) with
(i=0) or (1). Set (a=j-1) and (b=n-j). The exact counts are

```
unchanged: N = 3·2^n−1
flip at b₀: N += 2^(b+1)(3^a−1)+2^a+[a=0]−2[b=0]
flip at b₁: N += 2^b(2·3^a−3)−[b=0].
```

This construction varies the actual seam of the tournament, not
merely its residue. In particular the flip at (b₀) with (a=0)
and (b≥1) gives a positive exact change of two at every width.

## Proof by counting gained and lost paths

The unchanged count is the singleton-cycle formula with profiles
((1,4)) and ((1,2^n)). A path affected by the edge reversal uses
one of the two possible directed orientations consecutively.
Contract that ordered pair to one marked position (w). Arcs into
(w) follow the old first endpoint, and arcs out of (w) follow
the old second endpoint. Double or missing arcs are allowed in this
auxiliary directed graph. Each permutation counted lifts uniquely.

Let (z) be the other vertex of the two-vertex block. The ordinary
vertices split into left and right groups of sizes (a,b) in their
transitive order. Every ordinary vertex points to (p) and is
beaten by (z), and (p→z). For lost paths (w) points to the
right group only, has no incoming arcs from ordinary vertices,
and (w↔p). For gained paths (w) points to all ordinary vertices,
receives from the left group only, and has no arc with (p).
With (i=0), (w,z) have no connecting arc in the lost graph and
are bidirected in the gained graph. With (i=1), only (z→w)
holds in both.

For each order of the three special vertices, distribute ordinary
vertices into allowed gaps. Within any gap their path order is
fixed because they are transitive. The following table lists all
four counts; (A=[a=0]), (B=[b=0]):

| Special order | lost, (i=0) | gained, (i=0) | lost, (i=1) | gained, (i=1) |
|---|---:|---:|---:|---:|
| p,z,w | 0 | (3^a2^b) | (2^b) | (3^a2^b) |
| p,w,z | 0 | 0 | 0 | 0 |
| z,p,w | (2^b-A) | 0 | (2^b-A) | 0 |
| z,w,p | 0 | (2^a-B) | (A) | (2^a-B) |
| w,p,z | (2^b) | (3^a2^b-2^a) | (2^b) | (3^a2^b-2^a) |
| w,z,p | 0 | (2^a-B) | 0 | 0 |

For example, for gained order `w,p,z`, the left vertices may
occur in the starting gap (they point to (w)), in the gap
between (w,p), or at the end. Right vertices may occur only
in the middle gap or at the end. Since (w→p) is absent, the
middle gap is mandatory, excluding exactly (2^a) assignments.
For lost order `z,p,w`, the middle gap between (z,p) is
mandatory, with every left vertex there. Each right vertex
chooses this gap or the end, excluding the empty middle gap
only when (a=0). Other rows follow the same gap conditions.

## A deterministic absolute inverse

For the (b₁) flip with (b≥1),

```
N+1 = (3·2^(a+1)+2·3^a−3)·2^b.
```

The parenthesized mantissa is positive, odd, and strictly
increasing in (a≥0). At a prescribed odd target, compute
(b=v_2(N+1)) and binary search (a) to match the entire odd
mantissa. On acceptance, construct the actual tournament at
that position. This covers prime 59, using (a=1,b=2), and
119=7·17, using (a=1,b=3). The (b₀) flip at (a=2,b=1)
gives the prime 83. These examples are not congruence hits.
Neither this family nor the whole transitive-spine one-apex
family suffices to select every prescribed odd count.

The replay tests both exact formulas for every (n≤11) and
each reversal position, by literal subset-path dynamic
programming on the flipped tournament.
