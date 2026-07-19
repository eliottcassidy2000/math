---
id: THM-1221
title: THE SEVEN-WALL STRICT-SPECTRUM HUNTER FLOOR — every seven-speed packet has a maximum spanning tree and global safe set of mass at least 15/154
status: PROVED (THM-1166 pair formula and analytic tails; globally complete finite ratio-channel/cut referee; Lean arithmetic and provider-boundary consumer; no bounded speed-box assumption)
source: codex-2026-07-19-S82
depends_on: [LEM-042, THM-965, THM-1166]
related: [THM-1203, THM-1218, THM-1219, HYP-7678]
script: 04-computation/lrc14_seven_wall_strict_spectrum_hunter_floor_referee_codex_S82.py
output: 05-knowledge/results/lrc14_seven_wall_strict_spectrum_hunter_floor_referee_codex_S82.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCSevenWallStrictSpectrum.lean
script_sha256: 71b3389f46fc0844c0048f116428469f269b027d2dafaece3b3e2266b669cfa1
output_sha256: 640e191746bce43ec1d451568673480eadfc1afac71ed487895f36dbcd5ee506
formalization_sha256: 3f00a8585e66e1619bc2540a80e2c484c98c6c51c96699ad38d34123be7528a7
---

# THM-1221 — the seven-wall strict-spectrum Hunter floor

## Statement

At radius `1/14`, put

```text
D_s={t in R/Z: ||st||<1/14},
rho(a,b)=mu(D_a intersect D_b).                             (1)
```

For every seven distinct positive integer speeds `s1,...,s7`, the complete
weighted graph with edge weights `rho(si,sj)` has a spanning tree `T` with

```text
sum_(ij in T) rho(si,sj) >= 15/154.                        (2)
```

Consequently the common safe set

```text
S=(R/Z) \ union_i D_si
```

satisfies

```text
mu(S)>=15/154.                                             (3)
```

Both assertions are global and uniform: no upper bound on the speeds, their
ratios, or their gcds occurs in the statement.  The constant is a proved
floor, not asserted to be sharp.

There is also a deliberately simpler subcertificate

```text
sum_(ij in T) rho(si,sj) >=265/2772
                           =2/21+1/2772,                   (4)
```

which uses only the first strict spectral value plus the elementary
singleton-cut ledger.  The stronger value in (2) comes from resolving every
remaining threshold-component branch on the complete projective ratio bank.

If a common integer `g` divides all seven speeds, every connected interval
covered by their danger combs has length `L` satisfying

```text
gL<=1-15/154=139/154.                                     (5)
```

When the protected-needle input gives `L>=1/(7m)`, with `m` the retained-core
maximum, this improves the THM-1166 common-dilate consequence to

```text
g/m<=139/22.                                              (6)
```

## 1. The exact spectrum around `1/63`

For a reduced coprime pair `a<b`, THM-1166 gives

```text
rho(a,b)
 =1/49+[F((a+b) mod 14)-F((b-a) mod 14)]/(196ab),
F(r)=r(14-r),                                             (7)

rho(a,b)>=1/49-1/(4ab).                                   (8)
```

The right side of (8) is strictly greater than `1/63` for `ab>=56`.
Therefore the following finite banks are complete, rather than samples.
The strict-low bank is

```text
(1,10):1/70, (1,11):1/77, (1,12):1/84, (1,13):1/91,
(2,11):1/77, (3,10):1/70, (3,11):1/77.                   (9)
```

Adjoining the equality channels gives the complete `rho<=1/63` bank:

```text
(1,9), (1,27), (2,9), (4,9), (5,9)      at 1/63,
together with the seven rows in (9).                      (10)
```

The exact product table for the twelve positive ratios in (10) has only two
triangle factorizations,

```text
(9/4)*12=27,                 12*(9/4)=27.                 (11)
```

Thus the graph `rho<=1/63` has clique number three, not two.  There are two
three-cliques up to scaling and permutation,

```text
(1,12,27),                    (4,9,108).                  (11a)
```

They are reciprocal-dual, not scale/permutation-equivalent.  Exact
closed-neighborhood intersection shows that neither has a fourth common
low-or-equal neighbor.  In particular the graph has no `K4`.

The same tail argument is again exact at the next spectral step:

```text
1/49-1/(4*60)>5/308.                                     (12)
```

Scanning the complete coprime core `ab<=59` shows that the first value
strictly above `1/63` is

```text
c=5/308,                                                  (13)
```

uniquely at reduced ratio `4:11`.

For reference, the THM-1166 sharp three-speed sum is replayed independently:

```text
rho(x,y)+rho(x,z)+rho(y,z)>=51/1183.                      (14)
```

Its finite reduction has `124` oriented ratios of reduced product at most
`41`, `184` at product at most `57`, and exactly `22,692` distinct incident
ratio pairs.  The unique equality shape is `(1,13,169)`.  This replay is used
in the transparent proof of (4); the stronger proof of (2) also recomputes
every surviving edge directly from (7).

## 2. The simpler `265/2772` certificate

Let

```text
G_ge = graph of edges rho>=1/63.                          (15)
```

The subscript is mnemonic for the high wall; the inequality in (15) is weak.
If `G_ge` is connected, the no-`K4` result implies that among seven vertices
some edge has weight strictly above `1/63`.  By (13) it has weight at least
`5/308`.  Extend this edge to a spanning tree of `G_ge`; its other five edges
have weight at least `1/63`.  Hence

```text
5/308+5/63=265/2772.                                     (16)
```

Suppose instead that `G_ge` is disconnected.  Every edge between different
components is strict-low.  The strict-low graph from (9) is triangle-free, so
there are exactly two high components.  The following exact common-neighbor
calculation forces their sizes to be `1+6`.

Orient the seven strict ratios in (9) both ways, producing a set `R_<` of
fourteen ratios.  Normalize two possible centers to `1,y`.  Every common
strict-low neighbor lies in

```text
R_< intersect y R_<.                                     (17)
```

As `y` runs over the complete quotient bank `R_</R_<`, the nonempty
intersection-size histogram is

```text
size 1: 14 rows,       size 2: 72 rows,       size 4: 6 rows. (18)
```

The six four-neighbor center ratios are

```text
3/110, 1/3, 10/11, 11/10, 3, 110/3.                     (19)
```

For every row in (19), the four displayed neighbors have exactly the original
two common centers.  Thus a `2+5` cut is impossible by the maximum four in
(18), while a `3+4` cut would require a forbidden third center for one of the
six four-neighbor rows.  Only `1+6` remains.

Normalize the singleton speed to one.  Its fourteen possible crossings have
the exact weight capacities

```text
1/70: four,       1/77: six,       1/84: two,       1/91: two. (20)
```

The other six vertices form a high clique, since two strict-low neighbors of
the singleton cannot themselves form a strict-low edge.  Let `m` be the
number of `1/70` crossings.  For `m<=3`, every other crossing has mass at most
`1/77`, and (14) gives the mutual-edge credit

```text
u=51/1183-2/77=223/13013>1/63.                           (21)
```

Connect the `6-m` non-`1/70` endpoints with `5-m` edges of weight at least
`u`, attach the `m` remaining endpoints with high edges, and use the largest
singleton crossing.  The four exact floors are

```text
m=0: 1284/13013,
m=1: 115601/1171170,
m=2: 16303/167310,
m=3: 37547/390390.                                      (22)
```

When `m=0`, a crossing of weight at least `1/77` exists because only four
oriented strict channels lie below `1/77`.  When `m=4`, the four ratios are
forced to be

```text
1/10, 3/10, 10/3, 10.                                   (23)
```

Every mutual edge among these four endpoints has mass at least `32/1575`.
Using one such edge, four further high edges, and a `1/70` crossing gives

```text
1/70+32/1575+4/63=103/1050.                              (24)
```

Every value in (22)--(24) is strictly above (16).  This proves the simpler
global certificate (4), including the formerly worst `m=4` branch whose
relaxed triple ledger alone gave only `111161/1171170`.

The same data give an analytic proof that this entire disconnected branch
also clears the stronger target (2).  Put

```text
s=32/1575.
```

Every edge between two `1/70` endpoints has weight at least `s`, every edge
between two other endpoints has weight at least `u`, and every endpoint pair
has weight at least `h=1/63`.  For `m>=1`, build a tree on the `m` special
endpoints using `(m-1)` `s`-edges, a tree on the other `6-m` endpoints using
`(5-m)` `u`-edges, one `h`-edge joining the two trees, and one `1/70` crossing
to the singleton.  Together with the `m=0` row, this gives

```text
L_0 = 1/77+5u                         =1284/13013,
L_m = 1/70+(m-1)s+(5-m)u+1/63,
L_1 =115601/1171170,
L_2 =28411/278850,
L_3 =615257/5855850,
L_4 =633883/5855850.                                  (24a)
```

All five values exceed `15/154`; the uniform minimum is the `m=0` row, with
exact margin

```text
1284/13013-15/154=3/2366.                               (24b)
```

Thus the later `3,003`-packet replay is an independent diagnostic, not a
logical dependency of the strong theorem.

## 3. The strict-component sharpening to `15/154`

Return to the connected `G_ge` branch and now define

```text
G_gt = graph of edges rho>1/63.                           (25)
```

This strict/weak distinction is the key sharpening.  Between two components
of `G_gt`, every edge has `rho<=1/63`.  Choosing one representative from each
component makes a clique in the closed-low graph (10), so (11) implies that
`G_gt` has at most three components.

Three components cannot contain seven vertices.  Their representatives must
be a scaled copy of one of the reciprocal-dual triangles in (11a).  For
`(1,12,27)`, the exact common closed-neighbor sets of the opposite
representative pairs are

```text
N(12) intersect N(27) = {1,324},
N(1)  intersect N(27) = {9/4,12},
N(1)  intersect N(12) = {4/9,27}.                        (26)
```

For `(4,9,108)`, they are

```text
N(9) intersect N(108) = {4,243},
N(4) intersect N(108) = {9,48},
N(4) intersect N(9)   = {1/3,108}.                       (26a)
```

An additional vertex in one component must be closed-low to the two other
representatives, so (26)--(26a) permit at most two vertices in each component.
The total is at most six in either shape, a contradiction.

If `G_gt` has one component, it has a six-edge spanning tree.  Every edge of
that tree is strictly above `1/63`, hence at least `5/308` by (13).  Therefore

```text
6*(5/308)=15/154.                                        (27)
```

It remains only to exclude a smaller tree in the two-component branch.  All
cross-component ratios lie in the twenty-four oriented closed-low ratios
`Q` from (10), and at least one crossing has equality `rho=1/63` because
`G_ge` is connected.

For a `1+6` cut, the six opposite vertices lie in `Q`.  Exact evaluation of
the complete `C(24,2)` internal pair bank gives

```text
min{rho(x,y):x,y in Q, rho(x,y)>1/63}=1/56.              (27a)
```

Their `G_gt` component has a five-edge tree, so adjoining an equality crossing
gives

```text
5/56+1/63=53/504>15/154.                                 (27b)
```

For a `2+5` cut, normalize the two centers to `1,r`.  Five common closed-low
neighbors force, by the complete quotient bank `Q/Q`,

```text
r in {1/3,1/2,2,3}.                                     (27c)
```

There are six common neighbors in each of these four rows, and the center
edge has mass at least `1/21`.  The five-vertex `G_gt` component has a
four-edge tree, each edge at least `5/308`, so

```text
1/21+4*(5/308)+1/63=89/693>15/154.                       (27d)
```

For a `3+4` cut, normalize one center to one.  The other two centers lie in
`Q/Q`.  Exhausting the connected normalized center triples gives the exact
common-neighbor histogram

```text
size 1: 5,532 rows,       size 2: 243 rows,       size 3: 12 rows. (27e)
```

Thus no connected three-center row has the four common neighbors required by
a `3+4` cut.  This finishes the analytic two-component proof.

The referee additionally performs the larger, globally complete projective
cut census below as an independent replay:

```text
cut      complete normalized survivors     minimum exact MST
1+6                 131,593                    2228/18711
2+5                      12                   15473/103950
3+4                       0                       --       . (28)
```

The filters in (28) are exactly: each side is connected in `G_gt`, every
crossing belongs to the closed-low bank, and at least one crossing has equality
`rho=1/63` so that `G_ge` is connected.  For `1+6`, this is the complete
`C(24,6)` ratio bank.  For `2+5`, the second center lies in the finite quotient
bank `R_<=/R_<=`, and the other five vertices lie in the exact intersection of
its two closed neighborhoods.  For `3+4`, choose the four neighbors of the
normalized first center and intersect their four inverse neighborhoods; no
connected row survives.  Thus (28) is not a bounded speed census: every
positive-integer packet in this branch scales to one of its rows.  The census
is diagnostic because (27a)--(27e) already prove the branch.

Both nonempty minima in (28) are greater than `15/154`.  Finally, in the
disconnected `G_ge` branch of Section 2, normalization gives exactly
`C(14,6)=3,003` singleton packets.  Direct evaluation of their maximum
spanning trees has minimum

```text
18497/156156>15/154,                                     (29)
```

at normalized ratio packet

```text
(1; 1/13,1/12,2/11,11/2,12,13).                         (30)
```

This `3,003`-row value independently checks the analytic floors (24a)--(24b);
it is not needed to prove them.

Equations (24a)--(24b) and (27)--(27e), together with the impossibility of
three strict-high components, prove (2).  Equations (28)--(30) are independent
full-census confirmations.

## 4. Hunter, period, and incidence consumers

For a spanning tree `T` and any nonempty active vertex set `A`, the forest
induced by `A` has at most `|A|-1` edges.  Pointwise,

```text
1_(union_i D_si)
 <=sum_i 1_(D_si)-sum_(ij in T)1_(D_si intersect D_sj). (31)
```

The seven one-comb masses sum to one.  Integrating (31) and using (2) proves
(3).  If all speeds are multiples of `g`, the safe set is `1/g`-periodic and
has mass at least `15/(154g)` in each fundamental period.  A connected covered
interval cannot cross that safe mass, proving (5); composing with the
protected-needle lower bound proves (6).

There is a useful exact interface with THM-1203 and THM-1218.  Let `B` be any
four-speed BAD event controlled by one of THM-1218's deletion triangles.
THM-1203 gives `mu(B)<=2/21`, so for the seven-speed safe set `S`,

```text
mu(S \ B) >=15/154-2/21=1/462.                           (32)
```

If the quartet is not the arithmetic-progression completion singled out by
THM-1218, its contrapositive gives `mu(B)<=60/637`, and hence

```text
mu(S \ B) >=15/154-60/637=45/14014.                      (33)
```

These are finite-transfer/erosion budgets.  More explicitly, if a later
localization or discretization step loses at most `epsilon` of the global safe
set before excluding `B`, positivity survives whenever

```text
epsilon<1/462                         in the arbitrary branch,
epsilon<45/14014                      in the non-AP branch. (34)
```

The weaker aligned floor (4) gives the historical margins `1/2772` and
`355/252252`; (32)--(34) are the stronger budgets actually available now.
The AP alternative is precisely where THM-1218's equality of four direct
beat masks remains useful.  What is still missing is a theorem that transports
the global safe mass, with erosion below (34), into the particular protected
slow gap or beat block of a hypothetical LRC(14) cover.  The numerical crossing
of `60/637` does not by itself identify a heavy quartet: the left side of (2)
is a sum of six pair overlaps, while THM-1218 thresholds one four-speed BAD
event.

## 5. Carrier, tournament, and assumption audit

The pairwise observable is `rho(si,sj)`.  Switching at `1/63` gives two nested
binary relations, `G_ge` and `G_gt`.  A tournament gauge can orient high edges
up the speed order and reverse strict-low edges; the referee reports its score
histogram, directed triangles, SCCs, edge flips, a Redei Hamiltonian path, and
the Hamiltonian-path count for the extremal disconnected diagnostic packet.

That tournament preserves the threshold cut and global Hunter credit but
destroys the equality/strict distinction unless the edge color is retained.
It also destroys interval phase, tooth address, gcd positioning error, and the
four-deletion incidence defining a THM-1218 BAD event.  Candidate vertex sets
challenged here include runners, reduced ratio channels, threshold components,
cut obligations, BAD quartets, and safe-set atoms.  The faithful carrier for
this theorem is the edge-weighted complete graph together with the two nested
threshold graphs.  For composition with THM-1218 it must be enlarged by the
four-deletion circuit hypergraph; a bare runner tournament is insufficient.

## 6. Provider boundaries and verification

The analytic providers are (7)--(8), already proved in THM-1166, and the
standard Hunter forest inequality (31).  The product bounds `ab<=55` and
`ab<=59` are derived from (8) before enumeration.  All later finite work is on
the exact scale-invariant ratio banks forced by those analytic tails.  No
finite height, diameter, or speed cutoff is assumed.

The dependency-free referee uses `fractions.Fraction` in every proof-facing
calculation.  It independently matches the folded and tent formulas on `7,140`
pairs, replays all `22,692` sharp triple candidates, proves the strict and
closed channel banks, common-neighbor tables, exact singleton multiplicities,
all `3,003` disconnected-high packets, and every complete strict-high
two-component cut in (28).  Normal and optimized Python runs are byte-identical.

`LRCSevenWallStrictSpectrum.lean` kernel-checks the exact rational constants,
the exhaustive surviving-branch consumer, Hunter and common-period arithmetic,
and both erosion budgets.  Its provider structure makes the boundary honest:
the identification of arbitrary speed packets with one of the exact ratio-cut
branches remains supplied by the analytic formula plus the finite referee; it
is not smuggled into Lean as an unproved equality of Haar measures.

Frozen hashes are

```text
source         71b3389f46fc0844c0048f116428469f269b027d2dafaece3b3e2266b669cfa1
output         640e191746bce43ec1d451568673480eadfc1afac71ed487895f36dbcd5ee506
formalization  3f00a8585e66e1619bc2540a80e2c484c98c6c51c96699ad38d34123be7528a7
```
