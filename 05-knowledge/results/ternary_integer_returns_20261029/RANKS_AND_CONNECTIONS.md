# What 21 measures, and two exact arithmetic bridges

Status labels: **CITED** for current rank records and their family provenance;
**PROVED / FINITE-EXACT** for the local-certificate optimization and finite-field
bridge; **HEURISTIC** for rank-distribution models. No upper bound on rational
elliptic ranks is proved here.

## 1. Rank 21 is not an established ceiling

The [PPVW model](https://math.mit.edu/~poonen/papers/bounded-ranks.pdf),
Theorem 1.1.1 and Remarks 1.1.2–1.1.4, predicts finitely many rational
isomorphism classes of rank above 21. Its theorem concerns the model, not
actual elliptic curves. The separate minimalist prediction puts density one
in ranks zero and one. Finite high-rank examples contradict neither prediction.
The threshold balances `H^(5/6)` curves against a rank-tail scale
`H^(-(r-1)/24)`: `21=1+24*(5/6)`. The authors explicitly flag high-generic-rank
families as a limitation requiring a refined model.

The linked [curve 273](https://elliptic-rank.icarm.cloud/curve/273) has certified
rank **at least 30**. [Curve 302](https://elliptic-rank.icarm.cloud/curve/302)
has rank at least 31, and [curve 4239](https://elliptic-rank.icarm.cloud/curve/4239)
at least 32. Their reported exact ranks use GRH and BSD. The last is a
specialization of a rank-17 elliptic fibration on a K3 surface of discriminant
1152. A specialization jump is a different statistic from generic rank.

The 2026-10-29 public [database](https://elliptic-rank.icarm.cloud/database.json)
snapshot has 6127 entries, 563 with stated lower bound at least 22, and 66
with lower bound at least 30. The factual [snapshot summary](rank_population_snapshot.json)
pins the raw download hash, attribution, and high-rank IDs. The
[entry rules](https://elliptic-rank.icarm.cloud/api) select rank/height records;
these counts are not frequencies in a random or height-complete population.

The incoming [rank-family audit](../session97_overnight_arithmetic_carriers_20261029/ranks/README.md)
already does the stronger work: exact binary-character independence
certificates for all 93 displayed points, rational isomorphisms from the three
parameterized models, and transport of each observed local root. It also
derives the conditional family threshold `r0+1+24 alpha`, where `alpha`
counts distinct specializations by height and a shifted rank-tail law is
assumed. We reuse those results rather than replacing them by point incidence
or an uncalibrated regression against leaderboard rows.

## 2. Extend the forgotten August optimal-observer result

The [August rank-30 certificate](../elliptic_rank30_optimal_local_observer_certificate_codex_20260821.out)
already achieved the information-theoretic minimum of fifteen good-prime
observation sites. The new calculation extends that optimization:

| Curve | Independent witness points | Optimal binary-character prime sites | Largest prime in this bank | Extra no-2-torsion prime |
|---|---:|---:|---:|---:|
| 273 | 30 | 15 | 1291 | 23 |
| 302 | 31 | 16 | 751 | 31 |
| 4239 | 32 | 16 | 593 | 7 |

For curve 273, the older fifteen-prime bank has a smaller largest prime;
this table is a new certificate, not an improvement to that secondary cost.

At a good odd prime, all quadratic-character observations factor through
`E(F_p)/2E(F_p)`, whose dimension is at most two. Thus proving `r` independent
binary columns needs at least `ceil(r/2)` such prime sites. We achieve that
bound in each case by retaining full local root data. At a split prime the
three root characters sum to zero: their product is the square obtained from
the completed-square cubic. There are two available bits, not three.

The checker recomputes every root, observation bit, discriminant, and point
incidence. A second elimination algorithm verifies the full ranks. Removing
any selected prime lowers the rank; duplicating a point-column is rejected.
A separate good prime with no cubic root excludes rational 2-torsion, so a
full binary matrix proves independence modulo torsion by repeatedly halving
any hypothetical integer relation. That additional torsion test is explicitly
not included in the optimal **witness-bank** count. The optimum is within this
local-character model, not among all conceivable arithmetic proof systems.

See [check_rank_banks.py](check_rank_banks.py) and
[optimal_rank_banks.json](optimal_rank_banks.json). Coordinate signs remain
invisible to these mod-two probes; the exact point labels and transported
roots remain necessary for subsequent operations.

## 3. A genuine 7-to-21 map for the cheaper controller

There is an exact projective interpretation of the controller's probabilities
`p=1/7` and slope-tilted `q=1/21`. It is a map of **distribution models**, not
an elliptic-rank theorem or a conjugacy of integer Collatz orbits.

Take independent uniform digits in a finite field `F_b`, with `b=2` or `4`.
Let `A` be the index of the first nonzero digit; then
`P(A=a)=(b-1)b^(-a)`. Group the digits into triples and keep the first nonzero
triple `V`. It is uniform on `F_b^3\{0}`. Its projective class is therefore
uniform on `PG(2,b)`, with `b^2+b+1` points. The event `3|A` is exactly the
single projective point `[0:0:1]`. Consequently

\[
 P(3\mid A)=\frac{b-1}{b^3-1}
 =\begin{cases}1/7,&b=2,\\1/21,&b=4.\end{cases}          \tag{17}
\]

These are precisely the odd-Haar letter law `2^(-a)` and normalized
slope-tilted law `3*4^(-a)` from the earlier moment packet.

There is also an exact three-way refinement of the seven-point geometry.
Identify the three-dimensional `F_4` vector space with `F_64`. Then

\[
 PG(2,4)=F_{64}^{\times}/F_4^{\times},\qquad
 [x]\longmapsto x^{21}\in F_4^{\times}                 \tag{18}
\]

is well-defined. Its three fibres each contain seven projective points.
The norm-one fibre is the embedded `PG(2,2)` arising from `F_8`; the other
two are multiplicative translates of it. Indeed the norm-one subgroup has
order 21 and equals `F_4^* F_8^*`, with trivial intersection of those two
groups. Choose an `F_2` basis of `F_8` as the `F_4` basis of `F_64`; the
degrees three and two are coprime. Linear independence of vectors in `F_8`
over `F_2` is unchanged after this base extension, so each seven-point fibre
is a Fano subplane, not merely a seven-element set. A line through two of
its points meets that subplane in three points, as the exact checker confirms.

Thus conditioning the uniform 21-point law on norm phase one recovers the
uniform seven-point law. The factor `1/3` in `1/21=(1/3)(1/7)` is an actual
three-valued norm sidecar. The quotient forgets the number of initial zero
triples and the nonzero scalar; to recover `A` and the digit history one must
retain those clocks and digits. This is why the projective picture does not
grant a free resampling or a source-independent Collatz descent step.

The [checker](check_projective_bridge.py) constructs `F_64` explicitly,
verifies its `F_4,F_8` subfields, all three seven-point fibres, and their line
incidence. This recovers a useful finite-field connection to the older F64,
Fano and Hamming work while declaring exactly what is retained and lost.

## 4. Where the tournament comparison stops

The repository's THM-4812
[all-order presentation-fibre theorem](../../../01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md)
excludes fibre cardinalities 7 and 21 through orbit-product bounds and a closed
finite compatibility window. That is an actual impossibility theorem, whereas
rank 21 is a heuristic population threshold. Neither supplies a map identifying
Hamiltonian presentations with Mordell–Weil generators. The projective map
(18) relates two geometric carriers, not those forbidden tournament fibres.

What transfers is the proof architecture: retain the source and its allowed
operations, derive an unbounded reduction, then close the exact finite window.
The cheaper controller and the fruit-circle cover now give two more rigorous
examples. For ranks, the missing ingredients are family counting and rank-jump
laws; for Collatz renewal, they are pointwise reset payment. Numerical labels
do not fill either gap.
