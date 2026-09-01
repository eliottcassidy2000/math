---
id: THM-4326
title: "Rank-two wall-graph completion of the fixed-pool arbitrary two-outsider chart"
status: >
  PROVED RELATIVE TO THM-4231/4150 + FINITE-EXACT + CLEAN-ROOM
  INDEPENDENTLY AUDITED. For the displayed thirty-label pool, every labelled
  nine-body has Haar safe mass at least 4/63 after adjoining any two distinct
  positive outsiders. The new rank-at-most-two wall graph closes all 181,194
  pairs in THM-4231's exact finite remainder: 181,087 by one weighted-degree
  bound and 107 by exact optimization of all C(30,9) bodies. The exact global
  normalized minimum is for the retained rank-two lower bound, not the full
  safe mass. THM-4150 supplies the corresponding doubled-body/two-odd-tail
  thirteen-speed families. No arbitrary-row entry or LRC(14) follows.
source: root + literal_exit_theory + endpoint588_scout + endpoint586_independent / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
related:
  - THM-4286-signature-response-nonfactorization-and-two-deck-surgeries
  - THM-4287-repaired-carrier-endpoint-637-descent
  - THM-4324-lrc14-endpoint-586-direct-literal-closure
artifact_root: 05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326
artifact_manifest: 05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/SHA256SUMS
artifact_manifest_sha256: b5d04c30a1169e904f9a66b4470cfd39dd13efd3a9ab9ba43eb482926dbae9d8
verifier_output_sha256: 6edb482db01c95b20d496bb041fbe27f37c0080e99be583c677ad4a788f959e4
independent_artifact_root: 05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/independent_cleanroom
independent_artifact_manifest_sha256: 1fdd6025bae13f9bca08eb4fb97a82f149d845522e6dd6d9ac7d70629af79302
primary_scripts:
  - 04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/export_thm4231_remainder.py
  - 04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/rank2_degree_screen.cpp
  - 04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/rank2_allbody_audit.cpp
  - 04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/rank2_branch_bound_independent.cpp
  - 04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/verify_rank2_arbitrary_pair_packet.py
independent_scripts:
  - 04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/independent_cleanroom/rank2_event_sweep_cleanroom.cpp
  - 04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/independent_cleanroom/rank2_event_sweep_wide_cleanroom.cpp
  - 04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/independent_cleanroom/rank2_selected_ratio_exact.cpp
  - 04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/independent_cleanroom/rawcell_winner_replay.py
  - 04-computation/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/independent_cleanroom/verify_packet.py
audit: >
  PASS / ACCEPT for the closed outer and nested manifests. Wide O2/O3 primary
  screens agree byte-for-byte on all 181,194 rows. Flat all-body enumeration
  and a separate branch-and-bound optimizer agree on all 107 exact minima.
  A clean-room event-state sweep matches every primary graph field on all
  181,194 rows and every minimizing body on the 107 exceptions; an
  unaggregated raw-cell replay matches 110 selected winners. Both closed
  verifiers pass normally and under optimized Python invocation.
---

# THM-4326 -- rank-two wall-graph completion of the fixed-pool pair chart

**PROVED RELATIVE TO THM-4231/4150 + FINITE-EXACT + CLEAN-ROOM
INDEPENDENTLY AUDITED. LRC(14) REMAINS OPEN.**

## 1. Statement

Let

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},                                       (1)
```

and, for each finite positive label set `A`, put

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14}.               (2)
```

> **Fixed-pool arbitrary-pair theorem.** For every two distinct positive
> integers `q,r notin P` and every labelled `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {q,r}))>=4/63.                         (3)
> ```

The new finite computation proves a strict lower bound on THM-4231's exact
`181,194`-pair remainder. THM-4231 proves `(3)` for every pair outside that
remainder. Thus `(3)` is an arbitrary-pair theorem over the fixed pool, not
an endpoint-prefix or carrier-survival statement.

If `H=B union {q,r}`, `c` is any positive integer, and `a,b` are distinct
positive odd integers, then there is an `x in R/Z` with

```text
min_(v in 2cH union {a,b})||vx||>=1/14.                (4)
```

Indeed multiplication by `c` preserves Haar measure on the circle, so
`mu(G_(cH))=mu(G_H)`, and THM-4150 applies. Equation `(4)` is a structured
thirteen-speed family; it is not arbitrary-row entry.

## 2. Finite-null-wall graph lemma

The mechanism is abstract. Let `P={p_0,...,p_(n-1)}` be a finite labelled
pool and let `q,r` be distinct positive outsiders. Put

```text
D=lcm(14v:v in P union {q,r}).                         (5)
```

In the integer coordinate `t=Dx`, insert `0,D` and the walls

```text
(14k+1)D/(14v), (14k+13)D/(14v),      0<=k<v,         (6)
```

for every speed `v`. These are integers. On an open cell with consecutive
wall endpoints `u<w`, evaluate its midpoint `x=(u+w)/(2D)`. If
`s=v(u+w) mod 2D` is represented in `[0,2D)`, then

```text
||vx||>=1/14  iff  D<=7s<=13D.                         (7)
```

Every safe/unsafe indicator is therefore constant on the open cell. The wall
set is finite and Haar-null, so open-cell width summation computes the exact
measure even though `G_A` is closed and can contain wall points. The set
decompositions below are measure identities modulo this finite null set.

Retain cells on which `q,r` are safe. For such a cell `C`, let

```text
F(C)={i:||p_i x||<1/14}                               (8)
```

be its pool-failure set, and let `w_S` be the sum of integer widths of cells
with `F(C)=S`. For every body `B`,

```text
D mu(G_(B union {q,r}))=sum_(S intersect B=empty)w_S. (9)
```

Keep only failure ranks at most two:

```text
L2(B)=sum_(S intersect B=empty, |S|<=2)w_S.           (10)
```

All omitted widths are nonnegative, hence

```text
0<=L2(B)<=D mu(G_(B union {q,r})).                    (11)
```

Write

```text
w0=w_empty,  wi=w_{i},  wij=w_{i,j},
W=w0+sum_i wi+sum_(i<j)wij,
di=wi+sum_(j!=i)wij.                                  (12)
```

Exact weighted inclusion-exclusion gives

```text
L2(B)=w0+sum_(i notin B)wi+sum_(i<j, i,j notin B)wij
     =W-sum_(i in B)di+sum_(i<j, i,j in B)wij.        (13)
```

The last term restores the internal edges double-counted by the degree sum.
For `|B|=9`, nonnegativity of all `wij` yields

```text
L2(B)>=W-sum_(i in B)di
     >=W-(sum of the nine largest weighted degrees).  (14)
```

Thus strict positivity of

```text
63[W-(sum of the nine largest di)]-4D                 (15)
```

certifies every nine-body on the pair simultaneously. A nonpositive value
of `(15)` is an inconclusive relaxation, not an unsafe row.

For an inconclusive row, minimizing `L2` is weighted maximum nine-vertex
coverage. Equivalently maximize

```text
R(B)=sum_(i in B)di-sum_(i<j, i,j in B)wij.           (16)
```

At a partial body, the independent branch search sums the largest required
current marginal gains and omits only nonnegative pair penalties among
future choices. This is an upper bound for every completion of `(16)`, so
pruning strictly below the incumbent is sound. Equality is retained for the
least-mask tie break.

The observable `wij` is symmetric. The intrinsic object is an undirected
weighted graph with singleton loops, not a tournament.

## 3. Exact THM-4231 remainder

The maintained THM-4231 postprocessor exports

```text
E_rem count       181,194
ordered-pair FNV  3874fecac4ecbd8a
SHA-256           9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1
maximum endpoint  769.                                (17)
```

The verifier pins the complete ordered set, distinctness, outsider typing,
and equality with the degree-ledger pair set. The later THM-4287 universe is
also pinned as a strict `22,647`-row subset, but it is only a provenance
control: this theorem consumes all of `E_rem` directly and does not depend on
the later carrier descent.

The signed-128-bit degree screen gives

```text
strict degree-bound rows  181,087
coarse nonpositive rows       107
highest coarse exception      (50,554), unique at endpoint 554
degree-ledger FNV              f8da82c5a6d732ed
degree-ledger SHA-256          5b008404482b6e23006c0c1cb97407c22ddf34eb7cf4e09fd5f604206d8a0356.
                                                               (18)
```

Equation `(14)` proves all `binom(30,9)=14,307,150` bodies on each of the
first `181,087` rows. For the remaining `107`, flat enumeration checks

```text
107 binom(30,9)=1,530,865,050                         (19)
```

labelled body cases. Every minimum is strict and there are no equalities.
The separate branch-and-bound implementation agrees on the exact minimum
and least minimizing mask for all `107`; its complete search visits `55,104`
nodes and prunes `49,891`.

## 4. Normalized hostile and overflow controls

Raw ticks from different grids are not intrinsically ordered (MISTAKE-532).
Set

```text
ticks=63L2(B)-4D.                                     (20)
```

The scale-free certified surplus is `ticks/(63D)`. Among the `107` exact
rows, `(50,70)` has the smallest `ticks/D`. Exactly three coarse-positive
rows have a degree-bound ratio no larger than that candidate:

```text
(50,212), (50,274), (100,110).                        (21)
```

Fresh O2/O3 all-body optimization proves that all three exact ratios are
larger. Hence the unique minimizing pair for `L2/D` over the entire finite
remainder is `(50,70)`. Its least minimizing body under the frozen mask order
is

```text
(q,r)=(50,70),
B={80,85,88,95,143,145,168,193,240},
least mask=031c7400,
D=91,205,797,082,400,
L2(B)=5,794,739,949,188,
ticks=245,428,469,244.                                (22)
```

Consequently

```text
min_(E_rem,B) L2(B)/D
 =1448684987297/22801449270600
 =4/63+973922497/22801449270600
 >4/63.                                               (23)
```

Equation `(23)` is the exact minimum of the retained rank-two lower bound; no
uniqueness of the minimizing body is asserted.
By `(11)`, the actual Haar surplus is **at least** the displayed certified
surplus. No equality with the full safe mass and no full-mass minimizer are
claimed; discarded failure-rank-at-least-three cells may contribute more.

Exactly one residual row needs a wall grid larger than signed 64-bit:

```text
(713,719), D=9,351,275,651,380,222,560.               (24)
```

The proof sources use signed `i128` for grids, widths, masses, degrees,
scores, and ticks. The obsolete signed-64 scanner is excluded from the
packet. On the frozen domain `q,r<=769`, the elementary bound
`D<=D0(14q)(14r)<2.2*10^21` keeps every operation far below `2^127`.

The wide replay of the `22,647`-row control reproduces its numerical CSV and
uses the current two-limb FNV `611b1ba5c25594dd`; the stale pre-wide digest
`fe6111d297a72a3d` is explicitly rejected.

## 5. Independent audit and proof-graph consequence

The primary constructs all walls, sorts them, and classifies midpoints. The
clean-room implementation instead merges exact enter/leave events into a
state sweep. On all `181,194` pairs it matches every graph field; on all `107`
exceptions its independent optimizer matches every minimum and least body.
An unaggregated integer raw-cell evaluator also matches `110` selected
winners: the `107` exact exceptions plus the three normalized challengers.

The primary O2/O3 degree CSVs agree byte-for-byte. O2/O3 flat enumeration,
branch search, three-row normalized sidecar, and clean-room outputs agree at
their frozen semantic digests. Both manifests are closed and both verifiers
pass normally and under optimized invocation.

Sections 2--5 prove strict safety on every pair in `E_rem`; THM-4231 proves
`(3)` off `E_rem`. This proves the fixed-pool arbitrary-pair theorem and,
through THM-4150, `(4)`. **QED.**

## 6. Scope firewall

```text
source:       pair-safe wall cells with labelled pool-failure masks
target:       every nine-body over the fixed pool
map:          retain ranks 0,1,2 and aggregate exact widths
preserved:    retained Haar mass and every singleton/pair overlap
destroyed:    cyclic address, owners, cell order, and all rank>=3 masks
sidecar:      exact grid, ordered pair ledger, and normalized ticks/D
hostile:      (50,70), body 031c7400
decisive test: 63L2(B)-4D>0 on every exact minimizing body. (25)
```

This theorem does not say that higher-rank cells vanish, does not identify a
minimum full safe mass, and does not turn a carrier failure into physical
danger. It supplies neither a parity normalization nor a map from an
arbitrary normalized thirteen-speed row to the fixed `9+2` chart. It proves
no owner, arrival, or termination theorem. LRC(14) remains open.

## 7. Exact replay

From the repository root:

```text
python -B 05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/verify_rank2_arbitrary_pair_packet.py
python -B -O 05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/verify_rank2_arbitrary_pair_packet.py

python -B 05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/independent_cleanroom/verify_packet.py
python -B -O 05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/independent_cleanroom/verify_packet.py
```

The optimized independent invocation relaunches a nonoptimized verifier so
that Python cannot erase its assertion checks. Full C++ rebuild commands and
the optimization-invariant comparisons are frozen in `REPRODUCTION.md` and
the clean-room `reproduce.ps1`.
