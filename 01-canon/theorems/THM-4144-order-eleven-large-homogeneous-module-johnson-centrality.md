---
id: THM-4144
title: "Order-eleven large homogeneous-module Johnson centrality"
status: >
  FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every strong
  order-eleven tournament with a proper homogeneous module of size three
  through nine has only central rational and exact-coset Johnson support-floor
  optimizers. The 4,055,870-row substitution atlas has uniform exact-coset
  central-over-outer margin at least 380. Pair modules and prime order-eleven
  tournaments remain open, and no actual-maximizer centrality is claimed.
source: codex-frontier-synthesis-creative-20260825aa
depends_on:
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4131-strong-tournament-centrality-through-order-eight
related:
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4137-strong-tournament-centrality-complete-order-ten
script: 04-computation/tournament_order_eleven_large_module_centrality_thm4144.py
engine: 04-computation/tournament_order_eleven_large_module_centrality_thm4144_engine.cpp
output: 05-knowledge/results/tournament_order_eleven_large_module_centrality_thm4144.out
independent_audit_script: 04-computation/tournament_order_eleven_large_module_centrality_thm4144_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_order_eleven_large_module_centrality_thm4144_independent_audit.out
script_sha256: e72c8710a5de2aad77e7710880857f1b47148adf82b6fb4278417a30a982edfd
engine_sha256: 03dfb092830313d71d41c35dcf7fb7d300db9850bbd52b1291ace7e674bf63bc
output_sha256: ae29c3253ae2e21c4119e5a4fc6720fd55117b082a6c285a58303fa4ce0bff2e
semantic_sha256: 4c7ca82098f748137cc99e52a10b47c5eb59b405746c9daf20b664701bed10c9
independent_audit_script_sha256: 06b53688edcc3f37f0c3f583a7af5fe5eda58ebad663da9302291342415868d6
independent_audit_output_sha256: 6020c78ab6c6fe1f52704453b466c6fe1eff08489fee07948a62653f2fd8d28f
independent_semantic_sha256: 6516491b7097782ed50e200bdde26b032109f7275e95c5d5cd8c506c3a292829
hash_basis: raw LF bytes
primary_audit: >
  PASS. Homebrew nauty 2.9.3 gentourng streams supply every required strong
  quotient and arbitrary block representative. A separate exact-integer C++
  engine checks all 4,055,870 marked substitutions, strongness, the rational
  gate, and every exact layer coset. The Python replay freezes all 17 shards
  and independently reconstructs the 34 shard extrema through the audited
  THM-4131 response implementation. Normal, optimized, hash-seeded, frozen,
  and full generator replays agree.
independent_audit: >
  ACCEPT. A separate standard-library C++ implementation reads and pins the
  quotient/block generator streams, reconstructs response values by an
  ordered exposed-gap DP, and checks all 4,055,870 rows. It finds no rational
  boundary tie, no rational failure, and no exact-coset failure, with the same
  counts, extrema, minimum margin, hostile controls, and residual strata.
  Its full-universe semantic digest is independent of the primary certificate;
  optimized and pinned full replays agree, with additional unoptimized and
  sanitizer controls on the sharp q=5 regime.
---

# THM-4144 -- order-eleven large-module centrality

**FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The order-twelve counterexample of THM-4133 is built by substitution. The
first unresolved order is eleven, so the faithful inverse question is not
whether random order-eleven tournaments resemble that witness. It is whether
its homogeneous-module mechanism can occur one order earlier while retaining
strongness and a noncentral support floor.

## 1. Statement and inheritance

Call a nonempty proper vertex set `M` of a tournament homogeneous if every
vertex outside `M` has the same orientation to every vertex of `M`.

> **Theorem.** Let `T` be a strong tournament of order eleven. If
> `T` has a homogeneous module `M` with
>
> ```text
> 3 <= |M| <= 9,                                         (1)
> ```
>
> then every rational THM-4128 Johnson support-floor optimizer and every
> exact-coset THM-4123 optimizer lies on a central layer `t in {-1,1}`.
> For every such `T`, the best central exact-coset floor exceeds the best
> outer floor by at least `380`.

The closest proved mechanism is THM-4137's complete positive order-ten
census. The canonical hostile is THM-4133's strong modular counterexample at
order twelve. The corrected near miss is deleting a vertex from that witness:
one deletion retains normalized tilt greater than one, but loses strongness.
The least-used sidecar is the distinguished quotient vertex; quotient and
block isomorphism classes without it do not specify a substitution row.

This is not a classification of order eleven and says nothing about actual
response maximizers.

## 2. Why the substitution atlas covers the claimed stratum

Contract `M` to one marked vertex `v`. Homogeneity makes the result a
tournament `Q` of order

```text
q=12-|M|.                                                (2)
```

The quotient is strong. Indeed, a nontrivial one-way cut of `Q` lifts after
replacing `v` by `M` to a one-way cut of `T`, contradicting strongness.

Conversely, let `Q` be strong, mark `v in V(Q)`, and substitute an arbitrary
tournament `B` for `v`. The resulting lexicographic substitution

```text
Q(v <- B)                                                (3)
```

is strong even when `B` is not. Paths in `Q` between vertices outside the
marked block lift literally. A path entering or leaving `v` lifts through
any chosen block vertex because all cross-block arcs have the quotient
orientation. Finally, a directed closed walk through `v` in the strong
quotient supplies paths between any two vertices inside `B`, in either
direction, without using an internal `B` path.

Therefore isomorphism representatives of all strong `Q`, all tournaments
`B`, and every distinguished `v` cover every isomorphism type satisfying
`(1)`. Different triples may give isomorphic targets; the following are
construction rows, not distinct target-class counts:

| `q` | `|M|` | strong `q` classes | all block classes | marked rows |
|---:|---:|---:|---:|---:|
| 3 | 9 | 1 | 191,536 | 574,608 |
| 4 | 8 | 1 | 6,880 | 27,520 |
| 5 | 7 | 6 | 456 | 13,680 |
| 6 | 6 | 35 | 56 | 11,760 |
| 7 | 5 | 353 | 12 | 29,652 |
| 8 | 4 | 6,008 | 4 | 192,256 |
| 9 | 3 | 178,133 | 2 | 3,206,394 |
| **total** | | | | **4,055,870** |

A module of size ten cannot occur in a strong tournament because its
two-vertex quotient is not strong. Thus after this atlas the exact modular
residual is size two; tournaments with no nontrivial module form the prime
residual.

## 3. Exact response and centrality gates

For each constructed row, the Start/End/exposed-gap subset DP of THM-4131
reconstructs all Hamilton and one-ear responses. Write its directed exposure
packet as

```text
(H,W,D_4,C_hd).                                          (4)
```

THM-4128 puts every rational layer floor on one strict quadratic with vertex

```text
theta=4C_hd/D_4                                          (5)
```

at order eleven. Its admissible imbalance grid is
`{-9,-7,...,7,9}`. Hence all rational optimizers are central precisely when

```text
rho=2|C_hd|/D_4 < 1.                                    (6)
```

The strict inequality excludes a central/outer boundary tie.

For the exact-coset sidecar, let `F(S)` be the literal response on one
cardinality layer, let `N` be the layer size, and put

```text
A=sum_S F(S),                 B=sum_S F(S)^2,
J=(B-HA)/(A-NH).                                        (7)
```

If `a` is one layer value and

```text
g=gcd_S(F(S)-a),                                        (8)
```

then THM-4123's exact layer-coset floor is

```text
a+g ceil((J-a)/g),                                      (9)
```

with the constant-layer convention when `g=0`. The audit reconstructs all
`2^11` values, all ten nonconstant layers, and `(7)--(9)`. It then compares
the best central value at cardinalities five and six with every outer layer.

The primary C++ kernel stores `2W,4D_4,4C_hd` as integers and normalizes only
for its transcript. Its rational failure test is the exact cross-product

```text
2|4C_hd| >= 4D_4,                                      (10)
```

not a floating-point comparison.

## 4. Complete finite-exact atlas

The seventeen shards give

```text
construction rows                         4,055,870
rows independently rechecked strong       4,055,870
rational centrality failures                       0
exact-coset centrality failures                    0
minimum central-minus-outer coset margin         380. (11)
```

The sharp normalized tilt is

```text
rho_max=3,201,028/5,617,819 < 1.                       (12)
```

It occurs in the `q=5, |M|=7` regime. One exact marked presentation is

```text
quotient label:       1101111101
distinguished vertex: 1
block label:          110100110101101110111

target label:
1111111101110100111110101111101111110111111111111111101. (13)
```

The block is the regular order-seven tournament. The exact packet and layer
sidecars are

```text
(H,W,D_4,C_hd)=(6615,102123,3185303373,-907491438),
rational optimizer=(-1),          exact-coset optimizer=(-1),
actual optimizer=(-5),            coset margin=1494.    (14)
```

Thus `(14)` is also a hostile against confusing support-floor centrality
with actual-maximizer centrality.

The minimum margin is attained, among other repeated presentations, at

```text
label=1110111111111111111111111111111111111101111111111111101,
(H,W,D_4,C_hd)=(243,11871,42303227,-878818),
rho=1757636/42303227,
rational optimizer=exact-coset optimizer=actual optimizer=(-1). (15)
```

The primary certificate freezes every shard count, both extrema of every
shard, and zero failure counts. The audited THM-4131 Python implementation
independently reconstructs the complete packet, all layer cosets, and all
three optimizer types on those `34` extrema.

## 5. The deletion hostile and the strongness barrier

Delete the quotient vertex `c` from THM-4133's order-twelve witness. The
resulting order-eleven label is

```text
1111111101111100011111100111111011111111111111111111111. (16)
```

It has a sink, is not strong, and satisfies

```text
(H,W,D_4,C_hd)=(9253,147483,5815378668,-4396377498),
rho=732729583/484614889 > 1,
rational optimizer=exact-coset optimizer=(-3).          (17)
```

Thus deleting a vertex does preserve the noncentral scalar mechanism, but
destroys the hypothesis. Exactly ten one-edge reversals of `(16)` restore
strongness. Every one is central; the largest repaired tilt is

```text
5121379056/12507930943 < 1,                              (18)
```

and its exact-coset central-over-outer margin is `4078`. Strongness is
load-bearing rather than a cosmetic filter.

## 6. Typed connection and loss ledger

```text
source:       (Q,v,B), Q strong, v distinguished, B arbitrary
target:       T=Q(v <- B), a strong order-eleven tournament
map:          lexicographic substitution / module expansion
preserved:    full labelled adjacency, Hamilton responses, exposure packet,
              rational floors and every layer coset
destroyed:    uniqueness of presentation and target-class multiplicity
sidecar:      marked quotient vertex and module support
hostile:      the nonstrong deletion (16) retains rho>1
decisive test: all seven quotient/block regimes plus exact layer replay.
```

The destroyed multiplicity is harmless for a universal property but forbids
calling `4,055,870` an isomorphism-class count.

## 7. Scope and replay

The coverage proof and the two accepted full-universe audits establish the
theorem. The surviving order-eleven problem is exactly

```text
homogeneous pairs + prime tournaments.                  (19)
```

This theorem blocks the THM-4133 substitution mechanism through every module
of size at least three. It does not classify the pair-module or prime
strata, prove actual-maximizer centrality, establish an interval law, or
extend to order twelve.

The default exact certificate replay is

```text
python3 -B 04-computation/tournament_order_eleven_large_module_centrality_thm4144.py
python3 -B -O 04-computation/tournament_order_eleven_large_module_centrality_thm4144.py
PYTHONHASHSEED=271828 python3 -B 04-computation/tournament_order_eleven_large_module_centrality_thm4144.py
```

The full nauty/C++ replay is

```text
python3 -B 04-computation/tournament_order_eleven_large_module_centrality_thm4144.py --full
```

The independent full-universe replay is

```text
clang++ -std=c++17 -O3 \
  04-computation/tournament_order_eleven_large_module_centrality_thm4144_independent_audit.cpp \
  -o /tmp/thm4144-independent-audit
/tmp/thm4144-independent-audit
```

Both full replays pass. The independent implementation also proves by exact
cross-product that no row lies on the rational boundary `rho=1`. **QED.**
