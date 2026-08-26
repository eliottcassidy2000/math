---
id: THM-4133
title: "Strong cyclic-substitution counterexample to all-order Johnson centrality"
status: >
  VERIFIED-EXACT + INDEPENDENTLY AUDITED EXPLICIT COUNTEREXAMPLE. A strong tournament of order twelve,
  obtained by substituting R_9 minus one vertex into one block of a fixed
  strong five-vertex quotient, has normalized THM-4128 exposure tilt
  53092739331/40435524866>1. Its unique rational and exact-coset support-floor
  optimizer is the noncentral layer t=-2; the central t=0 exact-coset floor
  loses by 2,224. Thus all-order strong rational and exact-coset centrality is
  REFUTED. THM-4137 later proves order ten positive, but order eleven remains
  unclassified, so no globally minimal counterexample order is claimed.
source: codex-frontier-synthesis-creative-20260825q
depends_on:
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4131-strong-tournament-centrality-through-order-eight
related:
  - THM-3235-tournament-substitution-blowup-scaling-and-decimation
  - THM-4137-strong-tournament-centrality-complete-order-ten
script: 04-computation/tournament_strong_cyclic_substitution_centrality_counterexample_thm4133.py
output: 05-knowledge/results/tournament_strong_cyclic_substitution_centrality_counterexample_thm4133.out
script_sha256: 7bd4c518464d4c48baf9cb9c1c8c2012a9f79f029c8e07141c0e51c338ffeeb2
output_sha256: 52d6c229b46ac1f38afb61d54073eac2400757205f38f26a0c646b7a8cdf5bf1
semantic_sha256: 2d41d1a1bb6f8a936c6f8104d149cba898ae51d17371981dbf2d1332643d2873
independent_audit_script: 04-computation/tournament_strong_cyclic_substitution_centrality_counterexample_thm4133_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_strong_cyclic_substitution_centrality_counterexample_thm4133_independent_audit.out
independent_audit_script_sha256: c312eab367d2ace57ecc87b56383ed534db3bf3f119107d99d4bb57945e9c201
independent_audit_output_sha256: 00b906d460c74f1906b52124ae8e8ce7d6f0678e4d94af2f2356294e830fbbf1
hash_basis: raw LF bytes
primary_audit: >
  PASS. The construction is generated from its quotient arcs and cyclic
  difference rule, not a searched bit string. The pinned THM-4131 exact
  evaluator verifies strongness, every response, exposure packet, rational
  floor, layer lattice, exact-coset floor and actual maximum. The orders
  6,8,10 family prefix supplies positive controls, while the frozen order-12
  adjacency and all layers supply a direct hostile. Normal, optimized and
  hash-seeded replays agree with the frozen output.
independent_audit: >
  ACCEPT. A clean-room C++/GMP contracted-good/reversed-block evaluator
  imports no primary code. It reconstructs the quotient and cyclic deletion,
  checks strongness, recomputes all 4,096 responses by both a mask recurrence
  and direct arc sums, and independently recovers every packet, layer, lattice
  and optimizer. It also audits the larger regular-block family and the three
  central proper-prefix controls. O3 and O0 builds have byte-identical output.
---

# THM-4133 -- strong cyclic-substitution centrality counterexample

**VERIFIED-EXACT + INDEPENDENTLY AUDITED EXPLICIT COUNTEREXAMPLE.**

THM-4131 proves strong rational and exact-coset Johnson centrality through
order eight. Its worst classes at orders seven through nine suggest a common
substitution object. Following that object past its finite census refutes the
all-order extension.

## 1. Construction

Let `Q` be the tournament on the five ordered blocks `(a,B,c,d,e)` whose arcs
are

```text
a -> B,c,e;       B -> c,d,e;       c -> d;
d -> a,e;         e -> c.                                  (1)
```

These ten arcs specify every pair. The quotient is strong. Let `R_9` be the
regular cyclic tournament on `Z/9Z`, oriented by

```text
i -> i+j,                    j=1,2,3,4.                    (2)
```

Delete vertex zero and substitute the remaining eight-vertex tournament for
block `B` in `(1)`, with every inter-block arc made complete in its displayed
direction. Call the resulting order-twelve tournament `T_12`.

The quotient is strong and contains a directed route from `B` back to `B`
through other blocks. Hence any two vertices of the substituted block, as
well as vertices in distinct blocks, are mutually reachable; `T_12` is
strong. The exact adjacency bitmasks, in construction order, are

```text
(3070,3644,3704,3824,4064,4032,3970,3846,3598,1024,2049,512). (3)
```

## 2. Exposure packet and failure

With THM-4128's directed exposure coordinates, exact Hamilton/exposed-gap
dynamic programming gives

```text
(H,W,D_4,C_hd)
 =(27759,506085,80871049732,-23596773036).                (4)
```

Since `n=12` is even, its normalized centrality load is

```text
rho=(n-3)|C_hd|/(2D_4)
   =53092739331/40435524866
   >1.                                                     (5)
```

By THM-4128, `(5)` makes the rational support-floor optimizer noncentral.
The complete exact row is sharper:

```text
argmax J_m={t=-2},           argmax L_m={t=-2},
argmax max_(|S|=m)F_T(S)={t=-6}.                           (6)
```

The central and winning exact-coset floors are

```text
L(0)=350727,                 L(-2)=352951,                 (7)
```

so the best central layer loses by exactly `2,224`. Therefore strong
rational centrality and strong exact-coset centrality both fail in general.
**QED.**

## 3. Why this hostile was visible

Replace `R_9` in the construction by `R_r`, delete zero, and take
`r=3,5,7,9`. The exact normalized tilts and optimizers are

| `r` | order | `rho` | rational / coset optimizer |
|---:|---:|---:|---:|
| 3 | 6 | `13/557` | `0 / 0` |
| 5 | 8 | `14809/29733` | `0 / 0` |
| 7 | 10 | `97613103/107670631` | `0 / 0` |
| 9 | 12 | `53092739331/40435524866` | `-2 / -2` |

Thus order twelve is the first failure **within this family**. The order-eight
member is exactly a THM-4131 worst-tilt packet, and the order-ten member uses
more than ninety percent of the centrality budget. Substitution retains the
large cyclic block's near-regular internal exposure while the asymmetric
quotient creates a coherent field-degree correlation; the denominator
`D_4` no longer absorbs that correlation at `r=9`.

This identifies the first failed implication in the tempting all-order claim:
strong connectivity constrains reachability but does not bound
`(n-3)|C_hd|` by `2D_4`. The strongest survivor is the exact THM-4128
criterion itself, now serving as a counterexample generator.

## 4. Scope and replay

The theorem refutes only all-order centrality of the THM-4128 rational and
THM-4123 exact-coset **support floors**. It does not itself classify orders ten
or eleven, prove that twelve is globally minimal, classify all substitution
failures, or settle actual response maxima or interval completeness. THM-4137
later classifies order ten positively; order eleven remains open.

Run

```text
python3 -B 04-computation/tournament_strong_cyclic_substitution_centrality_counterexample_thm4133.py
python3 -B -O 04-computation/tournament_strong_cyclic_substitution_centrality_counterexample_thm4133.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_strong_cyclic_substitution_centrality_counterexample_thm4133.py
```

All three streams agree with the frozen output. **QED.**

The independent replay is

```text
clang++ -O3 -std=c++17 \
  04-computation/tournament_strong_cyclic_substitution_centrality_counterexample_thm4133_independent_audit.cpp \
  -I/opt/homebrew/opt/gmp/include -L/opt/homebrew/opt/gmp/lib \
  -lgmpxx -lgmp -o /tmp/thm4133_independent
/tmp/thm4133_independent
```

An `-O0` build produces the same frozen output. Besides the exact packet, it
checks the cut formula against all `4,096` masks and fingerprints the response
vector as `df5f7d507a5f27e2` (FNV-64). **QED.**
