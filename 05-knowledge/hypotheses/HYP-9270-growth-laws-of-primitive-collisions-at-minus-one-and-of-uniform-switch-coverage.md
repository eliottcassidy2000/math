---
id: HYP-9270
title: "Growth of primitive collisions at the -1 anchor: the number N(E) of values K with increasing {2,3}-representations of two distinct last exponents, all exponents <= E, satisfies N(E) = 2^E E^(-beta + o(1)) with beta about 0.47, so N(E) = o(2^E). (The former part (b), a logarithmic growth law for the switch coverage, is retired: it contradicts the bound c(E) <= 3/8, and THM-4640 (8) proves c(E) -> 3/8.)"
status: >
  OPEN (part (a)). NUMERICAL evidence: the exact census of THM-4640 (7), extended independently by audit A to E = 29.
  N(E) sqrt(E)/2^E = 0.025641, 0.025642, 0.025661, 0.025689, 0.025715 at E = 25, ..., 29, and the local exponent
  ln(2 N(E-1)/N(E)) / ln(E/(E-1)) is 0.507, 0.499, 0.481, 0.470, 0.471. Part (b) RETIRED 2026-10-14 by audit A (M2):
  REFUTED as stated (increments ~ kappa/E would make c(E) unbounded, but c(E) <= 3/8), and its open question is answered
  by THM-4640 (8)(d): c(E) -> 3/8, and with the root switch the coverage tends to 1/2 (MISTAKE-606).
session: mac-mini-2026-10-14-wilson
source: 05-knowledge/results/repair_algebra_minus_one_and_wilson_dictionary_20261014.md
related:
  - THM-4640 (the repair algebra at -1; the census and the coverage theorem), THM-4555 (uniform switches), THM-4581 (Haar coalescence)
  - 04-computation/experiments/repair_algebra_20261014/audit_A/a2_census.py, a2b_census28.py, a2c_census29.py (+ .out)
---

# HYP-9270 — growth of primitive collisions at the −1 anchor

## Statement

* **(a) Collision count.** `N(E) = 2^E · E^(−β + o(1))` with `β ≈ 0.47`. In particular `N(E) = o(2^E)`: a vanishing but only polynomially small fraction of all `2^(E−1)` canonical representations meets another representation of the same value with a different last exponent.
* **(b) Retired.** The former statement, `c(E) − c(E − 2) ~ κ/E` with `κ ≈ 0.072`, is false as a law: it would make `c(E)` unbounded, while `c(E) ≤ 3/8` for every `E`. THM-4640 (8) now proves `c(E) → 3/8`, and the coverage with the root switch tends to `1/2` (audit A, M2).

## Data (FINITE-EXACT counts)

| `E` | 14 | 18 | 22 | 24 | 26 | 28 | 29 |
|---|---|---|---|---|---|---|---|
| `N(E)` | 132 | 1,666 | 23,037 | 87,836 | 337,478 | 1,303,189 | 2,563,623 |
| `N(E)·√E/2^E` | 0.0301 | 0.0270 | 0.0258 | 0.02565 | 0.02564 | 0.02569 | 0.02572 |
| local exponent `β(E)` | 1.360 | 0.881 | 0.603 | 0.531 | 0.499 | 0.470 | 0.471 |

* The local exponent fell steadily through `E = 26` and has hovered near 0.47 for `E = 28, 29`. The `√E` normalisation is not constant: it bottoms out near `E = 25` and rises after.

## Why

* Each primitive collision is a vanishing `{2,3}`-unit sum with at least three differing positions (THM-4640 (5)). Matching two representations needs one coincidence in the 2-adic, 3-adic and length bookkeeping, which suggests a local-limit factor `E^(−1/2)` on top of the `2^E` contexts.
* The data give a slightly smaller exponent than that heuristic, and about 40% of the pairs are degenerate unit sums (THM-4640 (6)). So the local-limit picture is at best approximate; this is a heuristic, not a derivation.

## What would refute it

* The local exponent leaving `[0.4, 0.55]` for good beyond `E = 29`.
* `N(E)/2^E` not tending to 0.
