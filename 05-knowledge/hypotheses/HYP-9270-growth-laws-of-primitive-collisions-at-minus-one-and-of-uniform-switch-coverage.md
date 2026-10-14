---
id: HYP-9270
title: "Growth laws at the -1 anchor: the number N(E) of primitive collisions (values K with increasing {2,3}-representations of two distinct last exponents, all exponents <= E) satisfies N(E) ~ C 2^E / sqrt(E) with C ~ 0.0255; and the Haar coverage c(E) of odd 2-adic sources by downward uniform switches from these collisions grows like 0.036 ln E (E times the increment per two steps stays near 0.072), so it increases without a visible plateau through E = 22; whether it saturates below 1/2 (the mass of sources with a trailing one) is left open"
status: >
  OPEN. NUMERICAL evidence: the exact census of THM-4640 (7)-(8): N(E) sqrt(E)/2^E = 0.0282, 0.0270, 0.0261, 0.0258, 0.0257 at
  E = 16, 18, 20, 22, 24; c(E) = 0.0117, 0.0205, 0.0273, 0.0328, 0.0375, 0.0417, 0.0454, 0.0489 at E = 8, ..., 22, with
  E (c(E) - c(E-2)) = 0.079, 0.075, 0.072, 0.071, 0.071, 0.072, 0.074 at the step midpoints. Not audited.
session: mac-mini-2026-10-14-wilson
source: 05-knowledge/results/repair_algebra_minus_one_and_wilson_dictionary_20261014.md
related:
  - THM-4640 (the repair algebra at -1; the census and coverage data), THM-4555 (uniform switches), THM-4601 (the +1 barrier)
---

# HYP-9270 — growth laws at the −1 anchor

## Statement

* **(a) Collision count.** `N(E) ~ C·2^E/√E` with `C ≈ 0.0255`. A positive but slowly vanishing fraction of all `2^E` representations meets another representation of the same value with a different last exponent.
* **(b) Switch coverage.** `c(E) − c(E − 2) ~ κ/E` with `κ ≈ 0.072`, so `c(E) ≈ 0.036 ln E + const` while overlaps are negligible.
* **Open question (not part of the hypothesis).** Does `c(E)` saturate below `1/2`, the Haar mass of sources with at least one trailing one? THM-4601 (ii) blocks deletion merges inside the 2-run of the residual first-reset-2 branch (`D ≤ 8000`). It does not exclude switches whose shorter word extends past the run, so it does not decide this.

## Why

* **(a)** Each primitive collision is a vanishing `{2,3}`-unit sum with at least three differing positions (THM-4640 (5)). Matching two representations needs one coincidence in the 2-adic, 3-adic and length bookkeeping, which suggests a local-limit factor `E^(−1/2)` on top of the `2^E` contexts. This is a heuristic, not a derivation.
* **(b)** The switch cylinder of a collision at exponent scale `E` has Haar mass about `2^(−E)`, and the collisions number about `2^E/√E`. The fitted `1/E` decay of the increments is steeper than `1/√E`, so short and long cylinders overlap strongly.

## What would refute it

* **(a)** `N(E)·√E/2^E` drifting away from `0.025` as `E` grows beyond 24.
* **(b)** `E·(c(E) − c(E − 2))` drifting away from `0.072`, either way, beyond `E = 22`.
