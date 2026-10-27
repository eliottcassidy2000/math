---
id: HYP-9421
title: "After THM-4811-4813: (1) the strong tournaments without a black descent are exactly C3, R5 and the doubled-stack family E(M; Z); (2) the black distance to the non-strong locus never exceeds lambda + 1 (outside C3, R5); (3) the presentation-count spectrum: every odd v other than 7 and 21 is tc(T) for some tournament on at most 10 vertices, for v <= 2001 (PRE-REGISTERED); (4) the native ROOT forcing of a Collatz orbit is a tail property: sup F(s) over primitive seeds s < 10^8 is below 0.1127503 (PRE-REGISTERED); (5) strong presentation fibres are finite with order O(log^2 v), and the black-core connectivity reduction"
status: "OPEN (1), (2), (5) second clause; PRE-REGISTERED (3), (4) before the runs (this commit); (5) first clause PROVED in THM-4812 (f)"
session: mac-mini-2026-10-27-potentials
source: 05-knowledge/results/dipole_calculus_black_distance_tc_gaps_collatz_numerator_20261027.md
related:
  - 01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md
  - 01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md
  - 01-canon/theorems/THM-4813-collatz-numerator-coordinate-forcing-cocycle-gauss-sum-intertwiner-and-sharp-spectral-exception-budgets.md
  - 05-knowledge/hypotheses/HYP-9348-endpoint-layers-black-core-connectivity-and-the-next-zero-heavy-layers.md
  - 05-knowledge/hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md
---

# HYP-9421 — black-descent exceptions, the presentation spectrum, the forcing tail

## (1) Classification of black-descent failures

**Conjecture.** A strong tournament has no black word starting in a minimum out-cut and ending outside it if and only if
it is `C3`, `R5`, or a member of `E(M; Z)` (THM-4811 (C2)). Evidence: complete through order 10 (exactly `C3`, `R5`,
`E(C3; K1)`); members built through order 16 behave as predicted; the variant whose `M` starts with `K1` is not
exceptional. The family has only even orders, so the conjecture predicts no exception at order 11 and exactly two
classes at order 12.

## (2) One extra move suffices

**Conjecture.** Every strong tournament other than `C3` and `R5` reaches the non-strong locus by at most `lambda + 1`
black moves, with `lambda + 1` needed exactly for the members of `E(M; Z)`. Evidence: order <= 10 complete; six family
members through order 14 have black distance exactly 2.

## (3) The presentation spectrum (PRE-REGISTERED)

THM-4812 proves that `tc` never takes the values 7 and 21. **Prediction:** for every odd `v <= 2001` with
`v` not in `{7, 21}` there is a tournament on at most 10 vertices with `tc = v` (strong components multiply, so this
concerns products of strong values). Test: a sweep of all classes of orders 1-10 (`gentourng`), `tc = H/|Aut|` by
Held-Karp and nauty. Refutation: any odd `v <= 2001`, `v` not 7 or 21, missing.

## (4) The forcing tail (PRE-REGISTERED)

For a primitive seed `s` let `F(s)` be the total native price `sum lambda` along its compressed orbit to 1 (THM-4813 (A)).
Below `2 10^6` the maximum is `0.112749293418` (at `s = 1056939`), shared to many digits by seeds whose orbits share a
tail. **Prediction:** `max F(s) < 0.1127503` over every primitive `s < 10^8`, i.e. the excess over the small-seed maximum
stays below `10^-6`, because the price of a step at numerator `X` is `2^a/(6X + 2^a)`. Refutation: any `s < 10^8` with
`F(s) >= 0.1127503`.

## (5) Fibres and connectivity

The first clause (for every odd `v` the strong classes with `tc = v` have order at most an explicit `N(v) =
O(log^2 v)`) is PROVED in THM-4812 (f). Open: compute every strong `tc`-fibre for `v <= 101`. Separately, THM-4811 (B)-(C)
reduce black-core connectivity (HYP-9348 (b)) to connecting the non-strong basins through strong classes; no claim.
