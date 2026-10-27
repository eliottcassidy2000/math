---
id: HYP-9421
title: "After THM-4811-4813: (1) the strong tournaments without a black descent are exactly C3, R5 and the doubled-stack family E(M; Z) (refuted as stated; resolved by THM-4811 Theorem A with the larger family E); (2) the black distance to the non-strong locus never exceeds lambda + 1 (outside C3, R5; proved); (3) the presentation-count spectrum: every odd v other than 7 and 21 is tc(T) for some tournament on at most 10 vertices, for v <= 2001 (PRE-REGISTERED); (4) the native ROOT forcing of a Collatz orbit is a tail property: sup F(s) over primitive seeds s < 10^8 is below 0.1127503 (PRE-REGISTERED); (5) strong presentation fibres are finite with order O(log^2 v), and the black-core connectivity reduction; (6) every strong tournament other than C3, R5 is exactly delta^0 black moves from a singleton plus a strong core"
status: "(1) REFUTED as stated and RESOLVED (THM-4811 Theorem A: the exceptions are R5 and the larger family E); (2) PROVED (THM-4811 Theorem B); (3), (4) PRE-REGISTERED (commit 8ab6b3e358) and CONFIRMED; (5) first clause PROVED in THM-4812 (f), second OPEN; (6) OPEN (FINITE-EXACT through order 9 and the order-10 classes with lambda < delta^0)"
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

**Conjecture (as first stated).** A strong tournament has no black word starting in a minimum out-cut and ending outside
it if and only if it is `C3`, `R5`, or a member of `E(M; Z)`. Evidence then: complete through order 10.

**Outcome: REFUTED as stated, RESOLVED.** The order-11 census (all 884,464,590 strong classes) has exactly one exception,
`E_2(S5)`, which has odd order and a middle vertex, so it is not in `E(M; Z)` (the conjecture predicted none at order 11).
THM-4811 Theorem A proves the corrected statement at every order: the exceptions (`n >= 4`) are exactly `R5` and the
family `E = {E_t(S, x0)}` of bottom layers with `t - 1` middle vertices; `E(M; Z)` is its subfamily with `t = 1` and an
initial-segment pattern. The prediction "exactly two classes at order 12" was also wrong: there are three
(`E_3(S5)` and the two of `E(M; Z)`).

## (2) One extra move suffices

**Conjecture.** Every strong tournament other than `C3` and `R5` reaches the non-strong locus by at most `lambda + 1`
black moves, with `lambda + 1` needed exactly for the members of `E(M; Z)`. Evidence: order <= 10 complete; six family
members through order 14 have black distance exactly 2.

**Outcome: PROVED (THM-4811 Theorem B)** with the corrected family: the black distance is `lambda` except `R5` and the
members of `E`, where it is `2 = lambda + 1`.

## (3) The presentation spectrum (PRE-REGISTERED)

THM-4812 proves that `tc` never takes the values 7 and 21. **Prediction:** for every odd `v <= 2001` with
`v` not in `{7, 21}` there is a tournament on at most 10 vertices with `tc = v` (strong components multiply, so this
concerns products of strong values). Test: a sweep of all classes of orders 1-10 (`gentourng`), `tc = H/|Aut|` by
Held-Karp and nauty. Refutation: any odd `v <= 2001`, `v` not 7 or 21, missing.

**Outcome (CONFIRMED, `tc_spectrum_summary.out`).** Over all 9,733,056 + 191,536 + ... classes of orders 1-10 the
attained odd values `<= 6001` are exactly the odd numbers other than 7 and 21 (no even value occurs). So, with THM-4812,
the presentation spectrum below 6002 is exactly `{odd} \ {7, 21}`. The `tc = 1` counts by order, `1, 1, 2, 3, 4, 6, 9,
13, 19, 28`, are the `K1`/`C3`-stack counts (THM-4741; A000930), and the `tc = 3` counts `1, 2, 3, 6, 11, 18` (orders
5-10) are the one-`R5` stacks of HYP-9352 (1).

## (4) The forcing tail (PRE-REGISTERED)

For a primitive seed `s` let `F(s)` be the total native price `sum lambda` along its compressed orbit to 1 (THM-4813 (A)).
Below `2 10^6` the maximum is `0.112749293418` (at `s = 1056939`), shared to many digits by seeds whose orbits share a
tail. **Prediction:** `max F(s) < 0.1127503` over every primitive `s < 10^8`, i.e. the excess over the small-seed maximum
stays below `10^-6`, because the price of a step at numerator `X` is `2^a/(6X + 2^a)`. Refutation: any `s < 10^8` with
`F(s) >= 0.1127503`.

**Outcome (CONFIRMED, `forcing_tail_1e8.out`).** `max F(s) = 0.112749570730` at `s = 78027147` (longest native orbit
356 edges), `2.8 10^-7` above the small-seed maximum.

## (5) Fibres and connectivity

The first clause (for every odd `v` the strong classes with `tc = v` have order at most an explicit `N(v) =
O(log^2 v)`) is PROVED in THM-4812 (f). Open: compute every strong `tc`-fibre for `v <= 101`. Separately, THM-4811 (B)-(C)
reduce black-core connectivity (HYP-9348 (b)) to connecting the non-strong basins through strong classes; no claim.

## (6) Singleton targets

**Conjecture.** Every strong tournament other than `C3` and `R5` is exactly `delta^0 = min_v min(s_v, n - 1 - s_v)` black
moves from a singleton plus a strong core (and from the singleton locus). The lower bound is the singleton cut
potential of THM-4811 (A) (each move changes every `|s_v - (n-1)/2|` by at most one). Evidence: every strong class of
orders 6-9 (184,529 classes; `singleton_distance.cpp`, exact depth-first search over the moves that raise the maximal
deviation); there `delta^0 = lambda`, so the statement is the singleton refinement of THM-4811. The regular case is the
concurrent basin theorem (aliquot_power_124) and klein's HYP-9422; order 10, where `lambda < delta^0` first occurs, is
the first informative test.
