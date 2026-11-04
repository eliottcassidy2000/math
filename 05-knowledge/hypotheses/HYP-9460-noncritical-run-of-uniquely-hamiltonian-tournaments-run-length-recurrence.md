---
id: HYP-9460
title: "The non-critical run of a uniquely Hamiltonian tournament: pre-registered exhaustive checks at orders 14-18 of the structure proved for THM-4868 (local roundness; one run from a sole entry to a sole exit; #cut arcs + #non-critical vertices = n + 1; the refined Fibonacci counts by numbers of sole exits and sole entries; asymmetry), and the conjectured run-length recurrence N(n, r) = N(n - 1, r) + 2 N(n - 1, r - 1) - N(n - 2, r - 1), equivalently T(n, k) = T(n - 1, k - 1) + 2 T(n - 1, k) - T(n - 2, k - 1) for HYP-9458's cut-arc table"
status: CLOSED (kind-pasteur-2026-11-04-S10) - P1-P8 CONFIRMED at orders 14-18 (pre-registered in a8db6218f9 before the run); P7 (run-length recurrence) then PROVED as THM-4868 (H)
session: kind-pasteur-2026-11-04-S10
date: 2026-11-04
related: [THM-4868, THM-4867, HYP-9458, HYP-9457, HYP-9459, THM-4866, THM-4862, THM-4846]
note: 05-knowledge/results/uniquely_hamiltonian_noncritical_run_20261104.md (to be written after the run)
script: 04-computation/experiments/uh_noncritical_run_20261104.c
---

# HYP-9460 - the non-critical run of a uniquely Hamiltonian tournament

## Notation

- A tournament T of order n is **uniquely Hamiltonian** (UH) if it has exactly one directed Hamiltonian cycle C. Write
  v+ and v- for the successor and predecessor of v on C.
- v is **non-critical** if T - v is strong; NC(T) is the set of non-critical vertices and r = |NC(T)|.
- A **sole exit** is a vertex of out-degree 1, a **sole entry** a vertex of in-degree 1.
- An arc of C is a **cut arc** if deleting it leaves T non-strong (HYP-9458's notion; by THM-4867 these are the single
  backward arcs of one-arc cuts).
- N(n, r) is the number of UH classes of order n with r non-critical vertices; T(n, k) (HYP-9458, opus-S74) the number
  with exactly k cut arcs on C.

## What the session has proved (to be recorded in THM-4868 after the run and the audit)

For UH T of order n >= 4 (proofs in the note):
- (A) every non-critical v beats exactly the d+(v) vertices following it on C (local roundness), and if d+(v) >= 2 then
  v+ is non-critical (dually, d-(v) >= 2 makes v- non-critical);
- (B) NC(T) is one run of C, of length 2 <= r <= n - 1, starting at a sole entry and ending at a sole exit; the vertex
  after a sole exit and the vertex before a sole entry are critical;
- (C) an arc of C is a cut arc iff it is not internal to the run, so #cut arcs + r = n + 1, i.e. T(n, k) = N(n, n + 1 - k);
- (D) with Garey's count F(2n - 6) (CITED): the classes with (#sole exits, #sole entries) = (1,1), (1,2), (2,1), (2,2)
  number F(2n - 10) + 1, F(2n - 9) - 1, F(2n - 9) - 1, F(2n - 8) + 1 (n >= 5);
- (E) T has no non-trivial automorphism.
(A) with a "round" second cycle gives a proof of Shelor's theorem (THM-4867) that does not use Shelor's grammar.

## Exploratory evidence (before this commit)

`uh_noncritical_run_20261104 13 12 12` (orders 4-13, all F(2n - 6) classes, 0.3 s), output
`04-computation/experiments/uh_noncritical_run_20261104_explore.out`:
- (A)-(E) hold in every class; one-arc cuts exist; at orders <= 12 an exact Hamiltonian-cycle count (independent subset
  DP) gives 1 for every class, and a brute-force search over all vertex subsets agrees with the interval search;
- the cut-arc column reproduces HYP-9458's T(n, k) rows through order 13 from an independent program;
- the run-length rows N(n, r), r = 2..n - 1, satisfy the recurrence below at orders 6-13. The recurrence also
  reproduces HYP-9458's published order-16 row (a retrodiction: that row was read before the recurrence was found).

## Pre-registered predictions (confirmatory run: `uh_noncritical_run_20261104 18 14 14`)

- **(P1) Garey sanity check.** Class counts 17,711, 46,368, 121,393, 317,811, 832,040 at orders 14-18.
- **(P2) Local roundness (A).** `round_bad = 0` at every order.
- **(P3) One run (B).** `run_bad = 0` and `crit_bad = 0` at every order.
- **(P4) Cut arcs (C).** `cut_bad = 0` at every order, and the cut-arc rows at orders 14-16 equal HYP-9458's published
  T(n, k) rows (an independent recomputation).
- **(P5) Refined counts (D).** (#sole exits, #sole entries) = (1,1), (1,2), (2,1), (2,2):
  - n = 14: 2,585, 4,180, 4,180, 6,766;
  - n = 15: 6,766, 10,945, 10,945, 17,712;
  - n = 16: 17,712, 28,656, 28,656, 46,369;
  - n = 17: 46,369, 75,024, 75,024, 121,394;
  - n = 18: 121,394, 196,417, 196,417, 317,812.
- **(P6) Asymmetry (E).** `asym_bad = 0`.
- **(P7) The run-length recurrence (CONJECTURE).** For n >= 6,

      N(n, r) = N(n - 1, r) + 2 N(n - 1, r - 1) - N(n - 2, r - 1),   N(4, .) = (1, 0),  N(5, .) = (1, 1, 1),

  i.e. sum over n, r of N(n, r) t^r x^n = t^2 x^4 (1 + (t^2 - t) x) / (1 - (1 + 2t) x + t x^2) (at t = 1 this is
  Garey's x^4 / (1 - 3x + x^2)). Predicted rows r = 2..n - 1:
  - n = 14: 1, 10, 55, 210, 606, 1372, 2471, 3530, 3920, 3232, 1792, 512;
  - n = 15: 1, 11, 66, 275, 870, 2178, 4396, 7183, 9420, 9680, 7424, 3840, 1024;
  - n = 16: 1, 12, 78, 352, 1210, 3312, 7380, 13504, 20256, 24600, 23552, 16896, 8192, 2048;
  - n = 17: 1, 13, 91, 442, 1639, 4862, 11826, 23868, 40081, 55692, 63072, 56576, 38144, 17408, 4096;
  - n = 18: 1, 14, 105, 546, 2171, 6930, 18238, 40140, 74313, 115598, 149856, 159168, 134400, 85504, 36864, 8192.

  Orders 14-16 are not blind: HYP-9458 published T(n, k) there (its order-16 row was read). Orders 17 and 18 are blind.
  Through (C) the recurrence implies HYP-9458's fitted columns T(n, 2), T(n, 3), T(n, 4) and its top diagonals.
- **(P8) Sanity.** `no_one_arc_cut = 0`; at order 14, `hc_xcheck_bad = 0` and `bruteforce_cut_bad = 0`.

Every outcome is reported. A failure of P2-P6 would refute a proof step of THM-4868 and is reported with the class.
If order 18 does not finish, P1-P8 are reported at orders 14-17 and the order-18 items are left open.

## Results (confirmatory run `uh_noncritical_run_20261104 18 14 14`, started after the pre-registration commit a8db6218f9)

Output: `04-computation/experiments/uh_noncritical_run_20261104_confirm.out`.

- **(P1) CONFIRMED.** Class counts 17,711, 46,368, 121,393, 317,811 and 832,040 at orders 14-18.
- **(P2) CONFIRMED.** `round_bad = 0` at every order.
- **(P3) CONFIRMED.** `run_bad = 0` and `crit_bad = 0` at every order.
- **(P4) CONFIRMED.** `cut_bad = 0` at every order. The cut-arc rows at orders 14-16 equal HYP-9458's published T(n, k)
  rows; for example, order 16 reads 2048, 8192, 16896, 23552, 24600, 20256, 13504, 7380, 3312, 1210, 352, 78, 12, 1.
- **(P5) CONFIRMED.** The refined counts equal the predicted Fibonacci values at every order 14-18.
- **(P6) CONFIRMED.** `asym_bad = 0`.
- **(P7) CONFIRMED.** The rows N(n, r) equal the recurrence rows at orders 14-16 (non-blind) and at orders 17 and 18
  (blind).
- **(P8) CONFIRMED.** `no_one_arc_cut = 0` at every order. At order 14, `hc_xcheck_bad = 0` and
  `bruteforce_cut_bad = 0`.

**After the pre-registration: P7 is now PROVED.** During the run the recurrence was proved as THM-4868 (H). The proof
uses the admissible-parent structure THM-4868 (G):
- the vertices p with e_p(T) = 1 are y, y- and (when y- is a sole entry) y--;
- their children have run lengths r, r + 1 and r + 2.

(G) and (H) were first observed on all parents of order <= 12, then proved, and are checked on every parent of order
<= 14 by `04-computation/experiments/shelor_parents_20261104.c`. Through THM-4868 (D), every fitted formula of
opus-S74's HYP-9458 follows.

The blind order-17 confirmation preceded the proof. The order-18 confirmation ran in parallel with it and checks the
proof.

**Status.** Every item is CONFIRMED or PROVED. Nothing open remains in this hypothesis.
