---
source: kind-pasteur-2026-07-07-S63
status: REFUTATIONS-WITH-MECHANISM (two heresies died informatively) + one decisive-either-way
  LP experiment + a catalog of untested recombinations. Owner creativity directive.
tags:
  - lonely-runner
  - LRC14
  - creativity
  - schur
  - two-point
  - relation-lattice
  - catalog
---

# Three heresies and a catalog

**kind-pasteur-2026-07-07-S63 (HYP-4897).** Owner: *come up with even more creative possible
connections; recombine and break apart concepts; challenge all assumptions; explore.* Three
assumptions were put on trial with computations; two died fast and informatively, one opened
a genuinely new question. A catalog of further untested recombinations closes the note.

## Heresy A — "(A′) is secretly a convexity statement" — REFUTED, spectacularly

If `μ_{1/7}` were monotone under position-fixed Robin-Hood step transfers (move one unit
from a bigger step to a smaller one), iterating would reach the AP and (A′) would reduce to
a *local two-step exchange lemma*. Exhaustive exact test (order-cell engine) at diameters
13/14/15: **violation rates 0% / 77% / 60%**. Transfers toward equality can *raise* μ by
+0.34 (e.g. `(1¹¹,4) → (1⁶,2,1⁴,3)`: 0.560 → 0.903): moving mass into a mid-word position
creates a resonance-breaking irregularity whose μ-cost dwarfs the equalization gain — the
Wiener-weight physics of S62 (`l(14−l)`: middle positions are ~3.8× as expensive) made
quantitative. **The (A′) landscape is resonance-rugged, not convex; no rearrangement
inequality in step space will prove it.** This also explains, mechanically, why every naive
descent in the fleet's history kept getting trapped: the landscape is full of one-unit
uphill walls.

Bonus finding: at D=15 the exact minimizer orbit is `(1¹¹,4) ↔ (4,1¹¹)` — a **non-palindromic
mirror pair**. My S62 palindromic-extremizer conjecture is FALSE as an every-diameter
statement (no palindromic word attains the D=15 minimum) and survives only in the refined
form mac-mini's THM-639 independently found the same day: the deep (record-board) minimizers
are palindromic, but **symmetry-breaking mirror pairs occur** — the exact SC/NS coexistence
of the tournament metagraph, transplanted.

## Heresy B — "the sparse lane is governed by the shortest relation" — REFUTED with mechanism

The Poisson-summation frame says the deficit `iid − E[maxgap]` is a sum over the zero-sum
relation lattice. The tempting compression — "deficit controlled by the lattice **minimum**"
— fails immediately: essentially *every* 13-family (including all 8 random ones tested) has
an L1-norm-4 zero-sum relation (`e_a + e_b = e_c + e_d` coincidences are generic), yet the
deficit ranges over a factor of 100 within `minrel = 4` (records +0.048, random ±0.003).
Meanwhile the only families with *no* weight-≤3-coefficient relation in the search box
(Sidon-greedy, geometric) sit at deficit ≤ 0.006 in absolute value. **The law is the
counting/theta function of the lattice, not its minimum** — one short relation is noise,
what lowers the mean is *many, coherently aligned* short relations (the additive-energy
picture, now sharpened: it is the theta series of `L(E) ∩ {Σk=0}` weighted by the
functional's Fourier decay that any rigorous sparse-lane bound must sum).

## Heresy C — "is the mean floor a 2-POINT theorem?" — the new question (LP experiment)

Every integer family has **exactly uniform** unlabeled pair-distance statistics (S59). So:
does `E[maxgap] > 1/7` already follow from pair-uniformity alone, over *all* cyclic 13-point
processes — integer or not? If yes, a Cohn–Elkies-style 2-point LP proof of the entire mean
floor exists (and the fleet's arithmetic struggles were optional); if no, then **any proof
must consume weight-≥3 information** — a methodological theorem localizing the difficulty.

A first naive mixture optimization was inconclusive (constraint not enforced — reported
honestly). The decisive small version: the exact LP over all `C(25,12) = 5 200 300` 13-point
configurations on `ℤ/26` (deduped to **162 770** distinct `(maxgap, spectrum)` columns, 13
constraints). Note the logic: grid mixtures are valid continuum processes, so a *feasible
grid primal below 1/7 kills* the 2-point conjecture outright, while the grid dual only
bounds the grid version. **Result: the bracket straddles the threshold.** Converged dual
lower bound `0.126406`; best (infeasible, Linf 2.74/12) primal `0.14355`. The dual's
convergence plateau suggests the true grid-LP value sits near `0.127–0.135` — *leaning
below `1/7`*, i.e. leaning "pair data alone cannot force the mean floor, the weight-≥3
usage is necessary" — but without a feasible primal certificate the kill is not claimed.
My pre-experiment knife-edge estimate (`≈ 0.155`) appears to have over-credited the peak-
dilution constraint. **Handoff: the 162770×13 LP is trivial for any real solver
(HiGHS/scipy linprog) — one run decides the grid question exactly.** Either way it explains
a fleet-wide pattern: the 2-point-flavored bounds (CE ≈ 0.25, PZ ≈ 0.26, both above bar but
far from the truth 0.44) saturate near where the 2-point LP lives; the comfortable margin of
real integer families is a weight-≥3 phenomenon.

## The catalog (untested, ranked by how much I'd bet on them)

1. **CRT bi-polynomial vs composite-14** (the frontier heresy). The SOTA polynomial method
   needs `k+1` prime; `14 = 2·7` is THE obstruction (klein-S151). Challenge the assumption
   that the method's home is a prime field: `ℤ/14 ≅ ℤ/2 × ℤ/7`, and several fleet artifacts
   already CRT-split (saturation over `q ≤ 14`; monad's shift-sum prunes onto
   `Σm ≡ 0 mod 14` = a mod-2 AND a mod-7 condition). Concrete first question: does the
   S–T counting object factor as a convolution over the prime factorization of `k+1`? One
   session with the actual paper. Highest potential payoff in the whole catalog.
2. **The theta-series sparse-lane theorem** (Heresy B's constructive half): smooth the
   avoidance functional (choose the mollifier, control `|F̂|` decay explicitly), Poisson-sum
   over `L(E) ∩ {Σk = 0}`, and get: *families whose weight-3..W relation counts are below
   explicit thresholds have `μ ≥ iid − ε`*. The step-alphabet filtration (S62) then pinches
   from both ends: structured ⟹ grids/ledgers; unstructured ⟹ theta bound.
3. **The runner braid.** The 14 strands over one period form a *positive braid* whose
   pair-crossing numbers are `|vᵢ − vⱼ|` (closed form: the closure is a torus-link family);
   loneliness = a meridian disk avoiding all strands at some angle. Positive-braid geometry
   (Garside normal form, fractional Dehn twist ≈ `1/13`-ish) as a new invariant language
   for the movie/order-cell frame (opus-S136). Unknown payoff, genuinely fresh language.
4. **LP duality on the speed lattice.** LRC(14) ⟺ the Farey-indexed boxes
   `B_t = {v : ‖vᵢt‖ ≥ 1/14}` cover the nonzero speed lattice. The dual object — an optimal
   measure λ on times — is what every mean bound implicitly chooses (`λ = Lebesgue`). The
   ruler forces Lebesgue in Part A's counting, but the FLOOR half (`hlarge`) is free to
   reweight: `μ_λ` with λ concentrated near denominators ≤ 6 makes the roof windows heavy.
   Worth one session: is there a λ making the k=8..10 legs' G_P-intersected floors easy,
   while staying compatible with a λ-weighted witness count?
5. **The consensus coin-flip** (conceptual, free). The majority tournament of the runners is
   *exactly* fair on every pair — LRC families are the maximally intersubjective objects:
   zero pairwise signal, all structure collective. The owner's phrase lands precisely: the
   tournament project studies complete pairwise assignments; LRC studies systems whose
   pairwise data is *forced featureless* — two poles of one relational ontology, bridged by
   the order-cell walk (a loop in the permutohedron whose vertices are the transitive
   tournaments).
6. **p-adic/adelic framing** (bookkeeping only): loneliness = the archimedean place;
   covering/sieve certificates = finite places; the saturated reduction is a local-global
   statement. Might clean up write-ups; no new theorems expected.

## Ledger

- REFUTED: exchange-monotonicity of μ (59–77% violations, mechanism = Wiener-weight
  resonance cost); palindromic-extremizer at every diameter (D=15 mirror pair — but
  record-board palindromy + SC/NS refinement stands, converging with THM-639); shortest-
  relation compression of the deficit (counting function, not minimum).
- OPENED: the 2-point LP question (Heresy C) with the exact M=26 experiment; the catalog.
- Files: `lrc_three_heresies_kps_S63.py`, `lrc_twopoint_lp_M26_kps_S63.py` (+outs).
- Builds on: S59 (pair uniformity, deficit frame), S62 (step gauge, Wiener weights),
  opus-S136 (exact engine, order cells), mac-mini THM-639 (mirror pairs), klein-S157
  (voltage lifts).
- Does NOT prove LRC(14). Prunes three tempting routes and aims two new ones.
