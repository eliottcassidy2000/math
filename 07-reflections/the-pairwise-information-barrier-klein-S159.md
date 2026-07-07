---
source: klein-2026-07-07-S159 (HYP-4831)
status: BARRIER + delineation. The cherry tree deepened to its parent object — the
  partially-specified-moment LP — which is EXACTLY TIGHT at the AP (pairs alone determine
  W there) and exposes a PAIRWISE-INFORMATION BARRIER at spread shapes: no certificate
  using only THM-638 pair data can give a uniform k=8 floor above ~0.1233 (exhibit shape;
  Hunter's tree is near-dual-optimal there). Triples are LP-certified NECESSARY for the
  0.197 R-route and for the hard cores (conditional pairwise LP dips to 0.010-0.015 < m_P).
  Pairwise LP stays positive through k=12 at the AP where Hunter's base is negative.
  FKG hair: pairwise positive association holds (G >= 0) but the product bound fails at
  the AP by exactly -0.0052.
tags:
  - lonely-runner
  - LRC14
  - k8-leg
  - moment-LP
  - cherry-tree
  - barrier
  - hard-core
---

# The pairwise-information barrier

**klein-2026-07-07-S159.** Owner: explore cherry-tree concepts deeply, extend, work at the
hard core. The deepening: Hunter (spanning tree) and Bukszár–Prékopa (cherry) are
dual-feasible points of the **partially-specified-moment LP** — minimize the empty-atom
mass over all event-algebra measures consistent with the exact singles (`p_i = θ`, proved)
and the exact pairs (`m_ij`, THM-638, proved). Its value is the SHARP floor obtainable from
pairwise information; its dual optimal is the best possible "generalized cherry"
certificate. 128 atoms at k=8 — solved exactly per shape.

## 1. Two faces of the LP at k=8

- **At the AP the LP is exactly tight: floor = true W = 0.3347.** The AP's pair data
  (nested differences, maximal G-corrections) pins the entire joint distribution — the AP
  is *pairwise-determined*. Same tightness signature as the cherry (S158) and F₆ (opus),
  now at the information-theoretic level: nothing beyond pairs is needed AT the AP.
- **At spread shapes the LP collapses to Hunter: adversarial min 0.1233** at
  `{1,3,5,6,7,9,25,38}` (near-iid pair masses), a hair above `6/49 = 0.1224`. Since the LP
  is the SHARP pairwise bound, this shape is an **exhibit proving a barrier**: *no method
  consuming only pair masses can certify a uniform k=8 floor above ≈ 0.1233*. Hunter's
  spanning tree is near-dual-optimal in the pairwise regime — the S155 floor was already
  essentially the best pairs can do.

**Consequence (the honest floor hierarchy at k=8):**
`Hunter 6/49 (PROVED) ≤ LP ≈ 0.1233 (sharp pairwise BARRIER) ≪ cherry 0.197 (needs 5
chosen triple uppers) ≪ PZ-on-B 0.712 (needs two-moment bounds) ≪ truth ≈ 0.94.`
The 5-triple lemma (S158 handoff (a)) is not one route among several — it is **the gate**:
LP-certified necessary for anything beyond 0.123.

## 2. k = 9..12: disaggregated pairs rescue positivity where Hunter dies

Endpoint LP floors at the AP: **0.2776 (k=9), 0.1908 (k=10), 0.1266 (k=11), 0.0097 (k=12)**
— all positive despite Bonferroni bases of −1/7 … −4/7 (MISTAKE-122's bare zero at k=9 is
information-theoretically beatable: the disaggregated pair VALUES, not just their sum,
carry the floor). Adversarial k=9 min: 0.0026 — positive but razor (near-iid shapes;
solver-noise caveat, needs exact confirm). Pattern: pairs buy less as the supercriticality
(k−1)/7 grows.

## 3. The hard cores defeat pairwise information (the owner's "work at the hard core")

Conditional LP (atoms restricted to `G_P`, all singles/pairs measured within `G_P`) on
monad-S3's ledger-hard shapes: AP values 0.153/0.136, but **adversarial minima 0.0150
(P₈={9..13}) and 0.0104 (P₉={10..13}) — BELOW m_P = 0.0565.** The true ρ* at those shapes
remains ≈ 6× the bar (S158 dossier) — so this is not a threat to the leg but a proof that
**pairwise data is insufficient on the hard cores**: the G_P interaction hides mass that
only triples (or the R-route) can pin. My pre-run coupling heuristic
(`conditional ≈ meas(G_P)·unconditional + O(1/d)`) gives at best `0.4467·0.1233 ≈ 0.055`
— razor at m_P even before the adversary; the measured 0.015 kills it. **Coupling the
engine to the intersected ledger must carry triple information.** (Conditional cherry with
the 5-triple lemma, or THM-579-style R, are the two live carriers.)

## 4. The FKG hair (guardrail, pinned)

THM-638 gives all 21 same-sign pair correlations ≥ 0 — yet the 7-fold product bound FAILS
at the AP: `W(AP) = 0.334695 < (6/7)⁷ = 0.339917`, gap **−0.005222**. Pairwise positive
association does not extend to full association for these arithmetic events; any
"independence-domination" shortcut is dead by exactly half a percent. The LP/cherry layer
is genuinely necessary — now with the failure quantified.

## Honest status

The barrier (§1) is rigorous as an upper bound on pairwise-certificate power (exhibit
shape + LP duality + THM-638 values). The uniform lower "LP ≥ 0.1233 over all shapes" and
all adversarial minima are empirical (jump-move adversaries, two grids on measured data;
k=9's 0.0026 and the conditional minima deserve exact confirmation). Nothing closes any
leg. NEXT: (a) the 5-triple upper bounds for separated-scale cherries — now THE gate;
(b) extract the LP dual certificates at the adversarial minima (do the optimal duals have
a closed form = a "weighted cherry" theorem?); (c) exact-rational conditional pair/triple
masses on the hard cores (the intervals are all rational — mechanical); (d) R ≥ 0.75.

## Files
`04-computation/lrc14_pairwise_moment_lp_klein_S159.py` (+ .out). Pointers: THM-638,
HYP-4791/4801/4811/4821 (the k=8 program), monad-S3 (hard cores), mac-mini-S44 (palindrome
verdict + PZ side), Prékopa/Boros moment-LP lineage, MISTAKE-122 (whose bare-0 the k=9 LP
beats).
