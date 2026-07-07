---
source: klein-2026-07-07-S158 (HYP-4821)
status: COMPOSITE. The Bukszár–Prékopa CHERRY-TREE floor supersedes Hunter at the k=8
  criticality — EXACTLY TIGHT at the AP (0.3347 = true W), adversarial min 0.1969 at the
  conservative R-route bar; needs only 5 adaptively-chosen triple masses beyond THM-638.
  The smallest-frequency conditioning route to a generic triple bound is REFUTED (C ≈ 6.9,
  non-decaying for comparable frequencies). Mixed lifts do NOT erode the mean record
  (0/30597). Hard-core dossier: the two ledger-hard P-shapes carry true ρ* at 6.9×/8.3×
  the bar with R ≥ 0.93.
tags:
  - lonely-runner
  - LRC14
  - k8-leg
  - cherry-tree
  - bonferroni
  - voltage-lift
  - hard-core
---

# The cherry-tree floor is tight at the AP

**klein-2026-07-07-S158.** Owner worklist session (proofs before formalization). My lanes:
the triple-mass upper bound, mixed-lift erosion, the hard-core dossier, valuation trees.

## 1. The cherry-tree floor (the headline)

Hunter's spanning tree (7 events, 6 edges) gave the unconditional-but-weak `6/49 ≈ 0.122`.
The **t-cherry tree** (Bukszár–Prékopa: start with an edge, attach each new vertex to an
existing edge; 2n−3 pairs, n−2 triples; upper-bounds `P(∪)` by `S1 − Σ_pairs + Σ_triples`,
provable by the same leaf induction as Hunter) gives at the k=8 criticality (`S1 = 1`):

> `W_end ≥ Σ_{11 cherry pairs} m₂ − Σ_{5 cherry triples} m₃`, cherry structure OUR choice.

Greedy-optimized over cherry trees (numeric masses, two grids):

| shape | cherry floor | Hunter | true W |
|---|---|---|---|
| AP {1..8} | **0.3347 = true W exactly** | 0.1224 | 0.3347 |
| spread / geometric / two-cluster / mild | 0.218–0.276 | 0.1224 | 0.32–0.35 |
| adversarial min (jump moves) | **0.1969** at {3,4,5,7,8,9,14,22} | — | — |

**Tightness at the AP** is the correct-extremal-formulation signature (cf. opus-S134's F₆
anchors): the optimal cherry tree reproduces the AP's inclusion–exclusion exactly. The
adversarial min sits at the conservative R-route bar `0.197` (which assumed `R ≥ 0.75`;
measured `R ≥ 0.913` at k=8 lowers the true requirement to `≈ 0.162`, cleared with margin).
**Remaining lemma stack for the k=8 leg at the R-route bar:** Bukszár–Prékopa (classical,
Hunter-style induction), THM-638 (PROVED, the 11 pair masses), upper bounds on the **5
chosen** triple masses (adaptive choice = pick triples where the masses are provably small
— the entire remaining freedom), and `R ≥ 0.75`. The 7-adic G-bonus (S157) applies to the
cherry pairs identically.

## 2. The generic triple-conditioning route is DEAD (documented)

The S158 plan `m₁₂₃ ≤ θ·m₂₃ + C·θ·q₁/q₃` (condition on the smallest frequency) FAILS as an
assembly: measured sharp `C ≈ 6.9` (worst at comparable small frequencies), and the error
term does not decay when differences are comparable — summed over 35 triples the chain
gives floors of −12 to −18 (vs true Bonf3 ≈ +0.22). **Why:** the per-interval discrepancy
of `A₂∩A₃` is only small when `q₁ ≪ q₃`; spread shapes have many comparable differences.
The cherry route dodges this exactly because its 5 triples are CHOSEN — pick the ones with
separated scales (where the conditioning error is provably small), never the clustered ones.

## 3. Mixed lifts do not erode the record (0 / 30,597)

Extending S157: mixed lifts `c₁{1..a₁} ∪ c₂{1..a₂} ∪ B` (c ≤ 5) and 2-lifted perforated-AP
bases — **no candidate dips below the record** `12907/65520 = 0.196993` (numeric screen,
exact leaders re-confirmed). The mean sidecar's `+0.0057` margin now stands over pure
(S157: 8,726) + mixed + perforated (S158: 30,597) lift classes. Remaining erosion surface:
c ≥ 6, exotic bases, free search beyond the lift ansatz.

## 4. Hard-core dossier (for opus's avoidance-kernel / mac-mini's PZ-on-U)

monad-S3's two ledger-hard small parts, with adversarial-min TRUE objects (k=8/k=9):

| hard core | meas(G_P) | adv. min ρ* | at E | μ(E) | R | ρ*/m_P |
|---|---|---|---|---|---|---|
| P₈ = {9..13} (k=8) | 0.4467 | 0.3906 | AP {1..8} | 0.9402 | 0.930 | **6.9×** |
| P₉ = {10..13} (k=9) | 0.5251 | 0.4689 | {15..24}∖{16} | 0.9402 | 0.950 | **8.3×** |

The hard cores are hard **for the ILedger mechanism, not for the truth** — the true
intersected object never approaches the bar there, and quasi-independence is strongest
(R ≥ 0.93). Mean θ-degree: exactly 1 at k=8 (critical), 8/7 at k=9. These are the concrete
target numbers for the avoidance-kernel and PZ-on-U machines, and the first tournament-
bridge test cases (both cores sit at/above the criticality where the phase-digraph frame
is sharpest).

## Honest status

Cherry floor: verified-numeric (greedy over cherry trees, two grids, jump adversaries);
the inequality itself classical-provable; pairs proved (THM-638); the 5 triple bounds and
R ≥ 0.75 remain the open lemmas. Triple-conditioning: refuted as assembled (do not re-walk).
Mixed-lift: numeric screen + exact leaders. Nothing here closes the k=8 leg or LRC(14).

## Files
`lrc14_triple_upper_and_dossier_klein_S158.py`, `lrc14_cherrytree_floor_klein_S158.py`
(+ .outs). Pointers: THM-638, HYP-4791/4801/4811 (the k=8 program), monad-S3 (hard cores),
mac-mini-S44 (concurrent: hard-cores × PZ-on-U — complementary to this ρ*/R dossier),
opus-S136 (exact (A′) exhaustive + isolation gaps), THM-639 (runner tiling model).
