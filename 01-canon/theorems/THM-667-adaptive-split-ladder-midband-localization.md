---
id: THM-667
title: The adaptive-split ladder and the mid-band localization — re-splitting S = P̃ ∪ L̃ at the 9/14-threshold makes every covering 13-set's residual cluster WIDE (r̃ ≥ 14/5 = 2.8, the aliasing-window edge); the measure floor then closes at EVERY stratum (k̃ ≤ 7 by the LRC(≤13) LIPSCHITZ-FATTENING lemma — the apex-7/Fraenkel wall CANNOT occur for P̃; k̃ ≥ 8 by Bonferroni |P̃| = 13−k̃ conservation + the shape-universal THM-661/663 floors, margins +0.047…+0.309); the entire remaining content of LRC(14) localizes to REALIZATION FOR MID-BAND SPEEDS (Vmax/14, 9·Vmax/14) — speeds that can neither ride the snap window nor join the confined cluster; mid-band-free covering sets are closed end-to-end
status: PROVED (Lemma A fattening — 3 lines from the settled LRC(≤13) citation; Lemma B ladder arithmetic — exact rationals below; the zone-edge derivations) + VERIFIED (exact meas(G_P̃) on |P̃| ≤ 12 instances, always ≫ the Lipschitz floor; end-to-end composite on mid-band-free covering instances with exact ground truth M(S); census of the instrument firing map, HYP-5691). The REALIZATION on the mid-band remains OPEN — this theorem does not close LRC(14); it localizes what is left. Concurrent & complementary: boxeph LEM-014 (P-separated composed realization, wide regime) is the realization core this ladder FEEDS; both hit the same wall (mid-band members cannot be ε-eroded), which is the consistency check.
source: klein-2026-07-09-S208 (HYP-5691, HYP-5692)
depends_on:
  - THM-527   # the P∪L split, G_P, the lonely-density reformulation
  - THM-530   # k≤7 pigeonhole (μ_{1/7}=1 a.e.); the m_P floor
  - THM-661   # shape-universal moment floors μ_{1/7}(E) ≥ floor_k (min over ALL clusters)
  - THM-663   # the covering-case assembly; floor table k=8..13
  - LRC(≤13)  # SETTLED (owner directive 2026-07-02) — used as a citation, Lemma A
related:
  - LEM-014 (boxeph)  # the composed realization this ladder feeds (wide regime, P≠∅)
  - THM-665 (monad-explorer)  # aliasing window V₀ ≈ 2.8·spread — the 9/14 edge comes from it
  - THM-668 (mac-mini)        # pair-sum ruler — an independent realization axis for the mid-band
  - HYP-5691, HYP-5692, HYP-5707, HYP-5708, HYP-5710
  - OPEN-Q-108   # the apex-7 wall — Lemma A removes it for the P̃-positivity question
---

# THM-667 — the adaptive-split ladder and the mid-band localization

**Setting.** S a primitive covering 13-set (the open case of LRC(14) after THM-369/523),
Vmax = max S assumed > 13 (else S ⊆ [1,13] is the settled LRC(≤13) shape). THM-527 splits
S = P ∪ L at the fixed threshold 13. This theorem re-splits ADAPTIVELY at the ratio
threshold (9/14)·Vmax:

> P̃ = {v ∈ S : v < (9/14)·Vmax},  L̃ = {v ∈ S : v ≥ (9/14)·Vmax},  k̃ = |L̃| ≥ 1,
> Ẽ = {Vmax − u : u ∈ L̃} ⊆ [0, (5/14)·Vmax],  spread(Ẽ) ≤ (5/14)·Vmax.

So the residual cluster ALWAYS has split-ratio r̃ = Vmax/spread(Ẽ) ≥ 14/5 = 2.8 — the
aliasing-window edge (THM-665: V₀ ≈ 2.7–2.8·spread): **after the ladder, every covering
13-set's cluster is in the wide regime.** The cost is that P̃ may contain up to 12 members,
including mid-band speeds; the theorem is an exact account of that cost.

## Lemma A (LRC(≤13) Lipschitz fattening — the apex-7 wall cannot occur for P̃)

> For ANY set T of ≤ 12 distinct positive integer speeds,
> meas(G_T) := meas{x ∈ [0,1) : ‖t·x‖ ≥ 1/14 ∀t ∈ T} ≥ 1/(91·max T) > 0.

*Proof.* T ∪ {observer} is an LRC instance with ≤ 13 runners, so by LRC(≤13) (SETTLED)
there is τ₀ with ‖t·τ₀‖ ≥ 1/13 for all t ∈ T. Each x ↦ ‖t·x‖ is t-Lipschitz, so for
|x − τ₀| ≤ (1/13 − 1/14)/max T = 1/(182·max T) we keep ‖t·x‖ ≥ 1/14. That ball has
measure 2/(182·max T) = 1/(91·max T). ∎

*Consequence for OPEN-Q-108's apex-7 wall:* exact tiling of [0,1) by the 1/7-danger sets of
any ≤ 12 speeds would force meas(G_T) = 0, contradicting Lemma A. The tiling/Fraenkel
pathology is thereby EXCLUDED for the P̃-positivity question — at every |P̃| ≤ 12, uniformly,
with an explicit (if tiny) floor. (VERIFIED exact, 5 shapes incl. |P̃| = 12: actual
meas(G_P̃) = 0.16–0.28, three orders above the guarantee.)

## Lemma B (the ladder — the measure floor closes at every stratum)

|P̃| = 13 − k̃ exactly (the ladder conserves the total). Then ρ̃* := meas{x ∈ G_P̃ :
maxgap{frac(e·x) : e ∈ Ẽ} > 1/7} satisfies:

- **k̃ ≤ 7:** ≤ 7 circular gaps summing to 1 force maxgap ≥ 1/7 with equality only on the
  exact 1/7-grid (measure-zero in x), so μ_{1/7}(Ẽ) = 1 a.e. and
  ρ̃* ≥ meas(G_P̃) ≥ 1/(91·max P̃) > 0 (Lemma A).
- **k̃ ≥ 8:** |P̃| = 13 − k̃ ≤ 5, so Bonferroni gives meas(G_P̃) ≥ 1 − |P̃|/7, and the
  THM-661/663 floors (SHAPE-UNIVERSAL: the moment-LP bounds quantify over all integer
  clusters, no covering condition) give μ_{1/7}(Ẽ) ≥ floor_k̃. Union bound:
  ρ̃* ≥ (1 − |P̃|/7) + floor_k̃ − 1. Exact ledger:

  | k̃ | |P̃| ≤ | Bonferroni | need > 1 − floor_k̃ | margin |
  |---|---|---|---|---|
  | 8 | 5 | 2/7 = .2857 | .2389 | **+.0468** |
  | 9 | 4 | 3/7 = .4286 | .3554 | **+.0732** |
  | 10 | 3 | 4/7 = .5714 | .4469 | **+.1245** |
  | 11 | 2 | 5/7 = .7143 | .5952 | **+.1190** |
  | 12 | 1 | 6/7 = .8571 | .6441 | **+.2130** |
  | 13 | 0 | 1 | .6912 | **+.3088** |

**So the measure floor of the ladder-split system is positive for EVERY covering 13-set —
the floor is NEVER the obstruction, at any stratum.** (This supersedes the worry that the
adaptive split trades cluster members for an apex-7-walled P̃.)

## The zone decomposition and the localization

At the realization step (turn ρ̃* > 0 into a real lonely τ; the snap/embed places the fast
phase φ = frac(Vmax·τ) at the gap midpoint, moving τ by ≤ 1/Vmax — klein-S205 LRCDriftEmbed
/ boxeph LEM-014), the members of S behave by SPEED ZONE:

- **P = [1,13]:** phase move ≤ 13/Vmax over the full snap window — rides with ε-erosion
  cost O(1/Vmax) (LEM-014's ε = 20/Vmax). CLOSED.
- **T = (13, Vmax/14]:** phase move ≤ 6/98 = 3/49 < 1/14 over the realized sub-interval
  (width ≤ (6/7)/Vmax) — rides with slack; Bonferroni budget: each such member costs
  ≤ 13/49 of measure. CLOSED for the budgeted counts (exact ledger in the S208 script;
  |P| = 5 admits one rider, |P| ≤ 3 admits two).
- **MID-BAND (Vmax/14, 9·Vmax/14):** phase move over the snap window is between 1/14 and
  9/14 — too fast to ride (slack ≥ 1/14 unavailable on top of the 1/14 threshold), too slow
  to join the confined cluster (co-offset > (5/14)Vmax breaks r̃ ≥ 2.8). **No current
  instrument covers these members at the realization step.** Both edges are forced: the
  lower by the ride-slack arithmetic, the upper by THM-665's window edge.
- **CLUSTER [9·Vmax/14, Vmax]:** wide regime; aliasing existence (THM-665) + drift embed
  (klein-S205) / composed realization (LEM-014). CLOSED (modulo LEM-014's (H1) interface
  bookkeeping).

> **LOCALIZATION.** The open content of LRC(14) is exactly: REALIZATION FOR COVERING SETS
> WITH A SPEED IN (Vmax/14, 9·Vmax/14). Covering sets with empty mid-band (and budgeted
> T-zone) are closed end-to-end by [Lemma B floor] + [aliasing existence] + [drift-embed/
> LEM-014 realization] — VERIFIED exactly on constructed k̃ = 8, |P̃| = 5 instances at
> Vmax = 280, 420 (measure .44–.45, aliasing margins ≈ 1.7×, embed fires with gaps
> .25–.32, ground truth M(S) = 10/47, 73/349).

## What this changes

1. THE APEX-7 WALL IS GONE from the P̃-positivity question (Lemma A) — the wall that
   OPEN-Q-108 treats as the union-bound death is real for *danger-measure* bookkeeping but
   cannot produce meas(G_P̃) = 0: LRC(≤13) forbids it.
2. The three hrefl fronts (HYP-5707 grid / HYP-5710 pure-cluster / HYP-5708 = LEM-014
   P≠∅-wide) get a FOURTH: the proportional regime, where the ladder + LEM-014 compose and
   the only survivor is the mid-band. All four fronts are the same node (HYP-5692): the
   mid-band member's safety at the snap τ IS the Kronecker/equidistribution residue.
3. The residual family is now a RATIO-WINDOW statement: speeds in (V/14, 9V/14). A
   candidate independent axis for it: THM-668's pair-sum rulers (the witness lives at
   p/(v_i+v_j) — mid-band pairs give pair-sums in (V·15/14, V·9/7), a NEW ruler family
   nobody has attacked; see HYP-5720).

## Files

- `04-computation/lrc14_split_ratio_census_klein_S208.py` (+ `.out`) — HYP-5691 census:
  strata population, instrument firing map (aliasing 20/20 & embed 20/20 on confined;
  0/13 & 1/34 on proportional), ground truth margins 0.18–0.23 on proportional instances.
- `04-computation/lrc14_adaptive_ladder_klein_S208.py` (+ `.out`) — Lemma A exact
  verification, Lemma B ledger, end-to-end composites, mid-band population statistics
  (random k ≥ 8 covering sets: ~95% have ≥ 2 mid-band members — the residual family is
  the bulk by volume; the value is the localization, not a volume reduction).
