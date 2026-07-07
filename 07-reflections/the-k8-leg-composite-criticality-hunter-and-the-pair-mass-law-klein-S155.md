---
source: klein-2026-07-07-S155 (HYP-4791)
status: COMPOSITE ADVANCE on the binding k=8 witness-floor leg. PROVED-shape pieces (the
  pair-mass law — table-exact, proof an elementary lattice count; the Hunter-endpoint
  diameter-free floor 6/49 at the k=8 criticality; the union-route D_k* table), an
  exhaustive verified band (diam 12..26, margin +0.30), explicit-constant far-case
  skeleton, and the PZ-on-B two-moment reduction (adversarial min 0.712 >= bar 0.675).
  Nothing closes the leg; it now mirrors the k=13 composite's epistemic state.
tags:
  - lonely-runner
  - LRC14
  - route-1
  - density-floor
  - k8-leg
  - hunter
  - criticality
  - tournament-analysis
  - pair-correlation
---

# The k=8 leg: criticality, Hunter, and the pair-mass law

**klein-2026-07-07-S155.** Owner: attack target 1 (the k=8 window-avoidance inequality) at
the highest-leverage angles, integrating incoming work, connecting to the tournament side.
Target 1 (from my S154 lead): the binding THM-530 leg `μ_{1/7}(E) ≥ need_8 = thr_8 + m_P =
1702763/2522520 ≈ 0.67502` for every 8-element integer set E — equivalently (bisection
identity, k=7 base proved) `Bis ≤ 0.325`.

Concurrent landings integrated mid-session: **opus-S135 proved the Farey roof outright**
(every "modulo roof" caveat in kps-S59/monad-S2/klein-S154 is discharged — my S154 LOO≤36
composite is now unconditional); **kps-S60's intersection ledger** extended the diameter
floor to the G_P-intersected legs (k=8 bite: diam ≤ 11), superseding my union-route
crossings mid-run — my independently computed union bites (9/11/11 at k=8/9/10) exactly
match their MISTAKE-121-corrected values, a clean cross-validation.

## 1. The union-route D_k* table (exact; cross-check + waste quantification)

From the roof, `D_k* = max{D : μ_{1/7}(AP_{D+1}) ≥ need_k}`:

| k | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|----|----|----|----|
| need_k | .67502 | .56223 | .45209 | .33121 | .19934 | .05649 |
| union D_k* | 9 | 11 | 11 | 15 | 23 | 75 |
| kps-S60 intersected | 11 | 11 | 17 | 21 | 34 | 75 |

The intersected ledger dominates everywhere (ties k=9, k=13): intersecting with G_P beats
the union bound by 2–11 diameters per leg. k=13 row reproduces monad-S2 exactly.

## 2. Finite band (VERIFIED, exhaustive): diam 12..26 all clear at margin +0.30

All ~229k affine-normalized 8-shapes (0 = min, canonical vs mirror) of diameter 10..26,
two-stage grid screen + refinement: **min μ = 0.9719** (at `{0,1,3,4,5,6,7,10}`, D=10);
per-diameter minima 0.972–0.989, flat in D. Extends codex HYP-3530's span ≤ 13 scan. With
kps-S60 (diam ≤ 11 PROVED) the leg is now open only for **diam ≥ 27** — where every
observed μ is ≥ 0.97 vs the 0.675 bar.

## 3. The far case is nearly free at k=8 (and why)

For E = E₇ ∪ {f} (core diam T₀, f far): measured `Bis ≤ 0.0044` at every tested core and
ratio — μ ≥ 0.9956, far above the bar. Mechanism: before the far point can even matter,
the 7-point core orbit must be in a **near-net state** (unique big gap ≤ 2/7 — rare).
Rigorous skeleton with explicit constants:
`Bis ≤ Ind + Δ`, `Ind ≤ (1/7)·(1 − μ_{2/7}(E₇))`, `Δ ≤ C·T₀/dist(f,E₇)` (crossing count,
crude C ≈ 21: piece count Σ_{a<b}|e_a−e_b| ≤ 21·T₀, sweep rate ≥ dist − T₀) ⟹
**dist(f, E₇) ≥ ~96·T₀ ⟹ μ_{1/7}(E) ≥ 0.675.** The formal write-up of the crossing count
is the remaining (mechanical) step.

## 4. Criticality: the tournament/phase-clock frame

The θ-successor digraph `D_θ(x)` on the orbit (i→j iff `frac((e_j−e_i)x) ∈ (0,θ]`) has
**`E[#edges] = k(k−1)·θ` exactly for every E** (each pair distance is exactly uniform —
the kps-S59 forgotten factoid; verified 7.9997 ≈ 8). At `(k,θ) = (8, 1/7)`: **mean
out-degree exactly 1 — the binding k=8 leg is the critical-branching case** of the
phase-clock comparator movie (THM-373; the movie's walls ARE the crossing-lemma pieces;
k ≤ 7 is subcritical AND deterministic via pigeonhole; k=13 is supercritical at 12/7).
`B = #gaps > θ` = the count of out-degree-0 vertices; `μ = P(B ≥ 1)` exactly. This is the
observer-source lens (THM-381) applied per-runner at threshold θ.

## 5. THE PAIR-MASS LAW (new; table-exact; proof = elementary lattice count)

For coprime `q1, q2 ≥ 1`, θ = 1/7, same-sign windows:
> **`meas{frac(q1x) ∈ (0,1/7], frac(q2x) ∈ (0,1/7]} = 1/49 + G(r1,r2)/(q1q2)`,
> `G(r1,r2) = min(r1,r2)·(7 − max(r1,r2))/49`, `r_i = q_i mod 7 ∈ {1..7}` (0 ↦ 7).**
Verified CONSTANT on every one of the 49 residue classes over all coprime pairs q ≤ 40
(and the extended q₂ ≤ 200 scan's min is exactly 1/49). Hence **m ≥ θ² always** (equality
iff 7 | q1q2 — the apex-7 invisibility of opus-S134, reappearing as the vanishing rows of
G). MIXED-sign masses can vanish exactly ((1,−1), (−3,4), (2,−5) = 0) — so trees must
live on a single sign class. Proof path (half page): the joint residue offsets
`(jq2 − lq1)/(q1q2)` hit each 1/(q1q2)-grid point once (Bézout); the overlap sum is a
sampled cross-correlation of two 1/7-intervals = θ² + a triangle lattice count = G.

## 6. The Hunter-endpoint floor: diameter-free uniform μ-floor at the criticality

At k=8 the Bonferroni base `1 − 7θ = 0`, so Hunter's inequality on the 7 hit events of an
ENDPOINT vertex (all differences same-sign) gives, with any spanning tree T:
> **`μ_{1/7}(E₈) ≥ meas W_endpoint ≥ Σ_{(a,b)∈T} m(d_a,d_b) ≥ 6·θ² = 6/49 ≈ 0.1224`**
— the **first diameter-free rigorous uniform floor at the binding leg** (modulo the
half-page pair-mass lemma). The k=9 variant survives: `≥ 8θ² − 1/7 = 1/49 > 0`; dies at
k ≥ 10 (base too negative). Adversarial min of the best-vertex signed Hunter floor:
0.1753 (consistent, > 6/49). It does NOT reach 0.675 — its role is the R-route below.

## 7. The true object: the union bound wastes 6× at k=8

Direct adversarial minimization of `ρ* = meas(G_P ∩ Good_E)` over joint (|P|=5, E):
**min ρ* = 0.3403 = 6.0× m_P** (at P={1,5,9,11,13}, E = AP), and quasi-independence
**R = ρ*/(mG·μ) never dips below 0.913** (THM-530-D's probe floor was 0.796; kps-S60's
ledger anatomy 0.6–1.06 concerns the ILedger bound, not true ρ*). If `R ≥ 0.75` were
proved uniformly, the k=8 demand on μ drops from 0.675 to **0.197** — within reach of a
2-vertex union of Hunter floors. The three-piece program: pair-mass law (elementary) +
2-vertex Hunter ≥ 0.197 (needs pair-intersection upper bounds) + R ≥ 0.75 (spectral
overlap, THM-579-shaped) would close the k=8 spread residual without the 0.675 bar.

## 8. PZ-on-B: the two-moment route to the full bar

`μ = P(B ≥ 1) ≥ E[B]²/E[B²]`: adversarial min over 8-sets = **0.7122 ≥ 0.675**, attained
at the AP (margin 5.5%, thin but real; spread shapes ≥ 0.92). So target 1 reduces
(adversarially-verified) to two moment bounds: `E[B]² ≥ 0.675·E[B²]`, both sides exact
difference-set integrals (`E[B] = Σ_i meas W_i`; `E[B²] = E[B] + 2Σ_{i<j} meas(W_i∩W_j)`)
— the k=8 mirror of monad's k=13 (S1, pairSum) program, at a high bar with thin margin.

## Honest status

- PROVED (or one elementary lemma away): pair-mass law; Hunter-endpoint 6/49 floor;
  union D_k* table; (with opus-S135) everything previously "modulo roof".
- VERIFIED: band diam 12..26 (min 0.9719); far-case Bis ≤ 0.0044; PZ-on-B ≥ 0.712;
  ρ* ≥ 0.34, R ≥ 0.913 (adversarial, k=8).
- OPEN: the k=8 leg above diam 26 at the full 0.675 bar (via PZ-on-B moments), or at the
  0.197 bar via R ≥ 0.75 + 2-vertex Hunter. Nothing here closes the leg or LRC(14).
- Scripts: `lrc14_k8_leg_composite_klein_S155.py`, `lrc14_k8_criticality_hunter_tournament
  _klein_S155.py` (part 3 of it has the SIGN BUG — superseded by), `lrc14_k8_hunter_signed
  _klein_S155.py` (+ .outs). The unsigned part-2 Hunter numbers (0.62 "floor" at AP interior
  vertices) are INVALID — mixed-sign hit events were conflated; only the signed rerun counts.

## Pointers
kps-S60 HYP-4847 (intersected ledger; names the bisection recursion as next tool — the
far-case skeleton here is that tool's first k=8 instance), opus-S134/S135 (roof, F₇
anchors, apex-7 invisibility = G's zero rows), monad-S1/S2 (bars; CE program), mac-mini-S41
(PZ-on-U; CE=PZ unification — PZ-on-B is the count sibling at k=8), THM-373/374/381
(phase-clock/observer-source frame), THM-579 (projection/CV shape for the R-route),
LRCHunterLedger (the repo's Hunter lineage), HYP-3530/3531 + LTT-431 (the dormant thr_k
scan this extends), klein-S154 HYP-4781 (bisection identity), MISTAKE-119 (jump adversaries
used throughout).
