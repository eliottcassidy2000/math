# HYP-3950: LRC(14) r=2 residual, multi-outlier 11-cores — the ledger + uniform arc-count union floor

**Status:** VERIFIED (census + elementary floor lemma + finite windows; r>=7 slack probes) — kind-pasteur-2026-07-01-S28.
**Numbering note:** adopting klein-S68's per-machine HYP-block proposal (kps = 3950+).

## Context
The r=2 residual of the multi-far moment relaxation reduces LRC(14)'s open covering branch to
**inf meas(L_C) >= 1/36 over all 11-speed cores C** (band 1/14). kps-S27 closed the *bounded* case
(min = pentagon 313/9702 = 0.032261, 1.161x). Single-outlier min = {1,2,3,5,7,8,9,11,12,13}+50 at
40753/1261260 = 0.032311 (EXACT, 1.1632x) — the near-minimal 10-core plus one far runner.
This session closes the **multi-outlier residual** (>= 2 far speeds) to the same evidence standard.

## Results (all in lrc14_multi_outlier_flat_ratio_peel_kps.out + _extension_kps.out)

### 1. FLAT-RATIO LAW REFUTED — in the favorable direction (monotone-in-r minima)
min_m (bounded m-cores) x (6/7)^(11-m) GROWS as m drops: 0.0323 (m=11) -> 0.0484 (10) -> 0.0679 (9)
-> 0.1032 (8) -> 0.1174 (7). Bounded speeds are strictly more efficient killers than outliers.
**Observed minima by outlier count r: r=0: 0.032261 | r=1: 0.032311 | r=2: 0.034870 | r=3: >= 0.056.
The multi-outlier case is NON-BINDING and monotone increasing in r.** (Measure-side analogue of S23's
base-size domination for p0.) 98399 two-outlier + 21491 three-outlier configs, **0 below 1/36**.
2-outlier min = 733/21021 = 0.034870 (EXACT, 1.2553x) at {1,2,5,6,7,8,9,11,13} + (19,20) — a
CONSECUTIVE pair at moderate scale, and NOT over the minimal 9-core (arc-alignment beats base-minimality).

### 2. THE UNIFORM ARC-COUNT UNION FLOOR (the elementary theorem, = "the uniform arc-count bound")
Per-speed lemma (PROVED, elementary): for any arc-union L with N arcs and any integer speed w,
    meas(L ∩ D_w) <= meas(L)/7 + N/(7w)
(D_w = w teeth of width 1/(7w), 1/w-periodic; an arc I meets <= |I| w + 1 teeth). Verified: worst
danger/bound = 0.90 over ledger cores x w <= 1000. Hence the **multi-outlier union floor**
    meas(L_{B u W}) >= meas(L_B) (1 - r/7) - (N_B/7) Σ_j 1/w_j     (r = |W| outliers)
with N_B = arc count of the bounded core's lonely set (pentagon: N=8; all ledger minima: N <= 26).
**Finite cutoffs** W*(r) = sup over bounded bases of (N_B r/7)/(meas_B (1-r/7) - 1/36):
W* = 112 / 181 / 290 / 513 for r = 1..4 (over the ledger pool); outliers > W* are PROVED safe;
outliers <= W* are a finite box (scanned: 0 violations; global minima found inside, all clear).
**Union floor + ledger clears r <= 6** (m >= 5 bases: floors +0.054, +0.080, ... all >= 1/36).

### 3. LEDGER = CAP ATLAS (cross-validation with canon)
min over bounded m-cores = the sector route's cap_{13-m}, at the SAME extremizers:
min_2 = 66/91 = cap_11, min_3 = 55/91 = cap_10, min_4 = 1979/4004 = cap_9 @ {1,11,12,13},
min_5 = 2243/5880 = cap_8 @ {1,5,7,8,9} (exact matches). **Two extremal phases:** {1} u top
(Dirichlet-type) for m <= 4; odd-unit-heavy near-AP for m >= 5 — the phase transition is THM-576's
cap-pattern break at j=5. Greedy-drop NESTING within the near-AP phase: argmin_11 ⊃ argmin_10 ⊃ ...
⊃ argmin_7 (drops 4, 12, 13, 2 in order); nesting BREAKS at the m=4/5 phase boundary.

### 4. r >= 7 (union floor vacuous): packed clusters are SLACK and scale-invariant
meas({W..W+r-1}) = 0.85-0.88 x (6/7)^r, identical at W = 200/1000/5000 (the 2-torus offset limit is
already exact — scale-invariance). Worst mixed config found: B={1,2,3} + 8-packed cluster = 0.0895
(3.22x above 1/36). Adversarial variants (resonant 84-step AP cluster, paired doublets, two-scale)
all >= 3.9x. **Residual analytic lemma for full rigor at r >= 7:** the multi-dimensional
equidistribution rate for clusters (HYP-3787's O(1/w) is proved only <= 6-far) — with 3.2x margin.

### 5. EXACT PAIR-OVERLAP LAW (Farey-collision threshold AT n=14)
**[UPDATE 2026-07-01-S29: mac-mini THM-594 PROVES the first branch (2r/q iff r(p+q) <= 1, Fourier)
and CORRECTS my "-> 1/49 exact independence" gloss: beyond the threshold the exact overlap is the
wrapped-branch formula 4r^2 - 2(1-2rp)(1-rq)/(pq), which FLUCTUATES (1/55, 5/231, 3/143, ...);
independence 1/49 holds exactly only on the curve 7(q+2p-14) = pq. My coprime spot checks
(0.02037-0.02045) were consistent with fluctuation but I over-summarized. Defer to THM-594.]**
meas(D_w ∩ D_w') depends only on the reduced ratio p/q (scale-invariance, verified):
    **F(p,q) = 1/(7 max(p,q)) when p + q <= 14** (only the 0-teeth meet: Farey separation
    |i/p - j/q| >= 1/pq > (p+q)/(14pq) iff p+q < 14), transitioning to **1/49 (independence)** for
    p + q > 14 (torus filling; coprime spot checks 0.02041 ~ 1/49).
Exact multiples: F(k,1) = 1/(7k) — verified k=2..10. So a pair of outliers can beat independent
removal by at most 1/49 - 1/91 ~ 0.0094 (worst at ratio 13), and resonant small ratios (p+q <= 14,
max <= 7) OVERLAP MORE than independent — resonance wastes the adversary's budget. The n=14 apex
appears as the Farey-collision threshold in the two-runner correlation constant.

## Why the multi-outlier case closes (the abstract mechanism)
**The adversary pays in arc count.** To push meas(L_B) down, the core must concentrate its lonely set
into FEW arcs (near-tight cores have N = 8-14; loose ones N = 18-26). But the outlier bonus above the
1/7 cap is exactly the arc-boundary tax N_B/(7w) — so minimizing measure SIMULTANEOUSLY caps how much
outliers can help. The trade-off is self-defeating: meas-minimal bases give outliers nothing to bite.
That is the reason "a uniform arc-count bound" suffices — N_B is not an enemy to bound but a resource
the adversary must spend.

## CONVERGENCE with opus-S32 (same-prompt collision; opus pushed first — priority theirs on the lemma)
Discovered mid-session: **opus-2026-07-01-S32's HYP-3834 (simultaneous-peel lemma) IS this union floor**
— meas(L_C) >= (1-j/7)meas(L_low) - (2 c_low/7) Σ_{w∈F} 1/w, proved independently hours earlier, with
the Λ-gap-cut chaining + tower M_k and the j=7 covering-threshold boundary ("union bound dies at 7").
Their MISTAKE-090 retires the FALSE version of "uniform arc-count" (a bound on far-containing PREFIX
arc counts — impossible, c ~ w·meas grows; my sec-3 GROW verification N' <= 2N + w·meas exhibits the
same growth). **Terminology alignment: this file's "uniform arc-count" = the bound on the COMPACT
base's arc count N_B (= opus's c_low) entering the one-shot union floor — the CORRECT object, exactly
MISTAKE-090's fix, not the retired prefix version.**
My independent layer beyond opus-S32: (a) the sharper per-arc constant **N_B/7 vs their 2c_low/7**
(factor 2, via the ⌊ℓw⌋+1 tooth count); (b) explicit per-r finite cutoffs W*(r) = sup over bases
(112/181/290/513) — the finite-window quantification complementing their Λ=10^4 gap-cut; (c) the
ACTUAL minima with exact rationals + the monotone-in-r law (their guard table bounds, my census
locates: 733/21021 at (19,20), 40753/1261260 at +50); (d) LEDGER = CAP ATLAS (min_m = cap_{13-m});
(e) the pair-overlap Farey law (sec 5). Their unique layer beyond mine: the deep-cluster
renormalization program (HYP-3835: difference-core flow, AP = fixed point, boundary minimax 1.97x,
Fraenkel/exact-cover lead) — strictly stronger than my sec-4 packed-cluster probes (which match their
{N..N+10} limit 0.158 exactly).

## Honest status
- r <= 6: VERIFIED-CLOSED (ledger + PROVED per-speed lemma + finite windows; windows scanned at the
  binding bases, exhaustible in principle; ledger pool systematic-but-finite like S27's census).
  Convergent with opus-S32 HYP-3834 (their guard table, min margin 0.0130, same conclusion).
- r >= 7: 3.2x+ slack, scale-invariant limits; rigorous closure needs the r-far equidistribution rate
  (only <= 6-far is proved) = opus HYP-3835's tower floor. NOT a proof of the r=2 residual — the same
  evidence standard as S27's census, now covering ALL outlier counts. LRC(14) still NOT proved.

## Artifacts
- 04-computation/lrc14_multi_outlier_flat_ratio_peel_kps.py (+ .out)
- 04-computation/lrc14_multi_outlier_extension_kps.py (+ .out)

## Depends on / relates to
THM-522, THM-523, THM-576 (caps = ledger), HYP-3787, kps-S27 census, THM-563/564, HYP-3824,
klein HYP-3800 (difference-set/phase-residue: the p+q<=14 Farey law is its overlap avatar), OPEN-Q-108.
