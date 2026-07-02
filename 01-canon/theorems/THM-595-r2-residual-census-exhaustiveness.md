# THM-595: The r=2 residual census-exhaustiveness theorem (assembled written proof)

**Claim:** every 11-core `C` (a set of 11 distinct positive speeds arising as the bounded core of the r=2 multi-far residual of LRC(14)) satisfies `meas(L_C(1/14)) ≥ 1/36`, with minimum `313/9702 ≈ 0.03226` attained by the pentagon core `{1,2,3,4,5,7,8,9,11,12,13}`.

**Status:** ASSEMBLED PROOF — complete modulo two explicitly-scoped legs (G1 census standard, G2 compact-deep floor), each finite/structural and stated precisely below. This document is the written assembly requested by kind-pasteur-S27 (next-step a) and opus-S32 (handoff 1).
**Author:** mac-mini-2026-07-01-S94 (HYP-3850(b))

---

## 0. Why this statement matters

The LRC(14) covering program (OPEN-Q-108 critical path) discharges the r=2 residual — covering 13-sets with two "far" elements removed — by a lower bound on the lonely measure of the remaining 11-core at the critical radius. The constant `1/36` is the downstream requirement (kps-S25/S26: PoA ≤ 6.61 target; kps-S28: `1/36 = (1/φ(7))²` is the 6×6 sector-pair grid quantum, not an extremal accident).

## 1. Verified inputs (machine-checked, exact rational arithmetic)

- **(V1) The census** [kps-S27, `near_tight_11core_census_completion_kps.py`]: over 15,472 systematically enumerated 11-cores — all k-drops of `{1..V}` for the dilated-AP + Goddyn–Wong tight locus (`2/13, 3/14, 4/15, 5/16, 6/17`-sample), their dilations (2×, 3×, 5×), and 8,000 random primitive cores — the minimum of `meas(L_C(1/14))` is `313/9702` at the pentagon core; **zero cores below 1/36**; margin 1.161×.
- **(V2) Multi-outlier monotonicity** [kps-S28, HYP-3950]: the minimum over 11-cores with `r` outliers *increases* in `r`: `0.032261 (r=0) / 0.032311 (r=1, exact 40753/1261260) / 0.034870 (r=2, exact 733/21021) / 0.056+ (r=3)`; 120k configs, 0 below 1/36. Bounded speeds are strictly better killers than outliers.
- **(V3) Window cutoffs** [kps-S28]: finite far-thresholds `W*(r) = 112, 181, 290, 513` (`r = 1..4`) — beyond `W*`, the peel lemma (L1) discharges the far element with the ledger guard; the windows `[V_census, W*]` were scanned and all minima found inside, clear of 1/36.
- **(V4) Tower guards** [opus-S32, HYP-3834]: `(1 − j/7)·M_{11−j} > 1/36` for all `j = 1..6` with minimum margin `0.0130` (`Λ = 10^4` gap-cut chains), where `M_k` is the k-core floor from the census/tower recursion `M_k = min(census_k, min_j (1−j/7)M_{k−j})`, `M_k = 1 − k/7` unconditionally for `k ≤ 6`.

## 2. Proved lemmas (cited or this session)

- **(L1) Simultaneous-peel lemma** [opus-S32 HYP-3834; kps-S28 sharper constant]: for `C = B ∪ F` with `B` compact and `F` far (`j = |F|`),
  `meas(L_C) ≥ (1 − j/7)·meas(L_B) − (N_B/7)·Σ_{w∈F} 1/w`,
  where `N_B` is the arc count of the *compact base's* lonely set (the correct object per opus-MISTAKE-090: the far-prefix arc-count version is false, `c ~ w·meas`). At a multiplicative gap cut the error is scale-free (`≤ 22j/(7Λ)`).
- **(L2) Scale invariance + quantization** [THM-522]: `L(cS) = L(S)` (measure-preserving `t ↦ ct`), so cores reduce to `gcd = 1`; and `L(S) > 0 ⟹ L(S) ≥ 1/(14·lcm S)` — every census value is an exact rational, so "≥ 1/36" is a finite integer comparison, not a float claim.
- **(L3) Equidistribution of a huge speed** [HYP-3787; kps-S27 verification]: appending `w → ∞` to a fixed core multiplies the lonely measure by `(6/7) + O(1/w)`; verified `0.0854 → 0.0733–0.0737` for `w ∈ {200..5000}`. Hence a far element can never *create* a violation that its peeled core did not already have — the minimum lives at bounded speeds.
- **(L4) Critical-mass floor** [THM-594(E), this session]: any `j = 7` far sub-cluster (mass `A = 1`, exactly where the union-bound guard of (V4) dies) leaves uncovered measure `≥ 2sin²(π/7)/(7π²) ≈ 0.0055` unconditionally — the boundary brick between the peel regime and the renormalization regime.
- **(L5) No finite exact tiling** [THM-594(C)]: no finite distinct-speed danger system tiles the circle; exact-tiling limits (`D₇(k/7) = 0`) occur only along infinite divisor chains — the fixed locus of the scale action. Hence the renormalized far-cluster density `D_F` cannot vanish identically on any interval of `L_B`; the *quantitative* version is leg G2.

## 3. The assembled proof

Let `C` be any 11-core, `gcd(C) = 1` (L2). Let `F ⊂ C` be its far part above the relevant cutoff (`W*(|F|)` from (V3), or a `Λ`-gap cut), `B = C \ F`, `j = |F|`.

**Case j = 0** (all speeds bounded): `C` lies in the enumerated window; by (V1)+(V3), `meas(L_C) ≥ 313/9702 > 1/36`. *(Standard of this step = leg G1 below.)*

**Case 1 ≤ j ≤ 6:** By (L1), `meas(L_C) ≥ (1−j/7)·meas(L_B) − (N_B/7)Σ_{w∈F}1/w`. The base `B` is a bounded `(11−j)`-core, so `meas(L_B) ≥ M_{11−j}` by the tower recursion; the far error is controlled by the gap cut / `W*`; by (V4) the guard clears `1/36` with margin ≥ 0.0130 for every `j`. (V2) confirms the direction independently: adding outliers only raises the minimum.

**Case j ≥ 7:** the far part alone reaches critical mass (`2rj ≥ 1`), and the union-bound guard is vacuous. Two sub-cases:
- *Gapped far part:* if `F` itself splits at a gap cut, recurse — the tower loses ~6 runners per level and terminates in ≤ 2 levels (depth-1 residual ≈ shifted LRC(8), opus-S32).
- *Compact deep cluster* (bounded ratios, height → ∞): by renormalization (opus HYP-3835) `meas(L_C)` converges to `∫_{L_B} D_F^∞`, the difference-core limit; the AP is the fixed point of the difference flow, so the limit extremizers are the census extremizers at height ∞. The ~7,500-config boundary minimax probe (all 715 compact 4-cores × heights) found global minimum `1.97×` pentagon — nothing below the height-1 extremizer. (L4) gives an unconditional positive floor at exact criticality; (L5) says the limit density cannot vanish identically. *(Quantitative closure of this sub-case = leg G2.)*

Combining the cases: `meas(L_C(1/14)) ≥ 1/36` for every 11-core, minimum at the pentagon. ∎ *(modulo G1, G2)*

## 4. The two open legs, scoped exactly

- **(G1) Census standard.** (V1)+(V3) enumerate the tight-locus drop families + dilations + randoms + the `W*` windows — a *structured* enumeration, not the full `C(V, 11)` space. What is missing for symbolic end-to-end rigor: either (i) an exhaustive bounded enumeration below `W*(0)` (finite, large — engineering: the fast C/numpy interval kernel, opus handoff 2, extend `V = 19 → 25+`), or (ii) a structural theorem that the minimum over bounded cores is attained on the drop-family ledger (kps-S28's "ledger = cap atlas" — min over bounded m-cores = `cap_{13−m}` exactly, with known extremizers — is precisely this statement for `m ≤ 5`; extending it to `m = 11` bounded cores closes G1 structurally).
- **(G2) Compact-deep floor.** Prove `∫_{L_B} D_F ≥ c > 0` for compact deep clusters — the renormalized tower floor. Current state: probe evidence (1.97× pentagon, no descent direction), the fixed-point conjecture survived falsification, (L4) positivity at exact criticality, (L5) impossibility of identical vanishing. The missing piece is a quantitative rate: how fast can `∫ D_F` decay along the only escape route (infinite divisor chains)? By (L5) the route requires unbounded divisor towers, which bounded-ratio compact clusters do not contain — a formalization of *that* sentence is the shortest path to closing G2.

## 5. Cross-validation

The assembly is consistent with three independent computations: klein-S87's distribution-function machinery reproduces THM-523's single-perturbation infimum `1/1260` exactly; mac-mini-S93/THM-593 reproduces the tight-locus dichotomy `{AP, GW}` used in the census families; THM-594(B) turns kps-S28's empirical pair-overlap law (used in the sector-pair moment reading of 1/36) into exact arithmetic.

-> THM-522, THM-523, THM-594, HYP-3834, HYP-3835, HYP-3950, HYP-3787, kps-S27/S28, opus-S32, OPEN-Q-108.
