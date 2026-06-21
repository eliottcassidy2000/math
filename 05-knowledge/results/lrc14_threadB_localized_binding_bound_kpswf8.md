# LRC(14) Thread B — the localized err+ bound on the BINDING family via THM-546

**kind-pasteur-2026-06-21-S26-wf8.** Engines: `lrc14_threadB_*_kpswf8.py` (exact rationals);
CAP/Q from `lrc14_wide_branch_ridge_codex_s47.py`, moment-dual `boundary_value_direct` from
`lrc14_signed_multifar_boundary_hierarchy_codex_s51.py`. THM-546 `|Δ_w| ≤ (6/49)V/w` is the
PROVED comb bound (signed-Abel form).

## Setup (exact)

For `E = B ∪ F`, `B` a bounded base (`0∈B`, `max B ≤ 14`), `F` the far runners (`>14`):
- **Single far** `p0(B∪{w}) = Φ(B) + Δ_w`, `Φ(B)=p0(B)+(1/7)p1(B)` (HYP-2642). VERIFIED exact:
  `w·Δ_w = wDelta_signed(B,w)` matches on every tested config.
- **THM-546 comb bound** `|Δ_w| ≤ (6/49)V(B)/w`, `V(B)=Σ_{j=1..6} #(periodic arcs of B_j)`.
  VERIFIED: 0 violations over 2600 (B,w) pairs (w=15..79); worst ratio `|Δ_w|/[(6/49)V/w] = 0.149`
  (the bound is ~6.7× loose, more at small w).

caps / plateaus / margins:

| k | cap_k | Q(k-1) | margin |
|---|-------|--------|--------|
| 8 | 0.38146 | 0.19660 | 0.18486 |
| 9 | 0.49426 | 0.36210 | 0.13216 |
| 10 | 0.60440 | 0.44789 | 0.15651 |
| 11 | 0.72527 | 0.53125 | 0.19402 |
| 12 | 0.85714 | 0.60224 | 0.25490 |

## (1) SINGLE FAR — PROVED for the tail, finite-checked for the band

**Crude binding RHS does NOT close.** `Q(k-1) + (6/49)V(consec_{k-1})/15` = 0.548 / 0.729 / … —
exceeds cap for every k. At the smallest wide speed `w=15`, the THM-546 bound `(6/49)V/15 ≈ 0.35`
swamps the margin AND is ~26× loose. **THM-546 alone does NOT prove the binding (small-w) case.**

**Per-base combined cutoff (the right object).** `p0(B∪{w}) ≤ Φ(B) + (6/49)V(B)/w < cap_k` whenever
`w > W'(B) := (6/49)V(B)/(cap_k − Φ(B))`. Max over all bounded primitive bases (exact):

| k | #bases | uniform cutoff maxW' | argmax base (slack) | Φ | V |
|---|--------|----------------------|---------------------|---|---|
| 8 | 2996 | **78** | (0,8,9,10,11,12,13) | 0.116 | 167 |
| 9 | 3431 | **78** | (0,8,9,…,14) | 0.234 | 164 |
| 10 | 3003 | **63** | (0,7,8,…,14) | 0.339 | 136 |

So **single-far wide reduces to a finite band `15 ≤ w ≤ maxW' (≤78)`** over bounded bases.
- **TAIL `w > maxW'`: PROVED < cap** by THM-546 (`p0 ≤ Φ(B)+(6/49)V(B)/w < cap_k`).
- **BAND `15 ≤ w ≤ maxW'`: exact-checked.** Over all bounded primitive bases × wide w:
  - sup actual p0 = **0.2267 (k=8) / 0.3721 (k=9) / 0.4755 (k=10)**, all `< cap` (margin ≥ 0.122),
    **0 violations**. Argmax at small w (15–21): the binding regime is at small w, near consec.

**Binding mechanism (exact, k=9).** consec_8 ∪ {21}: `p0 = 0.37211`, `Φ(B)=Q(8)=0.36210`,
`Δ_w = +0.01001 > 0`. So **p0 > Q(8)** (this is MISTAKE-080: `wide ⟹ p0 ≤ Q` is FALSE), but
`p0 < cap_9 = 0.49426` with margin 0.122. Here `|Δ_w|=0.010` while `(6/49)V/w = 0.262`
(26× loose) — which is exactly why THM-546 cannot prove this config and the finite check must.

**VERDICT (single far): VERIFIED < cap (exact, 0 violations).** PROVED tail (w>78) + finite band
(checked). THM-546 supplies the tail only; the binding small-w part is the direct finite check.

## (2) DOUBLE FAR — simultaneous peel + the precise loose step

**Newton simultaneous-peel identity (EXACT, lhs−rhs = 0):**
```
p0(B∪{f1,f2}) = P_2(B) + [p0(B∪{f1})−Φ(B)] + [p0(B∪{f2})−Φ(B)] + [I_B(f1,f2) − Φ_2(B)]
```
`P_2(B)=boundary_value_direct(B,2)` (decorrelated 2-far limit), `Φ_2(B)=(2p2−p1)/49`,
`I_B = p0(B∪{f1,f2})−p0(B∪{f1})−p0(B∪{f2})+p0(B)` (joint Newton 2nd difference).
The two one-far residuals are over the **bounded** B ⟹ THM-546: `≤ (6/49)V(B)/f_i`.

**THE PRECISE LOOSE STEP (HYP-2776, confirmed).** Iterating THM-546 — peel f2 first, then f1 —
leaves the intermediate base `B∪{f1}`. For the **adjacent** pair `f2=f1+1`:

| f1 | V(B) | V(B∪{f1}) | (6/49)V(B∪{f1})/f2 | C = I_B−Φ_2 |
|----|------|-----------|--------------------|-------------|
| 15 | 45 | 63 | 0.4821 | 0.00212 |
| 80 | 45 | 163 | 0.2464 | 0.01116 |
| 320 | 45 | 546 | 0.2083 | 0.01152 |

`V(B∪{f1}) ~ 7·f1` (the far element fragments each B_j into ~7 arcs), so
`(6/49)V(B∪{f1})/f2 → (6/49)·7 = 6/7 ≈ 0.857`, a **CONSTANT, not →0**. The iterated single-peel
bound is loose by O(1) and **cannot** show the curvature → small. This is the exact failure point.

**The joint curvature C = I_B − Φ_2 is a genuine 2D object.**
- ADJACENT `f2=f1+1`: C **saturates** to a fixed nonzero value (~0.011–0.016) as f1→∞
  (the relation `f2−f1=1` never decorrelates). NOT bounded by `(6/49)V/f1 → 0`.
- DISSOCIATED (no small relation): C **decays** toward 0 (off-resonance Fourier).
- `sup|C| ≈ 0.026–0.029` over bounded bases (small + adjacent-saturated far pairs).

**Certified double-far majorant** `p0 ≤ P_2(B) + (6/49)V(B)(1/f1+1/f2) + sup|C|`:
- At w=15 the two peel terms (`2·(6/49)V/15 ≈ 0.70`) swamp everything ⟹ RHS ≈ 1.0–1.2 > cap.
- With the uniform `sup|C| ≈ 0.029` input, the RHS drops `< cap` for far speeds `≥ 99/166/116`
  (binding base consec_{k-1}, k=8/9/10); below that, a finite check closes it.

**VERDICT (double far): NOT closed by single-peel alone.** Identity exact; two residuals peel
(THM-546, bounded base); the joint curvature `C=I_B−Φ_2` is the SOLE open input — a uniform
`sup|C| = sup|I_B−Φ_2|` over all bounded B and far pairs (the joint 2D Erdős–Turán–Koksma
constant = **OPEN-Q-108**). All measured `|C| ≤ 0.029 ≪ margin`; 0 cap-threats observed.

## Reconciliation with concurrent THM-563 (mac-mini-S6) + HYP-2788

THM-563 (pushed concurrently) proves `w·Δ_w` is **exactly periodic** in w (period `7·lcm(B)`, a
generalized Dedekind sum), so the **signed period-max** replaces THM-546's lossy `(6/49)V`. For
consec_{k-1} the period-max is `1, 43/49, 1007/980 ≪ 15·margin` ⟹ single-far closes for **all
w≥15 with no w-window** — strictly sharper than my THM-546 tail+band. My finding "THM-546 is 26×
loose at small w, cannot prove the binding config" is exactly the HYP-2784 "125× lossy" wall that
THM-563's signed periodicity resolves.

HYP-2788's dichotomy + THM-563 give the actual closure that **sidesteps the joint 2D curvature**:
near-cap (`p0 > Q(k-1)`) ⟹ single-perturbation ⟹ single-far ⟹ THM-563 closes; genuine multi-far
⟹ `p0 < Q(k-1)` (slack floor) ≪ cap. **Confirmed exactly:** every double-far `consec_7∪{f1,f2}`
(k=9) is genuine-wide, `p0 ∈ [0.24, 0.29] < Q(8)=0.362`, margin ≥ 0.20. The joint curvature
`C=I_B−Φ_2` peaks at +0.029 (pair (18,19)) but **there p0=0.258** — deep in the slack regime. So
the joint curvature is **never cap-threatening**; the open `sup|C|` bound (OPEN-Q-108) is not on
the critical path for the cap-level wide bound once HYP-2788+THM-563 are used.

## Net

- **Single-far wide: VERIFIED < cap (exact, 0 violations); PROVED tail (w>78) + finite band.**
  THM-546 proves only the gapped tail; the binding (small-w, near-consec) case is the direct
  finite check, because THM-546 is 26× loose there.
- **Double-far wide: the loose step is pinned down** — iterated peeling gives a 6/7 (O(1)) bound
  in the adjacent slope-1 case, never →0. The residual is a uniform bound on the joint curvature
  `sup|I_B−Φ_2| ≈ 0.029` (OPEN-Q-108), which would close double-far with comfortable margin.

Scripts: `04-computation/lrc14_threadB_{localized_binding_bound,single_far_cutoff,combined_cutoff,
double_far_loose_step,cutoffs_fast,supC}_kpswf8.py`. Outputs in `05-knowledge/results/`.
