---
id: THM-548
title: The boundary-value decomposition for LRC(14) true-wide clusters — p0(B∪F) = P_r(B) + resonance corrections, where P_r(B) is the fully-decorrelated Fatou limit (≤cap with margin growing in r) and the corrections vanish off the resonance set (small additive relations among the far runners = Freiman structure = boundary-function "ambiguous points")
status: PARTIAL (the decomposition is exact and finite; the limits Φ_t and P_r(B) are derived and VERIFIED; P_r(B)≤cap−margin verified on sample bounded B with margin growing in r; the two-far curvature saturates/bounded. REMAINS: exhaustive P_r(B)≤cap over bounded B, and the signed resonance-correction bound + its Freiman reduction.) Complements codex HYP-2679 (exact two-far atlas).
source: mac-mini-2026-06-20-S3
depends_on:
  - THM-546   # one-far signed (6/49) bound
  - THM-547   # boundary-collar closure (bounded base, finite maxima)
  - THM-531   # scale-invariance — resonant far tuples reduce to bounded
  - HYP-2679  # codex's two-far curvature atlas (this is the analytic backbone)
related:
  - HYP-2678  # true-wide Ruzsa/Plünnecke program
  - HYP-2657  # QR reality (signedness)
  - HYP-2606  # covolume / relation-height (AP minimizes; resonance structure)
  - HYP-2637  # Freiman dimension penalty
  - OPEN-Q-108
external: Lonely Runner Conjecture; Fatou nontangential limits (bounded harmonic functions); Bagemihl ambiguous points; Kaczynski boundary functions & curvilinear convergence.
---

# THM-548 — The boundary-value decomposition for true-wide clusters

Region III of the LRC(14) sector crux is the **true-wide** case: `E = B ∪ F`, `B ⊆ {0,…,14}`
bounded, `F = {f_1<…<f_r}` the far runners (`>14`), `r ≥ 2`. THM-547 closed `r=1` (boundary
collar). This theorem handles `r ≥ 2` by the multivariate far-element recursion.

## 1. The exact finite Newton expansion

Since the runners in any `S ⊆ F` can newly complete a missed-sector configuration only when `B`
misses **exactly** the sectors they supply, and `B` misses at most 6 sectors, the mixed-difference
(Newton forward-difference / inclusion–exclusion) expansion **terminates at `|S|=6`**:
> `p0(B∪F) = Σ_{S⊆F, |S|≤6} Δ_S(B)`,  `Δ_S(B) := Σ_{T⊆S}(−1)^{|S|−|T|} p0(B∪T)`.
`Δ_∅=p0(B)`; `Δ_{f}` = one-far increment (THM-546/547); `Δ_{f,f'}=I_B(f,f')` = the **two-far
curvature** (codex HYP-2679); `Δ_S` for `|S|≥3` = higher curvatures.

## 2. The decorrelated limits (Fatou boundary values) — DERIVED + VERIFIED

As the runners in `S` go to `∞` independently (coprime), each `Δ_S(B)` tends to a limit `Φ_{|S|}(B)`:
- **`Φ_1(B) = p1(B)/7`** (one-far plateau, the THM-546 `Φ`).
- **`Φ_2(B) = (2 p2(B) − p1(B))/49`** — VERIFIED exact (`B=consec_8`: `47/24010≈0.00196`;
  `B=(0,1,2,4,8)`: `1/98`). `p_t(B)=meas{B misses exactly t sectors}`.
- In general `Φ_t` is the inclusion–exclusion of the miss-profile, weight `7^{-t}`.

Summing all limits gives the **fully-decorrelated boundary value** (all `r` far runners independent):
> **`P_r(B) := lim p0(B∪F) = Σ_{t=0}^{6} prof_t(B)·c_t(r)`,  `c_t(r)=Σ_{i=0}^{t}(−1)^i C(t,i)(1−i/7)^r`**,
`prof_t(B)=meas{B misses exactly t sectors}`, `c_t(r)=P(r iid uniform runners hit all t given sectors)`.
This is the **Fatou-type boundary limit** of `p0` as `F → ∞` (the user's bounded-harmonic-function lead).

**VERIFIED:** `P_r(B) ≤ cap_{|B|+r}` with margin **growing in `r`** — for `B=consec_8`,
`margin(P_r) = 0.132, 0.172, 0.277, 0.380, 0.482` at `r=1..5`; other bounded `B` have larger margin;
and the actual `p0(B∪F)` for a dissociated (Sidon-like) far set `F` tracks `P_r(B)` to `~0.01`.

## 3. The resonance corrections (Bagemihl ambiguous points) — the remaining content

> `p0(B∪F) = P_r(B) + R(B,F)`,  `R(B,F) := Σ_{S⊆F, 2≤|S|≤6} [Δ_S(B) − Φ_{|S|}(B)]`
> `         + Σ_f [Δ_f(B) − Φ_1(B)]`  (the one-far residuals, bounded by THM-547 since `B` is bounded).

Each curvature deviation has the Fourier form (two-far case shown):
> `Δ_{f,f'}(B) − Φ_2(B) = Σ_{(m,n)≠(0,0)} ŝ_j(m)ŝ_{j'}(n)·1̂_{A}(−(mf+nf'))`,
controlled by the frequency **`mf+nf'`**. Off resonance (`|mf+nf'|` large for all small `(m,n)`,
i.e. `f,f'` share no small additive relation) the deviation `→0`. The deviation is large only at
**resonances** `mf+nf'≈0`, e.g. `f'=f+1` gives `(m,n)=(1,−1)`, `mf+nf'=−1`; VERIFIED this saturates
to a **bounded** value `≈0.0139` as `f→∞` (not divergent), `≪` margin.

**Resonance = additive structure.** A small relation among far runners ⟺ they lie in a common
low-dimensional GAP (Freiman) ⟺ scale-invariance (THM-531) reduces the resonant tuple to a bounded
model. This is the LRC incarnation of the boundary-function dichotomy: a function has a clean
nontangential (Fatou) boundary value except at a structured set of **ambiguous points** (Bagemihl),
which here is the resonance set, and which is exactly where the additive-combinatorial (Ruzsa/
Plünnecke/Freiman) reduction applies.

## 3b. The two-far constant C₂ and the apex-prime hierarchy (NEW, verified)

A **second** Abel summation (antiderivative in each frequency) gives the parabolic analogue of the
one-far sawtooth `F_j`: `G_j` = second centered antiderivative of `1_{sector_j}−1/7` (piecewise
quadratic). Computed exactly:
> **`C₂ := max_j sup_y |G_j(y)| = 13/1372 = 13/(2²·7³)`**  (all `j` equal, by apex symmetry),
versus the one-far `sup|F_j| = 3/49 = 3/7²`. **Each Abel order adds one power of the apex prime `7`
to the denominator** (`1/7² → 1/7³`), so the `t`-fold curvature constants are suppressed
geometrically `~1/7^{t+1}` — the Newton expansion of §1 converges fast, governed entirely by `7`.

**QR-reality of the product kernel — VERIFIED (was conjectural).** `ŝ_j(m)ŝ_{j'}(n) +
ŝ_j(−m)ŝ_{j'}(−n)` is **real** (imag part `=0` exactly, all `j,j',1≤m,n≤7`) — the joint
`(m,n)↔(−m,−n)` pairing (the `6=−1` non-residue mod 7 argument, HYP-2657, extended to the product)
kills the imaginary part, so the two-far deviation is real and the signed/Abel form is licensed.

**Resonance-gated bound — VERIFIED.** `|I_B(u,v) − Φ₂(B)|·resdist(u,v)` is bounded `(~0.01)` where
`resdist = min_{small (m,n)} |mu+nv|`; so `|I_B − Φ₂| ≤ C·V(B)/resdist`. The worst case (resonant,
`resdist=1`, e.g. consecutive `u,u+1`) gives `|I_B−Φ₂| ≤ ~0.013 ≪` margin `0.25`. **The two-far
curvature is never cap-threatening.** (Consistent with the synthesis correction: the actual k=9
leader has *negative* curvature `I_B(15,16)=−13/1470`; positive curvature is sub-critical.)

## 4. The closure (target)

`p0(B∪F) ≤ cap_k` for all true-wide `E` follows from:
- **(i) `P_r(B) ≤ cap_{|B|+r} − margin_r`** for all bounded `B` (a finite check; margin grows in `r`).
  VERIFIED on samples; exhaustive sweep pending.
- **(ii) `|R(B,F)| ≤ margin_r`** — the resonance corrections: small for dissociated `F` (off-resonance
  Fourier bound), Freiman-reducible for structured `F`. The signed per-curvature bound (extending
  THM-546's `(6/49)` to `mf+nf'`) + a bound on the number/size of simultaneous resonances is the
  remaining analytic content. The growing margin in `r` is the budget.

**Net:** true-wide is reorganized as a **Fatou boundary value (safe, finite-base) + a resonance
correction (Bagemihl-structured, Freiman-reducible)**. The boundary-function lead is not a metaphor
here — it is the exact shape of the `F→∞` limit and its exceptional set. LRC(14) not proved; this is
the architecture and the verified main term.
