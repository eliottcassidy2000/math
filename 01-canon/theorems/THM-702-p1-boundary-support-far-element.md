---
id: THM-702
title: "The p1-boundary support of the far-element seven-sector correction (the missed-sector-phase collapse) — PROVED: the entire 127-atom signed far-element correction Δ_b = meas(S7(E'∪{b})) − Plat(E') equals the discrepancy of the far runner in the SINGLE missed sector, integrated over the p1-region ONLY; by one integration by parts it is a signed sum of the missed-sector antiderivative G_j over the p1-region BOUNDARY, giving |Δ_b| ≤ (6/49)·|∂P1(E')|/b. This is the rigorous form of the missed-sector-phase cancellation (opus-S218's 8× nut): the support collapses to |A|=1 by the linearity H_T = −Σ_{j∈T}G_j and the inclusion–exclusion identity Σ_{T⊆A}(−1)^{|T|}H_T = G_j·1{A={j}}."
status: PROVED (elementary: inclusion–exclusion collapse + one integration by parts; both the P1-integral form and the boundary form verified EXACT-rational with zero mismatch on consec_8, multiscale, and spread cores at b=50/137/200). SCOPE: this is the SUPPORT half of the sharp constant (THM-699/700 item 3) — it proves the correction sees only the p1-region and only the missed-sector phase. The remaining sharpness (the ~17× beyond |∂P1|/b) is the equidistribution of {frac(b·x_i)} over the p1-boundary points, a further discrepancy estimate. The tight margin still lives in the finite consec_k check.
source: opus-2026-07-11-S219
depends_on:
  - THM-699   # mac-mini's TV far-element contraction (this sharpens its constant to the p1-boundary)
  - THM-700   # kps's plateau decorrelation lemma (this is its "sharp constant", item 3, structural half)
  - THM-534   # the seven-sector object meas(S7) and the plateau Plat = p0 + (1/7)p1
related:
  - THM-538   # support concentrates on near-full coverage — here made exact: the correction lives on |A|=1
  - THM-699 (kps) # the zero-mean permanent weight Σ_c D7 = 0 — the WEIGHT side of the same collapse
  - HYP-2655/2664 # the joint plateau/Δ entanglement + the p1-tax — this replaces the crude V by |∂P1|
  - opus-S218 reflection # the missed-sector-phase Fourier nut this proves the support half of
external: bounded-variation Fourier decay; Koksma–Hlawka discrepancy.
---

# THM-702 — the p1-boundary support of the far-element correction

## Setup

`0 ∈ E'`, `b` a far integer offset, `p0(F) = meas(S7(F))` the seven-sector cover measure (all 7 sectors
`[s/7,(s+1)/7)` hit by `{frac(e·x)}_{e∈F}`), `p1(F) = meas{F misses exactly one sector}`,
`Plat(E') = p0(E') + (1/7)p1(E')`. The far-element correction is `Δ_b := p0(E'∪{b}) − Plat(E')`.
Since `0 ∈ E'`, sector 0 is always hit, so the free-sector set `A(x) = {0..6} ∖ {sectors of E' at x}`
satisfies `A(x) ⊆ {1..6}`; write `P1(E') = {x : |A(x)| = 1}` (E' misses exactly one sector `j(x)`).

## The theorem

> **(a) P1-restriction (exact).**
> `Δ_b = ∫_{P1(E')} [ 1{frac(bx) ∈ sector j(x)} − 1/7 ] dx.`
>
> **(b) Boundary form (one IBP).** With `G_j(y) = ∫_0^y (1_{sector j} − 1/7)` (the mean-zero
> missed-sector antiderivative, `G_j(0)=G_j(1)=0`, `‖G_j‖_∞ = 6/49`),
> `Δ_b = −(1/b) Σ_{x ∈ ∂P1(E')} [ P(A_after(x), frac(bx)) − P(A_before(x), frac(bx)) ]`,
> where `P(A, y) = G_j(y)` if `A = {j}` and `0` otherwise. **The sum is supported only on the boundary
> of the p1-region.** Hence
> `|Δ_b| ≤ (6/49) · |∂P1(E')| / b.`

The content is that the full inclusion–exclusion over the `2^6` sector atoms **collapses to a single
missed-sector antiderivative**, supported on `|A|=1`.

## Proof

**(a)** Adding `b` to `E'` fills at most the one sector containing `frac(bx)`. So `E'∪{b}` covers all 7
sectors iff `A(x) = ∅` (E' already covers) or `A(x) = {frac(bx)'s sector}` (E' misses exactly one, `b`
fills it). Thus `1_cover(E'∪{b}) = 1{|A|=0} + 1{|A|=1 and frac(bx) ∈ sector j(x)}`, and integrating,
`p0(E'∪{b}) = p0(E') + ∫_{P1} 1{frac(bx) ∈ sector j(x)} dx`. Subtract `Plat(E') = p0(E') + (1/7)p1(E')`
and use `∫_{P1} dx = p1(E')`. ∎

**(b)** Write `p0 = Σ_{T⊆{1..6}} (−1)^{|T|} A_T`, `A_T(F) = ∫∏_{e∈F} 1{frac(ex) ∉ U_T}`,
`U_T = ∪_{j∈T}[j/7,(j+1)/7)`. For each atom, `A_T(E'∪{b}) − (1−|T|/7)A_T(E') = ∫ F_T(x)(h_T(bx) − h̄_T)dx`
with `F_T` the E'-product (a step function, jumps `J_i^T` at the E'-breakpoints `x_i`), `h_T(y) =
1{frac(y)∉U_T}`, `h̄_T = 1−|T|/7`. The centered antiderivative is **linear in T**:
`h_T − h̄_T = −Σ_{j∈T} g_j` (`g_j = 1_{sector j} − 1/7`), so `H_T := ∫_0^·(h_T−h̄_T) = −Σ_{j∈T} G_j`.
Integration by parts (periodic `H_T`, integer `b`): `∫F_T(h_T(bx)−h̄_T)dx = −(1/b)Σ_i J_i^T H_T(frac(bx_i))`.
Summing over atoms with signs, `Δ_b = −(1/b) Σ_i H⃗_i(frac(bx_i))`, `H⃗_i(y) = Σ_T (−1)^{|T|} J_i^T H_T(y)`.
Now `J_i^T = 1{T⊆A_after} − 1{T⊆A_before}` (each `F_T(mid) = 1{T ⊆ A(mid)}`), so
`H⃗_i(y) = P(A_after, y) − P(A_before, y)`, `P(A,y) := Σ_{T⊆A}(−1)^{|T|}H_T(y)`. By the linearity,
`P(A,y) = −Σ_{j∈A} G_j(y) Σ_{T⊆A, j∈T}(−1)^{|T|} = −Σ_{j∈A} G_j(y)·(−1)(1−1)^{|A|−1} = G_j(y)·1{A={j}}`,
using `Σ_{T⊆A, j∈T}(−1)^{|T|} = (−1)(1−1)^{|A|−1}`, which is `−1` for `|A|=1` and `0` for `|A|≥2`.
So `H⃗_i` vanishes unless `A_before` or `A_after` is a single sector — the p1-boundary. ∎

## Why it matters

- **It is the "missed-sector-phase cancellation" of opus-S218, proved.** The 8×/17× looseness of the crude
  `V(E')/(6b)` bound was the inclusion–exclusion signs among the 127 atoms; THM-702 sums them: they collapse
  to `G_j` on `|A|=1`. The correction literally is *the far runner's discrepancy in the one missed sector,
  over the p1-region only* — nothing else.
- **It replaces the crude `V(E')` (≈ 4·span) by `|∂P1(E')|` (the p1-region boundary).** This is the correct
  small quantity: the p1-region is thin for a well-covering core. It sharpens THM-699/700's constant to the
  right support.
- **It meets THM-699 (kps) on both sides of one object.** kps proved the *weight* `Σ_c D7(c)=0` (zero-mean);
  THM-702 proves the *oscillation/support* side: the correction is a mean-zero `G_j` supported on the
  p1-boundary. Weight-zero-mean × oscillation-on-p1-boundary = the decorrelation.

## Scope — what remains (honest)

THM-702 is the **support** half of the sharp constant. It gives the rigorous `|Δ_b| ≤ (6/49)|∂P1|/b`, but
`|∂P1| ~ span`, so ungapped (`b ~ span`) this is still `O(1)` — the further ~17× is the **equidistribution
of `{frac(b·x_i)}` over the p1-boundary points** (`G_j` mean-zero ⟹ the signed `G_j` sum cancels). That
discrepancy estimate — a Koksma/three-distance bound on the p1-boundary sample under `×b mod 1` — is the
remaining analytic step for the bounded-ratio (ungapped) core. The tight margin remains in the finite check.

## Files
`04-computation/lrc14_sharp_boundary_support_opus_S219.py` (+ `.out`): both forms verified EXACT-rational,
zero mismatch (consec_8, multiscale k7, spread; b = 50/137/200); `‖G_j‖_∞ = 6/49` exact; `|∂P1|/V ≈ 1/3`.
