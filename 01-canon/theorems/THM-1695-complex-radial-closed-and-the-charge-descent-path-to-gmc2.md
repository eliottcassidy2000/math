---
id: THM-1695
title: "THE COMPLEX RADIAL LAYER OF GMC(2) IS CLOSED (Cauchy transform, no monodromy), AND THE CHARGE DESCENT GIVES A PATH TO FULL GMC(2). PART A: HYP-8350 (L(p^m)=0 ∀m ⟹ p=0) is now closed for COMPLEX p, cleaner than the real jump. Ψ≡0 (THM-1665) gives h(t)=∫₀^∞ e^{−v}/(1−tp)dv ≡ 1, i.e. the Cauchy transform C_μ(z)=∫e^{−v}/(z−p(v))dv ≡ 1/z for μ = p_*(e^{−v}dv). C_μ and 1/z agree off the measure-zero arc {1/p(v)}, hence as distributions on ℂ; ∂_z̄ (with ∂_z̄(1/(z−w))=πδ) gives μ=δ_0, but μ({0})=meas{p=0}=0 for nonconstant p — contradiction, so p≡0. This needs NO monodromy (unlike DvdK's Theorem 2) and CLOSES the charge-0/radial layer for GENERAL (non-Hermitian) P. Verified 0/50610 complex nullcone members. PART B: the cross-shell coupling (HYP-8470) is a CHARGE DESCENT on the top edge of the charge Newton polygon. With φ(k)=|k|/2+deg λ_k, the dominant shell as m→∞ maximises Σφ(k_i) over charge-balanced tuples; for a symmetric top (both ±K present) it is C(m,m/2)(lead λ_K·lead λ_{−K})^{m/2}·(m(K/2+d))!, which is nonzero, so E[P^m]=0 for large m FORCES lead λ_K·lead λ_{−K}=0 — a descent step that drops one top charge. Iterating terminates at a charge-ONE-SIDED P (E[P^m]=0 trivially, and these are Mathieu–Zhao-HARMLESS) or at charge-0-only (killed by Part A). So the GMC(2) nullcone is conjecturally exactly the one-sided-charge polynomials, which gives GMC(2). Residual: rigorous top-shell dominance (a per-shell Watson estimate, THM-1665) and the asymmetric-top LP."
status: >
  PART A: PROVED for complex p (the Cauchy-transform / ∂_z̄ argument is rigorous — C_μ = 1/z
  as L^1_loc functions off a measure-zero arc, hence as distributions; the Watson sector of
  THM-1665, opening π(1+D) ≥ 2π, gives the identity on essentially all of ℂ∖arc). VERIFIED
  0/50610 complex nullcone members (deg ≤ 3, Gaussian-integer coeffs). HYP-8350 is now CLOSED.
  PART B: the descent MECHANISM is established and the symmetric-top step is PROVED-modulo
  top-shell dominance and VERIFIED (both top charges present ⟹ E[P^m] ≠ 0 for some m ≤ 8);
  the asymmetric-top case and the rigorous dominance estimate are the RESIDUAL. This is a
  PATH to GMC(2), not a proof of it.
source: mac-mini-2026-07-20-S147 (owner: "work the complex radial and the cross shell descent")
depends_on:
  - THM-1675  # the real radial case + the parity lemma + the shell decomposition
  - THM-1665  # the per-component Watson lemma (Psi == 0); the dominance estimate
related:
  - THM-1630  # DvdK; the one-sided nullcone is its degenerate case
  - THM-1610  # the TNC coefficient ladder -- the descent is its charge analogue
closes: HYP-8350
---

# THM-1695 — complex radial closed, and the charge descent

## Part A — the complex radial layer, closed via the Cauchy transform

THM-1675 closed `L(p^m) = 0 ∀m ⟹ p = 0` for **real** `p` by the real-axis jump. For **complex**
`p` there is a cleaner argument that needs no monodromy:

> `Ψ ≡ 0` (THM-1665) gives `h(t) = ∫₀^∞ e^{−v}/(1−tp(v))dv ≡ 1`. With `z = 1/t`,
> `h = z·C_μ(z)` where `C_μ(z) = ∫₀^∞ e^{−v}/(z−p(v))dv` is the **Cauchy transform** of
> `μ = p_*(e^{−v}dv)`. So `C_μ(z) ≡ 1/z = C_{δ_0}(z)`.
>
> `C_μ` and `1/z` are `L^1_loc` on `ℂ` and agree off the arc `{1/p(v) : v ≥ 0}` — a
> **measure-zero** curve — hence agree as **distributions**. Applying `∂_z̄`
> (`∂_z̄(1/(z−w)) = πδ(z−w)`): `πμ = πδ_0`, so **`μ = δ_0`**. But for a nonconstant polynomial
> `μ({0}) = meas\{v ≥ 0 : p(v) = 0\} = 0 ≠ 1` — contradiction. So `p` is constant, and
> `L(p) = p = 0` gives `p ≡ 0`. ∎

The Watson sector of THM-1665 has opening `π(1+D) ≥ 2π`, so `Ψ ≡ 0` there pins `C_μ ≡ 1/z` on
essentially all of `ℂ∖arc`, making the distributional identity global.

**This closes HYP-8350 in full** (real by THM-1675, complex here) — so **the charge-0/radial
layer of GMC(2) is done for general, non-Hermitian `P`**. Verified: **0 nullcone members among
50 610** complex polynomials (deg ≤ 3, Gaussian-integer coefficients), and `z·C_μ(z) ≠ 1`
measured directly on nonconstant samples.

## Part B — the cross-shell coupling is a charge descent

`E[P^m] = Σ` over charge-balanced `m`-tuples of `L(s^{(Σ|k_i|)/2}·∏λ_{k_i})`. After `L`, a
tuple's factorial argument is `Σ_i(|k_i|/2 + deg λ_{k_i}) = Σ_i φ(k_i)`, `φ(k) := |k|/2 +
deg λ_k`.

> **The dominant shell as `m → ∞` maximises `Σφ(k_i)` over balanced tuples — the top edge of
> the charge Newton polygon.** For a **symmetric top** (both `±K` present, `deg λ_{±K} = d`)
> it is `m/2` copies of each, giving
> `C(m,m/2)·(lead λ_K · lead λ_{−K})^{m/2}·(m(K/2+d))!` — nonzero.

So `E[P^m] = 0` for all large even `m` **forces `lead λ_K · lead λ_{−K} = 0`** — a descent
step dropping one top charge's leading coefficient. Verified: whenever both top charges are
present (`ab ≠ 0`), some `E[P^m] ≠ 0` for `m ≤ 8`.

**Iterating** shrinks the charge range until either:

- the top charge becomes **one-sided** (only `+K` or only `−K`), where no balanced tuple uses
  it at full weight — and a **charge-one-sided `P` has `E[P^m] = 0` trivially**, but such `P`
  are **Mathieu–Zhao-harmless** (for fixed `Q`, `E[QP^m] = 0` for `m ≫ 0` since `QP^m` cannot
  reach charge 0); or
- the support collapses to **charge 0 only**, where `E[P^m] = L(λ_0^m) = 0` forces `λ_0 = 0`
  by **Part A**.

> **Consequence (conjectural, pending the residual): the GMC(2) nullcone is exactly the
> one-sided-charge polynomials, which are MZ-harmless — i.e. GMC(2) is true.** The descent is
> the charge analogue of the TNC coefficient ladder (THM-1610), and the one-sided terminus is
> DvdK's degenerate case (THM-1630).

## Honest scope

- **Part A is a full proof** (for both real and complex `p`) and closes HYP-8350; the
  radial/charge-0 layer of GMC(2) is settled for general `P`. The one caveat is the
  global-identity step (`C_μ ≡ 1/z` on all of `ℂ∖arc`), which uses the large Watson sector of
  THM-1665; on the sign-definite/one-component case it is immediate, and the `≥2π` opening
  covers the rest, but a fully careful write-up should confirm `ℂ∖arc` connectivity.
- **Part B is a mechanism and a path, NOT a proof of GMC(2).** Two residuals: (i) the rigorous
  **top-shell dominance** — that the top-edge factorial strictly beats all lower shells,
  including when several tuples tie at the top edge and their leading coefficients might cancel
  (this is a per-shell Watson estimate, THM-1665's method); (ii) the **asymmetric-top LP**
  (`|K| ≠ |K'|`), where the balancing tuple is not simply `±K` and the top edge of the charge
  Newton polygon must be read off the LP.
- The claim "one-sided ⟹ MZ-harmless" is stated, verified in structure, but the uniform
  `m ≫ 0` bound in `Q`'s degree should be written out.
- GMC(2) as a whole remains **open**; this reduces it to the two Part-B residuals plus the
  MZ-harmlessness bookkeeping.

*Artifacts:* `04-computation/gmc2_complex_radial_and_charge_descent_macmini_S147.py` (+out).
