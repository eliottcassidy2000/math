---
id: THM-1675
title: "THE RADIAL GMC(2) RESIDUAL IS CLOSED FOR REAL p (so the charge-0 layer closes for Hermitian P), AND THE CROSS-SHELL RESIDUAL IS ALL-INTEGER (the half-shell worry is a red herring). PIECE 1 (HYP-8350): THM-1665 left Ψ(t)=∫₀^∞ e^{−v}[1/(1−tp(v))−1]dv ≡ 0 ⟹ p ≡ 0. For REAL p of degree D ≥ 1 this is now PROVED by the discontinuity argument (DvdK's Theorem-2 mechanism transplanted from the constant-term functional to L): Ψ≡0 makes h(t)=∫e^{−v}/(1−tp)dv ≡ 1, analytic on the connected domain ℂ minus the curve {1/p(v):v≥0}, so its boundary-value jump across that curve, (2πi/t)·Σ_{p(v)=1/t} e^{−v}/|p'(v)|, must vanish — but every term is strictly positive, and a nonconstant polynomial is unbounded, so p is constant and L(p)=p=0 forces p≡0. Exhaustive check: 0 nullcone members among 16800 nonzero real p (deg ≤ 4, coeffs [−3,3], m = 1..12). So the charge-0/radial layer of GMC(2) closes whenever the effective radial polynomial is real, i.e. for Hermitian P. PIECE 2 (HYP-8470, cross-shell): E[P^m] = L(CT_w(Λ_s^m)) = Σ_j L(s^j·c_j(m)), where j = (Σ|k_i|)/2 over charge-balanced m-tuples. THE PARITY LEMMA — Σ|k_i| ≡ Σk_i ≡ 0 (mod 2) — makes j ALWAYS an integer, so E[P^m] never sees a half-integer moment and the 'half-shell' concern is void. Span-2 decouples (THM-1600); the genuine coupling is the fixed-m mixing Σ_j L(s^j c_j) with no per-shell separation at fixed m, and the top shell's (m·k_max/2)! growth is the descent handle. Framed precisely, not closed."
status: >
  PIECE 1: PROVED for real p (the jump argument is rigorous — h ≡ 1 on the connected
  domain ℂ∖curve by the identity theorem from the Watson open set, forcing the positive
  jump density to zero) and VERIFIED exhaustively (0/16800). The complex-p case is DvdK's
  actual Theorem 2 (same argument, curve in ℂ), cited not re-proved.
  PIECE 2: the parity lemma is PROVED (one line) and verified for three charge sets; the
  shell decomposition and span-2 decoupling are exact. HYP-8470 is FRAMED, NOT CLOSED.
source: mac-mini-2026-07-20-S146 (owner: "work the 2 GMC(2) residual pieces")
depends_on:
  - THM-1665  # the per-component Watson lemma; the reduction to Psi == 0 => p == 0
  - THM-1600  # span-2 elimination (the decoupled base case of piece 2)
related:
  - THM-1630  # DvdK; the Theorem-2 mechanism transplanted in piece 1
  - THM-1500  # the GMC master theorem; why n=2 needs U correlated with (Z,W)
closes: HYP-8350   # for real p
---

# THM-1675 — the two GMC(2) residuals

## Setup

For 2 real Gaussians write `Z = rw`, `W = Z̄ = rw^{−1}` (`w = e^{iθ}`, `s = r²`); the Gaussian
measure is `(1/2π)e^{−s}ds dθ`. Every `P ∈ ℂ[Z,W]` is a Laurent polynomial in `w` with
`s`-dependent coefficients, `Λ_s(w) = Σ_k w^k s^{|k|/2}λ_k(s)`, and

`E[P^m] = L(CT_w(Λ_s^m))`, `L(f) = ∫₀^∞ f e^{−s}ds`, `L(s^j) = j!`.

THM-1665 left two residual pieces.

## Piece 1 (HYP-8350): the radial Liouville step — CLOSED for real `p`

THM-1665 reduced the charge-0/radial layer to `Ψ(t) = ∫₀^∞ e^{−v}[1/(1−tp(v))−1]dv ≡ 0 ⟹ p ≡ 0`.

> **Theorem (real `p`).** For a real polynomial `p` of degree `D ≥ 1`, `Ψ ≡ 0 ⟹ p ≡ 0`.
>
> *Proof.* `Ψ ≡ 0` means `h(t) := ∫₀^∞ e^{−v}/(1−tp(v))dv ≡ 1` on the Watson open set
> (THM-1665). `h` is analytic on `ℂ` minus the curve `{1/p(v) : v ≥ 0} ⊂ ℝ`, which is
> connected, so by the identity theorem `h ≡ 1` on all of it. The boundary-value jump of `h`
> across the curve at `t = 1/p(v₀)` is
> `h(t+i0) − h(t−i0) = (2πi/t)·Σ_{v:\,p(v)=1/t} e^{−v}/|p'(v)|`,
> and `h ≡ 1` (analytic across) forces this to vanish. But every summand is **strictly
> positive**, so there is no `v` with `p(v) = 1/t` and `p'(v) ≠ 0` for any large `1/t` — and a
> nonconstant polynomial is unbounded on `[0,∞)`. Hence `p` is constant, and `L(p) = p = 0`
> gives `p ≡ 0`. ∎

This is **DvdK's Theorem-2 mechanism transplanted** from the constant-term functional to `L`:
a single-valued resolvent whose jump density is positive must be constant. Verified
exhaustively — **0 nullcone members among 16 800 nonzero real `p`** (deg ≤ 4, coeffs `[−3,3]`,
`m = 1..12`), and the jump density `Σ e^{−v}/|p'(v)|` computed strictly positive on samples.

> **Consequence.** The charge-0/radial layer of GMC(2) closes whenever the effective radial
> polynomial is real — i.e. **for Hermitian `P`** (`P(Z,Z̄)` real), whose charge-0 part `λ_0(s)`
> is a real polynomial.

The **complex-`p`** case is DvdK's actual Theorem 2 (the same jump argument with the curve
`{1/p(v)}` in `ℂ`), cited not re-proved.

## Piece 2 (HYP-8470): cross-shell coupling — all-integer, framed

`CT_w(Λ_s^m) = Σ` over charge-balanced `m`-tuples `(k_1,…,k_m)`, `Σk_i = 0`, of
`s^{(Σ|k_i|)/2}·∏λ_{k_i}(s)`. Hence `E[P^m] = Σ_j L(s^j·c_j(m))`, `j = (Σ|k_i|)/2` the
**charge-degree**.

> **Parity lemma.** `Σ|k_i| ≡ Σk_i ≡ 0 (mod 2)`, so `j` is **always an integer**.

So `E[P^m]` never involves a half-integer moment `Γ(j+½)` — the "half-shell" worry (from the
`s^{|k|/2}` in `Λ_s` itself) is a **red herring**: it cancels in the constant term by parity.
Verified: zero odd-`Σ|k|` balanced tuples for charge sets `{−1,0,1}`, `{−2,…,2}`, `{−3,0,3}`,
`m ≤ 5`.

**Where the coupling genuinely bites.** For **span 2** (`{−1,0,1}`) the decomposition is
`E[P^m] = Σ_j [m!/(j!²(m−2j)!)]·L(s^j(λ_1λ_{−1})^j λ_0^{m−2j})`, which THM-1600 **decouples**
(m=1 ⟹ λ_0 = 0, then λ_1λ_{−1} = 0). The genuine obstruction is that at fixed `m`, `E[P^m]`
is a **single number mixing all charge-degrees `j`** — there is no per-shell separation at
fixed `m`; separation must come from varying `m`. The **top shell** (all charges `±k_max`) has
the fastest-growing moment `(m·k_max/2)!`, dominating for large `m` — the shell-descent handle.

## Honest scope

- **HYP-8350 is CLOSED for real `p`.** That is the deliverable, and it closes the radial
  charge-0 layer of GMC(2) for Hermitian `P`. It does **not** cover complex `p` by a new
  argument — that is DvdK's theorem, cited.
- The jump argument's rigor rests on `h ≡ 1` on a connected domain (the identity theorem from
  the Watson open set of THM-1665) and on `p'(v) ≠ 0` on the relevant fibers (generic; the
  finitely many critical fibers contribute measure zero and do not save a nonconstant `p`).
- **HYP-8470 is FRAMED, NOT CLOSED.** The parity lemma and span-2 decoupling are exact; the
  fixed-`m` mixing is the real open content, and the top-shell descent is a *strategy*, not a
  proof.
- **GMC(2) as a whole remains open.** It now needs: the complex-`p` radial case (DvdK-style),
  and the cross-shell descent for span ≥ 3 / non-constant coefficients (HYP-8430/8470). The
  span-2 constant case (THM-1600) plus the real radial layer (here) are what is settled.

*Artifacts:* `04-computation/gmc2_two_residuals_macmini_S146.py` (+out).
