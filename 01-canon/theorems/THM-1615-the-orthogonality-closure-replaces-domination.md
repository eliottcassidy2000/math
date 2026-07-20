---
id: THM-1615
title: "THE ORTHOGONALITY CLOSURE — the GMC(2) radial layer needs no domination estimate at all, because the Gamma average turns the nullcone condition into a COMMON-ROOT question for a classical polynomial sequence, and no such sequence has a common root. (0) THE PRINCIPLE: if the radial average of the m-th moment equals (nonzero factor)·p_m(ζ) for a FIXED point ζ and a polynomial sequence (p_m) with p_0 = 1 obeying a recurrence that expresses p_{m−1} through its neighbours, then 'all moments vanish' would force ζ to be a common root of every p_m — impossible by descent to p_0 = 1. The nullcone therefore consists exactly of the elements where the prefactor degenerates, which is the ONE-SIDED locus. No asymptotics, no ℓ¹-mass comparison, no ESV saddle, no Eulerian-numbers bound. (I) DEGREE-2 INSTANCE, mine, on the {−1,0,1} stratum at M = 1: Lagrange–Bürmann on u = tφ(u) with H = log(φ/φ(0)) collapses ψ_m = (1/m)[uᵐ]φ(u)ᵐ, φ = ρa + bu + ρcu², which forces #(ρa) = #(ρc) so every term carries ρ^{2k} = r^k — ρ-FREE — giving ψ_m = (1/m)Σ_k m!/(k!²(m−2k)!)·(rac)^k·b^{m−2k}; then E_r[r^k] = k! cancels one k! and m·E_r[ψ_m] = s^m·He_m(b/s) with s = √(−2ac), He = PROBABILISTS' HERMITE. Verified exactly against the log-series at m ≤ 16 on six (a,b,c). (II) DEGREE-1 INSTANCE, mac-mini-S140's THM-1600, pushed the same day and reached independently: L((av+b)^m) = m!·aᵐ·e_m(b/a), e_m the TRUNCATED EXPONENTIAL. Same shape, different family. (III) BOTH CLOSE BY THE SAME TWO-LINE ARGUMENT: Hermite via He_{m+1} = x·He_m − m·He_{m−1}, so a common root of consecutive members descends to He_0 = 1 ≠ 0; truncated exponential via e_{n+1} − e_n = z^{n+1}/(n+1)!, so a common root forces z = 0 where e_n(0) = 1 ≠ 0. (IV) FORMALIZED: `GMC2HermiteNoCommonRoot.lean`, 12 theorems, sorry-free, no native_decide, wired into the root module and building clean under Mathlib v4.30.0. (V) CONSEQUENCE: on the CONSTANT-coefficient {−1,0,1} stratum at M = 1 the nullcone is exactly {ac = 0} = the one-sided locus — the one-sided conjecture, PROVED there, by orthogonality rather than by estimate. (VI) AND IT EXPLAINS THE DISPUTE OF THM-1585: domination was an analytic strategy for a fact whose real content is algebraic, which is why the top term's share could fall to 0.04% while the conclusion stayed true"
status: >
  (0) The principle is a statement pattern, not a theorem; it is discharged case by case
  and both discharges are below.
  (I) The Lagrange-Bürmann closed form is PROVED (the computation is three lines and is
  reproduced in the file) and VERIFIED-EXACT against the independent log-series at
  m <= 16 on six (a,b,c).  The identity m*E_r[psi_m] = T_m(b,ac) and T_m(b,-1/2) = He_m(b)
  are VERIFIED-EXACT in rational arithmetic; they are elementary binomial identities and
  are NOT formalized.
  (II) mac-mini-S140's, cited not reproved; the reading of it as the same shape is mine.
  (III) PROVED, and (IV) MACHINE-CHECKED in Lean 4 / Mathlib v4.30.0, sorry-free.
  (V) PROVED for CONSTANT a, b, c.  This is the honest scope and it is narrow: the step
  E_r[W^k B^{m-2k}] = k!*w^k*b^{m-2k} needs constants, so non-constant a(r), b(r), c(r)
  are NOT covered by the Hermite closure and remain open here.  death-star-S61h claims
  the non-constant case via THM-1515 by a different (domination) route; that claim is not
  adjudicated here and is not relied on.
  (VI) Commentary.
  What is NOT claimed: GMC(2).  It is open.  This closes one stratum and supplies a
  mechanism; it does not compose to the general case.
source: kind-pasteur-2026-07-20-S128c120 (owner: finish the GMC(2) math, then formalize; extend incoming ideas)
renumbered: "claimed as THM-1605; renumbered to THM-1615 by first-pusher rule -- boxeph-S175 and opus-S415 both pushed THM-1605 earlier the same day. Three-way collision, mine was last."
depends_on:
  - THM-1585    # the domination step is false -- why a new mechanism was needed
related: [THM-1600, THM-1540, THM-1530, THM-1550, THM-1515, THM-1580, THM-1590]
court: 02-court/active/CASE-gamma-bridge-domination-step.md
lean: 04-computation/lean/TournamentH7/TournamentH7/GMC2HermiteNoCommonRoot.lean
script: 04-computation/hermite_closure_gmc2_kps_S128c120.py (+ .out)
---

# THM-1615 — the orthogonality closure

THM-1585 showed klein-S351's Gamma-domination step is false: the top term's share of
`E_r[ψ_m]` falls to `0.04%`, so the mass sits at an **interior** index. The instinct is to
reach for a sharper estimate — ESV/Eulerian asymptotics, `ℓ¹` bounds, saddle analysis.

**That is the wrong move. The fact is not analytic.**

## 0. The principle

> If the Gamma/Laplace radial average of the `m`-th moment equals
> `(nonzero prefactor) · p_m(ζ)` for a **fixed** `ζ` and a polynomial sequence `(p_m)`
> with `p_0 = 1` obeying a recurrence expressing `p_{m−1}` through its neighbours, then
> `E[P^m] = 0 ∀m` would make `ζ` a common root of **every** `p_m` — impossible by descent
> to `p_0 = 1`. So the nullcone is exactly where the **prefactor** degenerates.

The radial average does not need to be *dominated*. It needs to be *recognised*.

## I. Degree 2 — Hermite

On the `{−1,0,1}` stratum at `M = 1`, Lagrange–Bürmann applied to `u = tφ(u)` with
`H = log(φ/φ(0))` gives `H'φ^m = φ'φ^{m−1} = (φ^m)'/m`, hence

> `ψ_m = (1/m)·[uᵐ] φ(u)ᵐ`,  `φ(u) = ρa + bu + ρcu²`.

Extracting `[uᵐ]` forces `#(ρa) = #(ρc) = k`, so every surviving term carries `ρ^{2k} = r^k`
— **the `ρ`-freeness is not an accident, it is the `u`-degree bookkeeping** — and

> `ψ_m = (1/m) Σ_k [m!/(k!²(m−2k)!)] · (r·a·c)^k · b^{m−2k}`.

Applying `E_r[r^k] = k!` cancels one `k!`:

> `m·E_r[ψ_m] = Σ_k [m!/(k!(m−2k)!)]·w^k·b^{m−2k} = s^m·He_m(b/s)`,  `w = ac`, `s = √(−2w)`,

since `He_m(x) = Σ_k m!/(k!(m−2k)!)·(−½)^k·x^{m−2k}` exactly. **The `−½` in the Hermite
coefficient is the same `½` as everywhere else in this corpus** (THM-1555).

Verified against the independent log-series at `m ≤ 16` on six `(a,b,c)`; and
`T_m(b,−½) = He_m(b)` checked exactly on four `b`.

## II. Degree 1 — truncated exponential (mac-mini-S140, THM-1600)

Pushed the same day, reached independently, and **the same shape**:

> `L((av + b)^m) = m!·aᵐ·e_m(b/a)`,  `e_m(z) = Σ_{k≤m} z^k/k!`.

Two different classical families, one phenomenon: *the Gamma average of an `m`-th power is
a classical polynomial sequence evaluated at one fixed point.*

## III. Both close in two lines

- **Hermite.** `He_{m+1} = x·He_m − m·He_{m−1}`. A common root of `He_m, He_{m+1}` forces
  `m·He_{m−1}(x) = 0`, so (characteristic zero) `He_{m−1}(x) = 0`; descend to `He_0 = 1`. ∎
- **Truncated exponential.** `e_{n+1}(z) − e_n(z) = z^{n+1}/(n+1)!`. A common root forces
  `z = 0`, but `e_n(0) = 1`. ∎

Neither uses an inequality.

## IV. Formalized

`04-computation/lean/TournamentH7/TournamentH7/GMC2HermiteNoCommonRoot.lean` —
**12 theorems, sorry-free, no `native_decide`**, building clean under Lean 4.30.0 /
Mathlib v4.30.0 and wired into the root module. Contents: `He` and `trExp` by recurrence,
`no_common_root`, `trExp_no_common_root`, and the corollaries `exists_nonvanishing`,
`trExp_exists_nonvanishing`, `nonvanishing_of_consecutive`.

**What is formalized is exactly the step that replaced the false one.** The
Lagrange–Bürmann collapse and the binomial identity `T_m(b,−½) = He_m(b)` are verified
exactly in rational arithmetic but are *not* in Lean; so the chain is machine-checked at
its new joint, not end to end. Stated so no one cites this as "GMC(2) in Lean".

## V. Consequence, with its scope

On the **constant-coefficient** `{−1,0,1}` stratum at `M = 1`:

> `P` is in the nullcone `⟺ ac = 0 ⟺ P` is one-sided.

*Proof.* If `ac ≠ 0` then `s ≠ 0` and `m·E_r[ψ_m] = s^m·He_m(b/s)`; by III some `m` has
`He_m(b/s) ≠ 0`, so `E_r[ψ_m] ≠ 0` and `P ∉ N`. (For `ac > 0`, `s` is imaginary and `b/s`
purely imaginary, while every root of `He_m` is real — so it fails at once, and at `b = 0`
already for even `m`.) Conversely `ac = 0` kills a extreme charge, so `P` is one-sided and
lies in `N`. ∎ Confirmed on all 252 constant `(a,b,c)` with `ac ≠ 0` in `[−3,3]³`: the
first nonvanishing `m` is always `≤ 2`.

**Constants only.** `E_r[W^k B^{m−2k}] = k!·w^k·b^{m−2k}` needs them. Non-constant
`a(r), b(r), c(r)` are open here — death-star-S61h claims them via THM-1515 by the
domination route, which THM-1585 disputes; that is not adjudicated here and nothing above
relies on it.

## VI. Why the dispute happened

klein and death-star were trying to prove a **nonvanishing** by showing one term outweighs
the rest. That is an analytic strategy for an algebraic fact, and it is why the top term's
share could collapse to `0.04%` (THM-1585) while the conclusion stayed true: the
conclusion never depended on the top term. Recognising the sum as `He_m(b/s)` makes the
question "is `b/s` a root of every Hermite polynomial", and orthogonal polynomial families
are precisely the objects that cannot do that.

> **Named next.** Push the recognition past constants. The obstruction is that
> `E_r[W(r)^k B(r)^{m−2k}]` is no longer a product; but it is still a *moment functional*
> applied to a product of powers, so the natural target is a **two-variable Appell/Sheffer
> family** with `ζ` replaced by a curve. If the degree-1 and degree-2 layers are both
> Sheffer, the general stratum should be too, and that — not a sharper estimate — is where
> the remaining GMC(2) content lives.
