---
id: THM-1610
title: "THE TNC INDUCTION EXISTS, IT IS THE COEFFICIENT LADDER, AND ITS FAILURE LOCUS IS EXACTLY {j : M ∤ j}. (A) RESTATEMENT: since Λ = u^{−M}R(u), CT(Λᵐ) = [u^{Mm}]R(u)ᵐ — the whole infinite family is a coefficient condition on powers of ONE polynomial, purely combinatorial and fully computable. (B) THE DEGENERATE CASE IS EXACTLY deg R < M: then Rᵐ has degree m·deg R < Mm, so the coefficient vanishes automatically. Verified — deg R < M gives identical vanishing, deg R ≥ M never does. That is what 'degenerate' means in this language, and it is why the hypothesis deg R = d = M+N (i.e. r_d ≠ 0) is load-bearing. (C) THE LADDER: normalise R(0)=1; then m=1 gives r_M = 0 outright, and each successive m peels one further coefficient LINEARLY — but only at multiples of M. For M=1 the peel is COMPLETE, forcing every r_j and proving TNC(1,N) outright, matching THM-1530. For M=2 it peels r₂ then r₄ and BREAKS at r₃ (degree 2); for M=3 it peels r₃ then r₅ and breaks at r₄. So the triangular induction handles exactly the r_j with M | j and stalls exactly at M ∤ j — INDEPENDENTLY REPRODUCING THM-1550's Σ_i ζ^{(k+1)i} arithmetic from a completely different route (coefficient extraction rather than Wiener–Hopf). (D) A CORRECTION TO MY OWN TEST: an initial brute search reported nonzero hits at M ≥ 2, but every one had r_d = 0, i.e. deg R < M — the automatic degenerate case, not counterexamples. With deg R = d enforced there are NO hits on the box for any of the ten bidegrees tested, including (2,4) and (3,3). (E) SEPARATELY, for the Laplace layer feeding GMC(2): contour rotation puts F(t) = ∫e^{tp(v)−v}dv in a sector of opening π(1+D), which exceeds Watson's threshold π for EVERY D ≥ 1 — so the missing input to HYP-8350 is now a named classical Gevrey-1 bound rather than an ad hoc estimate."
status: >
  (A) PROVED — one line from Λ = u^{−M}R(u).
  (B) PROVED (degree count) and VERIFIED for M = 2,3,4 across deg R = 1..M+1.
  (C) VERIFIED-EXACT by symbolic peeling at eight bidegrees. The M=1 completeness is a
  genuine proof of TNC(1,N) by this route. The M ≥ 2 breakdown is a VERIFIED OBSTRUCTION,
  not a proof of anything — it says this induction cannot close M ≥ 2 without a new idea.
  (D) VERIFIED: no nondegenerate solution with coefficients in [−4,4], m = 1..9, at
  (M,N) ∈ {(1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3)}. A bounded-box
  negative, not a proof.
  (E) The sector computation is elementary and done. THE GEVREY-1 BOUND IS NOT VERIFIED —
  it remains exactly the missing estimate, now correctly named.
source: mac-mini-2026-07-20-S141 (owner: "work induction to prove GMC(2) and TNC")
depends_on:
  - THM-1550  # the exact Wiener-Hopf criterion whose arithmetic this independently reproduces
  - THM-1530  # the M=1 theorem, re-proved here by the ladder
related:
  - THM-1595  # boxeph's ladder closures at (2,3), (2,4), (3,3)
  - THM-1600  # the Laplace layer at degree 1; HYP-8350
script: 04-computation/tnc_induction_macmini_S141.py (+ .out),
        05-knowledge/results/tnc_nondegenerate_macmini_S141.out
---

# THM-1610 — the TNC coefficient ladder

## (A) The restatement that makes induction possible

Since `Λ = u^{−M}R(u)`,

> **`CT(Λᵐ) = CT(u^{−Mm}R(u)ᵐ) = [u^{Mm}]R(u)ᵐ`.**

The whole infinite family is a **coefficient condition on powers of one polynomial** — purely
combinatorial, no Wiener–Hopf, no Puiseux branches, and completely computable. Normalising
`R(0) = 1` (scaling `R` rescales `t`), TNC(M,N) reads:

`[u^{Mm}]R(u)ᵐ = 0` for every `m ≥ 1`, `deg R = d = M+N` ⟹ `R` degenerate.

## (B) The degenerate case is exactly `deg R < M`

If `deg R < M` then `deg Rᵐ = m·deg R < Mm`, so `[u^{Mm}]Rᵐ = 0` **automatically**, for every
`m`, whatever the coefficients. Verified:

| `M` | `deg R` | vanishes for all `m`? |
|---|---|---|
| 2 | 1 | **yes** (`deg R < M`) |
| 2 | 2, 3 | no |
| 3 | 1, 2 | **yes** |
| 3 | 3, 4 | no |
| 4 | 1, 2, 3 | **yes** |
| 4 | 4, 5 | no |

So "degenerate" here means precisely `deg R < M`, and the hypothesis `deg R = d` (i.e.
`r_d ≠ 0`) is load-bearing — see (D).

## (C) The ladder, and exactly where it breaks

`m = 1` gives `[u^M]R = r_M = 0` outright. Each successive `m` can peel one more coefficient
**if the equation is linear in its highest unknown**. Peeling symbolically:

| `(M,N)` | `d` | peeled | ladder |
|---|---|---|---|
| (1,1) | 2 | **2/2** | `r₁`, `r₂` — then vacuous |
| (1,2) | 3 | **3/3** | `r₁, r₂, r₃` — then vacuous |
| (1,3) | 4 | **4/4** | `r₁..r₄` — then vacuous |
| (2,1) | 3 | 2/3 | `r₂`, `r₃` |
| (2,2) | 4 | 2/4 | `r₂`, `r₄`, then **breaks at `r₃` (degree 2)** |
| (2,3) | 5 | 3/5 | `r₂, r₄, r₅`, then **breaks at `r₃` (degree 3)** |
| (3,1) | 4 | 2/4 | `r₃`, `r₄` |
| (3,2) | 5 | 2/5 | `r₃`, `r₅`, then **breaks at `r₄` (degree 2)** |

> **For `M = 1` the peel is complete — every `r_j` is forced, so TNC(1,N) is proved outright
> by the ladder**, recovering THM-1530.
> **For `M ≥ 2` the ladder peels exactly the `r_j` with `M | j` and stalls exactly at `M ∤ j`.**

This **independently reproduces THM-1550's arithmetic** — its criterion carries
`Σ_i ζ^{(k+1)i}`, nonzero exactly when `M | (k+1)` — but from a completely different route
(coefficient extraction rather than an exact Wiener–Hopf factorisation). Two derivations, one
divisibility condition.

**The honest reading: the induction the owner asked for exists, is canonical, and provably
cannot close `M ≥ 2` on its own.** The obstruction is not a gap in the bookkeeping; it is the
non-multiples of `M`, and it is the same obstruction THM-1550 identified.

## (D) A correction to my own test

An initial brute search reported nonzero solutions at `M ≥ 2` — e.g. `R = 1 + r₁u` at
`(M,N) = (2,1)` for every `r₁`. **Those are not counterexamples.** Every one had `r_d = 0`,
so `deg R < M`, and (B) makes their vanishing automatic. The search had failed to enforce
`deg R = d`. With that enforced:

> **No nondegenerate solution exists with coefficients in `[−4,4]`, `m = 1..9`, at any of
> `(1,1), (1,2), (1,3), (2,1), (2,2), (2,3), (2,4), (3,1), (3,2), (3,3)`.**

A bounded-box negative, consistent with boxeph's THM-1595 closures at `(2,3), (2,4), (3,3)`.

## (E) The Laplace layer: the missing estimate now has a name

For HYP-8350 (`L(pᵐ) = 0 ∀m ⟹ p = 0`), the exponential generating function is
`F(t) = ∫₀^∞ e^{t·p(v)}e^{−v}dv`. On the rotated ray `v = ρe^{iφ}`, `|φ| < π/2` (which keeps
`e^{−v}` decaying and leaves `L(v^k) = k!` unchanged, the integrand being entire), the integral
converges when `Re(t·a_D e^{iDφ}) < 0`. Sweeping `φ` rotates that half-plane through `Dπ`, so
`F` continues analytically to a sector of opening

> **`π(1+D)` — which exceeds Watson's threshold `π` for every `D ≥ 1`.**

If all `L(pᵐ) = 0` then `F`'s asymptotic series is identically `1`, and by Watson–Nevanlinna a
function analytic in a sector of opening `> π` with **Gevrey-1** asymptotics is Borel-determined
by that series — giving `F ≡ 1`, hence the pushforward of `e^{−v}dv` under `p` is `δ₀`, hence
`p ≡ 0`.

**The Gevrey-1 bound is not verified here.** Numerically `|L(pᵐ)|/(Dm)!` stays `O(1)`, which is
Gevrey-`D` in `t` — i.e. Gevrey-1 in `t^{1/D}`, the same variable in which the sector opening
was measured. That is the right normalisation, but it is an observation, not a proof.

## Honest scope

- (C)'s `M ≥ 2` result is a **verified obstruction**, not a theorem about TNC. It says this
  induction stalls; it does not say TNC is false, and it does not close any new bidegree.
- (D) is a **bounded-box** negative (`[−4,4]`, `m ≤ 9`). It adds confidence, not proof, and it
  is weaker than boxeph's THM-1595 elimination at the bidegrees they cover.
- (E) **reduces** HYP-8350 to a named classical criterion; it does not establish the criterion's
  hypothesis. Anyone citing this should cite it as a reduction.
- **GMC(2) itself is not advanced here.** The Laplace layer is one input to it; the other
  inputs (THM-1600's span-2 elimination, spans ≥ 3) are untouched by this file.

*Artifacts:* `04-computation/tnc_induction_macmini_S141.py`,
`05-knowledge/results/tnc_nondegenerate_macmini_S141.out`.
