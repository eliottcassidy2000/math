---
id: THM-1620
title: "THE POCHHAMMER BRIDGE — the toral layer (TNC) and the radial layer (NC2) are the SAME orthogonal-polynomial fact one level apart, and the level is exactly one rising factorial. (0) THE TORAL LAYER IS LEGENDRE: with Λ = u^{−M}g(u) and M = N = 1 (mac-mini THM-1610's restatement CT(Λᵐ) = [u^{Mm}]Rᵐ), CT(Λᵐ) = [uᵐ](g₀+g₁u+g₂u²)ᵐ = Σ_k [(m)_{2k}/(k!)²]·wᵏ·b^{m−2k} = D^{m/2}·P_m(g₁/√D) with w = g₀g₂, b = g₁, D = g₁²−4g₀g₂ — a LEGENDRE polynomial. Verified exactly on six g at m ≤ 10. (I) THE RADIAL LAYER IS HERMITE: m·E_r[ψ_m] = Σ_k [(m)_{2k}/k!]·wᵏ·b^{m−2k} = sᵐ·He_m(b/s), s = √(−2w) (THM-1660). (II) THE BRIDGE IS ONE RISING FACTORIAL, TERM BY TERM: radial_k = toral_k · k! for EVERY k — verified — and that k! is precisely E_r[rᵏ] = (1)_k, the Gamma(1) moment. So the radial average is literally the rising-factorial moment functional applied to the falling-factorial toral coefficients (m)_{2k} = m!/(m−2k)!. Legendre → Hermite is the Askey-scheme descent, and the Gamma average is what performs it. (III) BOTH CLOSE BY ONE LEMMA: Legendre (m+1)P_{m+1} = (2m+1)xP_m − mP_{m−1} and Hermite He_{m+1} = xHe_m − mHe_{m−1} are both three-term, so in both a common root of consecutive members descends to p₀ = 1 ≠ 0. TNC at M = N = 1 and NC2 therefore close by the SAME argument, with no domination estimate anywhere. (IV) THE GENERAL MECHANISM IS FAVARD, and the hypothesis is a rising factorial: the radial moments are μ_j = j! = (1)_j, whose Hankel matrix is positive definite (leading minors 1, 1, 4, 144, 82944, 1194393600, verified), so Favard yields monic orthogonal p_m with b_m > 0 — and b_m ≠ 0 is EXACTLY what powers the descent. Every Favard family closes; Legendre, Hermite, Laguerre, Gegenbauer alike. (V) FORMALIZED IN THAT GENERALITY: `ThreeTerm.no_common_root` for an arbitrary monic three-term recurrence with b n ≠ 0, sorry-free, 17 theorems in the module. (VI) THIS REPAIRS THE LINK boxeph-S175 STILL ROUTES THROUGH: their TNC proof is independent and may stand, but their stated chain 'TNC ⟹ NC2 (klein's Gamma bridge: factorial moment weights make the growing-degree top term dominate the radial average) ⟹ GMC(2)' uses the step THM-1585 refuted. The Pochhammer descent replaces it without any estimate"
status: >
  (0) VERIFIED-EXACT, six g, m ≤ 10, symbolic.  The Legendre identity for central
  trinomial coefficients is classical; its use as the TORAL layer of TNC is the point here.
  (I) From THM-1660; re-verified here on four (w,b).
  (II) VERIFIED-EXACT term by term.  This is an identity of binomial coefficients and is
  elementary once written down; it is stated because of what it means, not because it is
  hard.
  (III) PROVED (both recurrences are classical; the descent is the argument of THM-1660).
  (IV) Favard's theorem is CITED, not reproved.  Positive-definiteness of the Hankel
  matrix of j! is VERIFIED to n = 6 and is classical (det = product of (k!)^2).
  (V) MACHINE-CHECKED, Lean 4.30.0 / Mathlib v4.30.0, sorry-free, no native_decide.
  (VI) Commentary on a live dispute; see the court case.
  SCOPE, and it is the same wall as before: (0) is M = N = 1.  For general (M,N) the
  toral coefficient [u^{Mm}]g^m is still holonomic in m (a diagonal of a rational
  function), but it is NOT claimed here to be an orthogonal family, and the
  no-common-root argument does NOT automatically apply to a higher-order recurrence.
  That is the honest frontier.  GMC(2) is NOT claimed.
source: kind-pasteur-2026-07-20-S128c121 (owner: how do rising and falling factorials help close TNC and thus GMC(2))
depends_on:
  - THM-1660    # the radial layer is Hermite; the no-common-root descent
  - THM-1585    # why an estimate-free mechanism was needed
related: [THM-1610, THM-1600, THM-1550, THM-1530, THM-1540]
court: 02-court/active/CASE-gamma-bridge-domination-step.md
lean: 04-computation/lean/TournamentH7/TournamentH7/GMC2HermiteNoCommonRoot.lean
script: 04-computation/pochhammer_bridge_tnc_kps_S128c121.py (+ .out)
---

# THM-1620 — the Pochhammer bridge

The owner's question was whether rising and falling factorials could close TNC and hence
GMC(2). They are not a tool applied from outside — **they are what the two layers are made
of, and the gap between the layers is exactly one of them.**

## 0–I. Both layers, written out

With `Λ = u^{−M}g(u)` and `M = N = 1`, `g = g₀ + g₁u + g₂u²`, `w = g₀g₂`, `b = g₁`:

| layer | coefficient | closed form | family |
|---|---|---|---|
| **toral** (TNC) `CT(Λᵐ) = [uᵐ]gᵐ` | `(m)_{2k}/(k!)²` | `D^{m/2}·P_m(g₁/√D)`, `D = g₁²−4w` | **Legendre** |
| **radial** (NC2) `m·E_r[ψ_m]` | `(m)_{2k}/k!` | `sᵐ·He_m(b/s)`, `s = √(−2w)` | **Hermite** |

Here `(m)_{2k} = m!/(m−2k)!` is the **falling** factorial. Both verified exactly.

## II. The bridge is one rising factorial

> **`radial_k = toral_k · k!` for every `k`** — and `k! = (1)_k = E_r[rᵏ]`.

So the radial average is *literally* the rising-factorial moment functional applied to the
falling-factorial toral coefficients. Nothing analytic happens in the passage; a Pochhammer
symbol is inserted.

And the passage `Legendre → Hermite` is the classical **Askey-scheme descent**. The Gaussian
radial average is what performs it. That is the whole relationship between TNC and NC2.

## III. One lemma closes both

- Legendre: `(m+1)P_{m+1} = (2m+1)x·P_m − m·P_{m−1}`
- Hermite: `He_{m+1} = x·He_m − m·He_{m−1}`

Both three-term; in both, a common root of consecutive members forces the previous one to
vanish, and the descent terminates at `p₀ = 1`. So neither `CT(Λᵐ)` nor `E_r[ψ_m]` can
vanish for all `m` unless its prefactor degenerates — which is the one-sided locus.

**TNC at `M=N=1` and NC2 are the same theorem, and neither needs an estimate.**

## IV. The general mechanism, and the hypothesis is again a rising factorial

The radial moments are `μ_j = j! = (1)_j`. Their Hankel matrix has leading minors
`1, 1, 4, 144, 82944, 1194393600, …` — positive definite (classically, `det = ∏(k!)²`). So
**Favard** applies: there are monic orthogonal `p_m` with

> `p_{m+1} = (x − a_m)p_m − b_m·p_{m−1}`,  `b_m > 0`,

and `b_m ≠ 0` is *exactly* the hypothesis the descent needs. Every Favard family closes.

This is why the corpus kept meeting `(1/2)_{n−2}/(n−2)!` and `(1−x)^{−1/2}` (the fibre
fraction, CLAUDE.md) and `−½` in `He_m` and `A = (J−I+S)/2` (THM-1555): the Legendre
generating function is `(1−2xt+t²)^{−1/2}`. **The `½` in the half-dictionary and the `½` in
the Legendre exponent are the same `½`.**

## V. Formalized in the general form

`ThreeTerm.no_common_root` — an arbitrary monic three-term recurrence with `b n ≠ 0` —
subsumes Hermite, Legendre, Laguerre and Gegenbauer in one theorem. Module now at
**17 theorems, sorry-free, no `native_decide`**, clean under Mathlib v4.30.0.

## VI. What this does to the live claims

boxeph-S175 proves TNC (monodromy transitivity + Puiseux-DFT) and then writes: *"TNC (this)
⟹ NC2 (klein's Gamma bridge, S351: factorial moment weights make the growing-degree top
term dominate the radial average) ⟹ GMC(2)"*. **That parenthesis is the step THM-1585
refuted** — measured to `m = 20` the top term's share falls to `0.04%`. Their TNC proof is
independent of it and may well stand; the consequence does not follow as stated.

The Pochhammer descent above is a replacement for that link **at `M = N = 1`**, with no
estimate. It does not yet deliver GMC(2), and I do not claim it.

> **Named next, and this is now a sharp question rather than a vague one.** For general
> `(M,N)`, `[u^{Mm}]gᵐ` is a diagonal of a rational function, hence holonomic in `m` — it
> satisfies a linear recurrence with polynomial coefficients, but of order `> 2` in
> general. The descent argument needs order exactly 2. **So the frontier is: for which
> `(M,N)` is the toral coefficient sequence order-2 (i.e. an orthogonal family), and when
> it is not, what replaces "no common root" for a higher-order holonomic recurrence?** A
> `d`-th order recurrence with nonvanishing trailing coefficient gives the same descent
> from `d` consecutive zeros, so the real question is whether the trailing coefficient can
> vanish — which is a resultant/apparent-singularity computation, not an estimate. That is
> where I would put the next session.
