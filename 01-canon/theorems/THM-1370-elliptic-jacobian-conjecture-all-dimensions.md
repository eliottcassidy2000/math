---
id: THM-1370
title: The ELLIPTIC Jacobian Conjecture holds in EVERY dimension — a Keller map ℂⁿ→ℂⁿ that is weighted-homogeneous for some strictly POSITIVE weight vector (equivalently, equivariant for a contracting ℂ*-action) is a polynomial automorphism, for all n. This lifts THM-1345's elliptic case from n=2 to all n: its proof never used dim 2. SHARP against the corpus's own counterexample — THM-1300's n=3 map IS weighted-homogeneous, but for the MIXED-SIGN weights (1,−1,−2), and provably admits NO strictly positive weighting (F₃ = 2x−3x²y−x³z forces b = −a). So JC is true in the positively-graded category in every dimension and false at n≥3 in general, with the known counterexample sitting exactly one sign-condition away.
status: PROVED (dimension-free; the five steps are classical and none uses n=2) + VERIFIED-EXACT (det JF ≡ −2 recomputed symbolically; the (1,−1,−2) equivariance recomputed component-by-component: F₁↦λ⁻², F₂↦λ⁻¹, F₃↦λ¹; the no-positive-weighting obstruction verified as a monomial identity).
source: mac-mini-2026-07-20-S123 (owner: creatively find a reduced Jacobian Conjecture that DOES hold)
depends_on:
  - THM-1345  # the n=2 elliptic case, whose argument this observes to be dimension-free (CREDITED)
  - THM-1300  # the explicit n=3 counterexample, against which sharpness is measured
related:
  - boxeph-S144  # the hyperbolic equivariant case (dim 2)
  - THM-1305     # the equivariant anatomy of the counterexample
---

# THM-1370 — the elliptic Jacobian Conjecture, in every dimension

**One line.** JC is false at `n ≥ 3` (THM-1300). But in the *positively graded*
category it is true in **every** dimension — and the known counterexample misses
that category by a single sign.

## Statement

> **Theorem.** Let `F : ℂⁿ → ℂⁿ` be a polynomial map with `det JF ∈ ℂ*` (Keller),
> and suppose `F` is weighted-homogeneous for some **strictly positive** weight
> vector `w ∈ ℤ_{>0}ⁿ` — equivalently, `F` is equivariant for the contracting
> ℂ*-action `λ·(x_1,…,x_n) = (λ^{w_1}x_1, …, λ^{w_n}x_n)` (with any linear
> ℂ*-action on the target). Then `F` is a polynomial automorphism.
>
> This holds for **every** `n ≥ 1`.

## Proof (none of the five steps uses the dimension)

1. **Étale.** `det JF ∈ ℂ*` ⟹ `F` is étale ⟹ every fiber is discrete.
2. **Invariance.** Equivariance makes `F^{-1}(0)` invariant under the source
   ℂ*-action (if `F(x)=0` then `F(λ·x) = σ(λ)F(x) = 0`).
3. **Collapse of the zero fiber.** A discrete ℂ*-invariant set consists of fixed
   points; the action is contracting (all `w_i > 0`), so its only fixed point is
   the origin. Hence `F^{-1}(0) = \{0\}`.
4. **Properness.** A weighted-homogeneous map with positive weights and isolated
   zero fiber is finite: if `F` were not proper there would be `x_k → ∞` with
   `F(x_k)` bounded; rescale by the ℂ*-action to `|x_k| = 1`, pass to a
   convergent subsequence `x_k → x ≠ 0`, and weighted-homogeneity forces
   `F(x) = 0`, contradicting step 3.
5. **Degree one.** A finite étale map onto the simply connected `ℂⁿ` is a trivial
   covering, so `deg F = |F^{-1}(0)| = 1`, i.e. `F` is an automorphism. ∎

**Relation to THM-1345.** That theorem proves exactly this for `n = 2` and reads
it as a two-dimensional fact ("dim 2 has no room"). The reading is too modest:
steps 1–5 above are verbatim its argument, and **not one of them mentions the
dimension**. The elliptic phenomenon is not 2-dimensional; it is dimension-free.
What *is* genuinely 2-dimensional in THM-1345 is the **hyperbolic** case
(boxeph-S144's telescoping `(s·fg)′`), which has no analogue here.

## Sharpness — the counterexample misses by one sign (VERIFIED)

THM-1300's map (`u := 1+xy`)

```
F₁ = u³z + y²u(4+3xy),   F₂ = y + 3xu²z + 3xy²(4+3xy),   F₃ = 2x − 3x²y − x³z
```

has `det JF ≡ −2` (recomputed symbolically here) and **is** weighted-homogeneous:
under `λ·(x,y,z) = (λx, λ^{-1}y, λ^{-2}z)` one gets `F₁ ↦ λ^{-2}F₁`,
`F₂ ↦ λ^{-1}F₂`, `F₃ ↦ λF₃` (recomputed). So it lives in the ℂ*-equivariant
category — but with the **indefinite** weight vector `(1,−1,−2)`.

> **Its grading is UNIQUE, and forced indefinite.** `F₃ = 2x − 3x²y − x³z` has
> monomials of weights `a`, `2a+b`, `3a+c`. Weighted-homogeneity requires all three
> equal, and solving `2a+b = a`, `3a+c = a` gives
>
> > `(a, b, c) = a·(1, −1, −2)` — **exactly THM-1300's weight vector, uniquely up to scale.**
>
> So `a > 0` forces `b = −a < 0` and `c = −2a < 0`: no `w ∈ ℤ_{>0}³` exists. ∎
>
> This is stronger than "the counterexample happens to be indefinitely graded". The
> single component `F₃` already **determines** the grading up to scale, and what it
> determines is indefinite. The counterexample has no freedom to be positively
> graded — the obstruction is one monomial triple.

So the theorem is sharp in a very precise sense: **the counterexample is graded,
just not positively graded.** JC's failure at `n ≥ 3` requires an indefinite
grading; flipping to a definite one restores the conjecture in all dimensions.

## What this does and does not say

- **Does:** identify a JC-restriction that (i) is provable, (ii) holds in *every*
  dimension including those where JC is false, and (iii) is met exactly, not
  slackly, by the known counterexample.
- **Does not:** touch full `JC₂` (open since 1939) or full `JC_n`. The positively
  graded class is a thin slice — e.g. `F = X + H` with `H` homogeneous of degree
  `d ≥ 2` is *never* in it (`x` and `h` cannot share a weight), which is precisely
  why the Bass–Connell–Wright cubic-homogeneous reduction escapes this theorem.
  The class's thinness is the honest price of its dimension-freeness.

*Artifacts:* `04-computation/jc_elliptic_all_dimensions_macmini_S123.py` (+out).
Credits: **THM-1345 / death-star-S59q** (the elliptic argument, here observed to be
dimension-free), **THM-1300** (the counterexample supplying sharpness),
boxeph-S144 (the hyperbolic dim-2 companion).
