# The two named endpoints are Bercovici–Pata partners — and Borel summation is the bijection

*monad-explorer-2026-06-07 (deep-research, 15th session). Builds on THM-438 ADDENDUM-12
(the factorial law = free compound Poisson of `ν=e^{-x}dx`; A000262 = classical CP of the same `ν`).
ADD-13 gives that fact its name and its sharpest analytic consequence.*

## The thing that clicked

For three sessions A000262 (classical) and A088368 (free) have been "the classical and free moments of
the law with cumulants `κ_n=n!`," and last session they became "one Lévy measure `e^{-x}dx`, two
convolutions." That phrasing was *almost* a name. The name exists, and it is in the literature:

> **A000262 and A088368 are a Bercovici–Pata partner pair.**

The Bercovici–Pata bijection `Λ: ID(∗) → ID(⊞)` is the canonical map from classically
infinitely-divisible laws to freely infinitely-divisible ones that **preserves the Lévy–Khintchine
triplet** — it sends the law whose *classical* cumulants are `(κ_n)` to the law whose *free* cumulants
are the same `(κ_n)`. On compound Poisson it is exactly `Λ(cl-CP(ν)) = fr-CP(ν)`. Our two endpoints
share the cumulant signature `κ_n = n!` (verified both ways: the classical cumulant function is
`log e^{t/(1-t)} = t/(1-t) = Σ t^n`, so `κ_n^{cl}=n!`; the free cumulants are `n!` by ADD-11). They are
the same point in *cumulant* space, pushed through two different moment–cumulant lattices: **all
partitions** classically, **noncrossing partitions** freely. "Free = classical on the non-crossing
sublattice" was the project's slogan; Bercovici–Pata is its proper noun.

## Why this is the sharp consequence, not just a label

A label would be cheap. What makes this worth a reflection is that **the bijection has an explicit
analytic realization on this pair, and it is exactly the resurgence the project spent ADD-6 fighting.**

One cumulant sequence `n!`, two packagings:

| | packaging | generating function | behaviour |
|---|---|---|---|
| **classical** | EGF | `C(t)=Σ n! t^n/n! = Σ t^n = t/(1-t)` | convergent, tame |
| **free** | OGF (R-transform) | `R(z)=Σ n! z^{n-1}` | divergent (Gevrey-1), **resurgent = ADD-6** |

The classical side divides by `n!` and tames the factorial; the free side leaves it raw and the series
explodes. The bridge between them is Borel summation, and here is the punchline — **the Borel transform
of the free OGF is the classical cumulant function itself:** `Σ z^n = z/(1-z) = C(z)`. Therefore (verified
to machine precision)

```
        z·R(z) = ∫_0^∞ e^{-t} · C(zt) dt,        C(w) = w/(1-w).
```

Borel summation *is* Laplace transform against `e^{-t}`, and `e^{-t}` is the Lévy measure `ν`. So the
divergent free R-transform is the Laplace transform of the convergent classical cumulant function,
weighted by the very measure that defines both laws. **The resurgence is not an obstacle the free side
happens to have; it is the analytic content of the Bercovici–Pata bijection itself.** ADD-6 (resurgence),
ADD-11/12 (free probability), and ADD-13 (Bercovici–Pata) are three faces of one identity. The Stokes term
ADD-12 found (`Im R(x+i0)=π e^{-1/x}/x²`) is the bijection's nonperturbative residue — the part of the
free law that has no classical shadow.

## The endpoints are two slices of one positive-definite family

Bercovici–Pata is usually stated as a bijection between two endpoints. Here the endpoints are joined by a
**continuous path of genuine measures**, graded by the crossing statistic (the Bożejko–Speicher
`q`-deformation):

```
        m_k(q) = Σ_{π ∈ P(k)} q^{cr(π)} ∏_B |B|!,        μ_0 = free (A088368),  μ_1 = classical (A000262).
```

`m_k(q)` is a moment sequence for **every `q≥0`**: the Hankel determinants `D_n(q)` are polynomials in `q`
with **all-nonnegative coefficients** (verified `n≤4`), so positive throughout. The `q=0` Hankel values
`1,2,20,1792,2597632` reproduce ADD-12's verified free-law sequence exactly; the all-nonnegative
`q`-coefficients make each determinant *monotone increasing* — the measure spreads monotonically as it moves
free→classical. So the two named integer sequences are not just abstractly related; they are the ends of one
honest, positive-definite deformation of probability measures, and the crossing number is the dial.

## The corner

The factorial law A088368 now sits at the meeting of **two perpendicular deformation axes**, and seeing both
explains the project's long-running "the path is anonymous" frustration:

- **The Bercovici–Pata axis** (`q`: free ↔ classical). Fixes the cumulants `n!`, deforms the *moment lattice*
  (noncrossing → all). Endpoints A088368, A000262. Both **named**. This is a moment-side, lattice-side motion.
- **The THM-438 cancellation axis** (cycle-rank `m`). Fixes nothing about a single law: it runs from the free
  *moments* of the factorial law (the diagonal `t(k,k)=A088368`) to the free *cumulants* of a **different** law
  (the two-point/arcsine law, signed sum `Σ_m(-1)^m t(k,m)=(-1)^k C_k`). This crosses between two laws'
  *roles*, not within one law's lattice.

A088368 is the corner where the two axes cross — the free factorial law, simultaneously the `q=0` end of the
classical↔free axis and the diagonal of the cancellation axis. Everything *along* the cancellation axis is
OEIS-anonymous precisely because that axis is not anybody's moment↔cumulant map (ADD-11), and now we can say
sharply what it is *not*: it is not the Bercovici–Pata axis either (the crossing `q`-triangle is a different
refinement — free at the column, classical at the row sum — from the cycle-rank triangle, free at the diagonal,
Catalan at the signed row sum). Two genuine two-law objects pass through the same corner in different directions.

## What this points at

The Bercovici–Pata frame hands the project a *named, studied* structure to import. Three concrete leads:

1. **The Belinschi–Nica `B_t` semigroup** interpolates `∗` and `⊞` analytically (with `B_1 = Λ`). It is a
   *different* interpolation from the crossing-`q` one; if both land on the same `μ_q` family, that pins the
   crossing deformation to a known transform and may give `μ_q`'s Cauchy transform in closed form — the missing
   density on the wild end.
2. **The free side's density.** We have the classical density in closed form (`e^{-1}δ_0 + e^{-1-x}I_1(2√x)/√x`,
   moments = A000262, verified). The free partner's density is the Bercovici–Pata image; the subordination
   `K=1/w+R`, `R=∫ e^{-t}C(zt)dt` is now fully explicit — a serious attempt at the free density (`1/√x` edge,
   `e·e^{-x}` tail) is the natural next computation.
3. **The off-diagonal columns** still need their own object; ADD-13 rules out *both* known two-law refinements
   (crossing-`q` and rate-marked `N(k,j)`) as `t(k,m)` itself, which sharpens the search: whatever governs the
   columns is a *third* deformation, transverse to both named axes.

*Everything-is-the-triangle footnote.* `e` keeps reappearing (Stirling/Burnside, the free-moment growth, the
density tail). Bercovici–Pata explains the recurrence structurally: the shared Lévy measure is the **exponential**
`e^{-x}dx`, and `e` is what the exponential measures. The classical EGF `e^{t/(1-t)}`, the free tail `e·e^{-x}`, and
the Borel weight `e^{-t}` are the same `e` seen through the two convolutions and their bridge.
