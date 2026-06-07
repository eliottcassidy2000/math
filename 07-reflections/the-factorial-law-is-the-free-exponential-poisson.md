# The factorial law is the free exponential-jump Poisson — one Lévy measure under two convolutions

*monad-explorer-2026-06-07 (deep-research, 14th session). Builds on THM-438 ADDENDUM-11
(the diagonal is free moments of the law `κ_n=n!`) and closes its explicit "name the law" lead.*

## The thing that clicked

For two sessions the diagonal `t(k,k)=A088368` and the over-count `A000262` have been "the free
and classical moments of the law with cumulants `κ_n=n!`." That phrasing made `κ_n=n!` sound like
a formal stipulation — a generating function with no inhabitant. It has one. Read the cumulants
as an integral:

```
        κ_n = n! = ∫_0^∞ x^n e^{-x} dx.
```

That is not a coincidence to admire; it is the **definition of a compound Poisson law**. A
cumulant sequence that equals the moment sequence of a finite measure `ν` is *exactly* the
free/classical cumulant signature of the compound Poisson with Lévy measure `ν`. Here
`ν = e^{-x}dx`, the exponential. So:

- **`A000262` = the classical compound Poisson of `e^{-x}dx`** (rate 1, exponential jumps).
- **`A088368` = the free compound Poisson of `e^{-x}dx`** (the diagonal).

Same Lévy measure. The only difference between the two named endpoints is *which convolution* —
classical (all partitions) or free (non-crossing) — you run it through. "Free = classical on the
non-crossing sublattice" stops being a slogan about lattices and becomes one concrete sentence:
**one exponential measure, two convolutions.**

## Why this had to be true (the combinatorics knew first)

The block weight in both sums is `|B|! =` the number of **linear orders** ("lists") of the block.
And `s! = ∫_0^∞ x^s e^{-x}dx`. The exponential measure is the analytic shadow of "a list of `s`
things has `s!` orderings." So `A000262` = "sets of lists" (EGF `e^{x/(1-x)}`) is the
exponential formula applied to lists = classical CP of the exponential; and `A088368` =
"non-crossing sets of lists" is its free shadow. The species was telling us the Lévy measure all
along: **lists are the exponential.**

## The wild end is wild because the exponential is unbounded

Three apparently separate "wildnesses" the project logged are now one fact about `ν=e^{-x}dx`
having unbounded support and factorial moments:

1. **ADD-6's resurgence.** `R(z)=Σ n! z^{n-1}` diverges (Gevrey-1). But as the `R`-transform of a
   *compound Poisson* it is a Borel sum: `R(z)=∫_0^∞ t e^{-t}/(1-zt)dt`, Borel transform
   `ζ/(1-ζ)`, one singularity at `ζ=1`. The divergence is just the factorial growth of `ν`'s
   moments. The Stokes/instanton term `Im R(x+i0)=π e^{-1/x}/x²` (verified) — the exponentially
   small nonperturbative content with action `1/x` — is precisely what the Cauchy inversion turns
   into the free law's spectral density. *Resurgence is not an obstacle here; it is the spectral
   measure being born.*

2. **`m_k ~ e·k!` (MISTAKE-063).** The free CP inherits the Lévy measure's tail: `μ` has an
   exponential right tail `~e·e^{-x}`. The growth constant `e` is the tail constant.

3. **The OEIS-structurelessness of the interior.** The law `μ` is genuine (Hankel `D_0..D_8>0`,
   support `[0,∞)`), determinate (Carleman), but **non-classical**: its Jacobi continued-fraction
   parameters are non-integer rationals (`b_n: 1, 4, 26/5, 1969/280,…`; `λ_n: 2, 5, 224/25,…`)
   with no polynomial pattern. There is no Hermite/Laguerre/Meixner family here. The
   structurelessness ADD-10/11 saw along the triangle is the law's own continued fraction having
   no closed form.

## The asymmetry that remains (and why the path is still anonymous)

The two endpoints are NOT moments-and-cumulants of one law. The diagonal is **moments** of the
free exponential Poisson; the signed-sum Catalan is **free cumulants** of the two-point law
(ADD-4) — a *different* measure. ADD-11 already noted this; ADD-12 sharpens the diagonal half to a
named, genuine measure while leaving the asymmetry intact. So the Möbius transit across the
triangle's interior connects the moments of one law to the cumulants of another; it is no single
law's moment↔cumulant map, and that is *why* nothing in between is catalogued.

## What this points at

Two laws are now fully explicit objects, not formal stipulations: the **free exponential Poisson**
(`K(w)=1/w+∫ te^{-t}/(1-wt)dt`, density with a `1/√x` edge and `e^{1-x}` tail) on the wild end,
and the **two-point/arcsine law** on the tame end. The off-diagonal columns `t(k,m)` — the
anonymous interior — might be a *marked/rate-deformed* free Poisson: `U(x,y)` as the
free-CP family in the marking variable `y`. We now have the explicit Cauchy transform to deform.
The mathematics handed us a measure where it had only handed us a divergent series; the next move
is to deform that measure and see whether the interior columns fall out of its `y`-family.

*Everything-is-the-triangle footnote:* `e` now appears a third way — not via Stirling/Gamma
(Burnside), not only as the free-moment growth rate (ADD-11), but as the **tail constant** of an
honest probability density `~e·e^{-x}` living on the staircase's wild diagonal. The same constant,
three independent doors.
