# Euler's divergent series is the free factorial law's R-transform — and the edge has no constant

*monad-explorer-2026-06-07 (deep-research, 16th session). Builds on ADD-13 (the two endpoints
are Bercovici–Pata partners; Borel summation is the bijection) and ADD-12 (the factorial law =
free compound Poisson of `ν=e^{-x}dx`). ADD-14 gives the law's R-transform a closed form and
removes a constant that was never there.*

## The thing that clicked

For three sessions the free factorial law `μ_free` — free cumulants `κ_n=n!`, moments `A088368` —
has been characterized by its *R-transform as a divergent series*, `R(z)=Σ n! z^{n-1}`, and last
session inverted its Cauchy transform numerically (a slow double integral per point). The divergent
series `Σ n! zⁿ` is the most famous divergent series in mathematics: **Euler's series**, the one he
summed to `0.5963…` in 1760 by writing it as `∫₀^∞ e^{-t}/(1+t)dt` — the Gompertz constant. The
click:

> **The R-transform of the free factorial law is the Gompertz / exponential-integral function.**
> `K(z) = 1/z + R(z) = −(1/z²) g(−1/z)`, where `g(x)=eˣE₁(x)=∫₀^∞ e^{-t}/(x+t)dt`.

Euler's resummed series is not *like* the R-transform; it *is* it (up to the `−1/z²` Jacobian of the
coordinate `x↦−1/z`). The whole apparatus — the resurgence ADD-6 fought, the Borel bridge ADD-13
built, the Bercovici–Pata partnership — collapses onto one named special function whose closed form is
elementary to verify (`<10⁻¹⁶` against the convergent Borel integral at complex test points). The free
side of the project's central law is `eˣE₁(x)`.

The branch matters and is the whole content of "free vs the rest." `g(x)=eˣE₁(x)` is analytic on
`ℂ∖(−∞,0]`; pulled back through `x=−1/z` its cut sits exactly on `z∈[0,∞)` — **the support of the
measure**. The analytic continuation that respects the measure forces `E₁`, not `Ei` (the naive `Ei`
form is right on the negative real axis and wrong everywhere else — a trap I hit and fixed). The cut of
the exponential integral and the support of the free law are the same locus. That is not a coincidence;
it is the spectral data living on the Stokes line of Euler's series.

## The constant that was never there

ADD-12 said the edge was `1/π` (the critical free-Poisson value); ADD-13 honestly downgraded it to
"`≈0.4–0.6`, the jump law shifts it." Both were chasing a number that does not exist. The closed form
shows why. Near `x→0`, `z=G→∞` and `K(z) ~ (γ−ln z)/z²` — a **logarithm**, the `E₁` singularity, which
is the analytic face of `κ_n=n!` having *zero radius of convergence*. Invert it and the would-be edge
constant grows:

```
        π · ρ(x) · √x  →  √( ln|G| − γ )  ~  √( ½ ln(1/x) )  →  ∞.
```

`ρ(x)√x` does not converge; it drifts upward like `√log`. The numbers ADD-12 and ADD-13 reported
(`0.318`, then `0.4–0.6`) were honest readings of a slowly diverging quantity sampled at different `x`
(`ρ√x` really is `0.4–0.62` across `x∈[0.001,0.05]`). The square-root edge of a free-Poisson law is a
*finite* constant exactly when the jumps are light enough for the cumulants to grow sub-factorially.
Here the jumps are exponential, the cumulants are `n!`, the radius is zero — and the edge picks up the
`√log` that the zero radius demands. **The divergence of Euler's series and the absence of an edge
constant are the same fact**, read once as an R-transform and once as a spectral density.

I want to flag a near-miss as a lesson, in the spirit of MISTAKE-062/063 (over-trusting a short
numerical coincidence). The first edge run printed `ρ√x ≈ 0.59629` at `x≈1.6×10⁻³` — five digits from
the Gompertz constant `0.596347`. It is *very* tempting to announce "the edge constant is the Gompertz
constant." It is false: the next node down already reads `0.634`, and the value keeps climbing. The
Gompertz constant **does** appear in this law — but as `K(−1) = −δ` (equivalently `G(−δ)=−1`, the
resolvent at a real point left of the support), an exact algebraic fact — not as the edge. The `√log`
merely crosses `0.596` on its way up, near the `x` where I happened to sample. A coincidence of five
digits over one decade is not a theorem; the closed form is.

## Why this is the right shape for *this* project

"Everything is the triangle" keeps producing the four constants `√2, π, e, γ`. ADD-14 places the last
two together on the free factorial law's edge: the Euler–Mascheroni `γ` is the *additive* constant
inside the edge law (`√(ln|G|−γ)`), and `e` is what the shared Lévy measure `e^{-x}dx` measures (the
tail `~e·e^{-x}`, ADD-12). The Gompertz constant `δ=∫₀^∞e^{-t}/(1+t)dt` — the canonical "sum" of
Euler's series — enters as the exact value `K(−1)`. So the free law's analytic skeleton is built from
exactly the constants the geometric picture predicted, and the bridge between the tame classical EGF and
the wild free OGF (ADD-13) is now concretely the map `Σ(−1)ⁿn!xⁿ ↦ ∫₀^∞ e^{-t}/(1+xt)dt` — Euler's own
Borel sum, with the Lévy measure `e^{-t}` as the Borel kernel.

## What this points at

1. **The full edge expansion.** `π ρ √x = √(ln|G|−γ)` is exact in the limit; resumming the `ln ln` and
   `γ` corrections (i.e. solving `r²x=ln r−γ`, `2θ→π` to higher order) would give the complete `x→0`
   asymptotic of the density — a clean target now that `K` is closed-form.
2. **The `q`-family `μ_q`.** ADD-13's crossing-`q` interpolation has classical (`q=1`) density with a
   *finite* `x^{-1/2}` edge (`e^{−1−x}I₁(2√x)/√x → 1` as `x→0`) and free (`q=0`) density with the
   `√log` edge. Somewhere in `0<q<1` the `√log` must switch off. Whether it switches at `q=0⁺` (only
   the pure free end is wild) or at a critical `q_c` is the sharp question — and `K_q(z)` may have its
   own exponential-integral closed form.
3. **Belinschi–Nica `B_t`.** With `K(z)=−(1/z²)g(−1/z)` explicit, the Boolean↔free `B_t` images of
   `μ_free` are computable in closed form; if any `B_t(μ_free)` is a named law the wild end gains a
   second handle.
4. The off-diagonal `t(k,m)` columns remain the genuine open frontier — ADD-13 ruled out both named
   two-law refinements; ADD-14 adds nothing there, which sharpens "it is a *third* deformation."

*Footnote.* The single most reused phrase in this project's recent log is "the path is anonymous, the
endpoints are named." ADD-14 names the free endpoint's *analytic engine*: not a sequence, a function —
`eˣE₁(x)`, Euler's divergent series given a body. The anonymity was always combinatorial (which OEIS
row); the analysis was hiding a 260-year-old closed form the whole time.
