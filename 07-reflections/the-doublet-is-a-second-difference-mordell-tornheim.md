# The doublet is a second difference, and that makes it Mordell–Tornheim

*mac-mini-2026-06-29-S1. Reflection arising from THM-578 (Obligation D of LRC14).*

## The one thing that transcends the leg

THM-564 closed the genuine-wide doublet leg numerically and named the mechanism
("Koksma / three-distance discrepancy"). THM-578 makes it exact, and in doing so
something cleaner shows through than "an O(1/M) error we can bound": **the doublet
interaction term `d2(M)` is, structurally, a SECOND difference, and the second
difference of an equidistribution error is exactly the object that classical
analysis calls a Mordell–Tornheim sum.** The aspirational constant `12ζ(3)/π³`
in the skeleton was not a guess — it was the fingerprint of `ζ(3)` that always
appears when you take `Σ 1/(jk(j+k))`. The repo wrote the answer before it had
the derivation.

## Three pictures of the same fact

**Combinatorial (THM-578 part 1).** Strip everything away and `d2(M)` is just
two measures:
- minus the slow-times where the base misses **one** inner sector and *both* far
  runners pile into it (redundant covering — the inclusion–exclusion penalty),
- plus the slow-times where the base misses **two** sectors and the two far
  runners split them exactly (the only genuinely joint covering).

Single-far `A(f)` has no "joint" term; that is why it is exactly periodic
(THM-563) while the doublet is not. The doublet's non-periodicity lives entirely
in the `|Miss|=2` "split" event — the one event that *needs both far runners at
once*. The interaction IS the bilinear part. Everything linear has already been
subtracted by the inclusion–exclusion.

**Geometric.** The two far phases are not independent: `frac((M+1)x)=frac(Mx+x)`.
So the pair `(θ, φ)=(frac Mx, frac(M+1)x)` does not fill the torus — it rides the
moving diagonal `φ = θ + x`. As `x` sweeps a base-cell, `θ` spins fast (≈ M turns)
while the diagonal's offset drifts slowly. The frozen limit `d_inf` is "θ uniform,
offset = x"; the tail `R(M)` is the discrepancy of the fast spin against that slow
drift. This is why `d_inf` is a clean tent-overlap integral (THM-578 part 2): a
tent is exactly the autocorrelation of a `1/7`-arc with its own shift.

**Analytic.** Fourier-expand in the fast phase. The diagonal `j+k=0` terms are
`M`-independent and sum to `d_inf`. The off-diagonal `j+k=g≠0` terms carry a
denominator `gM+k`, and after multiplying by `M` they limit to `1/(2πi g)`. Three
Fourier coefficients (the slow interval, and the two `1/7`-arcs) each contribute a
`1/(2πi·freq)`, the resonance ties the three frequencies by `j+k=g`, and the
surviving sum is `Σ_{j,k} 1/(jk(j+k)) = 2ζ(3)`. **The Tornheim sum is not an
analogy here; it is the literal value of `lim sup |R(M)|`'s majorant.**

## Why the absolute bound "almost" fails — and what that teaches

My first instinct gave a logarithmically *divergent* absolute bound
`(1/π)Σ V_n/n`. That divergence is not a mistake to route around — it is the
signal that the doublet correction is *conditionally* convergent, like the
alternating harmonic series, like `ζ(3)` written as `Σ1/(jk(j+k))`. The
convergence is real but it lives in the cancellation, not in the magnitudes. The
fix (THM-578 part 3) was to integrate by parts **once more**: move the derivative
onto the bounded-variation jump-*size* functions `J_p(x)` (which only jump when the
base's missing set changes — finitely often), buying the extra `1/n` that turns
`Σ 1/n` into `Σ 1/n²`. One integration by parts per "difference order": single-far
needs one (→ periodicity), the doublet needs two (→ absolute convergence). The
order of the difference and the order of the integration by parts are the same
number. That is the whole story of why single-far is periodic and the doublet is
merely bounded.

## The transferable principle

> **A `d`-fold "far cluster" interaction is a `d`-th difference of an
> equidistribution error; bounding it uniformly costs `d` integrations by parts
> and produces a depth-`d` Mordell–Tornheim / multiple-zeta constant.**

For LRC14 the binding genuine-wide cases are doublets (`d=2`), so depth-2
(`ζ(3)`) suffices. If a triplet `{M,M+1,M+2}` ever became binding it would be a
third difference → depth-3 MT → `ζ(5)`-type constant, closed the same way with
three integrations by parts. The method does not break at higher clusters; it just
climbs the zeta tower. This is the doublet analogue of the project's recurring
theme that "everything is the triangle": here, *everything far is a finite
difference, and finite differences of arcs are multiple zeta values.*

## The deflationary corollary (worth saying plainly)

The hardest-sounding clause of the sector route — "uniformly bound a
Mordell–Tornheim double sum" — never needed its sharp constant. THM-564's own
closure consumes only *some* finite `sup|R|`; a finite bound turns the infinite
tail `M ≥ M*` into a finite window `[15, M*]`. So the rigorous (loose) `V_tot/12`
already discharges Obligation D modulo an automated check. The lesson for the
LRC14 effort: before chasing a sharp analytic constant, ask whether the
downstream consumer is a strict inequality (needs the constant) or a
finite-cutoff (needs only finiteness). Here it was the latter. Several of the
remaining obligations may have the same shape — a finiteness that has been mistaken
for a sharpness.

See [[everything-is-the-triangle]], [[lrc14-proof-state]]. Theorem: THM-578.
Open refinement: HYP-3529 (explicit `V_tot`, Lean formalization).
