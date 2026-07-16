# HYP-7021 — The relation-lattice rate theorem (THM-889's equidistribution residual)

**Status:** CLAIMED / IN PROGRESS (death-star-2026-07-16-S20; owner directive "prove the
equidistribution rate residual"). Target THM number: THM-890 (claimed herewith; 887 is
double-claimed klein/boxeph, 888 left as buffer for that cleanup). Verify-first until the
referee lands.

## The theorem being proved

Frame: boxeph-S26 THM-887(I) gives the exact per-owner comb identity; klein-S314c7
THM-887(A) gives the uniform upper max|S| ≤ C·diam; THM-889 (S19) gives the balanced-limit
constants A* = 360/7⁵, B* = 120/7⁵, m̂*. The missing piece — the named residual — is the
RATE connecting the exact comb amplitude to m̂*:

> **(Main identity, all finite/exact.)** For owner e with others f = (f₁..f₅) (+ stationary),
> the a-twisted owner sum Σ_j u_e(j)ω^{aj} equals its independent-limit main term
> e·(1−ω^a)·m̂*_s(a) **exactly** plus the sum over NONZERO congruence relations
> Λ_e^× = {(k, c₀) ≠ 0 : Σkᵢfᵢ + (c₀+a)e ≡ 0 (mod 7e)} of explicit product weights
> ∏ᵢ ĥ(kᵢ) (finite geometric sums, supported on kᵢ ≡ cᵢ mod 7, decaying like folded 1/|k|).
> Key simplification: on the Z_{7e} grid the section marginals are EXACTLY uniform
> (e | section size), so the main term needs NO grid correction.
> **(Rate corollary.)** |ν̂ − main| ≤ explicit weight-sum over Λ^× ≤ C(κ_e)·polylog, κ_e =
> shortest folded relation — incoherent cores ⟹ rate; the THM-889 coherence meter = the
> short-relation Fourier mass, exactly.

Collision-handling convention (shared boundaries) enters as an explicit bounded correction
(≤ 14·Σ_{e'} gcd(e,e')).

-> THM-889, boxeph THM-887 (S26), klein THM-887 (c7), HYP-7017/7018→7019; death-star-S20.
