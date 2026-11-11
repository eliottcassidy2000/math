# Independent audit of the UV/IR separation note

Scope: the new argument in UV_IR_SEPARATION.md, not the trusted Yang–Mills proof. Reviewed the spectral integral, constants, small-radius estimates, interpretation of the examples, and the diagnostic script. No numerical execution or source rebuild was needed for this analytic audit.

## Verdict

**Proposition U1 is mathematically correct as a separated-point positive spectral covariance construction.** The coefficient 3D/π⁴ and all displayed integral constants are correct. Scoped corrections were made to the declaration name, the explicit dominated-convergence argument, and the contact/Euclidean-positivity interpretation. No change to the construction or its main conclusion was necessary.

## Checks

1. Under s=u²/r², ds=2u du/r², log s=2(L+log u). Multiplying the density c s²/(log s)² by sqrt(s)K₁(r sqrt(s))/(4π²r) produces c u⁶K₁(u)/(8π²r⁸(L+log u)²). The prefactor in the note is exact.
2. For r<e⁻² the splitting interval is ordered: er<sqrt(r). On u≥sqrt(r), L+log u≥L/2, hence the ratio factor is ≤4. The appropriate fixed-domain integrand is 1_{u≥sqrt(r)}u⁶K₁(u)L²/(L+log u)². It tends pointwise on (0,∞) to u⁶K₁(u) and is dominated by 4u⁶K₁(u).
3. On the lower interval, K₁(u)≤1/u and u=rv give L²r⁶∫_e^{r⁻¹ᐟ²}v⁵/(log v)²dv. Since log v≥1, this is ≤L²r⁶(r⁻³−e⁶)/6≤L²r³/6→0. Dividing by the full Bessel moment 384 gives exactly the diagnostic's lower omitted-integral upper bound L²e⁻³ᴸ/2304. This is an upper bound on the omitted lower interval, not a lower bound on its value.
4. The integral representation K₁(u)=(u/4)∫t⁻²exp(−t−u²/(4t))dt is correct. Its inner seventh moment is ∫u⁷exp(−u²/(4t))du=768t⁴. Tonelli then gives (768/4)Γ(3)=384. The same representation gives uK₁(u)=∫exp(−v−u²/(4v))dv≤1.
5. The final coefficient is 384c/(8π²)=48c/π²=3D/π⁴ for c=D/(16π²). Finite σ contributes at most σ([0,e²))/(4π²r²), which vanishes after the r⁸L² rescaling.
6. The spectral examples distinguish the intended properties. With σ=0 the one-particle invariant-mass support begins continuously at e; the Fock theory has a gap and no discrete mass shell. With σ=1_(0,1)ds the spectrum approaches zero with no atom. With σ=wδ_(m²), 0<m<e, there is an isolated one-particle shell and its stable generalized-free particle sector. In this specific example it also lies below min(e,2m), the competing ultraviolet and two-particle mass thresholds. The note's caution that a scalar atom should not generally be confused with isolation from all multiparticle spectrum is sound.
7. Polynomial spectral growth supplies a tempered positive Minkowski mass-shell superposition. However ∫ρ(ds)/(p²+s) diverges at fixed Euclidean p, so its unrenormalized Fourier covariance is not available as an ordinary positive Euclidean covariance on all tests. The original phrasing about a standard coincident Euclidean construction was too broad. It now distinguishes the generalized free Wightman field, its separated Euclidean continuation, and contact extensions, without claiming contact positivity or a continuum probability measure.
8. The numerical script correctly uses the normalized integral and omitted-region bounds for the listed L=10,20,40,80. Its breakpoints are ascending for these values. The numerical upper-tail calculation bounds the actual tail integrand because u≥100 implies L²/(L+log u)²≤1 for L>0, but evaluating that bound by ordinary quadrature does not produce a certified interval. The note already says this. Numerical monotonicity and ratios below one are diagnostics for the listed values, not premises of Proposition U1 or asserted universal inequalities.

## Source-name correction

There is no public declaration `ChallengeSpec.NonFree` in the pinned source. The exact structure is `YangMills.Hierarchy.TraceTest`; its fields include `centered`, `amplitudePositive`, `shortDistance`, and `excludes`. `YangMills.Hierarchy.Properties.nonfree` is the field carrying that structure. `YangMills.publicTrace` in TraceBridge.lean constructs it from `YM.Nonperturbative.Headline.TraceExclusions`.

Pinned declarations:
- https://github.com/yang-mills-group/OpenAI-YM-20261028/blob/3456e5e75cf8d4cd558955a5c7638fba2d68f464/ChallengeSpec.lean#L168-L181
- https://github.com/yang-mills-group/OpenAI-YM-20261028/blob/3456e5e75cf8d4cd558955a5c7638fba2d68f464/TraceBridge.lean#L36-L45

The example is not claimed to reproduce the source's physical stress tensor, Wilson origin, complete hierarchy, or Ward identities. Its valid conclusion is exactly that this one ultraviolet coefficient does not determine infrared spectral support.
