# Independent audit: tree moment measure and MRB reconstruction

Session31, 2026-10-22. Audited source: `../session31-tree-moments/TREE_MOMENT_BRIDGE.md` and its exact rational certificate. This packet changes no repository or Git state.

## Verdict

**PASS.** The density, Mellin transform, full right-half-plane reconstruction, and rational Gauss/modified-Gauss enclosure have valid proofs on their stated domains. No substantive mathematical defect was found. One initially ambiguous logarithm sentence was reported to the author for clarification: the condition is Re(s+1)>0, not Re Log(s+1)>0.

The mathematical content is appropriately labeled a synthesis of classical Lambert-W/Stieltjes, Nuttall–Bouwkamp, gamma, and Gaussian-quadrature mechanisms. It does not evaluate the MRB constant, establish a modular identity, or explain the small CM defect. The distinction between a full-measure reconstruction and a one-residue coincidence is essential and is preserved in the source.

The accompanying independent program checks the author's rational records without importing the author's program. It verifies the reversed-denominator orthogonality relations, all Padé matching equations, rational endpoint formulas, and outward decimal rounding. Its direct-density Mellin diagnostics use a different integrand from the author's integration-by-parts diagnostic. All floating-point integrations are explicitly diagnostics, not interval certificates.

## 1. Density, branch banks, and omitted contributions

Let T(z)=-W0(-z), R(z)=T(z)/z, and

F(w)=R(1/w)/w=-W0(-1/w), for w outside [0,e].

The parametrization x(v)=exp(v cot v) sin(v)/v has

x'(v)/x(v)=-[v²+(1-v cot v)²]/v<0, 0<v<π.

It therefore gives a genuine one-to-one coordinate from v in (0,π) to x in (0,e), with orientation reversed. Substitution into W exp(W)=-1/x gives W=-v cot v+iv on the upper bank. Since -1/(x+i0)=-1/x+i0, the Cauchy transform has

F(x+i0)=v cot v-iv.

Thus -Im F(x+i0)/π=v/π is positive, not negative. In the original z-plane the reciprocal change reverses the bank: R(t+i0)=(v cot v+iv)/t. The two formulas are consistent, and the source correctly distinguishes them.

The endpoint estimates needed for the keyhole proof are sufficient:

- At w=e, F is bounded with square-root variation, so an endpoint circle contributes O(ε).
- At w=0, F(w)=O(log(1/|w|)) uniformly on the slit neighborhood, so its endpoint contribution is O(ε log(1/ε)).
- At infinity, F(w)=1/w+O(1/w²); the outer Cauchy integral vanishes.

Consequently no endpoint residue, atom, or entire correction is missing from the reconstructed Cauchy transform. Comparing the coefficient of 1/w proves that the density has mass one. This reasoning uses both the discontinuity and the global analytic/growth conditions. It does not conflict with the separate MRB entire-perturbation obstruction, where those restrictive conditions are not supplied.

The endpoint expansions were independently differentiated/expanded:

log x(v)=1-v²/2-v⁴/36+O(v⁶),

x(π-ε)=(e ε/π) exp(-π/ε)(1+O(ε)).

They imply ρ(x)~π^(-1) sqrt(2(e-x)/e) at e and ρ(x)→1 at zero, with 1-ρ(x)~1/log(1/x). In particular, the integral for R(1/e) is locally of order (e-x)^(-1/2) and is finite. The value R(1/e)=e is legitimate, even though this is a branch point, and differentiating there is not legitimate as an ordinary finite derivative.

## 2. Independent verification of the Mellin proof

Write M(s)=∫₀ᵉ x^s ρ(x) dx, with real log x in the power. Since 0≤ρ≤1, every compact subset of Re s>-1 has an integrable majorant for x^s and each logarithmic derivative. Thus M is holomorphic on this whole connected half-plane.

For real -1<s<0, the positive integrand permits Tonelli:

∫₀∞ t^(-s-1)R(-t)dt
=∫₀ᵉ [∫₀∞ t^(-s-1)/(1+tx)dt] dμ(x)
=Γ(-s)Γ(1+s) M(s).

The x=0 endpoint carries no atom; substitution u=tx is applied for x>0. The beta integral is valid precisely in this strip. Near t=0 the integrand is O(t^(-s-1)), and at infinity it is O(t^(-s-2)log t), so there is no hidden divergence.

Independently substituting t=u exp(u) into R(-t)=W0(t)/t gives

∫₀∞ u^(-s-1)(1+u)exp(-(s+1)u)du
=Γ(-s)(s+1)^s+Γ(1-s)(s+1)^(s-1)
=Γ(-s)(s+1)^(s-1).

The last equality uses Γ(1-s)=-s Γ(-s). Cancelling the finite nonzero Γ(-s) yields M(s)=(s+1)^(s-1)/Γ(s+1) on a real interval. The proposed right-hand side is holomorphic on Re s>-1, using the principal log of s+1. The identity theorem therefore proves the entire claimed domain. An open complex strip of initial agreement is not necessary; a real interval with interior accumulation points suffices. No interpolation from integer values is being asserted.

Integration by parts gives

M(s)=[π(s+1)]^(-1)∫₀π x(v)^(s+1)dv.

Both boundary terms vanish: x^(s+1)ρ(x)→0 at zero for Re s>-1, and ρ(e)=0. This recovers the Nuttall–Bouwkamp identity without using it as a hidden assumption.

As an independent primary-reference cross-check, DLMF 5.9.2_5 states

1/Γ(p)=exp(p)p^(1-p)/(2π) ∫₋π^π exp(-p Φ(t))dt,

for Re p>0, where Φ(t)=1-t cot t+log(t/sin t). Substitution Φ=1-log x immediately reproduces the same integral. Bouwkamp's original solution also explicitly treats the real-positive parameter by a Hankel contour. These sources corroborate, rather than replace, the self-contained Mellin argument.

At p=s+1↓0 the expansion is

M(p-1)=1/p+log p+γ+(p/2)[(log p+γ)²-ζ(2)]
+O(p²(1+|log p|³)).

The source has the correct signs and coefficients. M(-1)=∞ because ρ(0)=1. The limiting Nuttall integral at p=0 is finite only because it equals p M(p-1); it must not be confused with a finite negative moment.

## 3. Gamma poles, branches, and full-kernel reconstruction

For Re w>0, p=1/w also has positive real part. Therefore:

- Γ(p) is finite and nonzero, with none of its nonpositive-integer poles in the domain.
- M(p-1)=p^(p-2)/Γ(p) is nonzero, by the proved formula.
- Principal logarithms obey Log p=-Log w, since both arguments lie in (-π/2,π/2).

Hence

w²/[Γ(1/w)M(1/w-1)]
=w² p^(2-p)=exp(Log(w)/w).

There is no branch multiplier, zero denominator, or gamma-pole exception inside this half-plane. Claims outside the half-plane would require analytic continuation and could no longer be interpreted as the same absolutely convergent Mellin integral. The source does not make that extension.

For the branch-kernel identity, 2p^p/[Γ(p)Γ(1-p)] agrees with (2/π)p^p sin(πp). At positive integers the reciprocal-gamma factor vanishes; it is wrong to substitute a literal infinite gamma value as an undefined numerical quotient. The source correctly specifies reciprocal-gamma continuation. The displayed real interval 0<p<1 is free of all such issues.

For integer k≥1, the transported MRB summand is exactly k^(1/k)-1. Its tail is log k/k+O((log k)²/k²), eventually decreases to zero, and is not absolutely summable. Consecutive paired differences are O(log k/k²), so pairing is allowed in the original order. No interchange of the alternating sum with an unbounded negative-moment integral is asserted or needed.

The full compactly supported measure is uniquely determined by its full integer moments, by polynomial approximation on [0,e]. That supplies all of its fractional moments and the full analytic profile. A finite coefficient prefix or a single branch-kernel/CM value does not supply this determination. In particular, the full reconstruction is nonlinear and does not furnish an MRB evaluation merely by its existence; substituting the gamma formula returns the original series.

## 4. Strict inequalities and coefficient bounds

For every real β>-1, integrating x^β|P(x)|² against the positive interior density proves strict positive definiteness of every finite shifted Hankel matrix. A nonzero polynomial cannot vanish on an interval. This also proves strict moment log-convexity and the required nonsingularity in the rational Padé construction.

The coefficient bound follows directly from

(log x(v))'=-v-(1-v cot v)²/v≤-v.

Thus x(v)≤exp(1-v²/2), strictly for v>0. In the integration-by-parts expression for a_m, extending the resulting Gaussian integral to infinity proves

a_m<e^(m+1)/[sqrt(2π)(m+1)^(3/2)].

The source's endpoint Taylor-tail bound then follows from the ordinary integral comparison for ∑n^(-3/2). Neither a floating-point Stirling approximation nor an unproved asymptotic error term is used.

For the normalized tilted measure x^(N+1)μ/a_(N+1), the mean is r_(N+1)=a_(N+2)/a_(N+1). The kernel 1/(1-zx) is strictly convex when z>0, so Jensen gives the claimed strict lower remainder bound. Its chord between 0 and e is

[1-(e-x)z]/(1-ez),

and integration gives exactly the stated upper bound. Strictness follows from the measure's interior support. On the negative real axis, 1/(1+tx) is also convex and lies strictly between zero and one. The finite geometric-remainder identity therefore yields the claimed sign and magnitude bounds even for t beyond the Taylor radius.

## 5. Rational quadrature certification and a necessary hostile

For a positive measure with infinitely many interior support points, the monic degree-N orthogonal polynomial has N distinct interior nodes and positive Gaussian weights. Its Gaussian formula is exact through degree 2N-1. Hermite interpolation of h(x)=1/(1-zx) at these nodes has error

h^(2N)(ξ)/(2N)! ×∏(x-x_j)²≥0,

strict except at the nodes, for 0<z<1/e. Therefore the Gaussian sum is strictly below R(z).

For ν_b=(b-x)μ with b>e, the same statement holds for I_b. The algebraic identity is

1-z I_b(z)=(1-bz)R(z).

If 0<z<1/b, replacing I_b by its lower Gaussian sum makes the numerator larger and preserves the denominator sign, giving a strict upper bound for R. Both z>0 and bz<1 matter. This is not a complex-domain inequality and is not an upper-bound statement on the negative real axis.

The submitted b=2719/1000 is validated by the exact enclosure

e<4697191/1728000<2719/1000.

The first inequality is the Taylor polynomial through degree ten plus the geometric majorant 1/(10·10!). For the exact rational z=1/3 the denominator is positive.

The independent program also includes a purely rational hostile showing why the denominator condition cannot be dropped. Put z=1839/5000. Then z·4697191/1728000<1, proving z<1/e, but bz>1. The N=1 modified-Gauss expression U is strictly negative, whereas R(z)>0. Thus retaining only z<1/e would make the advertised upper bound false.

The rational computation is a certificate of the displayed endpoints when coupled with these analytic theorems. Lambert-W floating-point values are not used to determine the certified inequality direction. Large rational numerators, high working precision, or a numerical midpoint comparison alone would not establish it.

## 6. Reproduction and source trail

Run:

    python session31-moment-audit/audit_tree_moments.py

Outputs are `independent_audit.json` and the captured `independent_audit.out`. The JSON records source SHA-256 digests so the replay is tied to an exact version. It separates exact arithmetic from numerical diagnostics and includes the denominator hostile.

Primary sources inspected independently on 2026-10-22:

- C. J. Bouwkamp, *A Conjectured Definite Integral (A. H. Nuttall)*, SIAM Review 28(4), 568–569 (1986), original article scan: https://pure.tue.nl/ws/files/3257614/Metis219248.pdf . The paper proves the positive-real-parameter identity through a Hankel contour and credits other contemporary solutions.
- NIST DLMF 5.9.2_5, https://dlmf.nist.gov/5.9.E2_5 . This records the equivalent gamma representation for Re p>0, including the required phase function.
- G. A. Kalugin, D. J. Jeffrey, R. M. Corless, and P. B. Borwein, *Stieltjes and other integral representations for functions of Lambert W*, author-hosted preprint: https://www.uwo.ca/apmaths/faculty/jeffrey/pdfs/StieltjesLamW.pdf . Theorem 2.1 and the boundary parametrization are classical antecedents of the compact reciprocal-coordinate representation.

No external source is needed to accept a computational output: the calculations are reproducible locally, and the proof audit above independently checks every analytic step used by the rational brackets.

## 7. Audit of the final mirror and endpoint-model additions

The later §§5.1, 6.1, and 6.2 were audited separately and also pass analytically.

### Actual mirror input

The positive atanh series for log 3 uses t=1/2. After K terms, replacing each omitted denominator by 2K+1 yields exactly the author's strictly positive geometric remainder bound. Division by three produces a rational enclosure for z3=log(3)/3, and the exact inequality b z+<1 keeps both endpoint evaluations within the permitted bracket domain. Monotonicity of R gives

G_N(z-)<R(z-)<R(z3)<R(z+)<U_b,N(z+).

No monotonicity assertion about the upper rational approximation itself is needed. The raw-series comparison uses precisely the stated moment set a0 through a2N, and the same rational input interval.

The added independent program `audit_mirror_endpoint.py` takes an entirely different route to the final emitted decimal intervals. For an endpoint r in (1,e), it rigorously encloses log r by the positive atanh series with t=(r-1)/(r+1), then compares log(r)/r with the rational interval for log(3)/3. Because log(r)/r increases strictly on (1,e), this certifies the lower and upper mirror bounds without using Lambert W, the author's Padé solver, or any floating-point comparison. If a wide upper endpoint already exceeds a rational upper bound for e, it is automatically above the smaller mirror root. Thus the decimal interval itself receives a second, independent exact certificate.

### Fixed-quadrature obstruction

For a fixed positive quadrature with positive nodes, the finite sum M_N(p-1) tends to B_N=∑ω_j/x_j in (0,∞). Therefore k²/[Γ(1/k) M_N(1/k-1)]~k/B_N, not one. The proposed summand fails the necessary term test. A positive node at zero instead gives infinite negative moments, so it is not a remedy. The qualification that N remains fixed is necessary; this argument says nothing about a controlled N=N(k) tending to infinity.

### Continuous endpoint model

The model is correctly defined on (0,1), zero elsewhere, and correctly called nonnormalized. With t=-log x it has positive density t/(1+t), hence

H(p)=∫₀∞ exp(-pt)t/(1+t)dt=1/p-exp(p)E1(p)>0.

From E1(p)=-γ-log p+p+O(p²), multiplication by exp(p) gives

H(p)=1/p+log p+γ+p(log p+γ-1)+O(p²|log p|).

Subtracting from the already-audited M expansion yields

M(p-1)-H(p)=(p/2)log²p+O(p|log p|).

Since Γ(p)~1/p and both H(p), M(p-1)~1/p,

[1/(p²Γ(p))][1/M(p-1)-1/H(p)]
=-(p²/2)log²p+O(p²|log p|).

This verifies both the sign and scale in the claimed true-minus-model kernel asymptotic. The correction is absolutely summable over k. Consequently the model alternating series converges by subtracting this absolutely convergent correction from the original convergent MRB series. No unsupported eventual monotonicity, acceleration claim, or nome-scale conclusion is required.

### Completed final replay

Both programs completed successfully. The exact benchmark replay covers N=1,2,4,8,12,16,24,32 at z=1/3 and independently verifies the 80-decimal outward rounding. Its narrowest exact bracket has width approximately 4.590587295233499953945335×10^(-18).

The independent rational-logarithm certificate covers every submitted actual-mirror interval, N=4,8,16,32,64. At N=64 it proves strictly:

2.478052666984991255845343341327140725699400658342168358878463
< R(log(3)/3) <
2.478052699570190519812004143717774468873226376897534097046115.

The width of the author's underlying exact bracket is approximately 3.258519926396666×10^(-8); the displayed decimal interval is separately proved by the independent root-comparison argument. The original-density Mellin diagnostics, complex reconstruction diagnostics, endpoint diagnostics, strict-remainder diagnostics, model integral diagnostics, and signed asymptotic diagnostics all pass. The floating-point outputs remain diagnostics throughout.
