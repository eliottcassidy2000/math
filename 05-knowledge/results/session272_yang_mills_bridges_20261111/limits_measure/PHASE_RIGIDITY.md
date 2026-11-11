# Full-support phase measures turn positive averaged identities into rigidity

Session272. **PROVED** abstract lemmas below; **CONDITIONAL/OPEN** their zero-defect or zero-variance premises for the trusted YM phase hierarchy. Full support, compact metric character space and continuity are inherited inputs, not new conjectures. No claim that the phase space is already a singleton.

Source: pinned `yang-mills-group/OpenAI-YM-20261028` commit `3456e5e75cf8d4cd558955a5c7638fba2d68f464`, `ChallengeSpec.lean:Hierarchy.Reconstructed` (normalized, fullSupport, characterContinuity, mixture) and the supplied compressed whole-word commutative C*-algebra interface. We use that algebra as a unital commutative C*-algebra A with its actual character space X, not an arbitrary auxiliary phase parameterization. The distinction matters for the singleton conclusion.

## P1. Positive defect rigidity and its quantitative form

Let X be compact metric and μ a full-support Borel probability measure. If f:X→[0,infinity) is continuous and integral f dμ=0, then f=0 everywhere.

**Proof.** If f(x)>0, continuity gives a nonempty open neighborhood on which f>f(x)/2. Full support gives positive measure there, contradicting the zero integral. Without nonnegativity, cancellation defeats the conclusion; without continuity, μ-null exceptional points can remain; without full support, points outside the support can remain.

For r>0 put b(r)=inf_x μ(B(x,r)). Then b(r)>0. Indeed choose finitely many centers x_j whose radius-r/3 balls cover X. For each x, one such center satisfies d(x,x_j)<r/3; B(x_j,r/3)⊂B(x,r). The finite minimum of these positive masses bounds b(r) below.

If f is nonnegative and L-Lipschitz, and integral f dμ<=ε, then for every r>0,

    ||f||_infinity <= L r + ε/b(r).

**Proof.** For y∈B(x,r), f(y)>=f(x)-Lr. If the right side is positive, integrate on that ball and use its mass lower bound. If it is nonpositive the desired inequality is immediate. Take the supremum in x.

If additionally b(r)>=c r^d for 0<r<=r_0, with c,d>0 and L>0, choosing r=(d ε/(c L))^(1/(d+1)) when it lies in (0,r_0] gives

    ||f||_infinity <= (d+1) d^(-d/(d+1))
                     L^(d/(d+1)) (ε/c)^(1/(d+1)).

For ε=0 use P1 directly; for L=0, f is constant and ||f||<=ε. For a general compact full-support measure there need not be a polynomial lower-mass bound, so only the b(r) version is automatic. Small averaged defect without any common modulus of continuity does not force small supremum: increasingly thin continuous peaks are the hostile.

**YM use.** For a continuous phase observable h, the defect |h-c|² is nonnegative. If its averaged value vanishes, h=c on every phase, including a selected phase of measure zero. This is stronger than an almost-everywhere conclusion, but the zero average is an additional obligation. The mixture identity for h alone does not provide it.

## P2. Conditional covariance and the mixture-clustering obstruction

Let ω_i be normalized states on a common star algebra, μ a probability measure on phases, and ω=integral ω_i dμ. Suppose α_t preserves each ω_i. For observables A,B write m_A(i)=ω_i(A), and define

    C_i(A,B;t)=ω_i(A* α_t(B))-conj(m_A(i))m_B(i).

Whenever these scalar functions are integrable, exact linearity gives

    C_ω(A,B;t)= integral C_i(A,B;t)dμ(i)
               + integral conj(m_A)m_B dμ
               - conj(integral m_A dμ) integral m_B dμ.

The second term is the phase covariance. For A=B it is Var_μ(m_A)>=0. No positivity of C_i at arbitrary t is needed for this identity.

If integral |C_i(A,A;t)|dμ→0, then

    lim_t C_ω(A,A;t)=Var_μ(m_A).

Sufficient hypotheses include a phase-uniform decaying bound, or pointwise decay with an integrable phase-independent dominating function. For bounded C*-observables, |C_i|<=2||A||² supplies domination whenever pointwise phase clustering holds. Unbounded fields need their own moment/domain bounds.

Thus phasewise clustering plus nonconstant continuous m_A and full support forces a strictly positive asymptotic mixture correlation. Conversely, mixture clustering plus integrated phasewise decay forces Var_μ(m_A)=0 and hence, by P1, m_A is constant on every phase. Cancellation among complex phasewise connected functions cannot evade this conclusion when their integral tends to zero.

The familiar total-variance identity is the t=0 case in a classical conditional-probability realization:

    Var(Y)=E[Var(Y|phase)]+Var(E[Y|phase]).

For states, ω(A*A)-|ω(A)|² has the same decomposition into averaged state variance and phase variance. All pieces are nonnegative by Cauchy–Schwarz for states. A single average first moment does not determine either piece.

**Conditional spectral link.** If the actual source's phasewise heat remainders give a common-family integrated decay for the observables in this identity, the only persistent mixture remainder is phase variance. This is a statement about the specified semigroup/observable realization; one must not replace heat evolution with spacelike translation clustering without a bridge theorem.

## P3. Zero variance on algebra generators collapses the actual character space

Let A be a nonzero unital commutative C*-algebra with compact metric character space X. Let μ have full support on X, and let S⊂A generate A as a unital C*-algebra. For a∈A write a-hat(χ)=χ(a). Assume for every a∈S,

    integral |a-hat - integral a-hat dμ|² dμ = 0.

Then A=C·1 and X has exactly one point.

**Proof.** Each Gelfand transform a-hat is continuous, so P1 makes it constant on X. Faithfulness of the Gelfand transform implies a=c_a 1. The unital C*-algebra generated by these scalar elements is C·1; by generation it is all A. A nonzero scalar unital algebra has exactly one character.

Equivalently, the full-support state φ(a)=integral χ(a)dμ is faithful on positive elements. The variance premise reads

    φ((a-φ(a)1)*(a-φ(a)1))=0.

Faithfulness forces a=φ(a)1. This proof also shows that generation, not a finite list of apparently useful observables, carries the singleton conclusion.

**Exact separation of three objects.**

1. Point observables χ↦χ(a) distinguish actual characters when a ranges over generators, because multiplicativity, star and continuity extend agreement to all A.
2. Measure observables μ↦integral χ(a)dμ need enough a's to identify μ. A point-separating list need not have measure-determining linear span.
3. Product observables χ(a*b)=conj(χ(a))χ(b) belong to the character algebra and supply variances. They are genuine algebra products. A product of two arbitrary phase correlation functions must not be asserted to be one of the original field correlations without an established map into this algebra.

This is the precise positive counterpart of Session267's reciprocal-charge obstruction: retaining the whole-word algebra restores products that scalar compression discards. It provides the *place* to formulate a zero-variance theorem, not the theorem that variances vanish.

**Scope.** If X is merely a compact phase label space mapping noninjectively to the character space, scalar algebra implies only that its algebraic character image is a singleton; duplicate phase labels can remain. Likewise scalar compressed algebra does not alone assert equality of every uncompressed field correlation unless the source's representation theorem connects those correlations to that algebra. No such additional identification is invented here.

## Research payoff and unresolved invoice

Trusted compactness, full support and continuity upgrade **proved averaged nonnegative identities** into every-phase identities. A quantitative upgrade needs a uniform continuity modulus and, for explicit powers, a lower ball-mass bound. To deduce phase collapse, prove zero variance for a generating family of the actual character algebra. To deduce it from clustering, first prove the relevant mixture clustering and integrated phasewise decay for those observables. Those additional zero-variance/clustering premises remain **OPEN in this bridge**, not consequences of normalization or the existence of a selected Wightman phase.
