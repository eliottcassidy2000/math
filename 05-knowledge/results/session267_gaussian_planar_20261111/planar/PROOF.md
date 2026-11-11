# Gaussian rotation, source-depth flags, and recursive boundary memory

Status: PROVED algebraic bridge and composition obstruction; OPEN chart-entry and termination hypotheses. This note does not prove GMC(2), JC(2), or any bound on the degree of a counterexample.

## Inheritance and the conceptual move

The source is Session254 `jacobian/ALL_ROW_DEPTH.md` and `PROOF.md` at commit `946fb3f88f1db2c20a1ad01de0284a0ae7a2ba88` in `eliottcassidy2000/math`. Its chart is u=x²t, p=t(1+x²t), y=xtp. Its all-row result gives terminal tangent degree d_m=floor((m−1)/2)+2, bracket operator

L_m θ = −xθ + (x²+6)θ'/(2m),

and bracket quotient of dimension floor(m/2)−1 for m≥5. Only one of those classes is payable by the same-row source packet. We inherit these facts in precisely their homogeneous scope, without re-auditing nonlinear prefixes.

Closest mechanism: the exact truncated Student sequence. Hostile: high monomials beyond d_m+1 cannot be repaired. Corrected near miss: identifying a one-dimensional source image with a one-dimensional obstruction space. Underused sidecar: homogeneous binary forms and rotation weights. The concept board is Student integration, Gaussian homogenization, SO(2) weights, depth flags, and recursive composition.

The new point is that the Student scalar is already a literal two-dimensional Gaussian observable, not merely asymptotically Gaussian. The other row debts come from a flag restriction. Recursive composition requires angular channels lost by taking only that scalar.

## 1. Exact Gaussian realization and rotation conjugacy

Work over C, with identities defined over Q after adjoining sqrt(6) and i. Let B_m=C[G,H]_(2m) be homogeneous binary forms of degree 2m. Set

T_m(F)=H^(2m) F(sqrt(6) G/H),  deg F≤2m,

and R=H∂_G−G∂_H. T_m is a linear isomorphism from C[x]_(≤2m) to B_m. Direct differentiation gives the all-polynomial identity

R T_m(θ) = (2m/sqrt(6)) T_m(L_m θ).               (1)

Although L_m can raise degree by one, its leading coefficient on x^(2m) vanishes, so it preserves C[x]_(≤2m). This top cancellation is essential.

Let G,H be independent standard real Gaussian variables, extending expectation C-linearly. Define

Λ_m(F)=E[T_m(F)]/E[H^(2m)].

Independence and Gaussian moments yield

Λ_m(x^(2j))=6^j (2j−1)!!(2m−2j−1)!!/(2m−1)!!,
Λ_m(x^(2j+1))=0, 0≤2j≤2m,                         (2)

where (−1)!!=1. These are exactly the inherited Student moments μ_j. Indeed successive even moments have ratio 6(2j−1)/(2m−2j+1). Gaussian rotation invariance, or integration by parts, gives E[R P]=0, recovering Λ_m L_m=0.

There is also a literal projective probability interpretation: tilt the joint Gaussian measure by H^(2m)/E[H^(2m)] and push it forward by x=sqrt(6)G/H. Integrating H in the change of variables gives density proportional to (x²+6)^(−m−1). This realizes all moments of order less than 2m+1. No claim is made about higher moments, which diverge.

## 2. Complete rotational quotient, and why depth creates more debts

Put U=H+iG and V=H−iG. Then R(U)=iU and R(V)=−iV. The basis

E_(m,k)=U^(m+k)V^(m−k), −m≤k≤m,

satisfies R E_(m,k)=2ik E_(m,k). Therefore

ker R=C(UV)^m,   im R=span{E_(m,k):k≠0}.

In particular the full quotient B_m/RB_m has dimension one, and its coordinate is the coefficient a_0 of (UV)^m. For P=Σ a_k E_(m,k),

E[P]=2^m m! a_0,     Λ_m(T_m^−1 P)=c_m a_0,
c_m=2^m m!/(2m−1)!!.                              (3)

Equation (3) gives an explicit projector π_m P=a_0(UV)^m and an explicit inverse to R on its image: divide channel k by 2ik.

This full inverse is usually NOT an admissible source-depth repair. The prior row problem restricts θ to degree≤d_m and F to degree≤m. Under T_m this is the flag

T_m(C[x]_(≤e))=H^(2m−e) C[G,H]_e.

The flag is not rotation invariant. Thus the prior exact quotient consists of the low Student scalar plus high coefficient debts, even though unrestricted Gaussian cohomology has only one class. More explicitly, for d<m and N≤2m, N≥d+1, its coordinates remain

(Σ_(j=0)^(d+1) μ_j F_j, F_(d+2),...,F_N).

The full Λ_m(F) equals the first coordinate plus Σ_(j=d+2)^N μ_j F_j. Only after the rigid high debts vanish do these two scalars coincide. This distinction prevents an invalid Gaussian shortcut.

Example: for the inherited row16 flag, F=x^11 has zero full Gaussian expectation by parity but is an unpaid rigid debt. Thus even exact Gaussian vanishing does not imply admissible row repair.

## 3. Exact recursive angular-channel calculus

If P∈B_m, Q∈B_n, then PQ∈B_(m+n), and

T_(m+n)(FG)=T_m(F)T_n(G).

Writing coefficient arrays a=(a_−m,...,a_m), b=(b_−n,...,b_n), the parent array is convolution:

(a*b)_r=Σ_(k+l=r) a_k b_l.                         (4)

The normalized rotational observable Ω_m(P)=E[P]/(2^m m!) is a_0. Hence

Ω_(m+n)(PQ)=Σ_k a_k b_−k.                         (5)

For the Student normalization, the corresponding formula is
Λ_(m+n)(FG)=c_(m+n) Σ_k a_k b_−k.

Equations (4)–(5) are an exact associative recursive carrier, valid for any binary multiplication tree. Parenthesization creates no anomaly; deleting channels creates lost information. The angular coefficients, radial grade m, and the source-depth flag are separate required sidecars.

A sharp hostile to scalar-only recursion is P=U² and Q=V². Both have zero Gaussian expectation, whereas PQ=(G²+H²)² has expectation 8. The two factors are even pure one-sided angular forms individually; opposing orientations create a neutral channel. P=U²+V² gives another hostile: Ω_1(P)=0 but Ω_2(P²)=2.

Positive cone theorem. Suppose every factor P_j has all nonzero angular channels k≥δ_j>0. Then every product has channels k≥Σδ_j and its Gaussian expectation is zero. If P has channels k≥δ>0 and Q is a fixed homogeneous form with smallest channel −M, then E[P^r Q]=0 whenever rδ>M. The same statement holds with all signs reversed. This is a proved homogeneous sector of the Gaussian eventual-multiplier phenomenon, with an explicit threshold. It is not a classification of all Gaussian moment-null polynomials.

## 4. No fixed-dimensional linear summary can handle all contexts

Theorem. Fix m≥0. Suppose σ:B_m→C^s is linear, and for every Q∈B_m the scalar Ω_(2m)(PQ) can be recovered from σ(P) (the recovery may depend arbitrarily on Q). Then s≥2m+1.

Proof. If σ(P)=0, its observable must agree with that of zero against every Q. Taking Q=E_(m,−k) in (5) extracts a_k. Thus every coefficient vanishes and σ is injective. Dimension gives the bound. No linearity assumption on the recovery map is needed. QED.

Consequently an exact linear finite-state recursion of fixed dimension cannot summarize unrestricted binary forms for all multiplication contexts. This is not an impossibility result for nonlinear summaries, bounded-degree families, or geometrically constrained actual Keller sources. It identifies what a proposed finite-state JC bridge must exploit: restrictions on realizable states, not Gaussian invariance alone. Even retaining the scalar plus a fixed number of extra linear moments cannot work unrestrictedly.

This gives a rigorous version of the discrete-structure analogy: recursive composition is lawful once the full boundary response is retained, but local zero observables do not compose. The carrier here is a graded convolution algebra, not a tournament; there is no intrinsic pairwise orientation that would justify manufacturing one.

## 5. A narrower open program, with its exact logical boundary

Define the actual-source set S_m to consist of bracket discrepancies arising from admissible normalized Keller prefixes in the inherited chart, before the row-m source payment. This definition does not assume that every Keller map enters the chart.

OPEN Rigid-Debt Vanishing: every F∈S_m has F_j=0 for j≥d_m+2, at every m≥5. This is a statement about actual nonlinear source equations, not arbitrary polynomials and not a consequence of (2).

OPEN Depth Lift: after the scalar source payment and unique terminal tangent repair, at least one next-row projected-depth solution exists for every actual prefix under consideration.

OPEN Recursive Compression: the actual admissible prefix family has a recursively closed boundary carrier that retains the source-depth flag and enough angular channels to determine all future rigid debts, with complexity bounded independently of raw polynomial degree. The carrier may be nonlinear; §4 rules out a fixed-dimensional linear carrier on unrestricted forms.

OPEN Effective Termination: there is a well-founded rank on actual repaired prefixes, unchanged by harmless normalizations and strictly decreasing under the chosen inverse-building descent. An infinite row-by-row formal lift does not supply this rank.

Quantified conditional consequence: for any specified class C of polynomial Keller maps for which (i) chart entry into these actual-source sets is proved, (ii) Rigid-Debt Vanishing and Depth Lift hold at every row, and (iii) an inverse-building algorithm starts from that entry, performs only these lawful repairs, has a well-founded decreasing rank, and on termination produces a polynomial inverse, every map in C is invertible. The proof is well-founded termination plus the stipulated correctness invariant. This implication deliberately exposes algorithm correctness and termination; it does not hide them inside the word 'recursive'. Extending C to all Keller maps is another required theorem.

GMC(2) alone is not shown here to imply Rigid-Debt Vanishing, Depth Lift, or chart entry. The concrete new link is (1)–(3), and the precise obstruction to a naive bridge is the non-invariant flag. A useful next theorem would show that actual nonlinear source formation preserves a restricted angular/flag carrier whose forbidden coefficients vanish identically. Such a theorem would target the growing rigid debts rather than repeatedly buying the single scalar.

## Reproduction and scope

Run `python3 planar/check.py`. The finite checks exercise exact conjugacy, Gaussian moments, channel convolution, hostile examples, and the inherited row16 parity obstruction. The all-parameter proofs are above; the script is not an extrapolation certificate. No inherited inhomogeneous lift is reconstructed and no literature-priority claim is made.
