# Rank-16 angular separator and all-temperature saddle certificate

Session29, 2026-10-22. Independently checked derivations; no novelty claim.

## Outcome

The lattices E8 ⊕ E8 and D16+ have identical radial theta series E4², yet different minimum-shell quartic operators. This provides a rigorous example of information invisible to every radial shell count. More strongly, E8 ⊕ E8 has Hessian inertia (71 positive, 64 negative, 0 zero) at every positive Gaussian temperature. D16+ has inertia (15,120,0) at the self-dual temperature and is a strict local minimum in both sufficiently remote temperature tails. The latter therefore must undergo at least one reciprocal pair of Hessian degeneracies.

The spectra and the modular Hessian method are classical. The elementary all-temperature proof below is supplied as a reusable certificate, not asserted to be new literature.

## Conventions

Set ΘL(t)=Σx∈L exp(−πt||x||²), t>0. Shapes vary by Lε=exp(εH/2)L, where H is real symmetric and tr H=0. The tangent inner product is tr(HK); infinitesimal rotations are already removed. Hessians below are second derivatives, not the ε² Taylor coefficient. Thus all signs and multiplicities are independent of a potential factor-of-two convention in other sources.

For either rank-16 lattice, define M by

⟨H,MH⟩=Σ||x||²=2 (xᵀHx)²,

and A=M−(40/3)I. The scalar is 2·2²·480/(16·18).

## Exact root-shell separation

E8 roots are the 112 vectors ±ei±ej in R8 and the 128 half-integral vectors (±1/2,…,±1/2) having even minus parity. Direct summation gives

ΣE8 (xᵀKx)²=12(tr K)²+24 tr(K²).

For H partitioned into two 8×8 blocks, the roots of E8⊕E8 occupy one block at a time. Therefore its M eigenspaces are:

- Two traceless diagonal blocks: eigenvalue 24, dimension 70.
- Relative scalar block diag(aI8,−aI8): eigenvalue 120, dimension 1.
- Off-block matrices [[0,K],[Kᵀ,0]]: eigenvalue 0, dimension 64.

The A spectrum is (32/3)^70, (320/3)^1, (−40/3)^64.

For D16+, the norm-2 shell is precisely ±ei±ej. Each pair i<j contributes

4(Hii+Hjj)²+16Hij².

Using tr H=0 gives Σi<j(Hii+Hjj)²=14Σi Hii². The M eigenvalues are 56 on trace-free diagonal matrices (dimension15), and 8 on off-diagonal symmetric matrices (dimension120). Its A spectrum is (128/3)^15, (−16/3)^120. Both A operators have trace zero.

The checker independently constructs all roots, a 135-dimensional trace-free basis, and every one of 18,225 entries of each bilinear Gram form. Doubled coordinates make every operation integral. The Gram forms agree with the displayed formulas exactly; no numerical eigensolver is used.

## Modular reduction and common scalar

Write θ(t)=E4(it)=ΘE8(t), F(t)=θ(t)². Every shell is a spherical 2-design. Degree-4 harmonic theta series have weight12 and are multiples of Δ, with their first coefficients prescribed by A. Harmonic decomposition consequently gives

HessL(t)=a(t)I+b(t)AL,

a(t)=t²F″(t)/144+tF′(t)/16,
b(t)=(πt)²Δ(it)>0.

The inequality b>0 follows from Δ(it)=q∏m≥1(1−q^m)^24, q=exp(−2πt). Both lattices therefore share the entire scalar temperature function r=a/b, but not the angular operator.

## A direct all-temperature negative subspace

For E8⊕E8 and H=[[0,K],[Kᵀ,0]], product factorization and the isotropic second moments of E8 give the Hessian eigenvalue

λcross(t)=t²θ′(t)²/32+tθ(t)θ′(t)/8
          =(tθ′(t)/32)(tθ′(t)+4θ(t)).

Indeed xᵀHx=2uᵀKv; the second-moment tensor in either factor is −θ′/(8π) times the identity, and ||H||²=2||K||F².

Put g(t)=tθ′(t)/θ(t). Nonzero lattice vectors imply θ′(t)<0. Poisson summation gives θ(t)=t^−4θ(1/t), hence

g(t)+g(1/t)=−4.

Both g values are negative, so −4<g(t)<0. It follows that λcross(t)<0 for every t>0. This proves a 64-dimensional strictly negative subspace without a numerical temperature scan.

## Strictly positive complementary directions

All E8 shells are spherical 4-designs: their degree-2 and degree-4 harmonic theta series vanish since the relevant cusp spaces of weights6 and8 vanish. For a traceless symmetric 8×8 deformation the E8 Hessian is the scalar

h8(t)=t²θ″(t)/40+tθ′(t)/8.

For t≥1, every nonzero norm is at least2, giving θ″≥2π(−θ′), and thus

h8(t) ≥ (−tθ′/8)(2πt/5−1)>0.

For t≤1, apply Poisson summation to the whole deformed determinant-one lattice. Its dual is exp(−εH/2)E8, so differentiating twice yields h8(t)=t^−4 h8(1/t)>0. No deformation of t is needed.

In E8⊕E8, the 70 traceless-block directions have eigenvalue θh8>0. The relative-scalar direction differs from this eigenvalue by (120−24)b=96b>0. Therefore all 71 block-diagonal directions are strictly positive. Together with the cross directions, this proves the exact inertia (71,64,0) for every t>0.

A reusable general principle emerges: an orthogonal square of a self-dual rank-d lattice with isotropic Gaussian second moment has cross-block eigenvalue

(2tθ′/d²)(tθ′+(d/2)θ).

Poisson duality traps tθ′/θ strictly between −d/2 and0, making every cross-block direction negative. The positive complementary directions need separate hypotheses; the E8 argument above supplies them here.

## Exact self-dual calculation

At t=1, E2(i)=3/π and E6(i)=0. Ramanujan differentiation gives

F′(1)=−4F(1),
F″(1)=18F(1)+(10π²/9)E4(i)³.

As Δ(i)=E4(i)³/1728, these imply

r(1)=40/3−216/y,  y=π²E4(i).

An elementary strict certificate is 9<y<27. For the lower bound, π>3 and E4(i)>1. For the upper bound, π²<10, exp(6)>400 and σ3(n)≤n⁴ give

E4(i) < 1+240 S(1/400),
S(q)=Σn≥1 n⁴qⁿ=q(1+11q+11q²+q³)/(1−q)^5.

Exact rational evaluation yields 10[1+240S(1/400)]<27. The exponential bound exp(6)>400 follows already from its positive Taylor terms through degree19, checked exactly.

Thus E8⊕E8 has signs r+32/3=24−216/y>0, r+320/3>0, r−40/3<0. D16+ has r+128/3=56−216/y>0 but r−16/3=8−216/y<0. Its self-dual inertia is (15,120,0), with no approximate evaluation needed.

## D16+ transitions and what is not proved

The minimal-shell M has smallest eigenvalue8>0. As t→∞ the Hessian divided by exp(−2πt) has leading operator π²t²M−πt(60)I, while higher shells are exponentially smaller. Hence it is positive definite for sufficiently large t. Poisson duality preserves its Hessian signs under t↔1/t, giving the other stable tail.

Moreover, r+128/3 exceeds the everywhere-positive E8-block coefficient r+32/3 by32, so the diagonal D16+ branch is positive at every t. Only its 120-dimensional off-diagonal branch can change sign. Continuity and the self-dual negative value force a zero above1 and its reciprocal below1. This does NOT prove uniqueness of that reciprocal pair, monotonicity, or locate a certified threshold.

## Interpretation and limits

Scalar theta data cannot reconstruct the angular response, even when every shell count agrees at every radius. The analogy to Dehn invariance is an obstruction to inference from an incomplete invariant, not an identification of the classical Dehn invariant with A. Nothing here proves a scissors-congruence theorem, a rank-88 angular spectrum, or a tournament/metagraph correspondence. The rank-16 example is an explicit test model for the same scalar-plus-angular organization used in rank88.

## Primary reference and literature scope

Heimendahl, Marafioti, Thiemeyer, Vallentin and Zimmermann, *Critical even unimodular lattices in the Gaussian core model*, arXiv:2105.07868v2 (2022):
https://arxiv.org/html/2105.07868v2

Sections2–4 develop harmonic theta Hessians and root-system operators. Section5.2 states equal theta series and the two M spectra, proves sufficiently-large-parameter behavior, and reports numerical evidence for an all-parameter E8⊕E8 saddle. Our calculations above independently rederive the spectra and provide an elementary uniform sign proof. A limited web search did not locate a later resolution; this is not a novelty certification.

## Reproduction

Run `python check_rank16.py` and `python -O check_rank16.py`. Both passed and produced byte-identical JSON. The exact finite checker certifies root moments and rational self-dual bounds; the all-temperature assertions are the analytic proofs above, not conclusions from a finite test.

## Spectral moments, tuple information, and the genus-four obstruction

The first three angular power traces coincide:

tr A=0, tr A²=30,720, tr A³=1,146,880.

But the fourth separates them:

tr A_E8²⁴=132,382,720,
tr A_D16⁴=49,807,360,
difference=82,575,360.

These exact integers are checked from the certified spectra. The notation A_E8² here denotes the angular operator for the lattice E8⊕E8, not the square of an operator.

This is related precisely, not just numerologically, to tuple theta invariants. For a root x let u_x=xxᵀ−(1/8)I16, the trace-free projection of xxᵀ. Then M=Σx u_x⊗u_x and

⟨u_x,u_y⟩=(x·y)²−1/4.

For each positive integer k,

tr M^k = Σ(x1,…,xk)∈R^k ∏j=1…k [(xj·xj+1)²−1/4],  xk+1=x1.

This rank-one trace identity makes tr M^k a linear functional of the ordered k-tuple Gram distribution. In the genus-k theta series convention

ΘL^(k)(Z)=ΣX∈L^k exp(πi tr(Gram(X)Z)),

the relevant coefficients are exactly those with diagonal entries2. The cyclic weight above is a polynomial in Gram entries; summing it against these finitely many coefficients gives tr M^k. The binomial identity A=M−(40/3)I likewise makes tr A^k a combination of genus-j coefficient functionals for j≤k, with the j=0 term interpreted as dimension135.

For these two lattices, the genus-1,2,3 theta series agree, while genus4 differs by a nonzero Schottky form. Therefore equality of the first three traces is forced by the stronger classical tuple theorem. The nonzero fourth-trace gap is an explicit separating linear functional of the genus-four difference: since all lower power traces of M agree, the fourth-trace gaps of A and M agree. This argument proves separation without enumerating 480⁴ ordered quadruples.

Important qualifications: power traces contain less information than full tuple Gram distributions; equality of the first three traces alone does not prove equality of the genus≤3 series. The genus-four difference vanishes upon restriction to the genus-four Jacobian locus; nonvanishing on the full Siegel upper half-space is the distinction meant here.

Primary reference: Oura, Poor, Salvati Manni and Yuen, *Modular Forms of Weight 8 for Γg(1,2)*, arXiv:0811.2259v2, p16 discussion of Witt's problem (PDF page16), identifies the equalities through genus3 and the genus4 Schottky difference, crediting Igusa and Kneser:
https://arxiv.org/pdf/0811.2259

This gives a rigorous bridge from the attachment's theme of incidence/tuple information to angular Hessian spectra. It does not establish a correspondence with tournament induced-profile invariants; that remains a separate proposed research direction.
