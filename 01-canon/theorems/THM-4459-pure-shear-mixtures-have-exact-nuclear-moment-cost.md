---
id: THM-4459
title: "Pure-shear mixtures: exact nuclear moment cost and spatial minor obstruction"
status: >
  PROVED ANALYTIC + FINITE-EXACT CONTROLS + INDEPENDENTLY AUDITED.
  Exact algebraic moment problem and affine-periodic W1,2 gradient
  obstruction. No Euler trajectory or PDE regularity claim.
source: cross-concepts-20260908
depends_on: []
related:
  - THM-4457-euler-sharp-transverse-shear-distance-budget
  - THM-4460-adverse-budget-continuation-state-and-shear-composition-collapse
  - THM-3163-universal-finite-prefix-markov-realization-and-physical-sidecar-boundary
  - THM-2352-q-adic-prefix-residue-collision-spectrum
script: 04-computation/cross_concepts_shear_mixture_20260908.py
output: 05-knowledge/results/cross_concepts_shear_mixture_20260908.out
script_sha256: 54c141d25f50130a4009515bd5fe1f7f48ad87c4b86afefb5369399a5457cf8f
output_sha256: cb82c2b7e4479d4b58d289e517b69e0cc20d687736a27521d2ee3add869944ff
hash_basis: raw LF bytes
spatial_script: 04-computation/cross_concepts_gradient_minors_20260908.py
spatial_output: 05-knowledge/results/cross_concepts_gradient_minors_20260908.out
spatial_script_sha256: 9bafc73232a03b19394cd484703cb0702440ef3d704e5db960e8f57c80c6451c
spatial_output_sha256: 88a378864b172f1efbdd0ed026768a180c9cf5dcef7055edb5a2f89b1db74231
audit: >
  Independent polar-factor and hollow-basis proof audit accepted.
  Exact rational polar certificates, rank-deficient boundaries, and
  two independent nonzero-stretching shear-composition hostiles.
  Two independent W1,2 minor/cofactor proof audits and exact Fourier
  convolution versus rational cubature; coefficient-one sharpness audited.
  Normal and optimized replays have all checks enabled. No Lean claim.
---

# THM-4459 -- pure-shear moments and spatial minor obstruction

**PROVED, independently audited; FINITE-EXACT controls.** Every real
trace-free matrix is the mean of pure shears. Their smallest possible
second moment is the square of the mean's nuclear norm. Thus the support
condition "every elementary gradient is a shear" imposes no restriction
on the trace-free mean. Its sharp surviving cost is fluctuation energy.
For an affine-periodic gradient field there is a further obstruction:
second minors must retain their means. Section 4 uses it to exclude such
shear mixtures at non-shear means, even in a quantitative approximation.
These are moment and spatial compatibility results, not a PDE construction
or a claim of literature priority.

## Inheritance and connection contract

The closest current mechanism is [THM-4457 / sharp transverse-shear
distance budget](THM-4457-euler-sharp-transverse-shear-distance-budget.md),
which retains the actual vorticity direction and full gradient. The older
[THM-3163 / finite-prefix Markov realization](THM-3163-universal-finite-prefix-markov-realization-and-physical-sidecar-boundary.md)
shows why an unconstrained realization can be automatic while the physical
transport identity remains missing. The support-versus-multiplicity
distinction of [THM-2352 / q-adic prefix collision spectrum](THM-2352-q-adic-prefix-residue-collision-spectrum.md),
and the continuation mechanism in THM-4460, suggest
asking what the proposed elementary pieces retain after aggregation.
These are methodological connections; the proof below is linear algebra.

The classical tool is a zero-diagonal orthogonal basis for a symmetric
trace-zero matrix. This is established mathematics, reviewed as Lemma 2 in
[Damm and Fassbender, *Simultaneous hollowisation, joint numerical range,
and stabilization by noise*](https://arxiv.org/html/1910.08813v2#S2).
We supply its short proof and combine it with the polar factorization.
No novelty claim is made for that basis lemma, polar decomposition, or
nuclear-norm convexity. A targeted repository search found no earlier
version of the exact shear-constrained moment formula below.

- **Source:** probability laws on trace-free rank-one real matrices.
- **Target/map:** take the mean matrix and the scalar second moment.
- **Preserved predicate:** the exact least second moment at that mean.
- **Destroyed information:** the shear labels, ordering, correlations,
  rank-one compatibility of spatial interfaces, and pressure evolution.
- **Needed sidecar:** the actual gradient field/trajectory, or a lawful
  spatial realization retaining those constraints and its fluctuations.
- **Cheapest hostile:** three equally weighted coordinate shears have
  zero individual self-stretching but their mean has positive stretching.

The canonical near miss is convexifying THM-4457's distance and assuming
that zero cost on each atom survives taking its mean. It does not. The
least-used sidecar is the second moment of the elementary gradients.

## 1. Rank-optimal shear decomposition with optimal total amplitude

Fix a positive integer n. On real n by n matrices use the Frobenius norm
`||.||_F` and nuclear norm `||M||_*`, the sum of the singular values. Let

```
N_n = {q p^T : p,q in R^n, p dot q=0},
```

including zero. Nonzero members are exactly rank-one trace-free matrices;
each satisfies `S^2=0` and `||S||_*=||S||_F`.

**Theorem 1 (PROVED).** If `tr M=0` and `rank M=r>0`, there exist exactly
r nonzero matrices `S_i in N_n` such that

```
M=sum_(i=1)^r S_i,        sum_(i=1)^r ||S_i||_F=||M||_* .       (1)
```

No decomposition into nonzero rank-one matrices uses fewer than r pieces,
and no decomposition into pure shears has a smaller total Frobenius norm.
For `M=0`, use the empty sum.

**Proof.** Write `M=QH`, where `H=(M^T M)^(1/2)` is positive semidefinite,
and extend the partial polar isometry to an orthogonal Q on `R^n`. Put
`V=ran H`; it has dimension r. The symmetric operator

```
K=H^(1/2) sym(Q) H^(1/2)
```

preserves V, vanishes on its orthogonal complement, and has trace
`tr K=tr(QH)=tr M=0`.

Every real symmetric trace-zero operator on V admits an orthonormal basis
`u_1,...,u_r` with `u_i^T K u_i=0`. Here is the elementary induction: if
the operator is zero, use any basis. Otherwise its smallest and largest
eigenvalues have opposite signs; continuity along their unit-circle arc
gives a unit vector of zero quadratic value. Compress to its orthogonal
complement. The compression is symmetric and still has trace zero, so
induction supplies the remaining vectors. The one-dimensional case is zero.

Set `p_i=H^(1/2)u_i`, `q_i=Qp_i`, and `S_i=q_i p_i^T`. Positive
definiteness of H on V makes `p_i` nonzero. Orthogonality of Q gives
`|q_i|=|p_i|`, while the hollow-basis condition gives `p_i dot q_i=0`.
Consequently `S_i` is a nonzero pure shear, and

```
sum S_i=Q H^(1/2) P_V H^(1/2)=QH=M,
sum ||S_i||_F=sum |p_i|^2=tr H=||M||_*.
```

The nuclear-norm triangle inequality gives the amplitude lower bound for
every competing decomposition. Rank subadditivity gives the count lower
bound. QED.

The polar weighting is essential to this proof of optimal amplitude. A
zero-diagonal basis for `sym M` alone gives a shear decomposition by
columns, but does not in general certify the optimal nuclear-norm sum.

## 2. Exact convex hull and least fluctuation energy

**Theorem 2 (PROVED).** For every real trace-free M, among all probability
laws on `N_n` with `E S=M`,

```
min E ||S||_F^2 = ||M||_*^2,
min E ||S-M||_F^2 = ||M||_*^2-||M||_F^2.                    (2)
```

For nonzero M, a minimizing law has exactly `rank M` atoms and every atom
has the same Frobenius norm `||M||_*`. Both minima cover arbitrary laws
with finite second moment, not just a preset finite set of directions.

**Proof.** For any admissible law, convexity and Cauchy--Schwarz give

```
||M||_* <= E ||S||_* = E ||S||_F <= (E ||S||_F^2)^(1/2).     (3)
```

Let `C=||M||_*>0` and take (1). Assign probability
`lambda_i=||S_i||_F/C` to `A_i=S_i/lambda_i`. The probabilities are
strictly positive and sum to one; `E A=M`; each `||A_i||_F=C`. Thus
the second moment attains `C^2`. The identity
`E ||S-M||_F^2=E ||S||_F^2-||M||_F^2` proves the variance formula.
At M=0 use the point mass at zero. Rank subadditivity shows no finite
law at a nonzero rank-r mean can use fewer than r nonzero atoms. QED.

For an arbitrary proposed law, equality requires constant Frobenius norm
almost surely in the Cauchy--Schwarz step and equality in the nuclear-norm
mean inequality. These are necessary and sufficient together; a constant
atom norm alone is not sufficient.

An equivalent exact bounded-amplitude statement is

```
conv {S in N_n : ||S||_F<=1}
  = {M : tr M=0, ||M||_*<=1}.                              (4)
```

One inclusion follows from trace linearity and nuclear-norm convexity.
For the other, the minimizing atoms above have norm `C<=1`; the zero
case is included. No closure operation is needed, and each nonzero mean
uses at most rank M atoms.

In particular `conv N_n` is the whole trace-free subspace. Any nonnegative
convex function on that subspace which vanishes on every pure shear must
vanish identically. The distance from the pure-shear cone cannot be replaced
by such a convex minorant while retaining its stretching obstruction.
[THM-4460 / adverse continuation and shear composition](THM-4460-adverse-budget-continuation-state-and-shear-composition-collapse.md)
proves the companion collapse for nonnegative subadditive minorants.

## 3. Hostiles with the actual vorticity direction

For the gradient convention `M_ij=d_j u_i`, write
`w(M)=(M_32-M_23,M_13-M_31,M_21-M_12)` and, when `w!=0`,
`alpha(M)=w^T M w/|w|^2` as in THM-4457. Every nonzero pure shear
satisfies `S w(S)=0`, so its own scalar stretching vanishes.

Take the law giving probability 1/3 to each of

```
3 E12, 3 E23, 3 E31.
```

Its mean and observables are

```
M=[[0,1,0],[0,0,1],[1,0,0]],
w=(-1,-1,-1), alpha(M)=1,
||M||_F^2=3, ||M||_*=3, min second moment=9, min variance=6.
```

All atoms have zero shear distance and zero self-stretching, but the
mean has `d_F(M)^2=2`. Indeed M is orthogonal and THM-4457's exact distance
formula has subtracted maximum
`max_(|p|=1)(1-(p^T M p)^2)=1`, attained at `p=e1`.
This contradicts the proposed Jensen step `d_F(E S)<=E d_F(S)` directly.

The smallest possible number of summands already suffices for a nonzero
stretching hostile. Put

```
S=E12, T=(1,1,1)^T(1,0,-1).
```

Both are pure shears, but `w(S+T)=(1,-2,0)` and
`alpha(S+T)=-3/5`. One summand can never be such a witness. This example
was independently recovered by the catalytic-localization lane.

## 4. Spatial compatibility restores a sharp quantitative obstruction

The arbitrary-law moment problem has now been solved, so the next test
must constrain its realization. Let `T^n` have normalized measure and let

```
A(x)=M+grad v(x),        v in W^(1,2)(T^n;R^n).
```

This is an **affine-plus-periodic cell gradient**. When M is nonzero it is
not the gradient of a globally periodic real-valued velocity; such a
velocity would have mean gradient zero. No Euler time law is assumed.

**Theorem 3 (PROVED).** For every pair of distinct rows and columns,

```
mean minor_(ik;jl)(A) = minor_(ik;jl)(M).                    (5)
```

Consequently, if `A(x) in N_n` almost everywhere, then `M in N_n`.
In dimension three there is the sharp quantitative bound

```
||cof M||_F^2
 <= (mean ||A||_F^2) (mean d_F(A,N_3)^2).                    (6)
```

In particular a sequence with fixed non-shear trace-free mean M and
bounded mean squared gradient cannot have its mean squared shear distance
tend to zero. No assumption about a finite number of shear directions is
needed. The coefficient one in (6) cannot be decreased.

**Proof of the mean identity.** Expand the minor of `M+grad v`. The
linear derivative terms have mean zero. For smooth v the quadratic term
is

```
(d_j v_i)(d_l v_k)-(d_l v_i)(d_j v_k)
 = d_j(v_i d_l v_k)-d_l(v_i d_j v_k),
```

whose integral vanishes by periodicity. Periodic smooth approximations
converging strongly in `W^(1,2)` have these products converging in `L^1`,
by Holder, so (5) holds at the stated regularity. If A has rank at most
one almost everywhere, every second minor vanishes and (5) forces all
second minors of M to vanish, which is equivalent to `rank M<=1`.
The mean trace is also zero when A is a pure shear. This proves the
qualitative conclusion, including M=0.

This mean identity is a classical null-Lagrangian property, not new
mathematics. For context see [Ball, Currie and Olver, *Null Lagrangians,
weak continuity, and variational problems of arbitrary order*](https://people.maths.ox.ac.uk/~ball/Papers/Ball%2C%20Currie%20%26%20Olver%201981.pdf).
The elementary proof here fixes exactly the boundary and regularity used.

**Proof of the quantitative form in dimension three.** Let
`r>=s>=t>=0` be the singular values of a real 3 by 3 matrix A. The
best unconstrained rank-one approximation has squared residual `s^2+t^2`
by the singular-value decomposition. Since `N_3` is a subset of that
rank-one set,

```
||cof A||_F^2 = r^2(s^2+t^2)+s^2 t^2
 <= (r^2+s^2+t^2)(s^2+t^2)
 <= ||A||_F^2 d_F(A,N_3)^2.
```

The slack in the first inequality is `s^4+s^2 t^2+t^4`. Equation (5)
says `mean cof A=cof M` in dimension three. Triangle inequality followed
by Cauchy--Schwarz gives (6). Both integrals are finite because
`d_F(A,N_3)<=||A||_F`.

For sharpness take the constant affine cell gradient
`M_epsilon=E12+epsilon E23`, with `0<epsilon<=1`. Its singular values
are `(1,epsilon,0)`, and `E12` is an admissible pure shear attaining the
unconstrained rank-one distance. Thus

```
||cof M_epsilon||_F^2=epsilon^2,
||M_epsilon||_F^2=1+epsilon^2, d_F(M_epsilon,N_3)^2=epsilon^2.
```

The ratio of the left side of (6) to its right side is
`1/(1+epsilon^2)`, tending to one. These matrices are trace-free, so
restricting to incompressible cell gradients does not improve the constant.
QED.

For the three-cycle mean in section 3, `||cof M||_F^2=3`. Thus every
affine-periodic cell gradient at that mean with `mean ||A||_F^2<=9`
satisfies

```
mean d_F(A,N_3)^2 >= 1/3.                                  (7)
```

The arbitrary three-atom law in section 3 has the same mean and second
moment 9, but has mean squared shear distance zero. It therefore fails
this precise spatial realization test. These two exact controls distinguish
the constraint on means and energy from the constraint of being a gradient.

The identity (5) concerns second minors in all dimensions. At `W^(1,2)`
we do not extend the cofactor assertion to dimensions greater than three,
where cofactors have higher degree. Nor do we infer (6) for a limiting
probability law merely from its support if gradient concentrations have
lost second-moment information.

## 5. Reproduction, boundaries, and next use

Run from the repository root:

```
python 04-computation/cross_concepts_shear_mixture_20260908.py
python -O 04-computation/cross_concepts_shear_mixture_20260908.py
```

The dependency-free exact verifier checks 324 positive weighted-cycle polar
certificates in dimensions 3 through 6, and 75 neutral-block certificates
in dimensions 2 through 4. The latter include rank-deficient means and
noncommuting polar factors. Rational orthogonal conjugations test basis
invariance. Positive semidefiniteness is supplied by explicit diagonal or
Gram factors; separate squared-singular-value, rank, nilpotence, moment,
and variance checks verify the certificates. Both stretching hostiles and
the zero and nonzero-trace boundaries are included. All **10,434 active
exact gates** pass in normal and optimized Python. The universal theorem
rests on the proof, not on extrapolation from these finite families.

The catalysis agent independently checked the hollow-basis induction,
orthogonal polar extension, strictly positive atom weights, exact lower
moment bound and attaining law. The rank-zero case is separated explicitly.

The additional [gradient-minor verifier](../../04-computation/cross_concepts_gradient_minors_20260908.py)
and [transcript](../../05-knowledge/results/cross_concepts_gradient_minors_20260908.out)
provide 444 active exact gates. They compare rational Gaussian Fourier
differentiation/convolution with exact 64-point cubature on four explicit
fields, check a compatible pure-shear wave and an aliasing hostile, and
test eight rational sharpness parameters. In the main divergence-free
three-wave field, the mean is the cycle matrix, gradient second moment is
9, and the squared Frobenius norm of its mean cofactor is 3. Two nonzero
quadratic minor terms both average to -1/2 and cancel. This is an
independent check of spatial incidence, not just a first-moment calculation.
The [supporting audit](../../03-artifacts/drafts/cross_concepts_gradient_minors_20260908.md)
gives the actual waves and the matched arbitrary-law control. The catalysis
and cyclic-repair agents independently accepted the regularity, gradient
type, quantitative inequality, and sharpness boundaries in section 4.

```
python 04-computation/cross_concepts_gradient_minors_20260908.py
python -O 04-computation/cross_concepts_gradient_minors_20260908.py
```

This result does not permit applying an Euler trajectory identity to an
abstract random matrix. A probability law with a specified mean need not
be a gradient Young measure with the required boundary data, and it need
not be generated by a solution satisfying pressure compatibility. Its
second moment is a moment of velocity gradients, not fluid kinetic energy.
No temporal ordering or spatial localization is supplied by (1) or (2).
Section 4 closes one specified affine-periodic spatial relaxation through
second minors; it does not remove the pressure or boundary obligations
of an actual proposed fluid construction.

The immediate research use is to test an aggregation of a proposed shear
cascade: compare its gradient fluctuation energy to (2), then test its
second minors and spatial boundary class against (5)--(6). Retain the
actual correlations, pressure and time evolution. If these data are
discarded, pure-shear support alone cannot obstruct any trace-free mean.
These results do not close the general smooth-fluid regularity or
singularity problems for three-dimensional Euler or Navier--Stokes.
