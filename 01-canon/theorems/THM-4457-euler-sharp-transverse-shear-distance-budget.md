---
id: THM-4457
title: "Sharp transverse-shear distance budget for vorticity stretching"
status: >
  PROVED ANALYTIC + FINITE-EXACT CONTROLS + INDEPENDENTLY AUDITED.
  Sharp trace-free matrix inequality and Euler trajectory budget.
  Euler continuation uses classical cited theory. Navier--Stokes
  continuation is a specialization of Miller's published criterion.
  No general regularity or singularity problem is settled.
source: euler-bridge-20260908
depends_on: []
external_inputs:
  - Beale--Kato--Majda continuation, CMP 94 (1984), 61--66
  - Tao H1 Navier--Stokes continuation, APDE 6 (2013), Corollary 5.8
related:
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
script: 04-computation/euler_shear_distance_budget_20260908.py
output: 05-knowledge/results/euler_shear_distance_budget_20260908.out
script_sha256: ed88cec1f7511450d0b6b5691dead7d4db46c32c7aa33bf27a533dc8eb926887
output_sha256: 92f1000ce9235fb2e7c985387f4c957223e35a2fc2381df3e81d93409777dc60
hash_basis: raw LF bytes
audit: >
  Independent analytic derivations; sharp constant checked by a rational
  limiting family. 275562 matrix/shear presentations, 20 rational rotations,
  and 20 explicit pressure-profile controls pass. Normal and optimized
  Python runs are identical with checks enabled. No Lean claim.
---

# THM-4457 -- sharp transverse-shear distance budget

**Status: PROVED elementary matrix and trajectory statements; CITED classical
Euler/Navier--Stokes continuation inputs; OPEN application to singularity
constructions.** The Navier--Stokes regularity consequence below is already
implied by the cited middle-strain-eigenvalue criterion, not a new general
regularity theorem.
Independently reviewed and promoted on 2026-09-08 after exact hostile tests.
This does not solve Euler, Navier--Stokes, or LRC(14). No literature-priority
claim is made for the inequality. There is no Lean formalization claim.

## Inheritance and connection contract

The closest classical mechanism is the vorticity equation
`D_t omega = (grad u) omega`, together with its scalar stretching rate.
The repo's relevant inherited discipline is lawful continuation: a product of
formally available operations need not be realized by a physical trajectory.
For a specific repo instance see
[THM-2680, dilation-reversed clock fibre products](../../01-canon/theorems/THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary.md).
That theorem is methodological context, not a mathematical dependency here.
The new proof below is elementary three-dimensional linear algebra.

- **Source:** Euler velocity-gradient trajectories, or a proposed composition
  of volume-preserving shear maps.
- **Target/map:** `M = grad u` goes to its distance from trace-free rank-one
  shear matrices, retaining the actual vorticity axis.
- **Preserved predicate:** the actual material logarithmic vorticity growth
  is bounded by this distance with the sharp Frobenius constant `sqrt(7/6)`.
- **Destroyed information:** a shear map alone forgets its generating time
  law, pressure compatibility, the full gradient, and the error/background.
- **Needed sidecar:** the pointwise full gradient or a certified error `E`,
  carried along the same physical trajectory and integrated in time.
- **Decisive test:** compare the claimed amplification `log(W1/W0)` with the
  error budget `sqrt(7/6) integral ||E||_F dt`.

The hostile example is an arbitrarily growing coefficient multiplying one
fixed rank-one shear: it preserves volume and remains nilpotent, but its
skew time derivative is incompatible with unforced Euler. The corrected
survivor is a constant-amplitude pure shear, or a growing shear with a
nonzero background/error that pays the bound below. The least-used relevant
sidecar is the distinction between material vorticity direction and the
instantaneous shear's vorticity direction.

## Definitions

For a real `3 x 3` matrix `M`, use the gradient convention `M_ij = d_j u_i`
and define

```text
w(M) = (M_32-M_23, M_13-M_31, M_21-M_12),
sym(M) = (M+M^T)/2,    skew(M) = (M-M^T)/2.
```

Then `skew(M) z = w(M) cross z / 2` and
`||skew(M)||_F = |w(M)|/sqrt(2)`.
Let the closed shear cone, including its zero element, be

```text
N = {a q p^T : a in R, |p|=|q|=1, p dot q=0},
d_F(M) = inf_{S in N} ||M-S||_F.
```

Every `S in N` has trace zero and square zero. Conversely a nonzero
rank-one trace-free real matrix is of this form. All the algebra below is
for trace-free `M`. When `w=w(M)` is nonzero, put

```text
xi = w/|w|,       alpha(M) = xi^T M xi.
```

At zero vorticity `alpha` is left undefined. The differential inequality for
`|omega|`, stated below, still holds there.

## Sharp algebraic stretching bound

**PROVED.** For every real trace-free `3 x 3` matrix with nonzero vorticity,

```text
|alpha(M)| <= sqrt(7/6) d_F(M).                         (1)
```

The constant is optimal. More precisely, for every decomposition `M=S+E`
with `S in N`,

```text
|alpha(M)| <= ||sym(E)||_op + |w(E)|/2
           <= sqrt(7/6) ||E||_F.                       (2)
```

**Proof.** Write `S=a q p^T`, `r=|a|=|w(S)|`, `eta=w(E)`. The case
`r=0` is immediate from the symmetric error term. Otherwise let `n=w(S)/r`
and let `theta` be the angle between `xi` and `n`. The two vectors `p,q`
span `n`'s perpendicular plane, so

```text
|xi^T S xi| = |a (xi dot q)(xi dot p)|
           <= (r/2) sin(theta)^2.
```

Since `w(M)=r n+eta` is parallel to `xi`, projection perpendicular to
`xi` gives `r sin(theta) <= |eta|`. Combining this with `sin(theta)<=1`
yields

```text
|xi^T S xi| <= min(r/2, |eta|^2/(2r)) <= |eta|/2.       (3)
```

The scalar quadratic form of `E` is its symmetric part, proving the first
inequality in (2). Put `D=sym(E)` and `K=skew(E)`. Because `tr(E)=0`, the
symmetric matrix `D` has trace zero. Its eigenvalues give

```text
||D||_op <= sqrt(2/3) ||D||_F,
|eta|/2 = ||K||_F/sqrt(2).
```

Cauchy--Schwarz and Frobenius orthogonality therefore give

```text
sqrt(2/3)||D||_F + ||K||_F/sqrt(2)
 <= sqrt(2/3+1/2) sqrt(||D||_F^2+||K||_F^2)
 = sqrt(7/6)||E||_F.
```

This holds for every `S`, so taking the pointwise infimum proves (1).
No differentiable or measurable choice of minimizing shear is needed.

An operator-norm alternative is `|alpha(M)| <= 2 d_op(M)`, since
`|w(E)|/2=||skew(E)||_op <= ||E||_op`. Optimality is claimed only for
the Frobenius constant above.

### Sharpness at the zero-vorticity boundary

For positive rational `epsilon`, set

```text
S = 1/2 [[ 1, 1,0],
         [-1,-1,0],
         [ 0, 0,0]],

D = diag(2/3,-1/3,-1/3),
K_epsilon = 1/2 [[0,0,0],[0,0,-epsilon],[0,epsilon,0]],
E_epsilon = D - skew(S) + K_epsilon,
M_epsilon = S + E_epsilon
          = diag(7/6,-5/6,-1/3) + K_epsilon.
```

Here `S` is a unit-Frobenius pure shear; `w(M_epsilon)=epsilon e1`,
`alpha(M_epsilon)=7/6`, and

```text
||E_epsilon||_F^2 = 7/6 + epsilon^2/2,
(7/6)/sqrt(7/6+epsilon^2/2) --> sqrt(7/6).
```

Consequently no smaller constant works even before minimizing the shear
distance; (1) together with `d_F(M_epsilon)<=||E_epsilon||_F` proves
optimality for the actual minimized distance as well. The nontrivial sharp
constant is approached as vorticity tends to zero, rather than by pretending
that `alpha` is defined at the limiting symmetric matrix.

## Exact shear-distance formula

**PROVED.** For any real `3 x 3` matrix, trace-free or otherwise,

```text
d_F(M)^2 = ||M||_F^2
 - max_{|p|=1} (|Mp|^2 - (p^T M p)^2).                 (4)
```

For fixed perpendicular unit `p,q`, expanding `||M-a q p^T||_F^2`
and minimizing in `a` gives `||M||_F^2-(q^T M p)^2`. Maximizing in
`q perpendicular p` gives the squared length of the projection of `Mp`
onto that plane, proving (4). Compactness of the unit sphere gives a
maximum and hence a nearest shear. The representation is invariant under
orthogonal conjugation and obeys `d_F(cM)=|c|d_F(M)`.

Formula (4) turns the diagnostic into a quartic optimization over a
two-dimensional compact sphere. A finite set of candidate `p` directions
provides an **upper** bound on the true distance, because it provides a
lower bound on the subtracted maximum. This is useful for certifying an
upper error budget. It must not be misreported as a lower distance bound.

## Euler trajectory consequences

**PROVED from classical Euler's equations.** Let a classical incompressible
Euler solution be given and let `X(t)` be a Lagrangian trajectory. Put
`M(t)=grad u(t,X(t))`, `omega(t)=curl u(t,X(t))`, `W(t)=|omega(t)|`.
The equations imply

```text
M' + M^2 = -H,       H = Hessian(p)(t,X(t)) = H^T,
omega' = M omega.
```

At `W>0`, `d log(W)/dt=alpha(M)`. Thus for `s<t`,

```text
|log(W(t)/W(s))| <= sqrt(7/6) integral_s^t d_F(M(tau)) d tau.    (5)
```

If `W(s)=0`, the homogeneous linear vorticity ODE keeps it zero along
this classical trajectory. Equivalently the inequality
`|W'| <= sqrt(7/6)d_F(M) W` covers both cases.

In particular, if `M(t)` is always a pure shear, then `omega'=0`, because
`S w(S)=0`. Its amplitude `||M(t)||_F=|omega(t)|` is constant, even when
the orthonormal shear frame moves. This statement permits some time-varying
frames with fixed amplitude; it does not assert that every such frame is a
global Euler field.

For a proposed amplification by a factor `A>1`, the necessary budget is

```text
integral_s^t ||E(tau)||_F d tau >= sqrt(6/7) log(A)       (6)
```

for **every** pointwise decomposition into a pure shear plus error along
that trajectory. Pressure compatibility and spatial realization are still
additional obligations. This is an obstruction to a proposed trajectory,
not a method for manufacturing one.

### Uniform continuation corollary

**PROVED using a CITED continuation theorem.** Consider an initially smooth
Euler solution in `H^s(R^3)` or `H^s(T^3)`, `s>5/2`, classical on `[0,T)`,
where `T` is finite. If

```text
integral_0^T ||d_F(grad u(t, .))||_L-infinity dt < infinity,
```

then (5) uniformly bounds `||omega(t)||_L-infinity` on `[0,T)`. The
classical Beale--Kato--Majda continuation theorem extends the solution
past `T`. Thus every finite-time singularity in this stated class must
have divergent integrated uniform shear distance. This is a sufficient
regularity criterion, not an a priori estimate establishing its premise.

There is also a direct full-gradient estimate:

```text
||M||_F <= |w(M)| + (1+sqrt(2)) d_F(M),                  (7)
```

because `||S||_F=|w(S)|<=|w(M)|+sqrt(2)||E||_F` and the triangle
inequality holds for every shear approximation. In particular the same
premise controls `integral ||grad u||_L-infinity dt` after (5).

## Navier--Stokes: scale-critical error budget

**PROVED with CITED continuation input; also a corollary of a prior geometric
criterion.** Fix `nu>0`, `0<T<infinity`, and a divergence-free initial datum
`u0 in H^s(R^3)`, `s>=3`. Let `u` be its unforced Navier--Stokes solution,
smooth in `C([0,T);H^s)` with the usual normalized pressure. In particular
it has finite energy and finite enstrophy on every compact time interval
before `T`. Put `d(t,x)=d_F(grad u(t,x))`. For

```text
3/2 < p < infinity,       q = 2p/(2p-3),
```

the condition `d in L^q((0,T);L^p(R^3))` guarantees continuation past
`T`. The endpoint pair `(p,q)=(infinity,1)` also suffices. This is a
sufficient condition on a given solution; its finiteness is not established
for arbitrary data. No claim is made for arbitrary weak solutions or for
the nonsmall endpoint `L^infinity_t L^(3/2)_x`.

### Independent enstrophy proof with viscosity tracked

Write `c0=sqrt(7/6)`, `omega=curl u`, `y=||omega||_2^2`,
`z=||grad omega||_2`, and `d_p=||d||_p`. The full viscous vorticity
equation is

```text
partial_t omega + (u dot grad)omega = M omega + nu Delta omega.
```

Multiplication by `omega`, spatial integration, incompressibility, and
integration by parts give

```text
(1/2)y' + nu z^2 = integral omega^T M omega dx
                <= c0 integral d |omega|^2 dx
                <= c0 d_p ||omega||_(2p/(p-1))^2.       (8)
```

The pointwise bound is used after multiplying by `|omega|^2`, so it remains
valid at vorticity zeros. The regularity assumptions justify these energy
calculations on every compact time interval; one then takes limits toward
`T`. Choose a homogeneous Sobolev constant `C_S` satisfying
`||f||_6 <= C_S ||grad f||_2` on `R^3`. The same bound for vector fields
follows by applying it to `|f|` and using `|grad|f||<=|grad f|`.
With `theta=3/(2p) in (0,1)`, interpolation gives

```text
||omega||_(2p/(p-1))^2
 <= C_S^(2 theta) y^(1-theta) z^(2 theta).
```

Thus the correctly normalized inequality is

```text
(1/2)y' + nu z^2
 <= B_p d_p y^(1-theta) z^(2 theta),
B_p = c0 C_S^(3/p).                                    (9)
```

The Sobolev factor is necessary; the sharp algebraic constant `c0` is not
the complete constant in (9). An explicit weighted Young inequality gives

```text
B_p d_p y^(1-theta) z^(2 theta)
 <= (nu/2) z^2 + A_p nu^(-3/(2p-3)) d_p^q y,
A_p = (1-theta) (2 theta)^(theta/(1-theta)) B_p^q.
```

Indeed this is obtained by maximizing
`B_p d_p y^(1-theta) X^theta - (nu/2)X` over `X>=0`.
Consequently

```text
y' + nu z^2 <= 2 A_p nu^(-3/(2p-3)) d_p^q y,
y(t) <= y(0) exp(2 A_p nu^(-3/(2p-3)) integral_0^t d_p(tau)^q d tau).
                                                               (10)
```

The energy equality bounds `||u(t)||_2`; Fourier analysis and `div u=0`
give `||grad u(t)||_2=||omega(t)||_2`. Therefore (10) uniformly bounds
`||u(t)||_H1` up to `T`. The classical local `H1` theory, for example
[Tao, Corollary 5.8](https://msp.org/apde/2013/6-1/apde-v6-n1-p02-s.pdf),
extends the solution past `T`. Equivalently, this gives a bounded
`L^infinity_t L^6_x` norm, which falls within classical velocity-based
continuation criteria. Only the stated `H1` continuation theorem is needed.

For `p=infinity`, (8) directly gives
`y'+2nu z^2<=2c0 d_infinity y`; no interpolation or viscosity loss is
needed. At `p=3/2`, interpolation instead yields
`(1/2)y'+nu z^2<=c0 C_S^2 d_(3/2) z^2`; dissipation can absorb this only
under a suitable smallness condition. This proof therefore does not
justify the nonsmall endpoint.

Under the Navier--Stokes scaling `u_lambda(t,x)=lambda u(lambda^2 t,lambda x)`,
the distance scales as
`d_lambda(t,x)=lambda^2 d(lambda^2 t,lambda x)`. Hence the mixed norm is
invariant exactly when `2/q+3/p=2`, as in the statement. The powers of
`nu` in (10) are consistent with the weighted Young step and are not
claimed optimal.

### Prior criterion that subsumes this sufficient condition

Let `lambda_1<=lambda_2<=lambda_3` be the eigenvalues of `sym(M)`.
For a pure shear `S` of amplitude `r`, the eigenvalues of `sym(S)` are
`(-r/2,0,r/2)`. The variational characterization of symmetric eigenvalues
therefore gives, for every `M=S+E`,

```text
|lambda_2(sym M)| <= ||sym E||_op <= sqrt(2/3)||E||_F,
lambda_2^+(sym M) <= sqrt(2/3) d_F(M).                   (11)
```

[Miller, arXiv:1710.05569v4, Theorem 1.1](https://arxiv.org/pdf/1710.05569)
already supplies Navier--Stokes continuation under scale-critical control
of `lambda_2^+` with these exponents. Equation (11) makes the criterion
above a direct sufficient specialization of that published result. The
independent proof (8)--(10) is useful for the explicit error interpretation
and the viscosity budget, not evidence of a new general regularity theorem.

The pointwise trajectory statement (1) is a different assertion. In the
sharpness family above, `sym M_epsilon` has eigenvalues
`(-5/6,-1/3,7/6)`, so `lambda_2^+=0` while `alpha(M_epsilon)=7/6>0`.
Thus the established middle-eigenvalue criterion does not provide a
pointwise bound of the scalar stretching rate by `lambda_2^+`; its
integrated strain identity uses the spatially compatible velocity field.
The shear-distance estimate retains the material vorticity axis and can
charge individual trajectories.

## Controls, failure boundary, and next use

- **Positive:** a constant shear `u=a x_2 e1` is Euler with constant
  pressure; `d_F=0` and vorticity stays constant.
- **Hostile:** `u=a(t)x_2 e1` with nonconstant `a` is divergence-free but
  its acceleration has nonzero curl `-a'(t)e3`; no scalar pressure repairs
  it into unforced Euler. It may serve as a prescribed transport field,
  which is a different role.
- **Direction:** the bound controls actual `omega`, not an arbitrary
  passive vector. A constant shear can amplify a passive vector even while
  its own vorticity is constant.
- **Sharp boundary:** the rational matrix family above checks the squared
  constant `7/6` and the fact that vorticity zeros need separate treatment.
- **Trace:** the `sqrt(2/3)` step uses `tr(E)=0`; dropping incompressibility
  invalidates that constant and changes the vorticity equation.
- **Spatial realization:** individual trace-free matrices and abstract ODE
  trajectories do not establish pressure nonlocality, localization, finite
  energy, boundary conditions, or compatibility with a proposed PDE ansatz.
- **Viscosity:** the material identity `omega'=M omega` belongs to Euler.
  The Navier--Stokes result above separately retains `nu Delta omega` and
  uses enstrophy dissipation; the Euler trajectory logarithm does not
  transfer unchanged to Navier--Stokes.

The next actionable calculation is to compute or bound the full-gradient
shear distance in the candidate shear-cascade construction, and compare
its integrated magnitude to the required logarithmic amplification.
Discarding a nominally lower-amplitude background without integrating its
time cost would erase the quantity that (6) requires.

## Classical source audit

- **CITED:** Beale, Kato, and Majda, *Remarks on the breakdown of smooth
  solutions for the 3-D Euler equations*, CMP 94 (1984), 61--66,
  [DOI](https://doi.org/10.1007/BF01212349),
  [author institution record](https://scholars.duke.edu/publication/759544).
  This supplies the continuation theorem, not the shear-distance bound.
- **CITED:** Gibbon and Holm, *Stretching and folding diagnostics in
  solutions of the three-dimensional Euler and Navier--Stokes equations*,
  [author-hosted paper](https://www.ma.ic.ac.uk/~jdg/jcrjrarx1.pdf),
  equations (1.2)--(1.6) and Theorem 1.1.1, records the classical vorticity
  law and periodic continuation setting. The present scalar stretching
  mechanism is part of that established framework, not a new vorticity law.
- **CITED:** Tao, *Localisation and compactness properties of the
  Navier--Stokes global regularity problem*, Analysis & PDE 6 (2013),
  25--107, [primary publisher PDF](https://msp.org/apde/2013/6-1/apde-v6-n1-p02-s.pdf),
  Corollary 5.8, printed page 56, gives the maximal `H1` continuation
  alternative used after (10). Remark 11.2 explicitly states that bounded
  `H1` norm permits continuation. The source treats viscosity one; a fixed
  positive viscosity is reduced to this by rescaling time and velocity.
- **CITED / PRIOR ART:** Miller, *A regularity criterion for the
  Navier--Stokes equation involving only the middle eigenvalue of the strain
  tensor*, [arXiv:1710.05569v4](https://arxiv.org/pdf/1710.05569),
  [published DOI](https://doi.org/10.1007/s00205-019-01419-z), ARMA 235
  (2020), 99--139, Theorem 1.1. Its time and space exponent letters are
  reversed relative to this note. Equation (11) establishes the precise
  inclusion of the present NS sufficient condition in that prior result.
- Targeted repository searches for `Navier`, `vorticity`, `Euler flow`,
  `Euler equation`, and `pure shear` found no directly overlapping PDE
  theorem. Searches for `Euler` alone mostly recover unrelated polynomial,
  graph, and arithmetic uses. This bounded search does not establish
  literature novelty.

## Reproducible exact controls

Run `python3 04-computation/euler_shear_distance_budget_20260908.py` from the
repo root. [Source](../../04-computation/euler_shear_distance_budget_20260908.py)
and [frozen output](../../05-knowledge/results/euler_shear_distance_budget_20260908.out)
travel together. The 6561 matrices have their first eight entries in
`{-1,0,1}` and the ninth fixed by trace zero. They are tested against six
coordinate shear directions and seven amplitudes `-3,...,3`; zero shears
are counted with their six presentations. The finite checks supplement
the all-matrix proof and do not compute a PDE solution or globally solve
the quartic optimization. Rational orthogonal conjugations, the sharp
zero-vorticity limit, and pure-shear hostiles audit covariance and boundaries.
