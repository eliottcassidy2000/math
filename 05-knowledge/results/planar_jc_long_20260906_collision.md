# All-order moving collision periods, arbitrary clocks, and the density boundary

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
In the fixed source coordinates below, every collision-compensated
constant-normal pencil in the inherited three-parameter Russell family,
with any nonzero formal clock, excludes a nonzero constant pullback density
for every regular or formal target two-form. The obstruction has exact
order between one and four times the clock order. The proof uses the whole
moving collision, rather than a truncation of the source map. An explicit
compatible formal unit density and source-coordinate change give a sharp
boundary: arbitrary source volume changes can escape this conclusion.

This is a theorem about a specified formal-local family. It proves no
all-order identification of target-lift debts, polynomial termination,
polynomial Keller pair, JC(2), or DC(2).

## 1. Inheritance, novelty check, and live board

The closest proved mechanisms are
[THM-4411 / first-order-collision-transgression-seminormal-tradeoff](../../01-canon/theorems/THM-4411-first-order-collision-transgression-seminormal-tradeoff.md),
[THM-4424 / russell-constant-normal-debt-discriminant-contact-correspondence](../../01-canon/theorems/THM-4424-russell-constant-normal-debt-discriminant-contact-correspondence.md),
and the independently audited
[ninth-source collision period](overnight8_20260906_jc_collision_period.md).
THM-4424 already constructs the unique formal compensator and identifies
its first live coefficient with each of four inherited retained debts.
Those finite formulas are dependencies here, not newly derived canon.
The extension proved below is an exact all-order period identity and its
consumer for arbitrary clocks throughout the full parameter family.

The canonical hostile is THM-4411's normal `h=x`, which fills the missing
period by splitting the collision. The corrected near miss is to treat
the finite four-level debt matching as an unproved all-order debt law, or
to forget the source volume form after a coordinate change. The least-used
sidecar is the weighted first moment of the moving labelled sections.
[THM-4404 / exceptional-quartic-descended-two-form-seminormal-cokernel](../../01-canon/theorems/THM-4404-exceptional-quartic-descended-two-form-seminormal-cokernel.md)
supplies the relevant differential, rather than mere equality of target
values. MISTAKE-494's graph-period warning is respected: no graph-period
completeness or equality of graph equalizers with target algebras is used.

The live board consists of five coupled objects:

| Object | Retained information | Decisive check |
|---|---|---|
| Common moving target | Every labelled branch has the same image | Differentiate the collision before taking a period |
| Full target two-form module | All three local coordinate slots and their common coefficients | Wedge with common motion |
| Transverse compensator | Its derivative and the section first moment | Exact identity `Pi+chi' Xi=0` |
| Source clock | Order, leading coefficient, and fixed coordinate `w=t` | Compose the period, not the source density |
| Source volume form | The Jacobian under a source change | Explicit compatible unit-density hostile |

Exact-statement, `Pi/chi/Xi`, moving-relation, clock, and source-unit
searches in the named canon and current JC results found only the finite
THM-4424 first-live identities. This is a repository inheritance check,
not a literature-priority claim. Incoming commit `d0208173b2` was also read
through `git show origin/main:05-knowledge/results/planar_jc48_sep06_collision.md`.
That independently proved result identifies a quadratic blind space for
scalar periods on tangents spanning a three-dimensional target, and gives
an all-order vertical graph test with fixed surface coordinates. Here the
tangents remain in a surface plane and carry a nontrivial relation. We
only use necessity of the period identity; we do not claim the period
test is sufficient in a three-dimensional target.

## 2. The all-order common-motion identity

Let `k` have characteristic zero. Let `f(x,t)` be a formal map near finitely
many labelled source points into one smooth target germ. Suppose formal
sections `z_i(t)` satisfy

```
f(z_i(t),t)=gamma(t)
```

for one common target section. Put `T_i(t)=f_x(z_i(t),t)`. If
`lambda_i(t) in k[[t]]` satisfy `sum_i lambda_i T_i=0`, then every formal
target two-form `Omega`, and hence every regular target two-form, obeys

```
f^*Omega=J_Omega(x,t) dx wedge dt,
sum_i lambda_i(t) J_Omega(z_i(t),t)=0.                 (1)
```

**Proof.** Differentiate the actual collision equation:

```
f_t(z_i,t)+z_i'(t)T_i=gamma'(t).
```

At the common target point, alternating bilinearity gives
`J_Omega(z_i,t)=Omega_gamma(T_i,gamma')`; the tangential term drops out.
Summing against the relation gives (1). Target coefficients are evaluated
at the same target section, which is essential. No closedness,
decomposability, dimension-two hypothesis, or period-completeness assertion
is required. This proves (1) at every formal order.

In particular a nonzero constant density `J_Omega=c in k^*` requires

```
Pi(t):=sum_i lambda_i(t)=0.                           (2)
```

For several labelled collisions, (1) holds separately for every collision
and every relation. The resulting vector of periods is retained; adding
its coordinates to one scalar can erase a nonzero obstruction over `C`.
This is a necessary system, not a sufficient descent criterion.

## 3. The exact compiler identity and its compensator consequence

Use the inherited compiler, with `x,q` independent in this calculation:

```
D=1+x^2 q,
B=(D-1)(D+2)^2,
C=xD(D+2),
E=q(D+3),                    C^2 E=B(B+4).
```

Direct differentiation gives the polynomial identity

```
det partial(C,E)/partial(x,q)=6(B+2).                 (3)
```

For example `C_x=5D^2+2D-4`, `C_q=2x^3(D+1)`,
`E_x=2xq^2`, and `E_q=2(D+1)`, where the `x` derivatives keep `q` fixed.
Their determinant is `6(D+1)(D^2+2D-2)=6(B+2)`.

Set

```
P=x^2(x^2-1)^2,
Q1=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4,
Q=Q1+P(a+bx+cx^2),                  (a,b,c) in k^3.
```

THM-4424 supplies a unique formal pencil and labelled sections

```
q_s=Q+s+chi_Q(s)x,
chi_Q(0)=chi_Q'(0)=0,
z_-(0)=-1, z_0(0)=0, z_+(0)=1,
(C(z_i,q_s(z_i)),E(z_i,q_s(z_i)))=(C_*(s),E_*(s)).
```

The six-equation implicit Jacobian has determinant `-288`, so this solution
exists uniquely over every characteristic-zero coefficient field. The
common target begins at `(B,C,E)=(0,0,-3)`. Since the surface relation has
derivative `-4` with respect to `B` there, the `B` values also agree. In
particular `6(B_*(s)+2)` is a unit with constant term `12`.

Let `T_i(s)` denote the derivative of `(C,E)` with respect to the source
`x`, including `q_s`'s `x` dependence, at `z_i(s)`. Use the inherited
normalization

```
u=(det(T_0,T_+), det(T_+,T_-), det(T_-,T_0)),
lambda=u/(-u_0),
lambda(0)=(5/18,-1,13/18).
```

The denominator has constant term `54`, so this is a formal relation at
every order. Define the additional first moment

```
Pi_Q(s)=sum_i lambda_i(s),
Xi_Q(s)=sum_i lambda_i(s) z_i(s),     Xi_Q(0)=4/9.
```

**Theorem.** The exact all-order identity is

```
Pi_Q(s)+chi_Q'(s) Xi_Q(s)=0.                          (4)
```

**Proof.** Apply (1) with the two-form `dC wedge dE`. Its pullback density,
by (3) and the chain rule, is
`6(B+2)(1+chi_Q'(s)x)`. At the common sections the first factor is the
same unit, so it cancels from the period identity. This leaves (4).

Since `Xi_Q` is a unit, (4) implies

```
Pi_Q identically zero  iff  chi_Q identically zero,
ord_s Pi_Q=ord_s chi_Q-1                              (chi_Q nonzero).
```

More precisely, if `chi_Q=chi_n s^n+O(s^(n+1))`, then

```
Pi_Q=-(4n/9)chi_n s^(n-1)+O(s^n).                    (5)
```

These statements use characteristic zero, including the implication
`chi_Q'=0 => chi_Q=0`. Equation (4) relates the actual period and
compensator; it makes no assertion identifying later target-lift debts.

## 4. Uniform nonvanishing throughout the three-parameter family

The following stratification is inherited from THM-4424:

```
chi_2=-(36a+16b+36c+259)/9.
```

If it vanishes, write

```
a=-259/36+p+4r, b=-9r, c=-p.
chi_3=16(729p-42120r^2+15192r+5717)/2187.
```

If this too vanishes, then

```
p=(520/9)r^2-(1688/81)r-5717/729,
chi_4=16F(r)/177147,
F=72783360r^4-77822208r^3-28419741r^2+7849770r-1276420.
```

On `F(r)=0`,

```
chi_5=64G(r)/3^17,
G=4049817787-19653044202r+87403912326r^2+104569906176r^3.
```

The polynomial `F` is irreducible over `Q` and `gcd(F,G)=1`. Thus `G`
is nonzero at every root of `F`, in every characteristic-zero extension.
Consequently every parameter triple has a finite first nonzero
compensator coefficient `chi_n`, with `n in {2,3,4,5}`. In particular
`chi_Q` and `Pi_Q` never vanish identically anywhere in this family.
This conclusion combines the inherited exhaustive finite stratification
with the new exact identity; it is not inferred from finitely many samples.

On the exceptional quartic, THM-4424 gives `chi_5=-9*kappa/20 != 0`.
Equation (5) therefore gives

```
Pi_Q(s)=kappa s^4+O(s^5).                            (6)
```

## 5. Arbitrary clocks, including earlier source terms

For any nonzero formal series `phi(t) in t k[[t]]`, form the specific map
into the smooth Russell surface times an affine line:

```
q_phi(x,t)=Q(x)+phi(t)+chi_Q(phi(t))x,
f_phi(x,t)=(B(x,q_phi),C(x,q_phi),E(x,q_phi),w=t).     (7)
```

The labelled sections `z_i(phi(t))` have one common target image, and their
`x` tangent vectors have last coordinate zero. The relation
`lambda_i(phi(t))` is therefore a relation in the full target as well.
By (1), every target two-form has weighted density period zero. A nonzero
constant pullback density would instead have period `c Pi_Q(phi(t))`.

Write `phi(t)=a_r t^r+O(t^(r+1))`, with `r>=1` and `a_r !=0`. If `chi_n`
is the first nonzero coefficient above, then composition of (5) gives

```
Pi_Q(phi(t))=-(4n/9)chi_n a_r^(n-1)t^(r(n-1))
             +O(t^(r(n-1)+1)).                       (8)
```

This is nonzero. Hence **no regular or formal target two-form pulls back
under (7) to `c dx wedge dt` with `c !=0`**, for every `Q_(a,b,c)` and every
nonzero clock. Formal target forms are interpreted in the completion of
the common target germ and pulled back to the labelled formal source
neighborhoods; a globally defined constant pullback would in particular
have this forbidden local property. The obstruction order is exactly
`r`, `2r`, `3r`, or `4r`, on the four strata.

This includes `phi=t`, which changes the first source term, and every
`phi=t^r`. On the exceptional stratum the order is `4r` and the leading
coefficient is `kappa*a_r^4`. These maps retain the full compensated
collision; they are not arbitrary earlier-jet perturbations. In particular
this is a different family from the uncorrected pencil `Q+t^2` in the
ninth-source theorem, and no inclusion of all of that theorem's sources
is asserted.

The construction (7) changes `q` while leaving `w=t`. It is not a silent
reparametrization of an old density. For example the exact normalized
surface form below pulls back under (7) with the additional factor
`phi'(t)`; that factor must not be discarded.

## 6. The compatible unit density and exact source-coordinate escape

Near the common target, the two-form

```
Omega_0=(dC wedge dE)/(6(B+2))
```

is regular. It is also a legitimate formal target two-form. By (3), its
pullback under (7) is exactly

```
f_phi^*Omega_0=phi'(t)(1+chi_Q'(phi(t))x) dx wedge dt. (9)
```

For `phi=t`, the density
`rho(x,t)=1+chi_Q'(t)x` belongs to `1+t k[x][[t]]`; it is a formal unit.
Equation (4) says exactly that its weighted period vanishes. Thus a
nonconstant unit density is compatible even though a constant one is
excluded in the fixed coordinates.

More sharply, the formal source change

```
X=x+(chi_Q'(t)/2)x^2,              T=t               (10)
```

has `X_x=rho` and reduces to the identity modulo `t`. It is consequently
an automorphism of the `t`-adically complete ring `k[x][[t]]`: the inverse
is constructed uniquely coefficient by coefficient and is polynomial in
`X` at each coefficient. It restricts to isomorphisms of the labelled
formal source neighborhoods. In these new coordinates,

```
f_t^*Omega_0=dX wedge dT.
```

Equivalently the reexpressed formal map `f_t(x(X,T),T)` pulls back the
same target form to constant density. It retains the collision after
transporting its sections. This is an exact source-coordinate escape
from the fixed-volume conclusion, not a polynomial Keller construction.
On the exceptional stratum the source change starts at order `t^4`.
No claim that (10), its inverse, the formal compensator, or a target
primitive pair is polynomial follows.

The covariance explains the escape. For a source change `X=psi(x,t)` with
unit `a_i=psi_x(z_i,t)` and fixed `t`, new tangents are `T_i/a_i`. A
transported relation has coefficients `lambda_i*a_i`, and the densities
obey `J_old=(J_new composed with psi)*psi_x`. Therefore the full weighted
period is unchanged, but the sum of relation coefficients becomes
`sum_i lambda_i*a_i`. Under (10) this is `Pi+chi' Xi=0`, rather than the
old nonzero `Pi`. A source change preserving `dx wedge dt` preserves the
constant-density question. Dropping its Jacobian is the first failed
implication in an unrestricted coordinate-invariance claim.

## 7. Connection contract and the paper's transferable mechanism

The source of the connection is the actual common-motion equation; its
target is the entire descended two-form measurement system. The map is
differentiation, wedge contraction, and relation summation. It preserves
necessary density identities, forgets normal information beyond those
identities, and needs the moving section moment `Xi` and the source
Jacobian as sidecars. The cheapest positive check is the compiler's exact
determinant; the decisive hostile is (9)--(10).

The explicitly requested paper, Yinjie Li's
[The S-matrix conjecture, v1, Sections 6.2--6.3](https://arxiv.org/html/2608.29750v1),
retains multiple constrained columns against one global defect budget;
when two measurements share a column, their joint constraint matters.
The transferable research move is to preserve the full measurement
operator and its shared data. Here common target evaluation couples all
two-form slots, and the common-motion identity couples the constant period
to a moving first moment. This proof is algebraic. No nonnegative energy,
real ordering, or numerical S-matrix bound is imported over `C`, and the
paper's theorem is not independently certified by this calculation.

The remaining question is precise: which allowed polynomial source
changes preserve the fixed descended differential while leaving this
compensated family, and can a polynomial target primitive pair realize
an escaping density? The present formal source escape gives neither
polynomial termination nor that primitive pair.

## 8. Exact controls and reproduction

The standalone [source](../../04-computation/planar_jc_long_20260906_collision.py)
imports no repository producer or output. Its symbolic engine checks
(3), the surface equation, irreducibility of `F`, and a Bezout certificate
for `gcd(F,G)=1`. Its separate truncated-series engine uses exact Fraction
arithmetic in `Q[alpha]/F`, with rational controls in the rational subfield.
It solves all six collision equations through order six for four named
coefficient strata and computes the normalized moving relation, not just
the first live raw numerator.

The four controls have `ord chi=2,3,4,5` respectively. For every control
the source checks (4) through degree five, every target two-form basis
slot through degree five, the common `B` image, and the compatible unit
density at all three sections. Differentiation loses one known order,
which is explicitly respected. Five clock controls are `t`, `t^2`, `t^3`,
`t+t^2`, and `2t^2-t^3`; their first-live orders and leading coefficients
are exact. The source-coordinate density and covariance identities are
also checked symbolically. These are FINITE-EXACT controls for the stated
universe; the all-order and all-parameter claims rest on the proofs above.

Run from the repository root:

```bash
python3 -B 04-computation/planar_jc_long_20260906_collision.py > /tmp/collision_normal.out
python3 -B -O 04-computation/planar_jc_long_20260906_collision.py > /tmp/collision_optimized.out
cmp /tmp/collision_normal.out /tmp/collision_optimized.out
cmp /tmp/collision_normal.out 05-knowledge/results/planar_jc_long_20260906_collision.out
```

Normal and optimized output agree byte for byte in **157 always-active
gates**. The [frozen output](planar_jc_long_20260906_collision.out) records
the coefficient strata and the source digest. Raw LF-byte SHA-256 hashes:

```
source 7320aecbc79c7c53c0118889a4c307f53102f6aa488a35641faac16873bf1e82
output 069811ba7e42bc8b83d0a0bb76c589b0f673406f12ecfde8fa4cd83b96d62d6c
```

The [independent analytic and exact audit](planar_jc_long_20260906_collision_audit.md)
accepts the result with no mathematical repair. Its separate sparse
polynomial-ring Hensel computation and native algebraic field arithmetic
pass 133 always-active gates, with byte-identical normal and optimized
output. It additionally reconstructs the source-coordinate inverse by
Catalan coefficients. No theorem identifier has been allocated in this note.
