---
id: THM-3759
title: "Central-product profile with constant flanks: smoothness and Keller classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over an
  algebraically closed field of characteristic zero, Q=aX+bT+chi(XT) is
  smooth exactly when chi is constant and (a,b) is nonzero, or chi is
  nonconstant, exactly one flank is nonzero, and chi'(0)=0.  It has a
  polynomial constant-Jacobian mate exactly in the constant-profile case,
  when Q is affine linear.  The smooth nonlinear one-flank boundary has an
  all-degree nonterminating Euler repair chain; THM-3551 independently
  strengthens this boundary to no rational mate.  This closes one radial
  slice, not arbitrary quartic components or JC(2).
source: root + jc_radial_escape_probe + lrc14_cover_defect_bridge / 2026-08-23
audit: >
  PASS after repair.  Independent hostile audit rederived every smoothness
  boundary, required a proof rather than an assumption of the Hamiltonian
  kernel k[Q], corrected source-fibre degree to source-foliation degree,
  and updated the live THM-3551/3755/3757 comparisons.  The two assertion-free
  probes check the torus eliminant, both axis controls, ten forced Euler
  layers, bounded centralizer equality through degree ten, the quartic
  foliation degree and its special-line hostile, and the complete degree-three
  coefficient cube.  Normal and optimized runs byte-match both frozen outputs.
depends_on:
  - THM-3551-one-ray-planar-jacobian-mate-no-go
related:
  - THM-2071-quadratic-fiber-square-parity-gate
  - THM-2118-all-degree-cubic-faber-boundary-flux-coprimality
  - THM-3741-radial-two-charge-keller-component-classification
  - THM-3754-affine-variable-euclidean-descent-classification
  - THM-3755-composite-monomial-generic-fibre-residue-obstruction
  - THM-3757-pell-chebyshev-three-charge-hyperelliptic-obstruction-tower
  - THM-3765-normalized-three-consecutive-charge-radial-keller-classification
script: 04-computation/jc2_central_product_constant_flank_classification_thm3759.py
output: 05-knowledge/results/jc2_central_product_constant_flank_classification_thm3759.out
script_sha256: 666ce636afcf2fd49b3dbfe987815f7bf2021075b2a50ee2e1af075a3abe71cf
output_sha256: d5cdfaa061dde959b2cf1500fd24a0afcf5851f4f8c290bb15dcb508ae118a11
independent_script: 04-computation/jc2_central_product_constant_flank_classification_independent_audit_thm3759.py
independent_output: 05-knowledge/results/jc2_central_product_constant_flank_classification_independent_audit_thm3759.out
independent_script_sha256: c7744f9264dfd0f60b3fea3fcc6e88b79eaae50afd094264620cbc56216d7e5a
independent_output_sha256: c55df29c1bb50bb21208d5608d9befca5e5a53c2df1552d626b55a30e116cabe
hash_basis: raw working-tree bytes
---

# THM-3759 -- the constant-flank central-product slice

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be
algebraically closed of characteristic zero.  For `a,b in k` and
`chi in k[z]` put

```text
z=XT,                  Q(X,T)=aX+bT+chi(z).             (1)
```

Then `Q` has no critical point if and only if exactly one of the following
holds:

1. `chi` is constant and `(a,b)!=(0,0)`;
2. `chi` is nonconstant, exactly one of `a,b` is nonzero, and
   `chi'(0)=0`.

Moreover, there are `P in k[X,T]` and `c in k*` with

```text
J(P,Q)=P_X Q_T-P_T Q_X=c                              (2)
```

if and only if the first case holds.  There `Q` is a nonzero affine linear
coordinate.  In the second case the direct proof below excludes every
polynomial mate, while THM-3551 supplies the stronger inherited rational
no-mate conclusion.

## 1. Smoothness, including all axes

The gradient is

```text
Q_X=a+T chi'(z),                  Q_T=b+X chi'(z).      (3)
```

For constant `chi` this is the constant vector `(a,b)`.  Suppose `chi` is
nonconstant and `ab!=0`.  The polynomial

```text
R(z)=z chi'(z)^2-ab                                      (4)
```

is nonconstant, so algebraic closure gives a root `z_0`.  It is nonzero and
`chi'(z_0)!=0`.  Consequently

```text
X_0=-b/chi'(z_0),                  T_0=-a/chi'(z_0)      (5)
```

has product `z_0` and kills both entries of (3).

If `a!=0,b=0`, a critical point must have `X chi'(z)=0`.  A zero of `chi'`
leaves `Q_X=a`, so only `X=0` can cause trouble.  There `z=0` and
`Q_X=a+T chi'(0)`.  It never vanishes precisely when `chi'(0)=0`.  Swapping
`X,T` proves the other axis.  When `a=b=0` the origin is critical.  These
cases are exhaustive.

The same calculation has a useful larger form.  For

```text
Q=X phi(z)+chi(z)+T psi(z),             w=z phi psi,
```

the torus chart `(X,T)->(X,z)` eliminates a critical point to

```text
X^2=z psi/phi,             X chi'=-w'/phi,
E(z)=w'(z)^2-w(z)chi'(z)^2.                              (6)
```

Every root of `E` off `w=0` reconstructs a torus critical point.  The axes
and the divisor `w=0` are a mandatory sidecar.  Constant nonzero flanks give
`E=ab[ab-z(chi')^2]`, recovering (4).  THM-3757 succeeds with a variable
flank precisely by forcing the eliminant roots onto the profile boundary;
a nonzero constant flank has no such absorption mechanism.

## 2. The Hamiltonian kernel must be proved first

By typed duality it suffices to treat `a!=0,b=0`.  A determinant-one source
scaling normalizes `a` if desired, but retaining it records every constant.
Let

```text
D=J(-,Q),                       Q=aX+chi(XT),            (7)
```

and give `X,T` weights `+1,-1`.  The `chi(z)` part preserves weight and the
`aX` part raises it by one:

```text
J(P_r,aX)=-a partial_T P_r.                              (8)
```

If a polynomial `P` has top weight `M>=1` and `D(P)` has no weight
`M+1` term, then `partial_T P_M=0`, so `P_M=lambda X^M`.  The exact target
shear

```text
P -> P-(lambda/a^M)Q^M                                  (9)
```

preserves `D(P)` and lowers the top weight.  Repetition removes all positive
weights.  The next equation makes the weight-zero part constant, which may
also be removed.

This proves, rather than assumes, the exact kernel.  If `D(P)=0` and a
reduced negative part remained, let `-k` be its closest-to-zero weight and
write it as `T^kF_k(z)`.  The weight `-(k-1)` equation is

```text
a[zF_k'(z)+kF_k(z)]=0.                                  (10)
```

The Euler operator `z d/dz+k` is invertible on `k[z]` in characteristic
zero.  Thus `F_k=0`, a contradiction, and

```text
ker_(k[X,T]) J(-,Q)=k[Q].                               (11)
```

This order of proof is load-bearing: quotienting by `k[Q]` before (11) would
be circular.

## 3. The forced mate tail never terminates

Assume `D(P)=c in k*`.  After the shears (9), a reduced mate would have a
finite expansion

```text
P=sum_(k=1)^K T^k F_k(z).                               (12)
```

Direct differentiation gives

```text
J(T^kF,aX)      =-aT^(k-1)(kF+zF'),
J(T^kF,chi(z)) =-kT^kF chi'.                            (13)
```

The weight-zero equation and successive negative-weight equations are

```text
F_1=-c/a,
a[zF_(k+1)' +(k+1)F_(k+1)]=-kF_k chi',     k>=1.       (14)
```

For every `k`, the operator `z d/dz+(k+1)` multiplies `z^j` by the nonzero
scalar `j+k+1`, so (14) has a unique polynomial successor.  If `chi` is
nonconstant, `F_1!=0` and induction makes every successor nonzero.  A finite
tail cannot terminate: its last term leaves the uncancelled debt
`-KT^K F_K chi'`.  Hence no polynomial mate exists.

If `chi` is constant, `P=-cT/a` is a mate when `a!=0`, and `P=cX/b` is a
mate when `b!=0`.  This proves the polynomial classification.

## 4. Source-foliation degree and the quartic control

Let `d=deg chi>=1`.  In an invertible affine source chart

```text
X=alpha u+beta v+x_0,             T=gamma u+delta v+t_0,
alpha delta-beta gamma !=0,
```

measure `deg_v Q` over the base ring `k[u]`.  If `beta delta!=0`, the leading
term of `(XT)^d` gives degree `2d`.  If `beta delta=0`, exactly one is
nonzero, and invertibility makes the coefficient of `v^d` a nonzero
polynomial in `u`.  Thus

```text
min over affine source foliations deg_v Q=d.            (15)
```

This is a foliation statement, not a claim about every special individual
line.  For the smooth control

```text
Q=X+(XT)^2+(XT)^4,                                     (16)
```

generic mixed directions have degree eight and an axis-parallel foliation
has degree four, but the special line `T=0` drops to degree one.  The missing
base coordinate is essential.

Every smooth nonlinear member is visibly a noncoordinate: for `b=0`,

```text
Q-chi(0)=X[a+T(chi(z)-chi(0))/z],                       (17)
```

and both factors are nonunits.  This also explains the comparison with the
degree-two and degree-three source-fibre closures.  THM-2071 and THM-2118,
combined with (17), independently reject the `d=2,3` cases.  The first
general source-foliation range not covered by those theorems is `d=4`;
THM-3759 closes only the special constant-flank central-product slice there.

## 5. Exact controls and scope

Reproduce with

```bash
python3 -B 04-computation/jc2_central_product_constant_flank_classification_thm3759.py
python3 -B -O 04-computation/jc2_central_product_constant_flank_classification_thm3759.py
python3 -B 04-computation/jc2_central_product_constant_flank_classification_independent_audit_thm3759.py
python3 -B -O 04-computation/jc2_central_product_constant_flank_classification_independent_audit_thm3759.py
```

The first companion derives the general torus eliminant, checks four
nonterminating repair chains, rejects bounded mates, and exhausts all `243`
degree-three coefficient tuples over `{-1,0,1}` with zero smoothness
mismatch.  The independent companion checks interior eliminants through
profile degree seven, both axes, ten repair layers, centralizer nullities
through total degree ten, foliation degrees `8,4`, and the hostile special
line of degree one.  The formulas above carry the arbitrary-degree proof.

THM-3741 has no central profile; THM-3755 treats the different exponent
boundary `X+F(X^mT)` for `m>=2`; THM-3765 closes the normalized variable-flank
three-charge ansatz.  The present theorem neither transports across an
arbitrary affine source change nor closes arbitrary quartic components,
nonradial supports, two interacting nonconstant Cohn factors, or JC(2).
**QED.**
