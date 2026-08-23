---
id: THM-3755
title: "Composite-monomial generic-fibre residue obstruction"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT AUDIT.  Let m>=2 and
  Q=X+F(X^mT), where F(0)=0.  Every such Q is nonsingular.  It is a
  coordinate exactly when F=0.  If deg F>=2, its generic-fibre response
  differential has a nonzero residue, so Q has no rational Jacobian mate.
  If deg F=1, an explicit rational mate exists but no polynomial mate does.
  Thus this missing q=1 thickened-ray boundary contains smooth reducible
  multi-charge components but no planar Keller pair.
source: root + jc_sparse_direct_search / 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The generic-fibre chart, residue sign,
  differential-field kernel, linear rational primitive, polynomial descent,
  smoothness identity, and reducible-fibre certificate have been rederived.
  The exact companion checks the residue formula, nonlinear invariant,
  boundary mates, a three-charge hostile, and bounded polynomial-mate
  systems.  Normal and optimized output agree with the frozen transcript;
  independent hostile audit remains due.
depends_on: []
related:
  - THM-3551-one-invariant-ray-and-residue-planar-jacobian-mate-no-go
  - THM-3598-danielewski-rational-exact-polar-graph-family-and-classification
  - THM-3716-monomial-broughton-hamiltonian-obstruction-family
  - THM-3741-radial-two-charge-keller-component-classification
  - THM-3754-affine-variable-euclidean-descent-classification
script: 04-computation/jc2_composite_monomial_generic_fibre_residue_thm3755.py
output: 05-knowledge/results/jc2_composite_monomial_generic_fibre_residue_thm3755.out
script_sha256: 0c842cb451a844ab2a96b58677cb34ef8c2acbb04ea8365c09ae28e52d12fbaa
output_sha256: 2d1bcb51d36f582f1c1f908b2b7fb84d44a8c97cbf3aaadd116c82767d2b8261
semantic_sha256: 400cd8344235d6e247270abbf175fceb0efd7444eb001ea023bf5f85a90dd6cb
hash_basis: raw LF bytes
---

# THM-3755 -- the missing thickened ray is residue-obstructed

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT AUDIT.**  THM-3551 closes
the additive one-invariant family `x+h(x^p y^q)` rationally at `p=q=1` and
polynomially for primitive `q>=2`.  The missing primitive exponent boundary
is `q=1,p=m>=2`.  It behaves differently at its linear edge, but it too is
empty as a polynomial Keller-component search space.

Let `k` be an algebraically closed field of characteristic zero, let `m>=2`,
take `F in k[W]` with `F(0)=0`, and put

```text
W=X^m T,                         Q=X+F(W).              (1)
```

Then:

1. `Q` has no critical point in `k^2` for every `F`;
2. `Q` is a coordinate polynomial if and only if `F=0`;
3. if `deg F>=2`, there is no `P in k(X,T)` with
   `J(P,Q)=c in k*`;
4. if `F=aW`, `a!=0`, rational mates exist, but no polynomial mate exists;
5. consequently a polynomial constant-Jacobian mate exists exactly on the
   coordinate boundary `F=0`.

The rational boundary in item 4 is explicit.  For any `c in k*`, all of its
rational mates are

```text
P=-c X^(1-m)/[a(m-1)]+H(Q),          H in k(Q).         (2)
```

For `F=0`, all rational mates are `P=-cT+H(X)`, and the polynomial ones are
obtained exactly when `H in k[X]`.

## 1. Every member is smooth, but every nonzero member is noncoordinate

Differentiate `(1)`:

```text
Q_X=1+mX^(m-1)T F'(W),             Q_T=X^m F'(W),
X(Q_X-1)=mT Q_T.                                      (3)
```

On `X=0`, the assumption `m>=2` gives `Q_X=1`.  If `X!=0` and `Q_T=0`,
then `F'(W)=0`, and the first equation again gives `Q_X=1`.  Thus the two
partial derivatives never vanish simultaneously.  The exponent boundary is
sharp for this axis argument: at `m=1`, the term discarded on `X=0` need not
vanish, and THM-3551 uses a different diagonal residue chart.

If `F!=0`, write `F(W)=WG(W)`.  Then

```text
Q=X[1+X^(m-1)T G(X^mT)].                              (4)
```

Both displayed factors are nonunits.  Hence the zero fibre is reducible and
`k[X,T]/(Q)` is not a domain.  A coordinate has quotient ring isomorphic to
`k[T]`, so `(4)` proves that `Q` is not a coordinate.  Conversely `F=0`
gives `Q=X`.

Thus this family supplies smooth noncoordinates before any mate equation is
solved.  Smoothness and irreducibility are genuinely different gates here.

## 2. The generic fibre turns the mate equation into one rational differential

The monomial change is birational because its second exponent is one:

```text
k(X,T)=k(X,W),                   T=W/X^m,
J_(X,T)(X,W)=X^m.                                      (5)
```

Let `L=Q=X+F(W)` and write a rational candidate in the `(X,W)` chart.  The
chain rule gives

```text
J_(X,T)(P,Q)=X^m(F'(W)P_X-P_W).                        (6)
```

Holding `L` fixed means `dX/dW=-F'(W)`.  Therefore `J(P,Q)=c` is equivalent
over `K=k(L)` to

```text
dP/dW |_L = -c/[L-F(W)]^m.                             (7)
```

The source-to-target map in this reduction is

```text
source       polynomial component Q=X+F(X^mT)
target       rational generic Q-fibre over K=k(L)
map          (X,T) -> (L,W)=(X+F(X^mT),X^mT)
preserved    existence of a rational, hence polynomial, Jacobian mate
destroyed    polynomial lattice and the special divisors X=0, infinity
sidecar      residues of the response differential and X-order of a primitive
cheap test   D_F^(m-1)(1/F') followed by the linear-profile boundary
```

The rational relaxation loses regularity information, but a nonzero residue
already obstructs exactness before that lost information can matter.

## 3. The higher-pole residue is a differential-field invariant

Suppose `F` is nonconstant and define on `k(W)`

```text
D_F=(1/F'(W)) d/dW,                 D_F(F)=1.           (8)
```

The polynomial `F(W)-L` is generically separable.  At one of its roots
`alpha`, use the local parameter `u=F(W)-L`.  Since

```text
dW=[1/F'(W(u))]du,                  L-F(W)=-u,          (9)
```

the residue of the response differential in `(7)` is

```text
res_alpha {-c dW/[L-F(W)]^m}
 =(-1)^(m+1)c/(m-1)!
   [D_F^(m-1)(1/F')](alpha).                           (10)
```

This is the same higher-pole classifier that appears in THM-3598, now on a
planar generic fibre.  A derivative of a rational function has zero residue
at every finite pole.  Hence a rational solution of `(7)` would force the
right side of `(10)` to vanish at every generic root.  A generic root is
transcendental over `k`; therefore the rational function itself must vanish:

```text
D_F^(m-1)(1/F')=0.                                    (11)
```

The kernel in `(11)` is exact, not merely bounded.  For every `r>=1`,

```text
ker_(k(W))(D_F^r)=span_k{1,F,...,F^(r-1)}.             (12)
```

Indeed, the constants of `D_F` are `k`.  If `D_F^r A=0`, then
`D_F^(r-1)A` is constant.  Subtract the corresponding multiple of
`F^(r-1)/(r-1)!` and induct.  Applying `(12)` to `(11)` gives

```text
1/F'=G(F),                         deg G<=m-2.          (13)
```

But then the polynomial identity

```text
F'(W)G(F(W))=1                                       (14)
```

forces both factors to be nonzero constants.  Thus `F` is linear.  We have
proved the stronger statement

```text
deg F>=2  =>  Q has no rational constant-Jacobian mate. (15)
```

This explains the uniform augmented-rank gap seen in finite mate systems:
it is the shadow of one nonzero de Rham residue class, not a delayed
coefficient accident.

## 4. The linear profile is rationally exact but polynomially impossible

Let `F(W)=aW`, `a!=0`.  Direct differentiation verifies

```text
J(-c X^(1-m)/[a(m-1)], X+aX^mT)=c.                   (16)
```

In `(7)`, two primitives differ by an element of the constant field `k(L)`,
which proves the complete rational formula `(2)`.  The pole at `X=0` cannot
be removed polynomially.  One can see this without guessing its order.
Write a hypothetical polynomial mate as

```text
P=sum_(j=0)^N p_j(X)T^j.                              (17)
```

If `N>=1`, the top `T` coefficient in `J(P,Q)=c` is

```text
aX^m p_N'-Na mX^(m-1)p_N=0,
(p_N/(aX^m)^N)'=0.                                   (18)
```

Thus `p_N=lambda(aX^m)^N`.  Subtracting `lambda Q^N` preserves the
Jacobian and lowers the `T` degree.  Descent ends with `P=p(X)`, where

```text
J(p(X),Q)=aX^m p'(X),                                 (19)
```

which cannot be a nonzero constant.  This is the exact `g=aX^m` boundary
of THM-3754 and the `n=1` monomial boundary of THM-3716.  It also shows why
the rational conclusion in `(15)` cannot be extended to linear `F`.

If `F=0`, then `Q=X` and `J(P,X)=-P_T`; the positive mate `P=-cT` and the
classification stated after `(2)` follow immediately.  Equations
`(15)--(19)` prove every mate claim.

## 5. A genuine three-charge hostile and the next boundary

Give `X,T` Euler weights `+1,-1`.  The smallest clean member of this family
that displays three distinct charges is

```text
m=3, F(W)=W+W^2,
Q=X+X^3T+X^6T^2,                  charges 1,2,4.        (20)
```

It is smooth, has the reducible zero fibre

```text
Q=X(1+X^2T+X^5T^2),                                  (21)
```

is nonlinear in both variables, is not affine in `T`, and has no rational
mate.  Thus merely adding a third Euler-charge sector crosses both the
radial wall of THM-3741 and the affine-variable wall of THM-3754, but it
still does not create an independent response channel: all nonlinear terms
remain powers of the single composite invariant `W`.

The lost-information ledger points to the next honest construction scale.
Replace `F(W)` by dependence on two monomials with nonzero exponent
determinant, or perturb the coefficient of one power so that `(X,T)->(L,W)`
no longer makes the generic fibre rational with a one-variable response
form.  A third charge is useful only when it carries a second invariant.

## 6. Exact controls

Reproduce with

```bash
python3 -B 04-computation/jc2_composite_monomial_generic_fibre_residue_thm3755.py
python3 -B -O 04-computation/jc2_composite_monomial_generic_fibre_residue_thm3755.py
```

The assertion-free companion checks the monomial-chart and smoothness
identities for `m=2,...,8`; 35 smooth and reducible-fibre controls; the
signed residue formula in 18 exact higher-pole cases; 390 nonlinear
`D_F`-invariant controls from a frozen coefficient cube; the rational linear
primitive for `m=2,...,9`; the coordinate boundary; and 55 full
constant-Jacobian linear systems through total mate degree ten.  Every
bounded obstruction has augmented-rank gap one.  These computations are
hostile controls; equations `(3)--(19)` prove the arbitrary-degree theorem.
**QED.**
