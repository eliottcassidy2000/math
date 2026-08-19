---
id: THM-3551
title: "One-invariant ray and residue planar Jacobian mate no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY SCOPE-AUDITED.  Over a
  characteristic-zero field, three all-degree one-invariant first-coordinate
  families have no polynomial Jacobian mate unless their nonlinear datum is
  constant: the multiplicative primitive ray x*phi(x^p y^q), the diagonal
  additive family x+h(xy), and the thickened additive ray x+h(x^p y^q) for
  primitive q>=2.  The first and third are killed by the unique
  constant-producing weight/residue sector; the diagonal family has the
  stronger obstruction that its generic response form dv/(z-h(v)) has
  nonzero residues and therefore no rational primitive.  This includes
  gradient-unimodular, positive-area and arbitrarily high-genus hostiles.  It
  is a classification of these stated families, not JC(2) or a general
  Newton-polygon theorem.
source: root-2026-08-18-planar-jacobian-counterexample-hostiles
depends_on: []
related:
  - THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate
  - THM-2118-all-degree-cubic-faber-boundary-flux-coprimality
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3544-planar-keller-target-pencil-total-degree-six-floor
script: 04-computation/jc_one_ray_mate_no_go_thm3551.py
output: 05-knowledge/results/jc_one_ray_mate_no_go_thm3551.out
script_sha256: 1e0fd05c85e3eb21bccf6dd85ea029277c5df255a059f420b1425f96b3b519f5
output_sha256: e14c9270a87f6438891e16b85f789882f93ce62803b1d448895ed3c725d30e28
semantic_sha256: cf050ca64d5054a65a01da346489de214f695d67d0f6b98b9470d61567a39db5
hash_basis: LF-normalized bytes
---

# THM-3551 -- one-invariant ray and residue mate no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY SCOPE-AUDITED.**

Let `k` be a field of characteristic zero and write

```text
D_P(Q)=Jac(P,Q)=P_x Q_y-P_y Q_x.                         (1)
```

This theorem closes three tempting counterexample-first-coordinate families.
The common failure is not low degree.  Each family retains only one toric
invariant, so the constant response lives in one exact sector and is forced
to run into either a terminal leading coefficient or a nonzero residue.

## 1. Three no-mate statements

Let `p,q` be coprime positive integers.

### A. Multiplicative primitive ray

For `phi in k[T]`, put

```text
P=x phi(x^p y^q).                                        (2)
```

There is a `Q in k[x,y]` and `kappa in k*` with `D_P(Q)=kappa` if and only if
`phi=a in k*` is constant.  In that case all mates are

```text
Q=(kappa/a)y+H(P),             H in k[T].                (3)
```

### B. Diagonal additive ray

For `h in k[T]`, put

```text
P=x+h(xy).                                                (4)
```

If `h` is nonconstant, `P` has no mate even in the rational function field:
there is no `Q in k(x,y)` with `D_P(Q) in k*`.  If `h` is constant, all
polynomial mates have the evident form

```text
Q=kappa y+H(P).                                           (5)
```

### C. Thickened primitive ray

Assume `q>=2` and put

```text
P=x+h(x^p y^q).                                          (6)
```

If `h` is nonconstant, `P` has no polynomial Jacobian mate.  If `h` is
constant, `(5)` again gives every mate.

The conclusions concern the displayed source chart.  An affine source
change is lawful, but replacing `x^p y^q` by a nonmonomial polynomial or by
two independent invariants leaves the theorem.

## 2. The multiplicative weight sector

Give `x,y` weights `q,-p`, and write `T=x^p y^q`.  Then `P=x phi(T)` has
weight `q`, and `D_P` sends weight `m` to weight `m+p`.  Hence only the
weight `-p` component of `Q` can contribute to a scalar.  Coprimality gives

```text
q i-p j=-p
iff (i,j)=(pn,qn+1),
```

so that component is exactly

```text
Q_(-p)=y f(T),                  f in k[T].               (7)
```

Direct differentiation gives

```text
D_P(yf(T))
 =phi f+pT phi' f+qT phi f'.                             (8)
```

If `d=deg(phi)>=1` and `N=deg(f)`, the coefficient of `T^(d+N)` in `(8)` is

```text
(1+pd+qN) lc(phi)lc(f),                                  (9)
```

which is nonzero in characteristic zero.  Therefore `(8)` cannot be a
nonzero constant.  No other weight sector can cancel `(9)`.  If `phi=a` is
constant, `(1)` is `a Q_y=kappa`, which gives `(3)`.  This proves A.

Polynomiality is load-bearing.  For example

```text
phi(T)=1+T+T^2,
P=x phi(xy),                 Q_rat=y/phi(xy)             (10)
```

has `Jac(P,Q_rat)=1` and the rational inverse

```text
x=u/(1+uv+(uv)^2),          y=v(1+uv+(uv)^2).            (11)
```

The construction is a symplectic automorphism only after deleting
`phi(xy)=0`; the pole cannot be cleared polynomially.

## 3. The diagonal residue obstruction

For `(4)`, introduce rational coordinates

```text
z=P,                         v=xy.                       (12)
```

Then

```text
x=z-h(v),                    y=v/(z-h(v)),
Jac_(x,y)(z,v)=x.                                        (13)
```

Writing a rational candidate as `Q_tilde(z,v)`, the equation
`D_P(Q)=kappa` becomes

```text
partial_v Q_tilde = kappa/(z-h(v)).                      (14)
```

If `h` is nonconstant, `z-h(v)` is separable over `k(z)`: a common divisor
with `h'(v)` would divide the transcendental coefficient `z`.  At every root
`r` of `z-h(v)`, the right side of `(14)` has the nonzero residue

```text
-kappa/h'(r).                                            (15)
```

The derivative of a rational function has residue zero at every finite pole.
Thus `(14)` has no rational solution.  This proves B and is stronger than a
failure of polynomial termination.

The submersion gate is strictly weaker.  The gradient of `(4)` is unimodular
exactly when `h'(0)=0`: if `h'(0)!=0`, then
`(x,y)=(0,-1/h'(0))` is critical; if `h'(0)=0`, the axes are safe and a torus
zero of `P_y=xh'(xy)` forces `P_x=1`.  Hence

```text
P=x+(xy)^2+(xy)^4                                         (16)
```

is a gradient-unimodular, degree-four-fibre hostile with no rational mate.

## 4. The thickened residue sector

Take `(6)` with `q>=2`.  Decompose `Q` by its `y` exponent modulo `q`.
Both pieces of `D_P` lower that residue by one, so the scalar response comes
entirely from the residue-one part `Q_1`.

Suppose `Q_1` is nonzero and let `A(x)y^N` be its term of largest `y`
exponent.  Then

```text
N=1 mod q.                                                (17)
```

If `a_d T^d` is the top term of `h`, the unique coefficient at the largest
output exponent `y^(qd+N-1)` is

```text
a_d d x^(pd-1) [p N A-q x A'].                          (18)
```

It must vanish.  A nonzero polynomial solution of

```text
q x A'=p N A                                             (19)
```

is a scalar monomial of degree `pN/q`, so `q|pN`.  Since `(p,q)=1`, this
forces `q|N`, contradicting `(17)`.  Thus no finite residue-one response
exists, proving C.

For example, with `T=xy^2`,

```text
P=x+T^2+T^3=x+x^2y^4+x^3y^6                             (20)
```

has positive Newton area and a unit gradient ideal.  Its forced formal
response cancels arbitrarily long finite prefixes, but `(18)` proves it can
never terminate.  A deep jet is therefore not evidence of a polynomial mate.

## 5. Failure of the apparent hyperelliptic escape

The family

```text
P=x+(xy)^a+lambda x^(2a-1)y^(2a),       a>=2            (21)
```

looks like a three-weight, positive-genus scaffold.  In fact, with

```text
T=x^(a-1)y^a,
```

it is exactly

```text
P=x(1+T+lambda T^2).                                    (22)
```

Statement A kills every value of `lambda`.  When `lambda!=1/4`, the quadratic
in `(22)` is squarefree and the gradient is unimodular; when `lambda=1/4`, it
has a repeated root and the corresponding torus curve is critical.  For
generic target value and `lambda!=0,1/4`, the fibre is a cyclic cover of
genus `a-1`.  Thus arbitrarily high generic-fibre genus and absence of affine
critical points still do not create a second Hamiltonian response channel.

The first failed implication is precise: more powers of one primitive
invariant are not more independent faces.

## 6. Boundary and next construction scale

THM-2045 supplies the closest binomial factorized wall; THM-3418 closes a
different one-top-monomial fibre family.  The present theorem permits
arbitrarily many terms but only one invariant.  It does not cover:

1. two nonparallel primitive monomials `T,S` with nonzero exponent
   determinant;
2. mixed coefficients depending on both `T` and `S`;
3. a non-invariant graph section of a higher-dimensional Keller map; or
4. a compactified fibre whose response one-form has a nontrivial but
   potentially cancellable pole divisor.

Those are the honest next construction lanes.  A proposed candidate should
be refactored into primitive invariants before degree, genus, or support size
is interpreted.

## 7. Exact verification

Reproduce with

```bash
python3 04-computation/jc_one_ray_mate_no_go_thm3551.py
python3 -O 04-computation/jc_one_ray_mate_no_go_thm3551.py
```

The companion checks nineteen primitive multiplicative sectors through mate
degree six, the rational symplectic cage and inverse, the diagonal residue
coordinate and bounded response systems through degree ten, twenty residue
top cells, the gradient hostiles, the hyperelliptic refactorization, and a
positive affine/shear control.  Ordinary and optimized transcripts agree
byte-for-byte.

**QED.**
