---
id: THM-3949
title: "Coprime one-variable internal factor splits are reducible or multi-ended"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Let p0=P and p1=P+A(Y)B(Y), with nonconstant coprime nonassociate
  one-variable factors assigned to the two complementary cube-difference
  rows.  On P=h^2 the common discriminant has the hidden monic cubic
  C=(B-A)AB+h^2[(1-omega^2)B-(1-omega)A]-2h^3.  Its norm is H, and after
  inverting P/h its curve is exactly the discriminant curve.  The Newton
  polygon has two distinct infinity slopes unless deg A=deg B with unequal
  leaders; in that last case its sole-edge residual cubic has no quadratic
  term and nonzero constant, hence cannot have one triple root.  Therefore
  the full discriminant is reducible or its normalization has at least two
  places on the standard line at infinity.  No arbitrary-line obstruction is
  claimed.
source: jc-zero-debt / all-degree one-factor Newton boundary analysis, 2026-08-24
depends_on: []
related:
  - THM-3942-affine-linear-double-torus-factor-split-one-place-obstruction
  - THM-3944-repeated-factor-double-torus-one-place-square-conductor-collapse
  - THM-3946-affine-internal-factor-split-two-end-conductor-collision-dichotomy
  - THM-3947-scalar-weighted-repeated-square-split-trichotomy
script: 04-computation/jc2_coprime_one_variable_internal_split_newton_thm3949.py
output: 05-knowledge/results/jc2_coprime_one_variable_internal_split_newton_thm3949.out
script_sha256: 41f35a90d19dad6af6dad807bc4ef8c6928d52e35677c46520336d00250c1c83
output_sha256: 6490b55d1a4b97b93421260fda896cd930201ebc31a433875414b36b7235a2ab
semantic_sha256: 4c7a0c10ed7397aacb60ddda5d0d63811bdf744688329777934a2e703e92e49b
hash_basis: raw LF bytes
---

# THM-3949 -- one-variable internal factors always leave at least two ends

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**
Work over an algebraically closed field `k` of characteristic zero.  Fix

```text
omega^2+omega+1=0.                                      (1)
```

Let `A,B in k[Y]` be nonconstant, coprime, and nonassociate.  Put

```text
p0=P,                         p1=P+A B,
L1=p1-omega p0,               L2=p1-omega^2 p0,          (2)
q1-q0=2 A L1,                 q1+q0=2 B L2.              (3)
```

Then the common discriminant

```text
H=q0^2-4p0^3=q1^2-4p1^3                              (4)
```

satisfies the following dichotomy:

```text
H is reducible in k[P,Y],
or
the normalization of H=0 has at least two places on Z=0
in the standard projective closure [P:Y:Z].             (5)
```

In particular, no curve in this one-variable one-factor grammar is an
irreducible polynomially parametrized one-place branch in the displayed
affine chart.  The conclusion is deliberately chart-specific: it does **not**
say that every projective line has at least two preimages, and hence is not an
arbitrary-line no-go.

The reciprocal scalar gauge is already included: replacing
`(A,B)` by `(lambda A,lambda^(-1)B)` preserves `AB` and gives the general
reciprocal weighting of `(3)`.

## 1. The hidden cubic and its exact localization

Set

```text
D=(1-omega^2)B-(1-omega)A,
E=(B-A)AB.                                               (6)
```

Solving `(3)` gives

```text
q0=B L2-A L1=E+D P.                                     (7)
```

Moreover,

```text
(q1-q0)(q1+q0)=4AB L1L2=4(p1^3-p0^3),                  (8)
```

which proves `(4)`.  On the square-root chart `P=h^2`, the branch
`q0=2h^3` is the cubic

```text
C(Y,h)=E+D h^2-2h^3=0.                                 (9)
```

It is not an auxiliary approximation.  There is an exact norm identity

```text
C(Y,h) C(Y,-h)=H(h^2,Y),                               (10)
```

and an exact localization isomorphism

```text
k[P,Y]/(H)[P^(-1)]  ~=  k[h,Y]/(C)[h^(-1)],            (11)
P |-> h^2,
h |-> q0/(2P).                                          (12)
```

Indeed `(4)` gives `(q0/(2P))^2=P`, while `(7)` and `(9)` give the
two inverse formulas.  Before quotienting, the exact residual identities are

```text
(q0/(2P))^2-P=H/(4P^2),
C(Y,q0/(2P))=-E H/(4P^3),
q0(h^2,Y)/(2h^2)-h=C/(2h^2).                           (13)
```

Since `A,B` are nonzero and nonassociate,
`E!=0`; hence neither `P` nor `h` is an entire component.  Consequently, if
`H` is irreducible, then `C` is irreducible and `(11)` identifies their
function fields.  A place has a pole in `(P,Y)` exactly when it has a pole in
`(h,Y)`, because `P=h^2`.  Thus their standard infinity-place ledgers agree.

One can also read irreducibility directly from `(10)`: a nontrivial factor of
`C` has a nonconstant invariant norm in `k[Y,h^2]`, and would factor `H`.

## 2. Newton polygon when the degrees differ

Write

```text
a=deg A,                  b=deg B.                     (14)
```

At `Y=infinity`, use the uniformizer `x=1/Y`.  As a cubic in `h`, `(9)` has
coefficient-valuation points

```text
(0,-deg E),                 (2,-deg D),                 (3,0). (15)
```

There is no `h` term.  After division by `-2`, equation `(9)` is monic in
`h`, so its projection to `A1_Y` is finite.  In particular, `h` cannot have a
pole over a finite value of `Y`; the Newton polygon over `Y=infinity`
accounts for every standard infinity place.

If `b>a`, then

```text
deg D=b,              deg E=a+2b,                       (16)
```

so the lower Newton polygon has two edges of horizontal widths two and one,
with pole slopes

```text
(a+b)/2,                       b.                       (17)
```

If `a>b`, the symmetric formulas are

```text
deg D=a,              deg E=2a+b,
slopes (a+b)/2,       a.                                (18)
```

In either case the slopes are distinct.  The width-two edge may represent
two unramified places or one ramified quadratic place, depending on parity
and its residual polynomial.  The width-one edge is a different local
factor because its pole slope is different.  Hence there are at least two
places over infinity; no parity assumption is being hidden.

## 3. Equal degree with equal leading coefficients

Now let `a=b=n` and write the leading coefficients as `alpha,beta`.  Suppose
first that `alpha=beta`.  Coprimality and nonassociateness imply `B-A!=0`; set

```text
r=deg(B-A)<n.                                           (19)
```

The leading coefficient of `D` is

```text
(omega-omega^2)alpha!=0,                               (20)
```

and therefore

```text
deg D=n,                 deg E=2n+r.                    (21)
```

The two Newton slopes are now

```text
(n+r)/2,                       n.                       (22)
```

They are again distinct because `r<n`, so the same two-place conclusion
follows.

## 4. Equal degree with unequal leading coefficients

It remains to take `a=b=n` and `alpha!=beta`.  Then

```text
deg E=3n,                         deg D<=n.              (23)
```

The Newton polygon has one edge of slope `n`; the `h^2` point is on or above
it.  A single slope is necessary for one infinity place, but not sufficient.
Put

```text
kappa=Y^n/h.                                            (24)
```

The edge residual polynomial is

```text
R(kappa)=alpha beta(beta-alpha) kappa^3
 +[(1-omega^2)beta-(1-omega)alpha] kappa-2.             (25)
```

Its leading coefficient and constant term are nonzero, and it has no
`kappa^2` coefficient.  Over `k`, a cubic with only one distinct root would
be a scalar multiple of `(kappa-c)^3`.  The missing quadratic coefficient
would force `c=0`, contradicting the constant `-2`.  Thus `(25)` has at least
two distinct residual roots.  Newton--Hensel factorization separates the
corresponding local factors, so this last case also has at least two infinity
places.

Together Sections 2--4 exhaust all degree and leading-coefficient cases and
prove `(5)`.

## 5. The smallest nonlinear collision still has two directions

The first nonlinear equal-degree seam already occurs in degree two.  Take

```text
A=Y^2+1,                     B=-omega^2 Y^2+1.           (26)
```

These polynomials are coprime because `B-A=omega Y^2` and
`gcd(Y^2+1,Y)=1`.  Here

```text
alpha=1,                  beta=-omega^2,
alpha beta(beta-alpha)=-1,
(1-omega^2)beta-(1-omega)alpha=3omega.                  (27)
```

The residual cubic is not generic: it deliberately has a double root,

```text
R(kappa)=-kappa^3+3omega kappa-2
        =-(kappa-omega^2)^2(kappa+2omega^2).             (28)
```

But the remaining simple root is distinct, so even this most economical
collision has at least two leading infinity directions.  It is a hostile
control against replacing the residual-root argument by a genericity claim.

## 6. Scope and reproduction

Run

```bash
python3 04-computation/jc2_coprime_one_variable_internal_split_newton_thm3949.py
python3 -O 04-computation/jc2_coprime_one_variable_internal_split_newton_thm3949.py
```

The companion verifies the two torus rows, `(7)-(13)`, every degree/leading
coefficient ledger, the coefficient-ideal impossibility of a triple residual
root, and the exact quadratic hostile `(26)-(28)`.

The next live internal-split coordinate must leave at least one hypothesis of
this theorem: factors depending on both target variables, a distribution
across more than one `Li`, a gcd/nonmaximal-order overlap, or a different
target chart.  The theorem does not rule out any of those, nor does it prove
`JC(2)`.
