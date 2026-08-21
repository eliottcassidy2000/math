---
id: THM-3567
title: "Separated rational Keller maximal-observable nodal-completion no-go"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.  Let f have squarefree
  derivative and pairwise distinct critical values over an algebraically
  closed characteristic-zero field.  The separated rational plane pair
  P=f(x), Q=y/f'(x) has unit Jacobian and degree-deg(f) collisions over
  regular values.  Its full-field polynomial-observable intersection
  k(P,Q) intersect k[x,y] is exactly
  k[P,A,B]/(A^2-Delta(P)B), with A=Delta(P)Q and
  B=Delta(P)Q^2.  This surface has one ordinary node per critical value.
  The induced polynomial map is finite etale of degree deg(f) over
  Delta!=0, is not globally etale, and for deg(f)>=3 contracts regular
  companion lines and is not quasi-finite.  This is distinct from
  THM-3561's polynomial-target-ring intersection and does not exclude its
  smooth nonseparated Danielewski completion.  No planar polynomial
  counterexample or conclusion about JC(2) is claimed.
source: codex-2026-08-20-frontier-overnight
depends_on: []
related:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3554-punctured-kummer-collision-surface-normal-form
  - THM-3559-affine-target-coordinate-pullback-no-go
  - THM-3560-jelonek-euler-gate-monomial-target-shear-no-go
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
script: 04-computation/jc_separated_rational_keller_observable_suspension_thm3567.py
output: 05-knowledge/results/jc_separated_rational_keller_observable_suspension_thm3567.out
script_sha256: 7c10c13ac53255643553ee3645b5a68ee916292b8feace4edb120642ff430f29
output_sha256: 8f290240535fe3ad5d82298a5110b5dd683a152dec64f31e49258eff4463380b
hash_basis: LF-normalized bytes
---

# THM-3567 -- separated rational Keller full-field observables form a nodal suspension

**PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.**

This theorem concerns a full rational target **field** intersection.  That
universe is deliberately different from the polynomial target **ring**
intersection in THM-3561.  The distinction is recorded explicitly in
Section 5.

## 1. Statement and rational Keller collision

Let `k` be an algebraically closed field of characteristic zero.  Let

```text
f in k[x],                         d=deg(f)>=2,          (1)
```

and assume:

1. `f'` is squarefree; and
2. the values `f(alpha)` at distinct roots `alpha` of `f'` are distinct.

Put

```text
P=f(x),                    Q=y/f'(x)                    (2)
```

inside `k(x,y)`.  On the open plane `f'(x)!=0`,

```text
Jac(P,Q)=f'(x) * 1/f'(x)=1.                             (3)
```

If `p0` is a regular value of `f`, then `f(x)=p0` has `d` distinct roots
`x1,...,xd`, and

```text
(x1,0),...,(xd,0)  ->  (p0,0).                         (4)
```

Thus `(2)` is a genuine rational Keller map with a `d`-fold collision.

Let the squarefree critical-value polynomial be

```text
Delta(T)=product_(f'(alpha)=0) (T-f(alpha)).            (5)
```

Multiplication of `Delta` by a nonzero scalar changes none of the claims.
Because every critical point is simple and the critical values are distinct,

```text
Delta(f(x))=f'(x)^2 H(x)                               (6)
```

for a polynomial `H in k[x]` coprime to `f'`.

Define the ambient target field and its full polynomial-observable
intersection by

```text
K=k(P,Q) subset k(x,y),
R_full=K intersect k[x,y].                             (7)
```

Then

```text
A=Delta(P)Q=f'(x)H(x)y,
B=Delta(P)Q^2=H(x)y^2                                 (8)
```

belong to `R_full`, and

```text
R_full=k[P,A,B]
      ~=k[P,A,B]/(A^2-Delta(P)B).                     (9)
```

The equality in `(9)` ranges over every rational function in `K` whose
pullback is polynomial on the source; it is not merely a presentation of a
chosen finitely generated subring.

## 2. Proof of the full-field intersection

Take `h in R_full`.  View it first as a reduced rational function of `Q`
over `k(P)`.  After the field extension `k(P) subset k(x)`, the substitution
`Q=y/f'(x)` is a nonzero linear rescaling of `y`.  Since `h` is in
`k[x,y]`, it has no finite pole in `y`, hence no finite pole in `Q`.
Coprime numerator and denominator over `k(P)` remain coprime after a field
extension, so the reduced denominator has `Q`-degree zero.  Therefore

```text
h=sum_(n=0)^N a_n(P)Q^n,              a_n in k(P).      (10)
```

After pullback, the coefficient of `y^n` is

```text
a_n(f(x))/f'(x)^n.                                      (11)
```

Suppose `a_n` had a finite pole at `t0`.  In a reduced presentation its
numerator is nonzero at `t0`, while its denominator vanishes there.  At any
root of `f(x)=t0`, composition creates a pole; the additional division by
`f'(x)^n` cannot cancel it.  This contradicts `(11) in k[x]`.  Hence `a_n`
has no finite pole and

```text
a_n in k[P].                                            (12)
```

At a critical point `alpha`, put

```text
m_alpha=ord_(T=f(alpha)) a_n(T).
```

The two local orders are

```text
ord_alpha(f')=1,
ord_alpha(f-f(alpha))=2.                                (13)
```

Thus `(11)` is regular at `alpha` exactly when

```text
2m_alpha>=n.                                            (14)
```

There are no other possible poles.  Distinctness of the critical values
turns `(14)` into the single sharp divisibility condition

```text
Delta(P)^ceil(n/2) divides a_n(P).                      (15)
```

Conversely `(6)` shows that `(15)` makes `(11)` polynomial.  For `n=2m`
and `n=2m+1`, respectively,

```text
Delta(P)^m Q^(2m)=B^m,
Delta(P)^(m+1) Q^(2m+1)=A B^m.                         (16)
```

Equations `(10)--(16)` prove `R_full=k[P,A,B]`.

The relation `A^2=Delta(P)B` follows from `(8)`.  Conversely,
`k[P,A,B]/(A^2-Delta(P)B)` embeds in `k(P,A)` by
`B=A^2/Delta(P)`, so it is a two-dimensional domain.  The displayed map to
`R_full` is surjective, and both fraction fields are `k(P,Q)` because
`Q=A/Delta(P)`.  Its prime kernel therefore has height zero and is zero.
This proves that the displayed relation is the full kernel.

The ceiling in `(15)` is exact, not an estimate.  One fewer copy of
`Delta` leaves a pole at every critical point at which the relevant factor
is missing.

## 3. Nodal geometry and failure of global etaleness

Let

```text
S_full=Spec(R_full)
      =V(A^2-Delta(P)B) subset A3_(P,A,B).              (17)
```

The three partial derivatives of the defining equation are

```text
-Delta'(P)B,                    2A,       -Delta(P).    (18)
```

Since `Delta` is squarefree, `(18)` vanishes exactly at

```text
(P,A,B)=(c,0,0),                 Delta(c)=0.            (19)
```

There is one such point for each critical value.  In local coordinates
`t=P-c`, the quadratic initial equation is

```text
A^2-Delta'(c)tB.                                      (20)
```

Its Hessian is nondegenerate, so every point in `(19)` is an ordinary
`A1` node.  In particular, `S_full` is not smooth and is not `A2`.

All full-field observables give the polynomial map

```text
Phi:A2_(x,y) -> S_full,
(x,y) |-> (f(x), f'(x)H(x)y, H(x)y^2).                 (21)
```

On `Delta!=0`, equation `(9)` recovers `Q=A/Delta`.  The source ring there
is

```text
k[x,y,1/Delta(f)]
 =k[P,Q,1/Delta][x]/(f(x)-P),             y=f'(x)Q.    (22)
```

The element `f'(x)` is a unit in `(22)` by `(6)`.  Hence `(21)` is finite
etale of degree `d` over this open subset.

Globally, the completion fails.  At every critical point `alpha`, the
source point `(alpha,0)` maps to the singular point `(f(alpha),0,0)`, so
`Phi` cannot be everywhere etale.  If `d>=3`, factor

```text
f(x)-f(alpha)=(x-alpha)^2 g_alpha(x),
deg(g_alpha)=d-2.                                      (23)
```

Every root `beta` of `g_alpha` is regular: a critical such root would give
two critical points with the same critical value.  Equation `(6)` then
gives `H(beta)=0`, and the entire affine line contracts:

```text
{x=beta, y arbitrary} -> (f(alpha),0,0).               (24)
```

Therefore `Phi` is not quasi-finite for `d>=3`.  For `d=2` the companion
line mechanism is absent; the theorem claims only the already proved global
failure of etaleness.  This degree-two boundary is important: the phrase
"etale over Delta!=0" is a proved sufficient open, not a claim that it is
the maximal etale locus.

A resolution or affine modification may replace `S_full` by a smooth
surface, but it changes `Spec(R_full)`.  It cannot result from adjoining one
more element of `K intersect k[x,y]`, because `(9)` already computes that
entire ring.  This does not exclude choosing a different subalgebra or a
nonseparated rational seed.

## 4. Rational cubic triple-collision control

Take

```text
f=x(x-3)(x-8)=x^3-11x^2+24x,
f'=(x-6)(3x-4).                                        (25)
```

Its critical values are

```text
f(6)=-36,                    f(4/3)=400/27,             (26)
```

and a convenient scaled critical-value polynomial is

```text
Delta(P)=(P+36)(27P-400).                              (27)
```

The exact pullback is

```text
Delta(f)
 =(x-6)^2(x+1)(3x-25)(3x-4)^2
 =f'(x)^2 (x+1)(3x-25).                                (28)
```

The pair

```text
P=x(x-3)(x-8),
Q=y/[(x-6)(3x-4)]                                      (29)
```

has the rational triple collision

```text
(0,0),(3,0),(8,0) -> (0,0),                            (30)
```

and all three denominators in `(30)` are nonzero.  Its full-field
observable surface is

```text
A^2=(P+36)(27P-400)B.                                  (31)
```

The two nodes are `(-36,0,0)` and `(400/27,0,0)`.  The factorizations

```text
f+36=(x-6)^2(x+1),
27f-400=(3x-25)(3x-4)^2                                (32)
```

show that the companion lines `x=-1` and `x=25/3` collapse to the nodes.
This gives rational collision points, rational critical values, and no
hidden field extension.

## 5. Exact distinction from THM-3561

The two nearby results use different universes and different rational maps:

| theorem | rational seed | intersection universe | resulting surface |
|---|---|---|---|
| THM-3561 | nonseparated `a=q/D^2`, `c=xD(D+2)` | `C[a,c] intersect C[x,q]`, only functions polynomial in the original target coordinates | smooth `c^2e=b(b+4)` |
| THM-3567 | separated `P=f(x)`, `Q=y/f'(x)` | `k(P,Q) intersect k[x,y]`, every rational target-field function with polynomial pullback | nodal `A^2=Delta(P)B` |

THM-3561 explicitly does not claim a full-field intersection.  THM-3567
does not apply to its mixed denominator `D=1+x^2q`.  Thus the smooth
Danielewski completion remains a positive near-counterexample, while the
separated all-degree family is now closed at its full-field completion.

## 6. Counterexample boundary and current passport

The audit yields the following exact status split:

```text
PROVED:       unit rational Jacobian and d-fold regular-value collisions;
              the full-field intersection (9);
              one node per critical value;
              degree-d finite etaleness over Delta!=0;
              non-quasi-finiteness for d>=3.

REFUTED:      smoothness of this full-field completion;
              global etaleness of this full-field completion;
              repair by adjoining another full-field polynomial observable.

OPEN:         nonseparated rational Keller seeds;
              smooth selected-subalgebra completions such as THM-3561;
              affine modifications using a changed target ring;
              JC(2) and DC(2).                         (33)
```

Increasing `deg(f)` cannot rescue this separated architecture: it creates
more critical values and therefore more nodes, while every critical fibre
of degree at least three supplies contracted companion lines.  A surviving
rational-completion variant must alter the local valuation semigroup
`(15)`, for example through a genuinely two-variable denominator.  The
canonical THM-3561 seed does exactly that, but its smooth target is not an
affine plane.  THM-3569 now forces any Darboux pair there to at least six
nonconstant weight pieces; the first live support shapes are `3x3` and
`2x4` (up to swapping outputs).

Any eventual planar polynomial counterexample must independently meet the
current cited degree passport.  Using the Guccione--Guccione--Horruitiner--
Valqui sub-`125` classification (arXiv:2204.14178) as a **CITED** input, not
as a theorem proved here, the only reduced pair below height `125`, up to
transposition, is `(72,108)`.  Its leading forms and the internally proved
first radical tax have the shape

```text
F_72=c H^2,                 G_108=d H^3,   deg(H)=36,
rad(H) | F_71,              H rad(H) | G_107.           (34)
```

Nothing in the rational family `(2)` supplies `(34)`.  The cheapest honest
next test is instead to combine the first surviving Danielewski `3x3` and
`2x4` support shapes with a squarefree two-factor model for `H` in `(34)`,
and then check the complete constant-Jacobian rows rather than only the
collision or leading form.

## 7. Exact verification

Reproduce with

```bash
python3 04-computation/jc_separated_rational_keller_observable_suspension_thm3567.py
python3 -O 04-computation/jc_separated_rational_keller_observable_suspension_thm3567.py
```

The ordinary and optimized transcripts agree byte-for-byte with the stored
output.  The companion checks the rational Jacobian, rational triple
collision, both critical values, cubic discriminant, `(28)`, polynomial
observables and relation, sharp parity-divisibility ledger through sixteen
independent exponents, both nodes and nondegenerate Hessians, both
critical-fibre factorizations, collapsed companion lines, and a regular
degree-three fibre.  It also checks the sharp quadratic hostile
`(x,y)->(x^2,2xy,y^2)`, which is etale at `Delta=0` away from the node and
prevents a false maximal-etale-locus claim.  The all-function intersection
equality and all-degree
geometry are the symbolic arguments in Sections 2--3, not extrapolations
from finite controls.

**QED.**
