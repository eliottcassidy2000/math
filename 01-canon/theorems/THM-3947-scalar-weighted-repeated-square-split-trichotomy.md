---
id: THM-3947
title: "Scalar-weighted repeated-square splits have a three-parabola trichotomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the repeated-factor row p1-p0=G^2, allow the two internal factors to
  carry arbitrary reciprocal scalars lambda and lambda^-1.  With t=lambda^2,
  the common discriminant is governed by one cubic h_t(P/G^2), whose exact
  discriminant is 48*delta*t*(t-1)^3*(t+omega)^3.  For generic t it is the
  union of three distinct smooth one-place parabolas.  At t=1 it has a
  doubled p0 component, and at t=-omega it has a doubled p1 component; the
  latter is the former after swapping the two torus rows and rescaling G.
  Neither endpoint is triple.  Thus no scalar weighting makes the full
  discriminant irreducible: the tempting one-place objects are individual
  reduced components, never the whole branch divisor.
source: jc-zero-debt / arbitrary-scalar completion of THM-3944, 2026-08-24
audit: >
  Independent hostile audit reconstructed the lambda-weighted reduction and
  exact -4 scalar, checked the root-cubic discriminant and both 2+1 endpoint
  factorizations, verified that t=0 is unavailable and no triple-root seam is
  omitted, replayed the full q-row endpoint swap (not merely the branch
  identity), and separated one-place components from the reducible full
  divisor. Normal/-O/frozen runs and all hashes match.
depends_on: []
related:
  - THM-3942-affine-linear-double-torus-factor-split-one-place-obstruction
  - THM-3944-repeated-factor-double-torus-one-place-square-conductor-collapse
  - THM-3946-affine-internal-factor-split-two-end-conductor-collision-dichotomy
script: 04-computation/jc2_scalar_weighted_repeated_square_split_thm3947.py
output: 05-knowledge/results/jc2_scalar_weighted_repeated_square_split_thm3947.out
script_sha256: a5a095054a92b8552d0f00194663e8f022b4e171efb4610f597d54cc1b553137
output_sha256: f030902fb4a3180c9c5169d7612e937ea923746e9527fcccb83fa253da5e2514
semantic_sha256: cd590de731560e5b91c4f81f1ee840247d500f459774880c801398b00b2125b6
hash_basis: raw LF bytes
---

# THM-3947 -- scalar imbalance never glues the repeated square

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero.  Fix

```text
omega^2+omega+1=0,       delta=omega-omega^2,       delta^2=-3,       (1)
```

and take `P,G` as independent affine coordinates.  Put

```text
p0=P,              p1=P+G^2,
L1=p1-omega*p0,    L2=p1-omega^2*p0.                         (2)
```

For an arbitrary `lambda in k^*`, make the reciprocal scalar-weighted
internal split

```text
q1-q0=2 lambda G L1,
q1+q0=2 lambda^(-1) G L2.                                (3)
```

Then

```text
H=q0^2-4p0^3=q1^2-4p1^3                                (4)
```

always has reducible support.  More precisely, set `t=lambda^2` and

```text
h_t(x)= [1+(1-omega^2)x-t(1+(1-omega)x)]^2-4t x^3.       (5)
```

If `r1,r2,r3` are its roots in `k`, counted with multiplicity, then

```text
H=-4 product_{j=1}^3 (P-rj G^2).                         (6)
```

The complete multiplicity trichotomy is

```text
t notin {1,-omega}:   three distinct factors;
t=1:                  H=-P^2(4P+3G^2);
t=-omega:             H=-(P+G^2)^2(4P+G^2).              (7)
```

There is no triple-root parameter.  Every distinct reduced factor in `(6)`
is a smooth `A1` with one normalization place at infinity, but the **full**
reduced branch has three components in the generic row and two at either
endpoint.  Hence an individual one-place component cannot be promoted to an
irreducible one-place discriminant.

The two exceptional rows are equivalent.  If `iota^2=-1`, swap the torus
presentations, replace `omega` by `omega^2`, and put

```text
P'=p1=P+G^2,     G'=iota G,     lambda'=-iota*omega*lambda. (8)
```

When `lambda^2=-omega`, one has `(lambda')^2=1`, and the transformed split is
exactly the `t=1` row.  On the discriminant this is the transparent identity

```text
-(P+G^2)^2(4P+G^2)=-P'^2(4P'+3G'^2).                    (9)
```

Thus THM-3944 is not hiding a scalar-imbalance escape: up to exchanging the
two presentations, it is the unique repeated-square conductor collision.
This theorem does **not** treat a genuinely unequal split
`p1-p0=F G` with coprime nonassociate factors, a repeated factor distributed
among several `Li`, or a non-coordinate pair `(P,G)`.

## 1. The common discriminant and the one-variable reduction

Equation `(3)` gives explicitly

```text
q0=G(lambda^(-1)L2-lambda L1),
q1=G(lambda^(-1)L2+lambda L1).                          (10)
```

Since `p1-p0=G^2`, subtraction of the two discriminant rows is

```text
q1^2-q0^2
 =(q1-q0)(q1+q0)
 =4G^2 L1 L2
 =4(p1^3-p0^3),                                        (11)
```

which proves `(4)`.  On putting `x=P/G^2`, equations `(5)` and `(10)` give

```text
H=lambda^(-2) G^6 h_t(P/G^2).                           (12)
```

The leading coefficient of `h_t` is `-4t`.  Because `t!=0`, all three roots
are finite.  Factoring `(12)` over the algebraically closed constant field
and using `lambda^(-2)t=1` gives `(6)` with the exact leading scalar `-4`.

This is the key structural point.  The scalar changes only a cubic in the
weighted ratio `P/G^2`; it cannot glue its constant roots into one polynomial
branch over `k`.

## 2. The discriminant of the root cubic is complete

Direct elimination in `k[t,delta]/(delta^2+3)` gives

```text
disc_x(h_t)=48 delta t (t-1)^3(t+omega)^3.               (13)
```

Since `delta!=0` and `t!=0`, a repeated root occurs only at `t=1` or
`t=-omega`.  At those values the exact factorizations are

```text
h_1(x)=-x^2(4x+3),
h_{-omega}(x)=omega (x+1)^2(4x+1).                      (14)
```

They have respective multiplicity patterns `2+1` at

```text
{0,-3/4},                         {-1,-1/4}.             (15)
```

Thus neither exceptional row has a triple root, and `(13)-(15)` exhaust all
parameters.  Substitution in `(12)` gives exactly the two endpoint formulas
in `(7)`.

The excluded value `t=0` is a useful hostile boundary: it removes the cubic
term in `(5)`, but it is unavailable because `(3)` contains
`lambda^(-1)`.  Algebraic closedness is also used essentially in `(6)`;
without it the constant cubic `h_t` need not split over the ground field.

## 3. Each factor is one-place, while their union is not

For a root `r`, the corresponding affine component is

```text
Dr: P=rG^2.                                              (16)
```

It is the graph of a polynomial in `G`, hence smooth and isomorphic to
`A1_G`.  If `r!=0`, its projective closure

```text
PZ-rG^2=0                                               (17)
```

has the unique infinity point `[1:0:0]`, where the `Z` derivative is nonzero.
For `r=0`, it is the line `P=0`, again with one infinity point.  Thus every
reduced component separately passes the one-place test.

But distinct roots give distinct polynomial factors.  In the generic row the
three parabolas meet at the affine origin and remain three irreducible
components.  At the exceptional rows, `(14)` leaves two distinct components
and doubles one of them.  This proves both the reducibility assertion and the
distinction between a one-place **component** and a one-place **full branch**.

## 4. The endpoint swap is an equivalence of presentations

At `t=-omega`, take `iota^2=-1` and perform `(8)`, together with

```text
(p0',q0')=(p1,q1),                  (p1',q1')=(p0,q0),
omega'=omega^2.                                            (18)
```

Then `p1'-p0'=(G')^2`.  The identities

```text
p1'-omega' p0'=-omega^2 L1,
p1'-(omega')^2 p0'=-omega L2                              (19)
```

and `lambda'=-iota*omega*lambda` turn `(3)` into

```text
q1'-q0'=2lambda'G'(p1'-omega' p0'),
q1'+q0'=2(lambda')^(-1)G'(p1'-(omega')^2p0').             (20)
```

Moreover `(lambda')^2=-omega^2t=1`.  Hence this is the first endpoint, not a
third geometry in disguise.  Formula `(9)` is its branch-level shadow.

## 5. Reproduction and next boundary

Run

```bash
python3 04-computation/jc2_scalar_weighted_repeated_square_split_thm3947.py
python3 -O 04-computation/jc2_scalar_weighted_repeated_square_split_thm3947.py
```

The companion verifies the two torus rows, the weighted reduction, the exact
root-cubic discriminant, both endpoint factorizations and multiplicities, the
absence of a triple-root seam, the projective one-place component ledger, and
the full endpoint change of variables.  It also includes the forbidden
`t=0` degree-drop and a generic distinct-root specialization as hostile
controls.

The next genuine internal-split coordinate is therefore not a scalar on the
same square.  It is the coprime unequal-factor row `p1-p0=F G` routed by
THM-3946, or a distribution that mixes more than one cube-difference factor.
