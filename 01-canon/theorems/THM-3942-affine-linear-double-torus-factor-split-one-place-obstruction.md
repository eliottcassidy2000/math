---
id: THM-3942
title: "Affine-linear double-torus factor splits cannot have a nonlinear one-place branch"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Two distinct affine-linear torus coefficients force, by the squarefree UFD
  factorization of their difference of cubes, exactly two complementary split
  types.  The 0|3 type is an irreducible sextic whose normalization is a
  Fermat elliptic curve with three infinity places.  The 1|2 type is an
  irreducible quartic whose normalization is Gm with two infinity places; all
  singleton choices, including l1|l0*l2, are linearly equivalent.  Both types
  really have two independent smooth-locus Cardano characters, so the failure
  is the one-place gate, not character scarcity.  Parallel coefficients give
  only a line, a reducible/nonreduced branch, or duplicate torus data.  More
  generally an A1-normalized branch with a nonconstant coefficient map cannot
  assign each whole factor p1-zeta*p0 to one side: a nonlinear escape must
  split a factor internally or exploit gcd/multiplicity overlap.
source: jc-degree6-place / post-THM-3940 universal double-torus split classification, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (double_torus_split_audit, 2026-08-24).
  The audit independently reconstructed both factor partitions, scalar gauges,
  normalization inverses and infinity ledgers; checked normality and the
  disjoint local A2 character witnesses; found and repaired the initially
  missing H|D=0 hypothesis in the nonlinear corollary; and exhausted the
  parallel and constant coefficient boundaries.  The 48-gate companion
  byte-matches in normal and optimized modes, frozen output and all hashes
  agree, and documentation/diff checks pass.
depends_on: []
related:
  - THM-3935-linear-conic-resolvent-class-group-unique-cubic-character
  - THM-3937-linear-conic-fold-three-family-uniform-resolvent-rigidity
  - THM-3939-two-boundary-elliptic-resolvent-three-character-rank-one-gate
  - THM-3940-i7-rank-two-linear-cross-term-resolvent-unique-character
  - THM-3944-repeated-factor-double-torus-one-place-square-conductor-collapse
  - THM-3947-scalar-weighted-repeated-square-split-trichotomy
script: 04-computation/jc2_affine_linear_double_torus_factor_split_thm3942.py
output: 05-knowledge/results/jc2_affine_linear_double_torus_factor_split_thm3942.out
script_sha256: 39a4019234caa6784652170056eea89a85c7442e41e409f7d29ae4adb78fb0f9
output_sha256: 6ad9a99de5fe37e956d02793665ed99846fe8754362abad3c54ef38ef5cd34cb
semantic_sha256: c0efbc25737af9d980219b0188e81cf60605c528ca2bded77e82ac0bb8609c41
hash_basis: raw LF bytes
---

# THM-3942 -- affine-linear double-torus splits pay at least two ends

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero.  Let

```text
H=q_0^2-4p_0^3=q_1^2-4p_1^3                         (1)
```

in `k[A,C]`, where `p_0,p_1` are affine-linear.  Suppose first that they are
independent affine coordinates.  After taking `x=p_0,y=p_1`, affine coordinate
changes, multiplication of `x` by a cube root of unity, nonzero scalar gauges,
and replacing `q_0` by `-q_0`, every such pair is exactly one of the following.

```text
type E (0|3):
q_0=y^3-x^3-1,                 q_1=y^3-x^3+1,
H_E=(y^3-x^3-1)^2-4x^3;

type C (1|2):
q_0=y^2+xy+x^2-y+x,            q_1=y^2+xy+x^2+y-x,
H_C=(y^2+xy+x^2-y+x)^2-4x^3.                    (2)
```

The two branches are absolutely irreducible, but their affine normalizations
are respectively

```text
H_E:  {v^3-h^3=1} minus its three points at infinity,
H_C:  P1 minus two points ~= Gm.                         (3)
```

Thus neither is `A1`.  This is sharp on the character side.  The normal
quadratic surfaces

```text
S_T=Spec k[x,y,W]/(W^2-H_T),             T in {E,C},     (4)
```

carry two linearly independent order-three divisor classes, and the two
Cardano radicands `q_0+W,q_1+W` give independent `C3` characters on
`(S_T)_reg`.

The requested ostensibly different split

```text
(y-omega*x) | (y-x)(y-omega^2*x),        omega^2+omega+1=0, (5)
```

is type C under `x -> omega*x`.  Its normalization is again `Gm`, with exactly
two infinity places.

If `p_0,p_1` are parallel rather than independent and are not related by
`p_1=omega^i p_0`, then `H` is univariate in an affine coordinate.  Its reduced
zero set is a union of parallel lines, and it is irreducible only when it is a
single line.  If `p_1=omega^i p_0`, equation `(1)` forces `q_1=+/-q_0`; the two
presentations are duplicates and cannot provide two independent characters.
Consequently there is no irreducible nonlinear `A1` branch with two genuinely
distinct affine-linear torus coefficients.

There is also a nonlinear consequence.  Let `D` be any irreducible affine
curve with normalization `A1`, and suppose the restrictions of polynomial
functions `p_i,q_i` to `D` satisfy `(1)`, the common branch equation

```text
H|_D=(q_0^2-4p_0^3)|_D=0,                              (6)
```

and a **whole-factor split**

```text
q_1-q_0=2 lambda product_{zeta in I}(p_1-zeta p_0),
q_1+q_0=2 lambda^(-1) product_{zeta notin I}(p_1-zeta p_0).  (7)
```

for one subset `I` of the cube roots of unity and `lambda in k^*`.  Then
`(p_0,p_1):D^nu=A1 -> A2` is constant.  Hence any nonconstant one-place
double-torus design must do something absent from the affine-linear atlas:
split the irreducible factors of at least one polynomial `p_1-zeta p_0`
between the two sides of `(7)`, or enter a common-factor/repeated-factor
overlap where the squarefree partition argument no longer applies.

The theorem classifies this two-affine-coefficient and whole-factor grammar.
It does not exclude internally split nonlinear coefficients, gcd overlap,
nonreduced factor packets, more than two torus presentations, higher-genus
coefficient carriers, or arbitrary one-place discriminant branches.

## 1. The UFD leaves only `0|3` and `1|2`

Put

```text
L_zeta=y-zeta*x,                    zeta^3=1.
```

Subtracting the two rows of `(1)` gives

```text
(q_1-q_0)(q_1+q_0)
     =4(y^3-x^3)=4 product_{zeta^3=1} L_zeta.           (8)
```

Because `x,y` are affine coordinates, the three `L_zeta` are nonassociate,
pairwise coprime affine primes.  The right side of `(8)` is squarefree in the
UFD `k[x,y]`.  Therefore for a unique subset `I` and a scalar `lambda!=0`,

```text
q_1-q_0=2 lambda product_{zeta in I}L_zeta,
q_1+q_0=2 lambda^(-1) product_{zeta notin I}L_zeta.     (9)
```

Replacing `I` by its complement changes `q_0` to `-q_0` and leaves `H`
unchanged.  Thus only `|I|=0,1` remain.  Multiplying `x` by a cube root of
unity permutes the three singleton choices without changing `x^3`.

The scalar in `(9)` carries no geometry.  For `|I|=0`, write `lambda=c^3`
and substitute `x=c^2X,y=c^2Y,W=c^3W_0`; this extracts the scalar `c^6`
from `(1)` and gives type E.  For `|I|=1`, the substitution

```text
x=lambda^2 X,              y=lambda^2 Y,
W=lambda^3 W_0                                            (10)
```

extracts `lambda^6` and gives type C.  This proves completeness of `(2)`.

## 2. The empty split is Fermat elliptic with three ends

For type E, the branch equation is

```text
H_E=(y^3-x^3-1)^2-4x^3=0.                              (11)
```

On the open set `xy!=0`, define

```text
h=(y^3-x^3-1)/(2x),                  v=(h^3+1)/y.       (12)
```

Equation `(11)` gives successively

```text
h^2=x,                 y^3=(h^3+1)^2,
v^2=y,                 v^3=h^3+1.                       (13)
```

Conversely `(x,y)=(h^2,v^2)` sends the smooth Fermat cubic

```text
E: v^3-h^3=1                                             (14)
```

to `(11)`, and `(12)` is its rational inverse.  In particular the localization
of `k[x,y]/(H_E)` at `xy` is a domain and is birational to `(14)`.  Neither
`x` nor `y` divides `H_E`, so this open meets every putative component;
therefore `H_E` is absolutely irreducible and `(14)` is its normalization.

Projectively the normalization map is

```text
[H:V:R] |-> [H^2:V^2:R^2],             V^3-H^3=R^3.    (15)
```

It is basepoint-free.  Target infinity pulls back to `R=0`, where

```text
[H:V:R]=[1:zeta:0],                    zeta^3=1.         (16)
```

These are three distinct points, and both affine target coordinates have a
pole at each.  The projective normalization has genus one and the affine
normalization has exactly three ends.

For later normality and character checks, differentiate `(11)`.  With
`q_0=y^3-x^3-1,q_1=y^3-x^3+1`,

```text
(H_E)_x=-6x^2q_1,                     (H_E)_y=6y^2q_0. (17)
```

Together with `H_E=0`, these equations give exactly six affine branch
singularities:

```text
x=0, y^3=1;                or                y=0, x^3=1. (18)
```

## 3. Every singleton split is the same two-ended conic

For type C, set `t=y-x`.  On the branch, define `h=q_0/(2x)`.  Then `h^2=x`
and `q_0=2h^3`, while `t` satisfies

```text
t^2+(3h^2-1)t+3h^4-2h^3=0,
disc_t=-(h-1)^3(3h+1).                                  (19)
```

Removing the square factor gives the smooth conic

```text
u^2=-(h-1)(3h+1),       u=(2t-1+3h^2)/(h-1).            (20)
```

The discriminant in `(19)` is not a square in `k(h)`, so this quadratic row
is irreducible.  Since `h=q_0/(2x)` is recovered from the branch and

```text
s=(u-1)/h,
```

the conic and the branch are birational.  Thus the branch is absolutely
irreducible and rational.  An exact normalization is

```text
x=4(s-1)^2/(s^2+3)^2,
y=(s^2-1)^2/(s^2+3)^2.                                  (21)
```

Its basepoint-free projective form is

```text
[S:R] |-> [
  4(S-R)^2R^2 : (S^2-R^2)^2 : (S^2+3R^2)^2
].                                                       (22)
```

Thus the target infinity line pulls back to `S^2+3R^2=0`.  It has two
distinct roots and neither cancels from the first two coordinates.  Hence the
affine normalization is `P1` minus two points, namely `Gm`.

For the split explicitly requested in `(5)`, put

```text
q_0=y^2+omega*x*y+omega^2*x^2-y+omega*x,
q_1=y^2+omega*x*y+omega^2*x^2+y-omega*x.                (23)
```

These are exactly the type-C rows after `x -> omega*x`.  Equivalently, with
`X=omega*x,t=y-X`, they become

```text
q_0=3X^2+3tX+t^2-t,             q_1=3X^2+3tX+t^2+t.    (24)
```

The normalization `(21)` simply acquires `omega^2` in its first coordinate.
It retains the same two infinity places.  Directly, the quadratic row in the
`x=h^2` chart has discriminant

```text
-omega^2(h-omega)^3(3h+omega),                          (25)
```

again a smooth conic after removing the square.

The exact singular ideal of type C has Groebner basis

```text
x^2+xy-x+y^3-2y^2+y,
y^2(2x+y-1),
y^2(y-1)^2,                                             (26)
```

so its affine singular support is exactly

```text
(0,0),                       (0,1),                      (1,0). (27)
```

## 4. The two characters are real, independent, and not the failure

Both branch polynomials are irreducible, and `(18)` and `(27)` show that each
surface `(4)` has only isolated singularities.  A hypersurface is `S2`; the
isolated singular locus has codimension two.  Hence both `S_T` are normal.

On either surface define the height-one primes

```text
D_0=(x,q_0+W),                         D_1=(y,q_1+W).    (28)
```

At the generic points of these primes the conjugate factors are units, so the
Cardano norm identities give

```text
div(q_0+W)=3D_0,                       div(q_1+W)=3D_1. (29)
```

For type E, inspect the surface points over `(x,y)=(0,1)` and `(1,0)`.
At the first, `q_0` is a transverse parameter and

```text
(W+q_0)(W-q_0)=-4x^3;                                  (30)
```

the completed local ring is `A2` and `D_0` is its nonzero class.  The divisor
`D_1` misses this point.  At the second point the same statement holds with
`q_1,y,D_1`, while `D_0` misses the point.  Type C has identical disjoint
`A2` witnesses at `(0,1)` and `(1,0)`.  Therefore a relation

```text
a[D_0]+b[D_1]=0                                         (31)
```

localizes first to `a=0 mod 3` and then to `b=0 mod 3`.  The classes have
exact order three by `(29)` and non-Cartierness, so they span `(Z/3)^2` inside
`Cl(S_T)[3]=Pic((S_T)_reg)[3]`.  Their Kummer covers are consequently two
independent smooth-locus `C3` characters.  The normalization ends, not a
collapse of the two Cardano formulas, are the exact obstruction.

## 5. Parallel affine coefficients cannot hide a nonlinear branch

Suppose now `p_0,p_1` are affine dependent.  Unless both are constant, choose
an affine coordinate `x` so that

```text
p_0=x,                       p_1=ax+b                    (32)
```

after possibly exchanging the indices.  If `p_1^3-p_0^3!=0`, equation `(8)`
has a nonzero right side in `k[x]`.  Regard its left side in `k[x,T]`, where
`T` is a complementary affine coordinate.  Polynomial degree in `T` is
additive under nonzero products, hence both `q_1-q_0` and `q_1+q_0` lie in
`k[x]`.  Thus `q_0,q_1,H` are univariate.  Over algebraically closed `k`, the
reduced zero set of `H` is a union of parallel lines; it is irreducible only
when it is one line.

If both `p_i` are constant, `(1)` makes `q_1^2-q_0^2` constant, so both
factors `q_1-q_0,q_1+q_0` are constant whenever their product is nonzero;
then `H` is constant, giving either the empty set or the whole plane rather
than a branch curve.  If their product is zero, the
duplicate row below applies.  Finally, if `p_1^3=p_0^3`, the domain property
gives `p_1=zeta p_0` for one cube root
of unity.  Equation `(8)` then gives `q_1=+/-q_0`.  This is the same torus
presentation, possibly conjugated, not a second cubic direction.

## 6. Whole-factor splitting is impossible on a nonconstant `A1` carrier

Let `D^nu=A1` and assume `(6)-(7)`.  The common branch hypothesis gives
`q_0^2=4p_0^3` on `D`.
After the scalar normalizations in Section 1, the coefficient map

```text
(p_0,p_1): A1 -> {H_E=0}             or             {H_C=0} (33)
```

is a morphism.  If nonconstant, it is dominant and the universal property of
normalization lifts it respectively to

```text
A1 -> E minus three points,                or          A1 -> Gm. (34)
```

The second map is constant because `k[t]^*=k^*`.  In the first case,
properness extends the map across infinity to `P1 -> E`; Riemann--Hurwitz
forbids a nonconstant map from a genus-zero curve to a genus-one curve.  Hence
both maps in `(34)`, and therefore `(33)`, are constant.

This proves the nonlinear escape statement.  Notice what is and is not lost:
the target of `(33)` remembers the full complementary factor partition but
forgets how a reducible nonlinear polynomial `p_1-zeta p_0` may distribute its
own irreducible factors.  That internal distribution, together with gcd and
multiplicity data, is exactly the next sidecar a one-place search must retain.

## 7. Reproduction and hostile boundary

Run

```bash
python3 04-computation/jc2_affine_linear_double_torus_factor_split_thm3942.py
python3 -O 04-computation/jc2_affine_linear_double_torus_factor_split_thm3942.py
```

The companion verifies the difference-of-cubes and both torus identities,
the two subset orbits, scalar gauges, Fermat and conic normalization maps,
projective basepoint freeness, infinity divisors, singular-support equations,
the requested `l_1` equivalence and direct discriminant, Cardano norms, and
the disjoint local `A2` witnesses.

Hostile boundaries are deliberate:

1. parallel coefficients can give a single affine line, but not a nonlinear
   branch or a second character;
2. `p_1=zeta p_0` makes `(8)` zero and duplicates the presentation;
3. the nonlinear corollary assumes both the branch equation `(6)` and the
   exact whole-factor formulas `(7)`;
4. common or repeated irreducible factors can invalidate the squarefree UFD
   partition and remain open;
5. two independent characters are necessary design data, not by themselves a
   Keller atlas or a nonmonogenic cubic order.
