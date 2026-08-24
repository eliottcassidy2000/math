---
id: THM-3944
title: "Repeated-factor double-torus one-place escape pays square-conductor and a rank-one boundary character"
status: >
  PROVED + VERIFIED-EXACT + CORRECTED AFTER INDEPENDENT HOSTILE RE-AUDIT.
  The cheapest internal factor split left open by THM-3942 sets
  p1-p0=G^2 and assigns one copy of G to each complementary cube factor.
  It produces a smooth one-place parabola, but the common discriminant is
  -P^2(4P+3G^2): the extra line is even and the quadratic order is nonnormal.
  Its integral closure is A2, with conductor (P,W), scalar units, trivial
  class group, and no cubic character extending across the conductor.  One
  Cardano radicand becomes an exact cube and its depressed cubic splits.  The
  other is a genuine nontrivial etale mu_3-torsor on the original order's
  regular locus Gm^2, with boundary vector (2,1), but ramifies along the two
  conductor-preimage lines on A2.  Thus the construction buys one infinity
  place at the price of nonnormality and drops the displayed Cardano rank
  from two to one, not to zero.
source: jc-degree6-place / nonlinear internal-split successor to THM-3942, 2026-08-24
audit: >
  CORRECTED INDEPENDENT HOSTILE RE-AUDIT (jc-cohn3709 and
  jc_zero_debt_lift, 2026-08-24). The initial promotion correctly checked
  independently checked the Gauss/UFD domain argument; the finite same-field
  polynomial normalization; the exact module quotient and conductor; the
  reduced one-place parabola versus full nonreduced branch scope; both cube
  factorizations and the split cubic; and the two regular-prime ramification
  residues of the noncube radicand. It incorrectly identified regularity on
  the normalization A2 with the smooth locus of the original nonnormal
  order. Re-audit proves B0_reg=D(P)=Gm^2 and recovers the nontrivial boundary
  class (2,1); see MISTAKE-471. Normal and optimized runs byte-match the
  repaired frozen output, and the theorem remains scoped to the displayed
  balanced split.
depends_on: []
related:
  - THM-3942-affine-linear-double-torus-factor-split-one-place-obstruction
script: 04-computation/jc2_repeated_factor_double_torus_square_conductor_thm3944.py
output: 05-knowledge/results/jc2_repeated_factor_double_torus_square_conductor_thm3944.out
script_sha256: 3698b325ff185a5d98e4edcafcadcad86977b1b4a41edcc6dc18c85146246585
output_sha256: 23d9f389f96dd5f2839c61b49bda5e2b6666f30329ab664b3d85c67ab3a62fdf
semantic_sha256: a93d0aa73a1c95e5e3a04ec1833ff337f3b52294bcebe5145ae62f184edb8d92
hash_basis: raw LF bytes
---

# THM-3944 -- the first internal split trades ends for conductor debt

**PROVED + VERIFIED-EXACT + CORRECTED AFTER INDEPENDENT HOSTILE RE-AUDIT.**
Work over an algebraically closed field `k` of characteristic zero.  Choose
`omega^2+omega+1=0`, put

```text
delta=omega-omega^2,                         delta^2=-3, (1)
```

and let `P,G` be independent affine coordinates.  Define

```text
p_0=P,                         p_1=P+G^2,
q_0=delta*P*G,                 q_1=G(3P+2G^2).          (2)
```

Then there are two exact torus presentations of the same discriminant:

```text
H=q_0^2-4p_0^3=q_1^2-4p_1^3
 =-P^2(4P+3G^2).                                      (3)
```

The reduced nonlinear component

```text
Gamma: 4P+3G^2=0                                      (4)
```

is a smooth `A1` with exactly one normalization place at infinity.  This is
the first whole construction after THM-3942's factor-split obstruction to
attain the desired one-place geometry.  It does so only because `(3)` also
contains the line `P=0` with even multiplicity.

Indeed the quadratic order

```text
B_0=k[P,G,W]/(W^2+P^2(4P+3G^2))                       (5)
```

is not normal.  Its integral closure is obtained by `V=W/P`:

```text
B=k[P,G,V]/(V^2+4P+3G^2)
 ~= k[G,V],                    P=-(V^2+3G^2)/4.         (6)
```

The exact conductor and index module are

```text
cond(B_0 subset B)=P B=(P,W),
B/B_0 ~= (k[P,G]/(P))*V.                               (7)
```

Thus `Spec(B)=A2`, so

```text
B^*=k^*,              Cl(B)=0,              H^1_et(Spec(B),mu_3)=0. (8)
```

This statement concerns the **full normalization**.  The original order has

```text
Sing(Spec(B_0))=V(P,W),
U_0=Spec(B_0)_reg=Spec(B_0[P^-1])
   ~=Spec(B[P^-1])~=Gm^2.                              (8a)
```

Indeed, with

```text
l_+=V+delta G,             l_-=V-delta G,
P=-l_+l_-/4,                                           (8b)
```

one has `B[P^-1]=k[l_+^{+-1},l_-^{+-1}]`.  Hence
`H^1_et(U_0,mu_3)~=(Z/3)^2`; conductor-supported Kummer characters can exist
on the nonnormal order's smooth locus even though none extends over `A2`.

The two apparent Cardano directions drop from rank two to rank one on `U_0`
and to rank zero after extension across the conductor.  After
`W=PV` and the substitution in `(6)`,

```text
q_0+W=-(V+delta G)^2(V-delta G)/4,
q_0-W= +(V+delta G)(V-delta G)^2/4,                     (9)

q_1+W=-(V+G)^3/4,
q_1-W= +(V-G)^3/4.                                     (10)
```

Because constants are cubes in `k`, `(10)` is the trivial Kummer class.  The
other radicand is not a cube.  On `U_0`, it is the unit

```text
q_0+W=-l_+^2 l_-/4,                                    (10a)
```

and therefore defines the nonzero etale Kummer class `(2,1)` in
`H^1_et(U_0,mu_3)`.  On the full smooth surface `A2_{G,V}`, its height-one
valuations are `2` and `1` along the two boundary lines

```text
V+delta G=0,                         V-delta G=0.        (11)
```

Consequently its cube-root cover ramifies when extended across the conductor
to `A2`.  It is not a class in `H^1_et(A2,mu_3)`, but it **is** a genuine
smooth-locus character of `B_0`.  This conductor-extension distinction is
essential and was missed in the theorem's first promoted version.

The result is a sharp near miss, not a counterexample construction.  It
classifies the balanced internal split of one repeated factor.  It does not
exclude unequal factors `p_1-p_0=FG`, internal splits in two or three distinct
cube factors, gcd packets not equivalent to a square, or nonnormal orders
whose integral closure has nontrivial class group.

## 1. The repeated factor is split internally

Write

```text
L_i=p_1-omega^i p_0.                                   (12)
```

Here `L_0=G^2`, and direct use of `(1)` gives

```text
q_1-q_0=2G L_1,
q_1+q_0=2G L_2.                                        (13)
```

Their product is

```text
(q_1-q_0)(q_1+q_0)=4G^2L_1L_2
 =4(p_1^3-p_0^3).                                      (14)
```

Thus one copy of the repeated prime `G` goes to each side of the factorization.
This is exactly the internal distribution that the whole-factor hypothesis of
THM-3942 deliberately leaves open.  Squaring the first row in `(2)` and using
`delta^2=-3` gives `(3)`; expansion of the second row gives the same identity.

For the associated depressed cubics

```text
F_i(T)=T^3-3p_iT-q_i,                                  (15)
```

one has `disc_T(F_i)=-27H` for both `i`.

## 2. The nonlinear reduced branch is exactly one-place

The affine curve `(4)` is the graph

```text
P=-3G^2/4,                                              (16)
```

so it is smooth and isomorphic to `A1_G`.  Its projective closure is

```text
4PZ+3G^2=0.                                             (17)
```

On `Z=0`, equation `(17)` forces `G=0`, hence there is exactly one infinity
point `[P:G:Z]=[1:0:0]`.  The derivative with respect to `Z` is `4P`, nonzero
there, so this is one smooth normalization place rather than several branches
coalesced at one target point.

The full divisor `(3)` is nevertheless not reduced: it is twice the line
`P=0` plus `Gamma`.  The line and parabola have no common divisor and meet only
at the affine origin.  Thus one-place geometry has been gained at the cost of
an even discriminant component, not within an irreducible reduced branch.

## 3. Normalization and conductor are exact

The factor `4P+3G^2` is irreducible in `k[P,G]`, so the square class of `H`
is nontrivial and `(5)` is a domain.  In its fraction field, set `V=W/P`.
Equation `(5)` makes `V` integral through

```text
V^2+4P+3G^2=0.                                         (18)
```

It is not already in `B_0`.  As a free `R=k[P,G]` module,

```text
B_0=R direct_sum R W.                                  (19)
```

If `V=a+bW` with `a,b in R`, multiplying by `P` and comparing the two basis
rows in `W=P(a+bW)` would give `a=0` and `Pb=1`, impossible.  Hence `B_0`
is nonnormal.

Adjoining `V` gives the first ring in `(6)`.  Solving `(18)` for `P` identifies
it with the polynomial ring `k[G,V]`; it is finite over `B_0`, has the same
fraction field, and is normal.  Therefore it is the full integral closure.

In the bases `{1,W}` and `{1,V}`, the inclusion `B_0 subset B` has coefficient
matrix

```text
diag(1,P).                                               (20)
```

This immediately gives the quotient in `(7)`.  For completeness, if
`a+bV in B` belongs to the conductor, then it and its product with `V` lie in
`B_0`.  The first condition gives `b in PR`; the second gives `a in PR`.
Thus the conductor is exactly `PB`, which contracts to `(P,W)` in `B_0`.

Geometrically, the normalized finite map is

```text
(G,V) |-> (P,G)=(-(V^2+3G^2)/4,G).                     (21)
```

Its Jacobian is `-V/2`; the ramification line `V=0` maps isomorphically to
the one-place parabola `(4)`.  The doubled conductor line `P=0` pulls back to
the two distinct regular lines `(11)`.

## 4. One character is trivial; one survives only off the conductor

Equations `(9)-(10)` are direct substitutions.  They also preserve the norm
identities

```text
(q_i+W)(q_i-W)=4p_i^3.                                 (22)
```

For `i=1`, choose cube roots of the scalar constants in `k`.  Both Cardano
factors are cubes, and the depressed cubic itself splits in `B[T]` as

```text
F_1(T)=(T+G)
       (T-(G-delta V)/2)
       (T-(G+delta V)/2).                               (23)
```

So this apparent second torus presentation supplies no cyclic field extension
at all after passage to the maximal quadratic order.

For `i=0`, equation `(9)` gives on `A2` the divisor

```text
div(q_0+W)=2[V+delta G]+[V-delta G].                    (24)
```

Both supporting lines lie in the smooth affine plane `Spec(B)`.  Their
coefficients are nonzero modulo three, so adjoining a cube root ramifies at
both corresponding divisorial valuations.  It therefore cannot define an
element of `H^1_et(Spec(B),mu_3)`.  Independently, `(8)` says this cohomology
group is zero: the units are `k^*=(k^*)^3` and `Pic(A2)=0` in the Kummer
sequence.

However, both lines are precisely the preimage of the singular conductor
removed in `(8a)`.  On `U_0~=Gm^2`, `(10a)` is a unit whose exponent vector
is `(2,1)` modulo three.  Thus it gives a nontrivial finite-etale `mu_3`
torsor there.  The two displayed Cardano directions have exact ranks

```text
rank on Spec(B_0)_reg = 1,       rank extending over Spec(B) = 0.  (24a)
```

This separates three notions that are easy to conflate:

```text
not a cube                         -- true for q_0+W;
class on the nonnormal smooth locus-- true, vector (2,1) by (10a);
extension over the normalization  -- false by (24) and (8).             (25)
```

## 5. Reproduction and next escape coordinate

Run

```bash
python3 04-computation/jc2_repeated_factor_double_torus_square_conductor_thm3944.py
python3 -O 04-computation/jc2_repeated_factor_double_torus_square_conductor_thm3944.py
```

The companion verifies both internal factor rows, both common-discriminant
and cubic-discriminant identities, the one-place projective branch ledger,
the normalization and module-index formulas, the conductor preimage and
ramification Jacobian, the exact singular locus, the `Gm^2` conductor
complement, all four radicand factorizations and valuation residues, the
surviving `(2,1)` boundary class, and the complete splitting of `F_1`.

The next live split is consequently not another balanced square.  It is an
**unequal internal factorization** `p_1-p_0=FG` with `F,G` coprime, or a split
distributed across more than one `L_i`, subject simultaneously to three gates:

1. the squarefree part of `H` must remain an irreducible one-place curve;
2. a conductor-supported smooth-locus class must extend across the normalized
   quadratic cover rather than ramify at its missing divisors; and
3. the two Cardano radicands must remain independent after conductor removal.
