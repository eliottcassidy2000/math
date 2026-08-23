---
id: THM-3844
title: "A two-cusp polynomial branch globalizes one anti-invariant three-class, but folds back to a monogenic cubic"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  An irreducible polynomially normalized quartic with two
  A2 cusps, one A1 node, and one place at infinity has quadratic-resolvent
  surface class group Z plus Z/3.  The three-class is primitive at both
  cusps, trivial at the node, and anti-invariant under the quadratic deck
  involution; its cyclic Kummer cover is explicit.  Nevertheless the branch
  is a polynomial pullback of the depressed-cubic cusp, the associated S3
  cubic is globally monogenic, and deleting its ramified sheet creates a
  nonconstant unit.  Thus this exact one-place design has no dominant plane
  atlas and is not a Jacobian counterexample.
source: jc_zero_debt_lift / one-place S3 completion design lane, 2026-08-23
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
related:
  - THM-3840-forced-cubic-two-arm-jelonek-passport
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3842-nonlinear-cubic-tower-trace-shift-eightfold-base-change
script: 04-computation/jc2_two_cusp_polynomial_branch_resolvent_thm3844.py
output: 05-knowledge/results/jc2_two_cusp_polynomial_branch_resolvent_thm3844.out
script_sha256: 23121e393bfcfc88ee1c0370636451291210ccf4cec5e8f4b4fbf916a5edff72
output_sha256: 2670b141ddbdd2b7a20b915ff2c9dee35086456f534fc1feae6e042c01e60ec2
semantic_sha256: 70efbb25eb7dfb096a1f582148a2c234586ec72dfc22f02976df5c133cc90e99
hash_basis: raw LF bytes
---

# THM-3844 -- two cusps produce one three-class, not two

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**  Work over an algebraically closed field `k` of
characteristic zero.  This theorem records a counterexample-design hostile,
not a Keller map.

Put

```text
delta=27X^4+16X^3+6X^2Y+48XY^2-16Y^3+27Y^2.                  (1)
```

Then the plane branch `C=V(delta)` is irreducible, has polynomial
normalization

```text
X=t^2(2t-3),                     Y=t^3(3t-4),                 (2)
```

and has exactly two `A2` cusps and one `A1` node in the affine plane.  Its
projective closure has one place at infinity, which is smooth.  The normal
quadratic-resolvent surface

```text
Q=Spec k[X,Y,W]/(W^2-delta)                                  (3)
```

satisfies

```text
Q*=k*,                            Cl(Q)=Z direct_sum Z/3.      (4)
```

The torsion generator is nonzero in the local class group at each cusp,
zero at the node, and is sent to its inverse by `W |-> -W`.  It therefore
does globalize the two local cusp directions to a cyclic cubic Kummer cover
of `Q`, with opposite local orientations.  However, that cover is the
quadratic-resolvent layer of the globally monogenic cubic

```text
F(T)=T^3+(X-Y)T+(X^2+Y)/2.                                   (5)
```

Deleting the ramified sheet of `(5)` makes `3T^2+X-Y` a nonconstant unit.
Consequently its maximal etale open admits no dominant morphism from the
polynomial plane.  The one-place branch geometry and the anti-invariant
three-class are therefore necessary-looking but insufficient design data.

## 1. Exact normalization and singularity ledger

Substitution of `(2)` in `(1)` gives zero.  Conversely, eliminating `t`
from

```text
2t^3-3t^2-X=0,                    3t^4-4t^3-Y=0                (6)
```

gives exactly `delta`, and

```text
t=(-X^2+4XY+3Y)/(6X^2+4X+2Y)                                  (7)
```

on a dense open of `C`.  Moreover

```text
t^3-(3/2)t^2-X/2=0.                                           (8)
```

Thus `t` is integral over `k[C]`, `(7)` identifies the fraction fields,
and the normal ring `k[t]` is the integral closure of `k[C]`.  The exact
factorization check in the companion shows that `delta` is irreducible.

The normalization-conductor packet is equally explicit.  Put

```text
c=t^2(t-1)^2(2t^2-2t-1).                                     (8a)
```

The two squared factors are the cusp conductors and the last factor gives
order one at each of the two node addresses.  Hence the conductor in `k[t]`
is `c k[t]`.  Directly,

```text
c    =(3X^2+2X+Y)/6,
tc   =-(X^2-4XY-3Y)/12,
t^2c =-(X-Y)(X+2Y)/9.                                        (8b)
```

Since `(8)` makes `1,t,t^2` generate `k[t]` over `k[C]`, `(8b)` verifies
`c k[t] subset k[C]` without a local-completion inference.  The complete
normalization-conductor sequence is

```text
0 -> k[C]/c k[t] -> k[t]/c k[t] -> k[t]/k[C] -> 0,            (8c)
```

with lengths `3,6,3`; the last length is the sum of the two cusp delta
invariants and the node delta invariant.

The singular ideal has reduced support encoded by the Groebner basis

```text
Y-8X^4-16X^3-7X^2,
X^2(X+1)^2(2X+1).                                              (9)
```

Hence the three and only three affine singular points are

```text
(0,0),                   (-1,-1),                   (-1/2,1/4). (10)
```

At `t=0`,

```text
X=-3t^2+2t^3,                    Y=-4t^3+3t^4,                 (11)
```

and at `t=1+s`,

```text
X+1=3s^2+2s^3,                  Y-2X-1=4s^3+3s^4.             (12)
```

These are the two `A2` cusps.  The node has the two distinct normalization
addresses

```text
2t^2-2t-1=0,                                                       (13)
```

whose tangent slopes are `2t`, hence are unequal.

Let `delta_h` be the degree-four homogenization.  On the line at infinity,

```text
delta_h(X,Y,0)=27X^4.                                          (14)
```

There is therefore only the point `[0:1:0]`; it is smooth because
`partial_Z delta_h=-16` there.  This proves the claimed one-place property,
not merely a one-point set-theoretic count.

## 2. The quadratic surface and its affine class group

The affine surface `(3)` is normal: it is a hypersurface and hence `S2`,
while the reduced branch makes it regular in codimension one.  Its only
singularities lie over the three codimension-two points in `(10)`.

The weighted-projective completion of `(3)` is the double plane

```text
Qbar: W^2=delta_h             in P(1,1,1,2).                  (15)
```

It has surface singularities `2A2+A1` over `(10)` and is smooth at infinity.
Its minimal resolution `S` is the standard rational weak del Pezzo surface
of degree two.  Thus

```text
Pic(S)=ZH direct_sum ZE1 direct_sum ... direct_sum ZE7,
(H^2,E1^2,...,E7^2)=(1,-1,...,-1),
-K_S=3H-E1-...-E7.                                            (16)
```

The inverse image of the line at infinity splits because `(14)` is a square:

```text
B_+ + B_-=-K_S,              B_+^2=B_-^2=-1,
B_+ B_-=2.                                                       (17)
```

The five exceptional roots over the affine ADE points form `2A2+A1` and
are disjoint from both boundary curves.  The class group of the affine
normal surface is therefore

```text
Cl(Q)=Pic(S)/<the five ADE roots,B_+,B_->.                      (18)
```

There is no marking assumption hidden in the answer.  The companion lists
all 126 `E7` roots, fixes one boundary `(-1)`-curve, obtains the 72
orthogonal `E6` roots and all 120 `A2` subsystems, and exhausts all 360
unordered compatible `2A2+A1` configurations.  For every configuration the
seven relation rows have rank seven and the gcd of their maximal minors is
three.  Hence `(18)` is uniformly `Z direct_sum Z/3`.

One marking makes the mechanism transparent.  Take

```text
r1=E1-E2,                 r2=E2-E3,
r3=E4-E5,                 r4=E5-E6,
r5=H-E1-E2-E3,
B_+=E7,
B_-=3H-E1-E2-E3-E4-E5-E6-2E7.                                (19)
```

The quotient relations reduce to

```text
H=3x,       E1=E2=E3=x,       E4=E5=E6=y,       E7=0,
3(2x-y)=0.                                                       (20)
```

Thus a torsion generator is

```text
tau=2E1-E4,                                                       (21)
```

and the exact integral relation is

```text
3tau=B_-+4r1+2r2-2r3-r4-3r5+2B_+.                            (22)
```

The Smith diagonal is `(1,1,1,1,1,1,3)`.  The companion additionally
constructs the saturated torsion lift for every one of the 360 markings,
not just `(19)`.

Finally, the divisor exact sequence for `S minus (B_+ union B_-)` shows

```text
Q*/k* = ker(ZB_+ direct_sum ZB_- -> Pic(S)).                   (23)
```

The two vectors in `(19)` are linearly independent in the torsion-free
lattice `Pic(S)`, so the kernel is zero.  This proves the unit assertion in
`(4)`.

## 3. Local addresses and the cyclic Kummer cover

For an ordered `A2` root basis, the local discriminant character sends an
intersection pair `(a,b)` to `a+2b mod 3`.  With `(r1,r2)` at the first cusp
and the oppositely oriented `(r4,r3)` at the second, `(21)` gives

```text
first cusp:  (-2,0) |-> 1 mod 3,
second cusp: (0,1)  |-> 2 mod 3,
node:         tau r5=2 |-> 0 mod 2.                            (24)
```

Thus the global class is primitive at both cusps and, after choosing the
opposite second orientation, is anti-diagonal.  The exhaustive saturation
test proves the marking-independent statement: the unique three-class is
nonzero at both `A2` points and zero at the `A1` point for all 360
configurations.

Let `sigma` denote the deck involution of `(15)`.  On the degree-two del
Pezzo lattice,

```text
sigma(D)=(D K_S)K_S-D.                                        (25)
```

For `(21)`, this gives

```text
sigma(tau)+tau=-K_S=B_++B_-.                                  (26)
```

Therefore `sigma(tau)=-tau` in `Cl(Q)`.  Since `tau` has exact order three,
the reflexive algebra

```text
O_Q direct_sum O_Q(-tau) direct_sum O_Q(-2tau)                (27)
```

with a chosen trivialization of `3tau` defines, after normalization, a
cyclic cubic cover of `Q` which is etale in codimension one.  Equation `(26)`
extends the quadratic involution by inversion, so the composite function-
field extension over `k(X,Y)` has group `C3 semidirect C2=S3`.

This answers the globalization question positively.  The next section shows
why it is nevertheless the wrong globalization for a Keller surface.

## 4. The two cusp addresses fold to one depressed cusp

Define polynomial coefficient functions and the normalization fold

```text
p=X-Y,                 q=(X^2+Y)/2,                 h=t(t-1). (28)
```

On `(2)` one has

```text
p=-3h^2,                         q=2h^3,                       (29)
```

and globally

```text
-4p^3-27q^2=-delta/4.                                            (30)
```

Thus the apparently two-cusp quartic is a polynomial pullback of the single
depressed-cubic cusp.  The map `t |-> h=t(t-1)` identifies the two cusp
addresses `t=0,1`; it also identifies the two branches of the node.  On the
coefficient plane this loss of address is visible in

```text
Jac_(X,Y)(p,q)=X+1/2,                                          (31)
```

which vanishes at the node.

Now take the cubic `(5)`.  Its discriminant is `-delta/4`.  It is irreducible:
at `X=0`, twice the polynomial is

```text
2T^3-2YT+Y,                                                     (32)
```

which is Eisenstein at `Y`; any factorization of the monic polynomial `(5)`
would specialize to a monic factorization of `(32)`.  Because `delta` is
irreducible and not a square, the generic Galois group is `S3` and its
quadratic resolvent is `(3)`, after a harmless constant rescaling of `W`.

The power-basis algebra

```text
A=k[X,Y,T]/(F)                                                   (33)
```

is already normal.  Indeed it is finite free over the regular ring `k[X,Y]`
and hence `S2`.  Away from `delta` it is etale.  At the height-one prime
`(delta)`, the power-basis discriminant has valuation one; the square-index
formula therefore forces index zero in the integral closure.  Hence `(33)`
is `R1`, so Serre's criterion proves normality.  This is a genuinely global
monogenic cubic, not merely a field presentation.

The cyclic layer in Section 3 is also explicit.  Over `(3)` put

```text
a=-q/2+W/(12 sqrt(3)).                                         (34)
```

The exact Cardano norm identity is

```text
a sigma(a)=-(p/3)^3.                                           (35)
```

Adjoin `u` with `u^3=a`; then

```text
T=u-p/(3u)                                                     (36)
```

satisfies `(5)`.  Irreducibility makes this a nontrivial cyclic cubic over
`k(Q)`.  It is unramified at every codimension-one point of `Q`: the only
codimension-one inertia over the coefficient plane is a transposition, and
the quadratic resolvent absorbs it.  Since `Q*=k*`, there is no unit Kummer
class, while `(4)` has only the two nonzero inverse three-classes.  Therefore
the cover `(34)` is exactly the cover of `(27)`, up to replacing `tau` by
`-tau`.

## 5. The different becomes a forbidden unit

In `(33)` let

```text
d=F_T=3T^2+p,                         c=-3T^2-4p.               (37)
```

Direct reduction modulo `F` gives the companion identity

```text
d^2 c=disc_T(F)=-delta/4.                                    (38)
```

Thus the ramified sheet is exactly `V(d)` and the simple companion is
`V(c)`.  The maximal etale open of the cubic map is

```text
U=Spec A[1/d].                                                  (39)
```

The element `d` is a nonconstant unit in `Gamma(U,O_U)`.  If a dominant
morphism `A2 -> U` existed, the induced map on coordinate rings would be
injective and would send `d` to a unit of `k[x,y]`, hence to a scalar.  This
contradicts injectivity because `d` is nonconstant.  In particular `(39)`
has no polynomial-plane atlas and cannot underlie a planar Keller
counterexample.

The conclusion is deliberately narrow.  It does not say that every
one-place branch, every nonmonogenic `S3` completion, or every
anti-invariant class fails.  It says that **one torsion direction shared by
two cusps can still be only one depressed-cusp direction in disguise**.
The next positive design should therefore seek a polynomially normalized
branch which defeats every polynomial factorization `(p,q)` through
`-4p^3-27q^2`, plausibly with torsion rank at least two so that local cusp
classes cannot all lie on one global line.

## 6. Exact companion and audit boundary

Reproduce with

```bash
python3 04-computation/jc2_two_cusp_polynomial_branch_resolvent_thm3844.py
python3 -O 04-computation/jc2_two_cusp_polynomial_branch_resolvent_thm3844.py
```

Both runs byte-match the frozen output.  The assertion-free companion checks
the eliminant, rational inverse, integral normalization, complete singular
locus, local analytic leading terms, node addresses, smooth unique infinity,
the normalization conductor, all 126 `E7` roots, all 120
boundary-orthogonal `A2` systems, all 360
compatible `2A2+A1` configurations, every torsion order and local saturation
pattern, the explicit Smith relation and deck action, the depressed-cusp
pullback, Cardano norm, irreducibility hostile, and different/companion
identity.  It reports 1,859 active exact gates.

The use of the standard weak degree-two del Pezzo resolution model in
Section 2 is the human geometric input; the finite lattice consequences are
exactly enumerated.  Independent hostile audit is still required before
promotion to proved canon.  **QED candidate.**
