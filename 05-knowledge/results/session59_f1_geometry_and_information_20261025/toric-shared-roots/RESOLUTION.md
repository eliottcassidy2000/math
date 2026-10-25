# Resolution restores the Frobenius sidecar that the singular count cancels

**Status: PROVED.** This is a continuation of `RESULT.md`, not an inference
from its point polynomial alone. The toric resolution construction is
classical; the explicit model-specific translation is new to this packet.

Write B=B_2 over a field k of characteristic different from two. Over k(i),

    B = {u v=c^2, w z=c^2} in P^4.

Its four nodes and conjugation action were proved in `RESULT.md`, Theorem 4.

## 1. The resolution, including its arithmetic form

### Theorem R1

Blowing up the reduced singular locus of B gives its smooth minimal
resolution pi:Y->B. Over an algebraic closure there are four exceptional
curves, each P^1 of self-intersection -2. A node with residue field L is
replaced by P^1_L, not by a possibly nonsplit conic.

**Proof.** Each node is defined over a field containing i: it is one of the
four points labelled u,v,w,z. Its local equation over that residue field is
`wz=c^2`. The blowup of the origin has a w-chart with `z=w Z, c=w C`, whose
strict-transform equation is `Z=C^2`. Thus this chart is smooth with
coordinates w,C. The z-chart is analogous. The c-chart has `w=c W,z=c Z`
and equation `WZ=1`, also smooth. The exceptional curve is the projective
conic `WZ=C^2` in P^2_L, with the L-point [1:0:0], and explicitly is P^1_L
via [s:t]->[s^2:t^2:st]. Its normal bundle is the tautological O(-1)
restricted from that conic in P^2, namely O_(P^1)(-2). Therefore its
self-intersection is -2. The centers are disjoint and every other point of
B is smooth, so these four local calculations give a global smooth
resolution over k. There is no exceptional (-1)-curve, so it is minimal.
QED.

This blowup is defined over k even when no individual geometric node is
k-rational: the entire reduced singular closed subscheme is defined over k.

## 2. A second construction from two AP conics

Let

    C={a^2+b^2=2c^2} in P^2.

The conic C is isomorphic over k to P^1, explicitly by

    [s:t] -> [c:a:b]
           =[s^2+t^2 : s^2-t^2+2st : s^2-t^2-2st].

The three displayed quadratics form a basis of the quadratic binary forms
when 2 is invertible, hence this is a Veronese parametrization of the smooth
conic. The points c=0 correspond to t=±i s.

### Theorem R2

The rational map C x C -> B given by

    ([c1:a:b],[c2:d:e])
       -> [c1c2 : a c2 : b c2 : d c1 : e c1]

has precisely the four geometric basepoints c1=c2=0. Blowing them up gives
the same Y as in Theorem R1. The four strict transforms of the fibers
`c1=0` or `c2=0` are the four (-2)-curves contracted to B's four nodes.

Over Q the four basepoints form two degree-two closed points. Consequently
Y is a Q-rational smooth surface, obtained by blowing up P^1 x P^1 in
those two degree-two points.

**Proof.** The displayed coordinates satisfy both defining equations of B.
They vanish simultaneously exactly when c1=c2=0. At one of these points
c1,c2 are local parameters and a,b,d,e are units, so the base ideal is
exactly `(c1,c2)`. Its blowup resolves the rational map. Away from the
boundary, normalize c1=c2=1; the image recovers the original two conic
points. Thus the resulting morphism is birational.

One can identify the contraction completely using the explicit torus charts.
The compactification C x C has the square fan with axis rays
`(±1,0),(0,±1)`. Blowing up its four torus fixed points adds the four
diagonal rays `(±1,±1)`. The cross-polytope model of B has the fan with
only those four diagonal rays. Thus the resolved map is the toric
contraction that removes the four axis rays. The fan description holds
after adjoining i and is equivariant under central inversion, so it
descends to the original k-form. Equivalently each original boundary fiber
has self-intersection zero and meets two basepoints; its strict transform
has self-intersection -2. These are exactly the four contracted curves.
The blowup from Theorem R1 is the same minimal resolution.

For clarity, the linear system has bidegree (2,2) under C x C≅P^1 x P^1,
not bidegree (1,1): each conic embedding has degree two. On the blowup its
class is `2H1+2H2-sum E_j=-K_Y`; its square is 4, agreeing with B's degree.

Conjugation simultaneously changes both boundary signs. It has two orbits
of size two on the four basepoints. The conics have the displayed rational
parametrization over Q, so Y is birational over Q to P^1 x P^1. QED.

The four blowup basepoints and the four singular nodes are different labelled
sets. Their Galois permutations have the same cycle type, but the map does
not identify the basepoints with the nodes: the strict transforms of boundary
fibers map to nodes, while the basepoint-exceptional (-1)-curves map to the
four boundary lines. Matching their counts is not a pointwise identification.

## 3. Extension counts and zeta function after resolution

### Theorem R3

Fix an odd finite field F_q and put epsilon=chi_q(-1). For every n>=1,

    #Y(F_(q^n))=q^(2n)+(4+2 epsilon^n)q^n+1.

Therefore

    Z(Y/F_q,t)=
      1/((1-t)(1-q^2t)(1-qt)^4(1-epsilon q t)^2).

In particular the split and nonsplit resolutions are distinguished by their
extension counts, even though B and its split form have identical extension
counts.

**Proof.** The node subscheme has `2(1+epsilon^n)` F_(q^n)-points. Each such
node is replaced by a P^1 over the same field, adding q^n to the count.
Theorem 4 of `RESULT.md` gives `#B=(q^n+1)^2`. Hence

    #Y-#B=2q^n(1+epsilon^n),

which is the formula. Exponentiating its power sums gives the zeta function.
No cohomological trace formula is needed for this derivation. Independently,
Theorem R2 gives the same count by blowing up P^1 x P^1 at the four
geometric basepoints, whose Galois permutation also consists of two
transpositions in the nonsplit case. QED.

For example at q=3, `#B=16` and `#Y=16`: neither node nor basepoint is
rational, so there is no rational-point count change. At q=9, `#B=100`
and `#Y=136`. These two values exhibit the extension-level sidecar that a
single-base-field resolution comparison might miss.

## 4. The Picard lattice explains the two Artin factors

### Proposition R4

Let tau be the nontrivial Galois action of k(i)/k in the nonsplit case.
The rational geometric Picard representation of Y is

    Pic(Y_kbar) tensor Q ≅ Q_triv^4 direct_sum Q_sign^2.

Thus tau^n has trace `4+2(-1)^n`. For a finite field the corresponding
Frobenius trace is `4+2 epsilon^n`, exactly the coefficient of q^n in
Theorem R3.

**Proof.** The smooth complete fan of Y has eight rays, cyclically

    (1,0),(1,1),(0,1),(-1,1),
    (-1,0),(-1,-1),(0,-1),(1,-1).

Consecutive primitive vectors have determinant one, so the fan is smooth.
Galois acts by central inversion, exchanging the rays in four pairs.

For this toric surface the boundary-divisor sequence is

    0 -> M -> Z^8 -> Pic(Y_kbar) -> 0.

Here M is the rank-two character lattice, whose map sends a character to
its orders along the eight boundary divisors. To see the sequence directly,
the dense torus has trivial divisor class group, so every divisor class has
a representative supported on the boundary. A principal divisor supported
there comes from a rational function invertible on the torus; its restriction
is a constant times a Laurent monomial, hence comes from M. Smoothness
identifies divisor class group with Picard group.

On Q^8 the involution has four trivial and four sign summands. On M tensor Q
it is -I, hence two sign summands. The quotient therefore has four trivial
and two sign summands, proving the claim. The invariant rational rank is
four in the nonsplit case and six in the split case. QED.

This is an explicit positive transfer: the boundary fan yields the Picard
representation, which predicts the resolution's two Artin sign factors.
The singular total zeta had canceled those factors. No analogous recovery
from the singular total count alone is possible without the boundary/node
sidecar.

## 5. Connection contract and limits

- **Source -> target:** marked singular toric form B -> its canonical
  node blowup Y. The map pi is a proper birational resolution.
- **Preserved:** dense torus, function field, Galois descent, and all square
  labels as rational functions pulled back along pi.
- **Changed:** four points become curves; total extension counts acquire
  a q^n-weighted node contribution.
- **Needed sidecar:** residue fields of the four nodes; geometric node
  number alone would incorrectly predict a uniform addition of 4q^n.
- **Decisive test:** q=3 then q=9, or the Picard involution trace 2 then 6.
- **Boundary:** rationality of B or Y does not solve the simultaneous added
  square-root conditions on H± or M. The interior-branch obstruction for
  the actual covers remains. The construction is over characteristic
  different from two; it is not a resolution of the nonreduced bad fiber.

`check_resolution.py` independently checks the blowup/fan arithmetic,
Galois traces, symbolic parametrization and finite extension count identities.
