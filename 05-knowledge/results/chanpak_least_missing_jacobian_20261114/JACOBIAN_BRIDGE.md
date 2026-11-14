# Jacobian comparison: rigidity of an affine sweep, and the missing coherence in an integer lift

Two distinct statements are involved. The classical **complex**
Jacobian conjecture demands a nonzero *constant* determinant.
The higher-dimensional premise is correct: Alpöge's 2026
three-dimensional example has constant determinant and three
points in a single fiber, and adjoining identity coordinates gives
counterexamples in larger dimensions. Gao's [research preprint,
especially §3](https://arxiv.org/html/2608.00222v2#S3) verifies
the construction in a tangent-sweep description and supplies
different higher-dimensional examples. The two-dimensional
constant-determinant problem is not disposed of by the example.
On the other hand the **strong real** conjecture asking for an
everywhere nonzero, not necessarily constant, determinant already
fails in the plane by Pinchuk's examples; see the explicitly
constructed planar examples in [Braun and Fernandes dos Santos,
The Pinchuk example revisited](https://doi.org/10.5802/crmath.570).
No dimensional analogy should merge these two hypotheses.

## A polynomial affine sweep cannot be a Keller counterexample, in any dimension

Let n≥2 and let F:C×C^(n−1)→C^n have the form
F(s,t)=a(s)+Σ_j t_j b_j(s), with polynomial vector functions.
If det DF is a nonzero constant, then F is a polynomial
automorphism. In particular a direct planar tangent-line sweep
with genuinely variable slope cannot have constant nonzero
Jacobian. The same conclusion is valid for an affine-moving-
hyperplane sweep in higher dimensions, so just adding another
linearly swept coordinate is not the counterexample mechanism.

Proof. The n−1 columns b_j(s) have full rank at every s.
The determinant is affine in t, with the coefficient of t_j
equal to det(b_j',b_1,...,b_(n−1)). Each of these is zero,
so each b_j' lies in the hyperplane they span over C(s).
Let ν(s) be the nonzero covector v↦det(v,b_1,...,b_(n−1)).
Differentiating its exterior product shows ν'=r(s)ν for
a rational function r. Thus the ratio of any two nonzero
coordinates of ν has derivative zero and is constant. Consequently
ν(s)=q(s)ν_0, for a polynomial q and a fixed nonzero covector ν_0.
Now det DF=q(s)ν_0(a'(s)) is a nonzero constant. Both factors
are nonzero constants. Hence ν_0(F(s,t)) is an affine-linear
function of s with nonzero slope and it recovers s polynomially.
After using this fixed target linear coordinate, the remaining
linear system in t has a polynomial coefficient matrix with
nonzero constant determinant. The adjugate gives a polynomial
inverse.

This is an elementary rigidity *subclass*, not a proof of JC in
dimension two. The planar one-fiber literature is already present
in [THM-2063](../../../01-canon/theorems/THM-2063-one-fiber-linear-planar-keller-pairs.md).
The useful difference here is to show that passive addition of more
affine fiber dimensions cannot do the job either. Gao and Speyer's
counterexample uses a genuine twist: the tangent-sweep description
allows a nonconstant factor in intermediate coordinates and
monomial changes cancel it in the final Jacobian. A plane chart
of this affine kind has no such compensation. The repo's
[non-slice theorem](../../../01-canon/theorems/THM-4562-slice-lemma-and-the-jacobian-counterexample-is-fully-non-slice.md)
expresses a stronger, different obstruction to merely slicing off
one coordinate, by comparing exceptional and generic fibers.

A compatible topological diagnostic is properness: a proper
local diffeomorphism of C^n to C^n is a covering, since it has open
and closed image and evenly covered neighborhoods. As C^n is
simply connected, a connected covering is one sheet. The known
counterexamples are nonproper, so local invertibility is valid
without a globally valid inverse. This does not assert that any
equation for a missing tournament count is a Keller map.

## Lattice separation: a one-corner determinant is not a covering test

Let F(k,t)=a_k+t b_k for k,t∈Z, with a_k,b_k∈Z², and assume
for every k,t that
$$
 det(F(k+1,t)-F(k,t),,F(k,t+1)-F(k,t))=1.
$$
Comparing coefficients of t forces b_(k+1) parallel to b_k.
The determinant condition also makes each b_k primitive.
Thus after an oriented integer linear coordinate change there
is ε_k∈{−1,1} with b_k=ε_k(0,1); the first coordinate x_k
of a_k satisfies x_(k+1)−x_k=ε_k. Conversely these data give
exactly the displayed one-corner determinant condition,
regardless of the second coordinates of a_k.

If the orientation sign is coherent (all ε_k equal), the map
is a bijection of Z²: x_k recovers k and the second coordinate
recovers t. If instead ε_k=(−1)^k, the map
$$
 F(k,t)=(kmod 2,,(-1)^k t)
$$
passes the same test but has infinite fibers and omits all
but two columns. This example is **not** a locally bijective
map of the square lattice graph: both horizontal neighbors
can fold onto one column. The test is too weak to imply that.
Demanding determinant +1 at the opposite corner too, i.e.
using b_(k+1) in the determinant with the same horizontal edge,
forces ε_k=ε_(k+1) and removes the fold.

## Precise scope of the Chan–Pak analogy

For H>1 a tournament can be trimmed by deleting its singleton
strong components without changing H, since its SCCs are ordinal
and H multiplies over them. In this reduced domain, Busch's
strong-component minimum makes the order O(log H): each nontrivial
component costs exponential growth in its order, so their sizes
sum to O(log H). For H=1 one can use the singleton. Thus even the
escape-to-infinity intuition has a boundary here: fixed bounded
targets admit only finitely many trimmed isomorphism types, after
forgetting any marks not preserved by trimming. Arithmetic holes
can remain despite this form of properness.

The useful transfer is an **oriented lift** question, not a theorem
that congruence or local switches suffice. If N is a least odd
unrealized count other than 7 and 21, with N>23, then N−2 is
realized. Every graph at N−2 fails to have *any* allowed forward
unit move into N, since such a move would construct N. Showing
merely that many states carry a ±2 switch, or that a switch exists
somewhere in a projected profile, does not eliminate a boundary
like this. The symmetry that matters is whether the same physical
state can continue in the **increasing** orientation and whether
the hypotheses are regained after gluing. The alternating-fold
example makes the distinction exact in a tractable model.

Even a numerical spectrum abstractly satisfying multiplicative
closure, every local residue condition compatible with oddness,
and arbitrarily long consecutive windows need not have only the
known two holes. For any arbitrarily large odd prime p≠7 the set
of positive odds with exactly 7,21,p removed is multiplicatively
closed: 21=3·7 cannot be a product of its elements and neither
missing prime can be a nontrivial product. It is cofinite among
odds so it passes all residue and window tests. This is not a
candidate actual tournament spectrum; it rules out a proof based
only on those abstract data. A constructive contradiction needs
another physical invariant or coherent gluing rule that
disallows such an isolated prime hole.
