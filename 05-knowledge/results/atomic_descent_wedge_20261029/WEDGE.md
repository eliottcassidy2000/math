# Three counting dimensions, four ambient coordinates, six bivectors

**PROVED elementary matrix identities and obstruction.** The tournament
grammar and its audited all-order observer criterion are inherited from
[actual grammar refinements](../session103_words_necklaces_20261029/languages/ACTUAL_GRAMMAR_REFINEMENTS.md).
This is not a new all-tournament enumeration formula.

## 1. The count at thirteen has an arithmetic structure

The productive count matrix, initial row and terminal column are

    C=[[1,0,1],[1,1,1],[1,1,0]], alpha=(0,1,0), beta=(1,1,1)^t.

For tournament order n=m+4 the count is v_m=alpha C^m beta. Thus

    v_0,v_1,v_2=1,3,7;
    v_m=2v_(m-1)+v_(m-2)-v_(m-3).

The proved grammar criterion makes v_9=2065 an exact observer count at
order 13; v_6=182 is the count at order 10. This supersedes the earlier
finite-only interpretation in operation_closed_synthesis, whose state
minima and matrix calculations remain dependencies.

The characteristic polynomial is `x^3-2x^2-x+1`, of discriminant 49.
Put u=x-1. The shifted polynomial is

    f(u)=u^3+u^2-2u-1.

The identity `z^3 f(z+z^-1)=1+z+...+z^6` identifies its roots with
`2 cos(2 pi j/7)`, j=1,2,3. Irreducibility follows from the rational-root
test. Therefore the counting spectrum lies in the real cubic seventh-root
field. Within this field

    sigma(u)=u^2-2,  sigma^3(u)=u,
    (1+u)(u^2-2)=-1.                                  (7)

This supplies an operation, not just a shared seven. Modulo seven,
`char_C=(x-3)^3`, so `(C-3I)^3=0` and `C^7=3I` modulo seven.
Using the actual boundary vectors gives

    v_m = 3^m (1-binomial(m,2)) mod 7.

Indeed, for N=C-3I, the three boundary moments are
`alpha beta=1`, `alpha N beta=0`, and `alpha N^2 beta=-2`.
The nilpotent binomial expansion proves the formula for every m, rather
than only the checked initial range.

Consequently the observer count at order n is divisible by seven exactly
when `n=3 or 6 mod 7` (n>=4). The divisibility of 2065 and 182 follows
from this law. It does not follow from the index 13 being prime. Separately
`C^12=I mod13`; reduction loses integer counts and cannot reconstruct them.

## 2. An exact exterior-square lift

In the ordered basis `(e1 wedge e2,e1 wedge e3,e2 wedge e3)`,

    W=Lambda^2 C=[[1,0,-1],[1,-1,-1],[0,-1,-1]].

The integral matrix

    S=[[-1,0,-1],[0,0,-1],[0,1,-1]], det S=-1,

satisfies `W S=S(C-I)`. This is integral conjugacy, including the lattice,
not just agreement of eigenvalues. In three dimensions the general identity
is `Lambda^2 C ~= det(C) C^-T` under volume duality. Here det C=-1;
equation (7) explains why its reciprocal spectrum is the shifted spectrum.

Now take the explicitly defined lift `G=diag(1,C)` on four coordinates.
Its exterior square has six coordinates, naturally indexed by K4's edges:

    Lambda^2(R e0 direct_sum R^3)
      = (e0 wedge R^3) direct_sum Lambda^2 R^3,
    Lambda^2 G = C direct_sum W ~= C direct_sum (C-I).   (8)

This is a faithful algebraic connection from the aggregate counting matrix
to a six-dimensional bivector carrier. The first block, with its original
alpha,beta, recovers every count. The lift does not identify the four
coordinates with the four acceptance states or the six coordinates with
the legality states.

## 3. Why the six-state machine is not this exterior square

Use acceptance states `(A,B0,B1,F)`. On letter e the images are
`(F,B1,F,F)`, so its linearized transition has rank two and its exterior
square has rank one. For the legality-aware states `(A,B0,B1,FA,FB,X)`,
the images under e are `(X,B1,FB,X,FB,X)`, of rank three. Conjugate linear
representations would have equal ranks for each letter. Thus no linear
isomorphism identifies these two natural six-dimensional constructions.

Counting, acceptance and legality have different future observations.
The inherited minima 3,4,6 remain valid, but the numerical pattern alone
does not commute with the named letter operations. The scalar lift (8)
is useful precisely with this boundary stated.

## 4. Edge signs, Hodge stars, and fields of bivectors

The six components of a decomposable bivector `p=x wedge y` satisfy

    p01 p23 - p02 p13 + p03 p12 = 0.                   (9)

This follows by expanding its two-by-two minors. Conversely a nonzero
six-tuple over a field satisfying (9) is decomposable: on a nonzero
coordinate chart solve for the two rows, and (9) supplies the remaining
minor. Thus the three opposite-edge matchings of K4 index the three terms.

In oriented Euclidean three-space, a cross product is a wedge followed by
metric/volume duality. In oriented Euclidean four-space, the Hodge star
instead splits bivectors into two three-dimensional spaces. This splitting
requires the chosen metric and orientation. Our G is not orthogonal and
does not preserve those Hodge eigenspaces. It preserves the wedge pairing
up to det G=-1; the checker verifies `B^t J B=-J` for B=Lambda^2 G.

Unsigned K4 edges and signed exterior coordinates also have different S4
actions. A vertex transposition has trace two on unsigned edges, but zero
on bivectors: one fixed edge reverses orientation. Hence their rational
representations are not isomorphic. The old
[THM-3045 unsigned-edge clutch](../../../01-canon/theorems/THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch.md)
has decomposition 1+2+3 and integral index 24. It cannot be silently replaced
by a signed Hodge 3+3 decomposition.

Finally, a bivector field adds differential obligations. For

    P=partial_x wedge (partial_y+x partial_z)

on R4, (9) holds at every point. Nevertheless the bracket has
`{x,y}=1,{x,z}=x`, and Jacobi(x,y,z)=1. Pointwise decomposability therefore
does not imply the Poisson condition. Likewise a sum of decomposable
bivectors need not be decomposable: `e0 wedge e1+e2 wedge e3` violates (9).
These are concrete tests before transferring a small algebraic gadget to
weighted composition or planar-Jacobian dynamics.

The remaining constructive question is to find an actual letter-compatible
geometric carrier for the six-state machine, or a controlled quotient of
one. The rank-one/rank-three e test rules out this specific naive choice.
