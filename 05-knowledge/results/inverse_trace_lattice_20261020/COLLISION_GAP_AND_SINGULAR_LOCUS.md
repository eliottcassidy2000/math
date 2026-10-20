# An exact collision-gap unit and the singular-locus projection gate

20 October 2026. **Scoped proof; independent review pending.** This is a separate bounded addendum. The three completed trace, finite-selector, and source-line packets remain unchanged. The collision identity was proposed by the coordinating researcher.

## 1. Setup and the concrete objects

Retain the actual complex polynomial Keller map F=(P,Q), J=P_xQ_y−P_yQ_x∈C*, primitive affine source coordinate x, and P monic in its companion y. For the boundary statements, retain the global low-trace hypothesis from `FINITE_SELECTOR_EULER_GATE.md`, and its finite graph

    G(x,y)=(P(x,y),Q(x,y),x),
    W=V(R)⊂A³_(u,v,Z),
    R=H Z^d+Hc₁Z^(d−1)+…+Hc_(d−2)Z²+A Z+B.

Thus G is finite birational and surjective, V(H,A,B)=∅, the nonproperness curve is C=V(H), and every F-fiber on C has at most one point. The reduced hypersurface W is the actual graph image, not an arbitrary inverse-polynomial model.

Introduce independent variables a,s,t and write Δ=s−t. Define polynomial divided differences

    U(a,s,t)=[P(a,s)−P(a,t)]/Δ,
    V(a,s,t)=[Q(a,s)−Q(a,t)]/Δ,
    I=(U,V)⊂C[a,s,t],    T=Spec(C[a,s,t]/I).

The coordinate a specifies the source line x=a; s,t specify an ordered pair on that line. The quotient notation is exact polynomial division, including its values on Δ=0.

## 2. The collision gap is an explicit unit

Set

    Ω=[P_x(a,s)V−Q_x(a,s)U−J]/Δ.                 (1)

This is a polynomial. Indeed, on t=s the values of U,V are P_y(a,s),Q_y(a,s); the numerator in (1) becomes J−J=0. The factor theorem gives divisibility by s−t.

Consequently the following identity holds in the polynomial ring:

    J=P_x(a,s)V−Q_x(a,s)U−ΔΩ.                  (2)

In the full quotient C[a,s,t]/I, even if it is nonreduced,

    Δ·(−Ω/J)=1.                                (3)

Therefore the gap Δ is a unit, I:Δ^∞=I, and T has no diagonal point. Its complex points are exactly the ordered distinct pairs on the same source line with equal F-image. There is no hidden diagonal component to remove by saturation.

This is an algebraic unit certificate, not an analytic lower bound on |s−t|. A unit can approach zero along a noncompact curve. The identity does not establish that T is empty.

## 3. A finite bridge to the actual singular graph locus

Let Σ=Sing(W), with its reduced structure. The map

    Φ:T_red→W,
    (a,s,t)↦(P(a,s),Q(a,s),a)                  (4)

is finite. If P has y-degree n, U is monic of degree n−1 in t. Thus T is finite over the source plane with coordinates (a,s); compose with finite G. If n=1, U=1 and T is empty, consistently with this argument.

The image of Φ is exactly Σ. To see both directions, use the following elementary property of the actual finite graph:

    w∈W is smooth ⇔ G^(-1)(w) has one point.     (5)

For the forward direction, a finite birational map from the normal source A² is the normalization of W and is an isomorphism over its smooth locus. For the reverse direction, a singleton fiber and finiteness put the full inverse image of a sufficiently small analytic neighborhood inside one source inverse-chart neighborhood. The Keller inverse-function theorem then makes the entire germ W_w a single smooth graph Z=x(F^(-1)(u,v)).

A point of T gives two distinct source preimages of the same w, so its image is singular. Conversely, every singular w has two such preimages by (5); their x-coordinates equal Z(w), so they give a point of T. Hence (4) is finite and surjective onto Σ after reduction. It preserves the scalar parameter exactly:

    a = Z∘Φ.                                   (6)

This is the promised concrete source-to-target map: ordered source-line collisions map finitely onto graph singularities, preserving nonemptiness of every parameter fiber.

There is also an exact equation for Σ. Differentiating R(P,Q,x)=0 gives

    R_u P_x+R_v Q_x+R_Z=0,
    R_u P_y+R_v Q_y=0,                          (7)

where the derivatives of R are evaluated at G(x,y). Since J is nonzero, (7) implies that R_Z=0 if and only if all three partial derivatives of R vanish. Surjectivity of G therefore gives the set equality

    Σ=V(R,R_Z)=V(R,R_u,R_v,R_Z).                (8)

No scheme-theoretic equality of these ideals is asserted. For every a, differentiation of R_a(u,v)=R(u,v,a), together with (7), similarly gives

    Sing(V(R_a)) = Σ∩{Z=a},                    (9)

identifying the plane Z=a with the target plane. The source-line packet proves that each R_a is reduced irreducible.

## 4. Surjective projection, entirely away from the boundary

For d>1, the one-line injectivity theorem used in the source-line packet says F cannot be injective on any source line. Thus every a∈C has at least one ordered collision:

    π_a:T(C)→C is surjective,
    π_Z:Σ(C)→C is surjective.                  (10)

The second statement follows equivalently from the finite bridge (4). These are surjectivity statements on every complex parameter, not merely dominance or Zariski-dense image. Each fiber is finite because the source-line parametrization is a finite birational normalization and has only finitely many singular-image collisions. Finiteness of each fiber does not assert that either projection in (10) is proper or finite.

At H=0, the graph equation and its Z derivative specialize to

    R=AZ+B,       R_Z=A.

If H=R=R_Z=0, then H=A=B=0, contrary to the frozen finite-selector gate. Therefore

    Σ⊂D(H),                                   (11)

and all collisions occur off the nonproperness curve. The function H(P(a,s),Q(a,s)) is consequently a unit on the collision quotient as well, by the Nullstellensatz. Unlike (3), this last argument supplies existence of an inverse without an explicit coefficient formula.

The joint necessary condition is now fully algebraic: a noninvertible map in this sector requires a collision curve with two units, Δ and H∘F, whose parameter projection is surjective onto all of A¹. Its finite image is exactly the singular graph locus cut out set-theoretically by R=R_Z=0.

## 5. Sharp controls and the remaining obstruction

The checker `verify_collision_gap.py` tests (2) on exact polynomial automorphisms, where the collision ideal is the unit ideal; tests a non-Keller graph where the gap is not invertible and V(R,R_Z) is larger than Sing(W); and tests the following two abstract collision-like curves. Neither abstract curve is claimed to arise from a Keller map.

**Dominance is insufficient.** The ideal

    (a(s−t)−1, t)

makes Δ a unit and has zero elimination ideal in C[a], but its projection omits a=0. Dominance or a zero elimination ideal does not prove (10).

**Even surjectivity plus a unit gap is insufficient for properness.** The ideal

    (t(s−t)−1, a−t²+t)

has quotient C[t,t^(-1)], with a=t²−t and Δ=t^(-1). Every complex a occurs: for a≠0 the roots of t²−t−a are nonzero; for a=0 the root t=1 survives. Nevertheless t→0 gives a→0 and s=t+t^(-1)→∞, so the projection is not proper. This is compatible with finite fibers, a unit gap, and surjectivity.

The next bounded target is to show, for a prescribed polynomial family or source compactification, that the actual equations U=V=0 have an empty fiber at some a, or that the actual singular locus Σ misses some scalar Z=a. Either would contradict the one-line theorem in the noninvertible case. Equivalently, one could prove that the required collision curve cannot exist with its two units and its finite source-to-graph structure. No such general exclusion, analytic bound, or properness statement follows from the present identities alone.

The only external input beyond the frozen proofs is the same explicitly versioned one-line theorem recorded as Theorem 2.7 in [Z. Jelonek, A note on the Jacobian Conjecture, arXiv:2011.03472v3](https://arxiv.org/html/2011.03472v3), citing Gwoździewicz. Its use and source scope were independently checked in the source-line audit. The unit identity and finite bridge above have direct proofs here.
