# Preserving the elliptic pencil already forces preservation of its parameter

2026-10-11. **PROVED, scoped and self-audited, with CITED elliptic-curve
inputs for P1:** P1–P3. **FINITE-EXACT:** rational identity controls. **OPEN:**
changing to a genuinely different source pencil and planar JC(2).

## Inheritance and the enlargement

[M5, the fixed-invariant classification](mobius_affine_source_return_20261011.md)
already closes *all* birational symplectic transformations preserving

    I=x^2t^2+(1+x^2)t=u+p.

It is stronger than the earlier fibre-affine calculation. We do not
reopen that solved restriction. Here the new question is whether a
birational symplectic F can **permute** this pencil's fibres:

    I composed with F=h(I),        h a nonconstant rational function. (1)

The hostile is a genuine fibre permutation which changes the area form.
The least-used sidecar is the relative elliptic differential, whose
multiplier must pay for the derivative of h. Matching j-invariants alone
forgets precisely this information.

The closest taxicab connection is the genus-one cubic carrier, but
the [taxicab note](golden_taxicab_boundary_20261011.md) keeps the distinction
between its constant j=0 curve and this nonisotrivial pencil.

## P1. The symplectic pencil normalizer is the same Klein four group

Over any characteristic-zero field K, a K-birational transformation F
of K(x,t), with F*(dx wedge dt)=dx wedge dt, satisfying(1), must have
h(I)=I. Consequently its complete list is the inherited four maps

    (x,t),  (−1/x,x^2t),  (1/x,−1−x^2(1+t)),
    (−x,−1−t−x^−2).                                      (2)

**Proof.** Put k=I, X=−t and Y=xt(t+1). The generic fibre is

    E_k: Y^2=X(X−1)(X+k).

It is geometrically integral and smooth over K(k). Thus K(k) is
relatively algebraically closed in K(x,t). A field automorphism F sends
it to another relatively algebraically closed subfield, K(h(k)). Since
k is algebraic over K(h(k)), this forces K(k)=K(h(k)). Therefore h is
a fractional-linear map. In particular, a higher-degree rational h
cannot escape the normalizer calculation by wrapping the base.

The two generic elliptic curves E_k and E_h(k) must be isomorphic as
curves; transporting the image of the origin back to the origin makes
this a pointed isomorphism. Its j-invariants agree:

    j(k)=256(k^2+k+1)^3/[k^2(k+1)^2].                       (3)

The poles of j are exactly0,−1,infinity. Hence h permutes these points
and belongs to the six-element list

    k, −1−k, 1/k, −k/(k+1), −1/(k+1), −(k+1)/k.            (4)

All six preserve(3), also checked directly. This uses only the explicit
rational function, not a numerical special fibre.

After shifting X to remove its quadratic term, the short Weierstrass
coefficients are

    A(k)=−(k^2+k+1)/3,
    B(k)=(k−1)(2k+1)(k+2)/27.

An origin-preserving isomorphism has coordinate scale
Z' = v(k)^2 Z, Y'=v(k)^3Y, with

    A(h)=v^4 A(k),     B(h)=v^6 B(k).                      (5)

Its pullback on the relative differential dX/Y is multiplication by
1/v. Any translation used to move the origin preserves that relative
differential. Translation depending on k may add a multiple of dk
to an absolute differential, but that term vanishes after wedging with dk.

Crucially the original area form is

    dx wedge dt=−dk wedge dX/(2Y).

Thus symplecticity requires h'/v=1, or v=h'. Substitution into(5) gives
the necessary conditions

    A(h)=h'^4 A(k),       B(h)=h'^6 B(k).                  (6)

If h in(4) is nonaffine, it has a simple finite pole at0 or−1. The
left side A(h) has pole order2. The right side has pole order8, since
h' has pole order2 and A is nonzero at both0 and−1. This is impossible.
The remaining nonidentity h=−1−k has h'=−1; A(h)=A(k), but B(h)=−B(k)
is not h'^6B(k)=B(k). Therefore only h=k survives.

The complete fixed-invariant classification M5 then gives(2). QED.

For the elementary Weierstrass isomorphism criterion and proof, see
Andrew Sutherland's [Lecture13, Theorem13.13](https://math.mit.edu/classes/18.783/2022/LectureNotes13.pdf).
The four-point section-group input is cited and scoped in M5, with its
primary source. The rational identity and pole-order arguments above are
our deduction; they do not follow merely from equality of j-invariants.

**Hostile and boundary.** Over a field containing i, the polynomial map

    (x,t) -> (ix,−1−t)

satisfies I' =−1−I. Its Jacobian is −i, so it is a real example of the
base permutation rejected by the area-form test. Rescaling its output
to repair the Jacobian changes the claimed pencil action and has to be
checked afresh. Maps between different pencils, nonbirational maps,
base extensions with new source data, and arbitrary Keller maps remain
outside this theorem.

## P2. The first linear change of polynomial correction changes genus

The original source family is G_H=−u/2+H(p,y), where

    v=xt,   u=v^2/t,   p=t+v^2,   y=v(t+v^2).

Already the linear corrections H=ap+by have a complete generic-genus
trichotomy:

| Coefficients | Genus of the generic fibre G_H=E |
|---|---:|
| a=b=0 | 0 |
| b=0, a nonzero | 1 |
| b nonzero | 2 |

Constants added to H just shift E. This classifies this family of
generic fibres, not the full source problem.

**Proof.** Multiplying the fibre equation by t gives a quadratic in t:

    (a+bv)t^2 + [v^2(a+bv)−E]t −v^2/2=0.

For (a,b) nonzero its discriminant is

    Delta(v)=(f(v)−E)^2+2f(v),        f(v)=v^2(a+bv).       (7)

Here E is transcendental. The quadratic in f has distinct roots because
its discriminant is4(1−2E). Neither root can equal a constant critical
value of f. If b nonzero the critical values are0 and4a^3/(27b^2),
from f'=v(2a+3bv); hence Delta is squarefree of degree6. If b=0 and
a nonzero, the only critical value is0 and Delta is squarefree of
degree4. The smooth projective hyperelliptic model therefore has genus2
or1, respectively. The quadratic transformation is birational on a
nonempty open set; clearing t does not declare a spurious extra component
part of the generic fibre. If a=b=0, the equation gives t=−v^2/(2E),+a rational generic fibre. QED.

This points to a productive change of target: compute genuine
cross-source compatibility on the genus-two branch rather than continue
enlarging the already closed Legendre normalizer. A birational
transformation carrying whole pencils cannot transfer the elliptic branch
to the genus-two branch, since generic genus is preserved. An isolated
return, a singular special fibre, or a map not preserving the source
pencil has different quantifiers and remains unexcluded.

## P3. A positive cross-source transfer; it becomes polynomial at a=0

Write G_(a,b)=−u/2+ap+by. For (a,b) not both zero define

    F_(a,b)(x,t) = (2(a+bxt)/x, −x^2t/[2(a+bxt)]).         (8)

It is birational and symplectic, and satisfies the exact source equation

    G_(a,−b) composed with F_(a,b)=G_(a,b),
    F_(a,−b) composed with F_(a,b)=identity.               (9)

Most usefully, for a=0 and b nonzero it is the polynomial linear map

    F_(0,b)(x,t)=(2bt,−x/(2b)),             det J=1.         (10)

This transfers two genuinely different genus-two sources. It does not
preserve the old invariant I, so P1 is no obstruction to it.

**Mechanism and proof.** In the (v,t) chart the area form is
dv wedge dt/t. The quadratic in P2 has root product
−v^2/[2(a+bv)]. Its root exchange is therefore

    (v,t) -> (v,−v^2/[2(a+bv)t]).

It preserves the source and reverses the area form. The sign map
(v,t)->(−v,t) also reverses the area form, while exchanging b with−b.
Their composition is(8), which preserves the area form and gives(9).
This independently explains the direct rational-identity checks. At a=0
the denominators cancel to give(10). Conversely, if a is nonzero, the
first coordinate2a/x+2bt has a pole along x=0, so this particular map is
not polynomial. This is the exact cancellation boundary.

At b=0 the construction fixes one elliptic source and agrees with the
previous inversion-type symmetry, after its parameter is matched. The
new useful case is b nonzero, especially the polynomial a=0 subfamily.
Transporting a Keller complement, if one is supplied, through(10) is
legitimate because the Jacobian is1. The map does not create such a
complement, prove it exists, or settle the general Jacobian conjecture.

This is the constructive reroute after P1: change the actual source
pencil, retain its affine offset, and compose two sign-reversing operations.
The source coefficients and the domain of pole cancellation are essential
parts of the receipt. Counting the two operations alone would lose them.

## Reproduction

[Checker](../../04-computation/experiments/golden_taxicab_boundary_20261011.py)
and [output](golden_taxicab_boundary_20261011.json): six exact rational base
transformations are checked against(3) and(6); the four existing maps and
the non-symplectic hostile are checked directly; six coefficient pairs
give independent squarefreeness/degree controls for(7). The cross-source
identity, Jacobian, inverse and polynomial cancellation are symbolic
checks with a,b left indeterminate. The all-parameter
proof is the critical-value argument, not the six-case computation.
