# The general fixed-phase response is controlled by degenerate boundary shapes

**Status: PROVED + FINITE-EXACT;
[independent analytic and exact audit passed](long_frontier2_sep06_boundary_audit.md).** This reduces the general B-only and two-interlacer
sign questions to explicit boundary families. It does not prove the signs
on those families or close the general Laurent model.

## 1. Inheritance and the precise objects

The [coefficient-prism theorem](long_frontier2_sep06_rectangle.md) gives
the complete original response and its degree-eight elimination H.
The earlier [ordinary-carrier audit](long_frontier_sep06_anchor_audit.md)
already identified an indefinite coefficient Hessian as an obstruction
to a convex positivity shortcut. Here that same curvature becomes a
boundary reduction. The [degree-eight decoder](continuing4_20260906_moments_packet.md)
supplies the exact weak two-interlacer geometry, including common roots.

The canonical hostile is the positive off-original-root response at the
valid centre. The corrected near miss is replacing the full geometry
by lower moments or assuming the whole coefficient domain is convex.
The least-used sidecar is the original affine phase plane and an explicit
direction of strict positive curvature within that plane.

Use the same polynomials B,C,D and complete P,Q as the prism note:

    B(v)=v^5-13v^4+55v^3-av^2+bv-f,
    C(v)=v^4-12v^3+45v^2-(2a/3)v+3b/7,
    D(v)=v^4-11v^3+36v^2-(5a/12)v+b/7,
    P(-s)=182-20020s+2002as^2-3432bs^3+2002fs^4.

For each fixed real s>0, let K_B(s) consist of coefficients (a,b,f)
for which B has five nonnegative real roots, with multiplicities, and
P(-s)=0. Let K_CD(s) be the subset where C and D both weakly interlace B.
The two anchors e1=13,e2=55 are built into B. Thus a,b,f are nonnegative
without a separate positivity relaxation.

Define the boundary sets

    E_B(s)={points of K_B(s): B has a zero or a repeated root},
    E_CD(s)={points of K_CD(s): B has a zero or a repeated root,
             or gcd(B,C) has positive degree,
             or gcd(B,D) has positive degree}.             (1)

All these conditions are on the original polynomials. A canceled common
root is retained, and weak boundaries are not silently discarded.

## 2. Exact boundary extremum theorem

**Theorem.** For either (K,E)=(K_B,E_B) or (K_CD,E_CD), if K(s) is
nonempty then E(s) is nonempty and

    max_(K(s)) Q(-s)=max_(E(s)) Q(-s),
    min_(K(s)) Q(-s)=min_(E(s)) Q(-s).                       (2)

Every maximizer and every minimizer belongs to E(s). Consequently,
at each fixed phase, strict negativity, nonpositivity, or the existence
of a nonnegative/positive response can be tested on the named boundary.
The same equivalences hold for all s>0 by applying the result separately
at each phase; no interchange of limits and extrema is used.

There is a constructive strengthening for the maximum. Every point of
K(s) outside E(s) can be moved along an explicitly given coefficient line
to a point of E(s) with a strictly larger response, preserving the
original phase and the entire relevant root/interlacer predicate.

## 3. Compactness and the genuinely open two-dimensional chart

The possible ordered B roots lie in the compact set

    0<=beta1<=...<=beta5<=13,
    sum beta_i=13, sum beta_i^2=59.

Their elementary-symmetric coefficient image is compact. Intersecting
with P(-s)=0 is closed, so K_B(s) is compact. Weak interlacing is a
closed condition on the ordered real roots of these fixed monic degrees,
so K_CD(s) is compact too. The named boundary subsets are closed.

If B has five distinct positive roots, sufficiently small perturbations
of a,b,f keep that property: five disjoint positive sign brackets persist
and exhaust its degree. At such a point a,b,f are strictly positive.
If C,D also weakly interlace and neither shares a B root, both interlacings
are strict. Their four disjoint interlacing root brackets persist under
small perturbations as well.

Since s>0, the original equation solves for f with nonzero coefficient:

    f=12b/(7s)-a/s^2+10/s^3-1/(11s^4).                      (3)

Therefore every point of K(s) outside E(s) has an open two-dimensional
neighborhood in the free coordinates a,b on this plane. No full-domain
convexity or unstated tangent-rank assumption is needed.

## 4. Indefinite curvature and an explicit boundary move

The complete carried identity at(3) is

    Qbar(a,b;s):=s Q(-s)=-(14/11)H(a,b,s),                   (4)

where H is the19-term polynomial in the prism theorem. For fixed s it
is quadratic in a,b. Its Hessian entries are

    h11=26558675s^6+(593856780/7)s^5,
    h12=-(845791650/49)s^7-(3563140680/49)s^6,
    h22=(286833690/49)s^8+(1972792800/49)s^7.

Its determinant is exactly

    -31194786150*s^12
      *(10966105s^2+72692884s+144097056)/2401 <0.            (5)

Thus both H and Qbar have indefinite Hessians at every s>0. An indefinite
quadratic has no local minimum or maximum on an open plane: along an
eigenvector with either curvature sign, the average of its two nearby
values lies strictly on that side of its value at the midpoint.
Compactness supplies extrema, and the open-chart argument excludes every
point outside E(s). This proves(2), including nonemptiness of E(s).

An explicit line gives the stronger maximum statement. Put

    da=-h12, db=h11,
    df=12db/(7s)-da/s^2.                                   (6)

The phase equation stays identically zero along (a,b,f)+t(da,db,df).
Moreover h11>0 and

    (da,db) Hess(H) (da,db)^T=h11 det Hess(H)<0,

so Qbar is strictly convex along this line. Start at a point outside E(s)
and take the connected component containing t=0 of the line's intersection
with K(s). It is a nontrivial closed bounded interval, with zero in its
interior. This uses only compactness and local openness; other components
of the line intersection may exist. At least one interval endpoint has
strictly larger Qbar by strict convexity. An endpoint outside E(s) would
have a further open neighborhood in K(s), contradicting maximality of
the component. The endpoint therefore lies in E(s). Since s>0 is fixed,
the same strict increase holds for Q(-s).

## 5. Consequence, controls and the remaining problem

The map is now a native coefficient operation: retain the original phase,
move in direction(6), and stop on the actual geometry boundary. It preserves
the two anchors and all required root/interlacer constraints along the
chosen interval, while increasing the response. It deliberately loses
the starting shape, but the boundary root multiplicity/common factor
records the endpoint obstruction. No representative is selected merely
because it shares a discriminant value.

For the B-only model, the next sign problem is therefore zero-root and
repeated-root configurations. For the full model, shared C or D roots
are additional branches. The common-root example in the exact degree-eight
decoder shows why such points cannot be discarded. The repeated/zero
control B(v)=v^2(v-3)(v-5)^2 also belongs to the B-only boundary.

The [standalone source](../../04-computation/long_frontier2_sep06_boundary.py)
reconstructs P,Q from ordinary carriers, eliminates f, and independently
checks every Hessian entry, determinant, lifted phase direction and strict
curvature factor. Four rational phase controls and the repeated-root
polynomial are exact checks. Compactness and the endpoint proof are
analytic; they are not inferred from those controls.

    python3 04-computation/long_frontier2_sep06_boundary.py
    python3 -O 04-computation/long_frontier2_sep06_boundary.py

All20 always-active gates pass in normal and optimized Python with raw
output identical to the frozen transcript. The independent audit also
reconstructs the full eliminated polynomial through the original Hadamard
convolutions, separately from this source's ordinary-carrier path.

    source SHA256 badbd9ba520c93b39c4386fdcc9360510181ce903938d950198a5bb74ab4af1e
    output SHA256 a15d43905603b0131d5f810ef56769c5ca017cee2828d82422fcf70d469a33a4

The all-boundary sign inequality is still OPEN. This theorem supplies
neither a finite boundary census nor a proof that every boundary branch
is safe, and it does not make the B-only and two-interlacer domains equal.
