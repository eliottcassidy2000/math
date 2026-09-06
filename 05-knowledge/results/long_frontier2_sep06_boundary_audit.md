# Independent audit of fixed-phase boundary extrema and the explicit line move

**Status: INDEPENDENT ANALYTIC + EXACT SOURCE AUDIT PASS.** The complete
[written proof](long_frontier2_sep06_boundary.md),
[source](../../04-computation/long_frontier2_sep06_boundary.py), and
[frozen output](long_frontier2_sep06_boundary.out) pass without a mathematical
correction. This audit covers minima as well as maxima, nonemptiness of
the named boundary, and the explicit strictly improving line operation.
The sign on the boundary families remains open.

## 1. Compact domains and the local chart

The two fixed elementary-symmetric anchors give sum beta_i=13 and
sum beta_i^2=13^2-2*55=59. The ordered nonnegative root domain is therefore
a closed subset of [0,13]^5. Its coefficient image is compact. The fixed
original-phase equation is closed, so K_B(s) is compact for every fixed
real s>0. The weak C,D interlacing subset is closed as well: all root
degrees are fixed and monic, and the intervening quartic roots remain in
[0,13]. Ordered weak root inequalities survive limits. Thus K_CD(s) and
both named discriminant/resultant boundary sets are compact.

A B polynomial with five distinct positive roots has disjoint positive
sign brackets. They persist under sufficiently small independent changes
of a,b,f, and degree exhaustion keeps all five roots distinct and positive.
All three coefficients a,b,f are then strictly positive automatically.
Consequently no extra coefficient-positivity face or f=0 obstruction is
missing from the argument: those cannot occur at such an interior point.

For the full model, weak interlacing with simple positive B and no common
B/C or B/D root is strict interlacing. The four root intervals for each
monic quartic are disjoint, and their strict brackets persist. In particular
there is no unlisted repeated C- or D-root phenomenon in this interior:
one root in each of the four disjoint open gaps is already distinct.

For s>0 the coefficient of f in P(-s) is 2002s^4, which is nonzero. Solving
for f therefore identifies the original-phase plane with the full free
(a,b) coordinate plane. Intersecting the open coefficient neighborhood
with this graph gives an actual open two-dimensional neighborhood in
K(s). No hidden rank condition, quotient of the phase, or global convexity
assumption is needed. Points with zero or repeated B roots, or shared
interlacer roots in the full model, are retained in E(s).

## 2. Independent complete identity and curvature reconstruction

The producer uses ordinary carriers built from (1+u)^14 and B,C,D.
An independent calculation instead used the original Hadamard rows and
explicit coefficient convolutions, retaining every Q index -1 through 8.
After the exact original-zero substitution it recovered all 19 H monomials
in the already independently audited rectangle certificate. Differentiating
that independently recovered polynomial gives precisely the three printed
Hessian entries.

Direct multiplication of the entries independently gives

```text
det Hess(H)
 =-31194786150*s^12*(10966105*s^2+72692884*s+144097056)/2401,

det Hess(Qbar)
 =-1031232600*s^12*(10966105*s^2+72692884*s+144097056)/49.
```

The second expression is the first multiplied by (14/11)^2, agreeing
with the earlier anchor audit. Both determinants are strictly negative
for every s>0. Thus each quadratic has both curvature signs. At any point
in an open plane, the average of its two sufficiently close values along
a direction of positive curvature is strictly larger; along a direction
of negative curvature it is strictly smaller. No gradient hypothesis is
required. There can be neither a local maximum nor a local minimum there.

Every compact nonempty K(s) has both extrema. The open-chart argument
therefore places every maximizer and every minimizer in E(s). In particular
E(s) is nonempty whenever K(s) is. Since E(s) is compact, its extrema
exist too and equal the corresponding extrema over K(s). This proves
both directions of every printed extremum equality.

The stated sign equivalences follow from these extrema, including strict
negativity on the entire compact fibre. Applying the argument separately
for each s>0 is legitimate; it asserts neither a phase-uniform margin nor
an interchange of limits and extrema. Empty fibres contribute no exception.

## 3. The explicit coefficient line preserves the full predicates

Write the Hessian entries as A=h11, B=h12, C=h22. Here A>0, and
v=(-B,A) is a nonzero fixed direction at each fixed phase. Independent
matrix multiplication yields

```text
v^T Hess(H) v = A*(A*C-B^2) = A*det Hess(H) <0.
```

Hence Qbar=-(14/11)H is strictly convex along the line. The printed
`df=12*db/(7s)-da/s^2` preserves the original affine phase equation
identically. Both anchors remain built into B throughout the move.

At a starting point outside E(s), local openness supplies a line interval
around parameter zero. The inverse image of K(s) on this nonconstant line
is compact; the containing connected component is a nontrivial compact
interval with zero in its interior. Connected subsets of the real line
are intervals. Other disconnected line sections may exist and are not
silently joined. Every point of the chosen interval retains the complete
B-only or weak-two-interlacer predicate by construction.

Strict convexity makes the starting value strictly smaller than a convex
combination of the two endpoint values, so at least one endpoint gives a
strict increase. An endpoint outside E(s) would have an open coefficient
neighborhood and extend that same interval, a contradiction. The improved
endpoint lies in E(s). Since s>0, the increase for Qbar is also an increase
for Q(-s). The operation preserves the original phase and all geometric
constraints; it does not evaluate the sign of its endpoint.

The repeated/zero polynomial v^2(v-3)(v-5)^2 has the stated fixed anchors
and a=75,b=f=0. It is an admissible degenerate B shape; membership in a
particular fixed-phase fibre still includes its original phase equation.
The finite control does not supply the compactness or endpoint argument.

## 4. Replay and frozen scope

The source has 20 always-active exact gates. Independent normal and
optimized replays agree byte-for-byte with the saved output. The source
has no numerical root inference or imported research producer. The complete
carrier, Hessian, direction and sign-factor identities were additionally
reconstructed by the independent path described above.

```text
source SHA256 badbd9ba520c93b39c4386fdcc9360510181ce903938d950198a5bb74ab4af1e
output SHA256 a15d43905603b0131d5f810ef56769c5ca017cee2828d82422fcf70d469a33a4
semantic SHA256 adff08236e7e7012141c279dd3870a6db56a263819450191af7536493d99248d
```

Targeted incoming searches located the earlier indefinite-Hessian
obstruction, but not the present fixed-phase extremum reduction or its
explicit line operation. This records a new use of that inherited
calculation without claiming external priority. The result is safe to
promote in its stated model scope. It gives neither a finite boundary
census nor a sign closure, and it does not identify the B-only and full
two-interlacer domains or supply general actual Laurent transport.
