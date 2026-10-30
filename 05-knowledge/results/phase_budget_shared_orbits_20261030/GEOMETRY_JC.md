# Tight palette obstructions and arithmetic determinant resonance

**PROVED elementary; FINITE-EXACT independent computations.** Neither result
resolves Hadwiger–Nelson or the planar Jacobian conjecture.

## 1. Exactly nine tight pair-palette covers

The inherited input is
[THM-4614, three-colourable contexts](../../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md).
Fix a proper three-colouring q of the exterior graph. Each exterior vertex
has at most two neighbours on the boundary, whose proper six-colouring is
fixed. A pair-palette option assigns disjoint two-colour palettes to the
three q-classes. There are15 perfect matchings of the six colours and six
orders of their pairs:90 options.

An exterior vertex blocks an option exactly when its two boundary-neighbour
colours are distinct and equal the pair assigned to its class. A blocker
type is(c,e), c one of three classes, e an unordered colour pair. There are
45 types, each blocking six options. Fewer than15 blocking vertices cannot
cover all90 options. If an option survives, each exterior vertex has at least
one available colour in its pair; disjoint palettes across adjacent classes
give a proper extension.

**Theorem.** Exactly nine sets of15 blocker types cover all90 options:

1. Three row covers: all15 pairs in one fixed exterior class.
2. Six star covers: for a fixed boundary colour alpha, all five pairs
   {alpha,beta} in each of the three exterior classes.

Any cover by15 blocking vertices uses distinct types and blocks every option
exactly once, since15*6=90. Write f_c(e) for its binary incidence. For every
ordered perfect matching (e1,e2,e3),

    f_1(e1)+f_2(e2)+f_3(e3)=1.                            (1)

Swapping the classes assigned to two disjoint edges gives
f_i(e)-f_j(e)=f_i(f)-f_j(f). The disjointness graph on pairs of six colours is
connected: overlapping pairs have a common disjoint pair. Hence f_i-f_j is
constant on all15 edges. If any such constant is1 or-1, one row is all1 and
one row all0; equation(1) forces the third row all0. These are the row covers.

Otherwise all three functions coincide with a binary edge function f such
that every perfect matching has f-sum1. Exchanging two edges inside a matching
gives, on any four distinct vertices,

    f_ab+f_cd=f_ac+f_bd=f_ad+f_bc.                         (2)

These relations imply f_ij=a_i+a_j. For completeness, f_ik-f_jk is independent
of k outside i,j by(2); the resulting differences satisfy the additive cocycle
law using a fourth vertex. Choose vertex1 as base. The numbers
f_1j-(a_j-a_1) are constant, fixing2a_1, and recover every f_ij. The matching
condition is sum a_i=1.

Binary f implies every difference a_i-a_j is0 or+/-1. Thus the a_i take at
most two values, differing by1. They cannot both occur at least twice:
the values2a,a+(a+1),2(a+1) would all need to lie in{0,1}. If they are all
equal, the matching sum is0 or3, not1. If one value occurs once, the possibilities
are a single1 among five0s (a star, sum a_i=1), or a single-1/2 among five1/2s
(the complementary star, sum a_i=2). Only the star satisfies(1). This proves
the six remaining cases and their completeness. QED.

`palette_audit.py` independently enumerates the90-option exact-cover problem
(124 search nodes, nine solutions) and all32768 binary graphs on six named
vertices (exactly six meet every perfect matching once).

## 2. A geometric refinement at the fifteen-blocker boundary

For a boundary colour alpha let n_alpha be the number of distinct boundary
points carrying it. Two distinct plane points have at most two common
unit-distance neighbours: intersections of their unit circles. Therefore
the number of distinct exterior vertices adjacent to a point of colour alpha
and a point of colour beta is at most

    2 n_alpha n_beta.                                    (3)

A tight star cover needs three distinct blockers for every pair{alpha,beta},
one in each exterior class. If both colours occur on singleton boundary
points, (3) permits at most two, so that star is impossible.

**Corollary.** Suppose six boundary points have six distinct colours, every
exterior vertex has at most two boundary neighbours, and the exterior has a
specified proper three-colouring. With at most15 blockers, the boundary
colouring extends unless one exterior colour class realizes all15 distinct
boundary-colour pairs. In particular, ruling out that row configuration is
enough even at equality15, improving the inherited strict threshold.

For eight boundary points using all six colours, at least four colours occur
once. A possible star centre must therefore be a repeated colour. This is a
necessary restriction, not an existence construction in Q(zeta55).

The blocking condition concerns one extension method. A hostile consists of
15 independent exterior vertices, one for each boundary-colour pair, all put
in class0. It blocks every pair-palette option, yet each vertex independently
has four legal colours, so the actual six-colouring extends. This abstract
graph is not asserted to have the required unit-distance realization.
Even a geometrically realized row obstruction would not alone force six
colours or shrink the complete boundary relation. Other extension methods
must still be excluded.

## 3. An all-parameter resonance criterion for projected Keller maps

The incoming Session113
[Jacobian symmetry analysis](../session113_symmetry_consumers_20261030/jacobian/README.md)
proves that Reynolds projection need not preserve a Keller determinant; its
quadratic two-shear example loses a balanced bracket channel. The following
explicit family determines exactly when that happens.

Work over a characteristic-zero field. Let r,s>=2, w>=1 be integers and
a,b nonzero scalars. Set

    P=x+a y^r,       Q=y+b(x+a y^r)^s.                    (4)

This is a composition of polynomial shears, so det J(P,Q)=1. Give x weight1
and y weight-w; project P to weight1 and Q to weight-w, obtaining P0,Q0.
Since -wr!=1, P0=x. A term indexed by j in Q has weight
s-(wr+1)j. Therefore

    Q0=y                                          if wr+1 does not divide s+w,
       y+b binom(s,j) a^j x^(s-j)y^(rj)           otherwise,
    j=(s+w)/(wr+1).                                    (5)

When the divisibility holds,1<=j<=s, and

    det J(P0,Q0)=1+b binom(s,j)a^j rj x^(s-j)y^(rj-1).    (6)

The correction is nonzero and nonconstant. Hence **the projected pair is
Keller if and only if wr+1 does not divide s+w**. This is an exact statement
about this family and this projection; the original maps are automorphisms
in all cases. The resonance also gives s-j=w(rj-1), so the retained pair has
the incoming normal form (x,y B(x^w y)). The incoming torus classification
already excludes nonlinear B from being Keller; formula(5) adds the precise
arithmetic test for when this projection produces such a B.

For w=1,r=2, the s=2 example has defect4abxy. At s=5, j=2 and the defect is
40a^2 b x^3 y^3. At s=3 there is no resonance and the projection is(x,y).
Thus projection can succeed or fail within the same elementary family; its
failure is an arithmetic compatibility of weights, not merely nonlinearity.

For a finite cyclic weight action of order M, assume M does not divide wr+1,
so P still projects to x. The exact criterion becomes:

    no j in{1,...,s} satisfies (wr+1)j=s+w modulo M.       (7)

If a solution exists, differentiation of its retained monomial contributes
a distinct nonzero monomial to Q0_y. Different j cannot cancel because their
exponent pairs differ. A retained j=0 term depends only on x and has zero
y derivative, so it is correctly omitted from the obstruction test. Dropping
the assumption on P would invalidate this simple criterion; both components
could then retain additional terms.

`jacobian_audit.py` checks140 weighted-torus cases and1897 finite-cyclic cases
by exact symbolic differentiation, including success and failure controls.
The proof is the expansion above, not the finite test range. This provides a
precise channel guard for symmetry reductions, without manufacturing equivariance
of an arbitrary Keller pair.
