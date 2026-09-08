# Continuing synthesis: the missing unit order, forced cuts, and independent repair

**CURRENT CHECKPOINT, 2026-09-08. Status: PROVED scoped theorems,
FINITE-EXACT certificates, and independent audits.** A concrete global
quadratic resolves the previous session's unit-order-two existence question.
General LRC(14), planar JC, and extremal no-three-in-line remain **OPEN**.
No external priority claim is made.

The [previous synthesis](continuing10_20260908_synthesis.md) contains the
source-linear classification, clock7200 closure and two-direction repair
bound. The [manifest](continuing11_20260908_manifest.json) pins this session's
sources, outputs, certificates, independent verification and incoming review.
Research began at47bf7e622. The first proof checkpoint is88652c5f3;
incoming mathematical review reaches840f1e47c, followed by the status-only
planar closeoutd8bd800bf.

## 1. The portfolio and the successful change of assumptions

The anchor was LRC's connected-complement entry, the niche was the missing
distinguished unit order in global source geometry, and the wildcard was a
third direction in no-three-in-line repair. Six concepts stayed live:
complete profile domains; forced connections between zero-credit components;
original-cell retention; the actual coefficient ring; every labelled fibre
component; and the next native operation on a response.

The inherited mechanisms were the complete126-profile minimum-tree consumer,
unit-cell-capacity bipartite flow, exact integration on quadratic pencils,
and the canonical principal-part connection. The hostiles were jointly
realizable zero-credit trees, overlapping directional deletion costs,
source-linear order-two examples with a boundary critical point, and scalar
pole orders that lose independent component directions. The corrected near
miss was treating a classification of source-linear functions as a global
obstruction. The underused sidecars were zero-edge graph cuts, the originally
safe cell set, translated repeated roots in an exact pencil, and a regular
third component whose principal part is zero.

The niche overtook the anchor after a direct two-chart calculation produced
a globally smooth quadratic with the missing order. Restoring its complete
component tuple then supplied a stronger theorem about the whole torsion.

## 2. A genuine global quadratic realizes unit order two

On the fixed surface W2=(P1 x P1) minus{Z=x^2}, use the actual charts
x=1/r, t=-r^2-r^4b. Set u=x-1, w=1+u^2t and L=(u+1)w-4. The new object is

    F=u w L
     =x(x-1)^5 t^2+2(x-1)^3(x-2)t+(x-1)(x-4),
    G=((1-2u)w+2)/(3u^2w^2).

**PROVED.** F is regular and has no critical point anywhere on W2, and
J_(x,t)(F,G)=1. In the ORIGINAL source response
C_F=C[x,t]/D_F(C[x,t]), the distinguished theta=[1] has exact annihilator
(F^2). Thus source-t degree two realizes the order excluded by the complete
source-linear classification on every W_m. Degree two is minimal in that
fixed source filtration; no minimal total-degree claim is made.

The [primary proof](continuing11_20260908_quadratic_unit_two.md),
[independent source audit](continuing11_20260908_quadratic_unit_two_audit.md),
and [independent Weyl audit](continuing11_20260908_quadratic_weyl_audit.md)
pay all source and boundary points, rational constants, and every affine
fibre. On the added boundary F=-4 and F_r=4. On the omitted source line
u=0, F_x=-3. Elsewhere y=uw gives valid coordinates with
F=(1+1/u)y^2-4y, proving submersivity directly.

The zero fibre has three reduced disjoint source components u=0, w=0,
and L=0. Every other original-source fibre is irreducible, including -4,
by the odd-order zero of its quadratic discriminant. Globally the added
boundary is another component over -4, so the source qualifier matters.
The rational field is C(x,t)=C(F)(y), with
u=y^2/(F-y^2+4y); its Hamiltonian constants are exactly C(F).

The complete scalar principal parts of G in the SAME parameter g=F are

    (9/g^2+4/g, (32/3)/g^2+4/g, 0).

In particular g^2G is polynomial. Cancelling a remaining pole of gG by
H(F) creates a pole on the regular L component, proving exact order two
and excluding polynomial mates in every degree.

The inherited component theorem identifies the full torsion with two
principal-part arms. Modulo their common diagonal write
theta=A/g^2+B/g, where A=(9,32/3), B=(4,4). Their determinant is -20/3.
The canonical connection nabla differentiates these scalar principal parts,
and the two lawful operations

    g theta=A/g,       (g nabla+2)theta=B/g

recover both first levels. Differentiation generates every higher level.
For the Weyl algebra D=C<g,nabla>/(nabla g-g nabla-1), the exact result is

    tors_(C[g])(C_F)=D theta ~= D/(D g^2).

The proof covers arbitrary finite operators: modulo the LEFT ideal Dg^2,
write P=a(nabla)+b(nabla)g. Independence of A and B first forces a=0,
then b=0. Each j-th canonical derivative of the unit has order j+2.

The regular third component is essential. Discarding its zero principal
part and quotienting only the two polar components would incorrectly erase
B as a diagonal vector. Scalar order two alone also does not prove cyclicity:
A/g^2 has that order but generates just one arm. The ring remains C[x,t];
D_F does not preserve O(W2), since D_F(b)=F_r/r^2 has a genuine boundary pole.

## 3. A complete pencil family has a wall invisible to unit order

The [family theorem](continuing11_20260908_quadratic_family.md) and its
[independent analytic audit](continuing11_20260908_quadratic_family_audit.md)
extend the example to an exact classification within the fixed whole
discriminant pencil (x-h)^5 span{1,x-h}. Its globally submersive members
are precisely

    u=x-h, w=1+u^2t,
    F=s+(a u^2+d u)w^2+k u w,
    k=-4ah,      a*d*h*(d-4ah)!=0.

Every member has a rational mate, three reduced disjoint source components
over s, no other reducible source fibre, and exact distinguished-unit
annihilator (g^2), g=F-s. The complete torsion always has two arms.
Nevertheless the unit generates both arms **if and only if d differs from
6ah**. On d=6ah it generates exactly one arm, with exact Weyl left annihilator

    D g^2 + D(g nabla+2-(6a/k^2)g).

Off that wall the left annihilator is Dg^2. The proof retains both scalar
principal coefficients, whose two component directions have determinant
-ak(2d+3k)/(6d^2). The wall lies inside the global-submersion locus. For
a=h=1,d=6, the source and boundary normal derivatives are2 and24; the unit
still has order two and the extra annihilator is g nabla+2-(3/8)g.
This is a real geometric counterexample to inferring cyclicity from scalar
order and the number of available arms.

The separate allowed values d=-2ah and d=2ah drop one leading component
pole and raise gcd(N,P) from u^3 to u^4, respectively. Neither is the
cyclicity wall. Recording these boundaries distinguishes three operations
that would look alike after retaining only a coarse valuation.

The reusable mechanism is exact in the inherited labelled torsion model:
if eta=sum_(k=1)^N A_k/g^k, then

    D eta=span_C{A_1,...,A_N} tensor g^-1 C[g^-1].

Euler interpolation separates the eigenvalues -k of g nabla, then
multiplication and differentiation generate every level in each coefficient
direction. If all N coefficient vectors are independent, the exact left
annihilator is Dg^N by triangular PBW reduction. Thus coefficient rank,
rather than highest pole order, determines the generated arms. Also
D(nabla^j eta)=D eta for every j>=0: differentiation shifts the pole
indices and rescales each vector by a nonzero scalar. Pole order can grow
while the generated module stays unchanged.

## 4. Three composite LRC clocks close by forced graph connections

**PROVED + FINITE-EXACT.** Exactly the declared clocks10080,10800,11520
are removed from the inherited connected-complement necessary array. The
[producer](continuing11_20260908_lrc_composite_clocks.md) and
[independent unpruned C++ audit](continuing11_20260908_lrc_composite_clocks_audit.md)
cover8,937,423 full multisets and212,562 valid word/clock evaluations.

| Clock | Valid words | Least guaranteed margin |
|---:|---:|---:|
| 10080 | 62,905 | 156 |
| 10800 | 79,399 | 48 |
| 11520 | 70,258 | 96 |

Every worst-case owner has five zero-credit tree edges and one forced
positive connection. The complete zero-edge graph has two components;
any actual connected complement must cross that cut. The cheapest bridges
cost156,156,180 at the actual quotient clocks2520,2700,2880. At11520,
retaining interval placement improves the separate-component bound by16.
The global smallest edge cost is zero, so a scalar edge floor loses the
entire gain. The full positional graph restores it.

The physical theorem concerns primitive rows of thirteen distinct positive speeds with a selected
six-label A, gcd(A) equal to one of the three clocks, and the complementary
seven labels' ACTUAL strict-atlas graph connected. It imposes no height box,
physical unit speed, connectivity of A, decoder equality or unit gcd of B.
Possible-edge minima are used as a lower bound; they are not asserted to be
jointly realizable. All126 proper-subset profiles are retained. No partial
minimum-tree monotonicity is assumed.

The necessary array now has **7,622 clocks, maximum11,935**, canonical SHA256
`832714c179ab4f425a76172c41d7d8eecf05dc15adfdd3c5e08bc7eb5a13e48b`.
The declared three-clock experiment is complete. Disconnected-complement
entry and the remaining clocks need further arguments.

## 5. A third direction forces separately chargeable repair

Let tau2 be the minimum number of original cells deleted to clear both
slopes+1 and-1, and tau3 also impose slope2. Let U consist of cells whose
two ORIGINAL old diagonals each have occupancy at most two. Then

    J3=sum_(slope2 lines l) (|U intersect l|-2)_+,
    tau3 >= tau2+J3.

**PROVED.** Deleting U cannot help an originally overfull old line, while
the third-direction lines partition U. The two costs therefore require
disjoint original cells. This survives every adaptive repair that adds new
points: any final no-three set must discard at least tau3 original cells.

The [three-direction theorem](continuing11_20260908_no3_three_directions.md)
also exactly compiles the Boolean line-selection dual, with value beta3:

    beta3=max_Qold [|N(Qold)|-2|Qold|
                    +sum_l (|l minus N(Qold)|-2)_+],
    tau3 >= beta3 >= tau2+J3.

The complete n<=5 census finds a genuine n=5 triangle with
beta3=1 < deletion LP=3/2 < tau3=2. The third direction destroys the
integrality of the two-direction flow. This is a precise obstruction to
an attractive extension, with the compiled lower bound as its survivor.
Among two-regular square boards, the first isolated extra triple occurs
at n=6. Its six old diagonals are
all singletons, so the additional deletion cannot repair an old violation.

Uniformly for EVERY fixed simple two-regular skeleton and EVERY fixed row
order under a random column permutation, the mean repair bound improves to

    liminf E tau3/n >= gamma2+exp(-416/45)/4
                        =0.3286222170740135... .

Here gamma2=0.3285980553011038... is the inherited two-direction constant.
The numerical increment is modest. The new structure is the independently
chargeable region and a uniform isolated-event argument: only O(n^3) of
the asymptotic n^4/32 target triples fail companion-row eligibility; all
eight column assignments avoid a bounded-degree forbidden permutation
matrix. Exact seven-line geometry and Jensen's inequality give416/45.

The [independent exhaustive/geometric audit](continuing11_20260908_no3_three_audit.md)
and [independent analytic referee](continuing11_20260908_no3_three_directions_analytic_audit.md)
verify the fixed-row uniformity, constants and scope. Column swaps change
tau3 by at most four, giving

    P(tau3<=k) <= exp(-(E tau3-k)_+^2/[8(n-1)]).

This improves the exponential obstruction to sublinear original-cell
repair. It does not pay a union bound over n! row orders or settle the
extremal no-three-in-line number.

## 6. Incoming work changes the next frontier

The [incoming review](continuing11_20260908_incoming_review.md) verifies the
new fixed-DG quartic closure, moving-index all-m obstruction, prime-degree
leading filter and four cross-concept theorems. The incoming prime filter
actually reuses the previous session's reciprocal-polynomial integration
lemma. On W2 a quadratic with a polynomial mate must have constant leading
coefficient or a single-root fourth, sixth or eighth power. Our leading
coefficient x(x-1)^5 fails that necessary filter, giving an independent
polynomial-mate exclusion while leaving its rational mate intact.

The new continuation-state theorem, THM-4460 in
[its exact file](../../01-canon/theorems/THM-4460-adverse-budget-continuation-state-and-shear-composition-collapse.md),
has a concrete analogue INSIDE our response module. The actual classes
theta and -nabla(g theta)=A/g^2 have the same scalar annihilator, but the
same future operator g nabla+2 separates them. This is a paid failure of
the scalar observable, without claiming a map between the two problems.
The incoming cyclic capacity and spatial shear results similarly require
their native incidence and realization constraints; neither supplies
additional forest credit or an integral three-direction flow.

Full fixed-W2 quartic closure ends the previous last-location search. For an
arbitrary quartic first coordinate, both coordinates must be global on W2;
the square-prefix subfamily already excludes arbitrary polynomial mates.
These hypotheses do not exclude our rational examples. Once the all-m constant-boundary gate holds,
further inverse infinity residues in that test vanish automatically. The
next observation must involve finite places or a different invariant.
The remaining quintic leading strata are necessary filters, not entry.

## 7. Next decisive tests

1. In the quadratic lane, extend coefficient-rank classification to another
   exact pencil only after matching its full global section space. The
   fifth-order pencil and its d=6ah wall are now closed; revisiting that
   determinant cannot provide a new case. Test whether other pencils admit
   more than two torsion arms with a distinguished cyclic unit.
2. For another LRC clock, compute the complete minimum-tree residual first,
   then recover zero-edge cuts and native endpoint labels only where needed.
   A cheap scalar floor of zero is a prompt to inspect the graph, not a
   proof of failure. Disconnected entry requires its own paid consumer.
3. For more grid directions, retain original-cell incidence and test an
   actual integral triangle before invoking flow duality. Extend the safe
   deletion charge or the exact compiled dual before optimizing constants.
4. Incoming fixed-DG quartic polynomial-pair closure ends the old remaining
   placement search. Higher-degree entry must use the current planar board;
   rational positive examples and full polynomial mates have different scopes.

All new transfers keep an explicit source, target, preserved predicate and
lost coordinate. None turns a successful scoped certificate into general
conjecture closure.

The final manifest contains37 frozen artifacts and nine reproducible engines.
Normal and optimized raw LF outputs agree across438,633 always-active gates:
435,458 algebraic/combinatorial controls and3,175 historical snapshot checks.
Seven certificates regenerate unchanged. Of872 inherited pins,868 retain
their baseline bytes and four have explicitly audited incoming status-only
footer migrations. Independent final integration review corrected the
minimal isolated-event statement to its two-regular-board universe.
