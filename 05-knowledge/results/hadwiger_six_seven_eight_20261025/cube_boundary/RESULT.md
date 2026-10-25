# Eight boundary assignments, seven Walsh modes, and a genuine six-state unit interface

2026-10-25. **PROVED in the stated scopes; [independent audit PASS](../audits/BOUNDARY_AUDIT.md).**
Elementary statements below have exhaustive FINITE-EXACT controls. The
prescribed-colour examples are exact induced planar unit graphs. They do
not supply an unprescribed six-colour obstruction or a seven-chromatic
unit graph. No historical-priority claim is made.

## 1. Inheritance and the changed question

The closest mechanism is the complete boundary relation in
[spindle_tournament_interfaces_20261011g, sections 3--4](../../spindle_tournament_interfaces_20261011g.md):
an induced triangle with three pendant ports has a three-colour
NOT-ALL-EQUAL relation, invisible to existential pair projections.
The [fixed-ternary parity compiler](../../plane_six_reframes_parity_20261011b.md)
retains exact binary boundary equations and witnesses. The
[Session61 palette-deck packet](../../session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md)
distinguishes three transported bits from their parity and exhibits
pairwise blindness among Fano kernels. We refine a different issue:
**which information is lost when the actual three-bit boundary table is
compressed, and what composition can recover it?**

The K3 carrier is the eight central branch conics, indexed by F2^3, with
intersection matrix

    M=-2I+2A(Q3).

The [root-deletion proof, D3--D4](../../six_seven_geometry_20261025/root_deletion/DELETION_AND_CUBE.md)
and [Gale/cube continuation](../../six_seven_gale_cube_20261025/BRANCH_CUBE_AND_DOUBLE_COVER.md)
give its Walsh eigenvalues `4,0,-4,-8` in weights `0,1,2,3`. The three
weight-one zero modes are actual divisor relations. The
[regular-action geometry packet](../../ports_regular_actions_hn_20261025/geometry/README.md)
identifies the same three modes with the centered coordinate space of
injective unit cube drawings. It does not retain arbitrary colouring
boundary tables.

The canonical hostile is forgetting whether a marginal means **existence**
or **number of extensions**. The corrected near miss is that pairwise
counting marginals are much stronger than pairwise feasibility shadows.
The least-used sidecar is the integer weight of every boundary assignment.
Board: cube adjacency / exact joint support / counting marginal / Fourier
mode / pinned unit interface / unprescribed colouring quantifier.

**Concurrent recovery:** live `origin/main` received `9807ae1678` while
this lane was running. Its
[Session62 palette packet](../../session62_agent_theme_bridges_20261025/palettes/RESULT.md)
develops weighted path/port response ranks and recovers the inherited
three-colour NAE example. Complete weighted boundary tables and joint
composition are shared mechanisms, not novelty claims here. The present
additions are the exact 35-relation cube classification, the unique Boolean
counting-marginal ambiguity, and the explicit prescribed six-palette pair of
unit graphs isolated in the three weight-two K3 modes.

## 2. All pairwise feasibility tests admit exactly 35 different joint relations

Let `B={0,1}^3`. A relation `R subset B` has full pair shadows if every
assignment to any two coordinates extends to an element of R.

**Theorem B1.** R has full pair shadows if and only if its complement is
an independent set of the cube graph Q3. Exactly 35 such R exist:

| Number of allowed assignments | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|
| Relations with all pair shadows full | 2 | 8 | 16 | 8 | 1 |

**Proof.** Fixing two coordinates leaves precisely one cube edge as the
fiber of the projection. The projection is full exactly when no such edge
has both endpoints outside R. This is the independent-set criterion.

For the count, the independent sets have sizes zero through four. There
are one empty set, eight singletons, and `binom(8,2)-12=16` nonadjacent
pairs. A pair from opposite parity classes must be antipodal if it is
nonadjacent; every other cube vertex is adjacent to one of that antipodal
pair. Thus an independent set of size at least three lies in one parity
class. There are eight triples and two quadruples, with no larger set.
Complementing gives the table. QED.

Consequently **six, seven and eight allowed states can have identical
complete pairwise feasibility data**. This statement has a precise map:
project R onto its three named coordinate pairs. Six states do not mean
six colours. A distinguished six-state example is

    R_NAE=B minus {000,111}.

Its induced cube graph is a six-cycle. The excluded antipodal pair is a
sidecar: of the sixteen six-state relations above, only four remove an
antipodal pair; the other twelve remove a distance-two pair.

## 3. Counting marginals leave just one cubic degree of freedom

For an integer or real boundary table `w:B->R`, define unnormalized Walsh
coefficients

    W_S = sum_(x in B) w(x) chi_S(x),
    chi_S(x)=(-1)^(sum_(j in S) x_j).

Fourier inversion is `w(x)=1/8 sum_S W_S chi_S(x)`.

**Theorem B2.** All numerical marginals on proper coordinate subsets
determine w up to adding a scalar multiple of `chi_{123}`. Equivalently,
they determine the mass, three weight-one coefficients, and three
weight-two coefficients; only the cubic coefficient remains.

**Proof.** Summing over an omitted coordinate kills exactly the characters
using that coordinate. Hence all proper marginals determine every character
except full support. Conversely that full-support character sums to zero
in each proper marginal. Fourier inversion proves both directions. QED.

For nonnegative integer tables, write one possible table as `w_0`. All
tables with the same proper marginals are exactly

    w_lambda(x)=w_0(x)+lambda chi_123(x),
    -min_(chi=+1) w_0 <= lambda <= min_(chi=-1) w_0,
    lambda integral.

This is an explicit one-dimensional interval, not a vague information-loss
claim. For **Boolean indicator tables**, there is a sharper rigidity:

**Corollary B3.** Two distinct relations have identical numerical proper
marginals if and only if they are the even and odd parity tetrahedra.

Indeed their indicator difference is `lambda chi_123`; any nonzero
coefficient compatible with values zero and one is `lambda=+1` or `-1`.
It forces exactly the two parity indicator functions. These two four-state
relations therefore share every counting marginal but are disjoint. A
six-state NAE relation and the eight-state full cube share existential pair
shadows, but not counting pair marginals. Keeping these statements separate
is necessary.

## 4. Exact gluing uses all seven nonconstant modes

Fix three named binary ports and a fixed palette and exterior prescription.
Let w_F(x) count proper extensions in a fragment F. For fragments with
disjoint interiors and exactly those shared ports, take their literal graph
union, with no additional edges joining the private interiors or introducing
any other undecomposed cross-constraint. Extension choices are independent
after x is fixed:

    w_(F union G)(x)=w_F(x)w_G(x).

Thus their Walsh data compose by xor convolution,

    W_(F union G)(S)=1/8 sum_T W_F(T) W_G(S xor T),
    total extensions=1/8 sum_S W_F(S)W_G(S).          (B4)

Forgetting a port is ordinary marginalization. The support is precisely
the set where w is positive; counting thus retains feasibility as well.
Extra edges between shared ports can be included as zero-one indicator
factors. Edges involving private interiors require enlarging the boundary
or retaining additional joint data; merely listing those edges does not
make the product formula valid.
The same formulas hold for Boolean feasibility tables if products are
taken before converting positive weights back to support.

No proper linear quotient of the eight-dimensional table space preserves
all counting compositions: the continuation pinning ports to any x is
the delta table at x and reads w(x). These eight legal boundary tests
separate every table. This is a minimal dimension statement for a linear
counting interface, not a minimal automaton-state or complexity theorem.

The K3 matrix M is not such a complete interface. It kills three singleton
Walsh modes, so the distinct tables `1+chi_j` and `1-chi_j` have the same
image under M and different pinned-port responses. Passing from eight
component labels to divisor classes must therefore not be used as a
colouring compression. Conversely the singleton modes used by planar
cube coordinate maps discard the weight-two modes needed by the next
explicit unit-distance relation.

## 5. Two exact prescribed six-palette unit graphs with matching coarse data

Use `alpha=(1+i sqrt(3))/2` and the six base points

    z0=0, z1=1, z2=alpha;
    p0=-alpha, p1=2, p2=2alpha.

Their induced unit graph is the triangle on z0,z1,z2 and the three spokes
`z_i p_i`. This is the inherited exact geometric hostile, now with an
explicit **six-colour prescribed interface**. Let the palette be
`{0,1,2,3,4,5}`. Attach fixed-colour leaves to each z_i forbidding
3,4,5, and to each p_i forbidding 2,3,4,5. The three port bits now specify
their actual colours 0 or 1. The inner triangle has exactly two extensions
when those bits are not all equal and none when all equal.

Add a separate vertex with allowed colours {0,1}, enforced by four
fixed-colour leaves. This gives a 32-vertex induced unit graph F with

    w_F(x)=4*1_NAE(x)=3-chi_12(x)-chi_13(x)-chi_23(x).

For comparison keep the same three port positions, restrict them to
{0,1} by leaves, and add a separate vertex restricted to {0,1,2} by
three leaves. This 19-vertex induced unit graph G has

    w_G(x)=3 for every x.

Every introduced leaf has one intended neighbour. Exact placements can
be chosen successively: on its parent's unit circle, collision with any
earlier point excludes finitely many choices, and unintended unit contact
with another earlier vertex excludes at most two. The checker goes further
and constructs all coordinates rationally in the alpha basis, using

    u(t)=((1-t^2)/(1+t+t^2), t(t+2)/(1+t+t^2)),

whose norm `a^2+ab+b^2` is one, and choosing the first nonnegative integer
t avoiding every forbidden coincidence/contact. It then recomputes every
pair distance independently of the intended edge list. The two graphs
have 31 and 15 edges respectively. Their prescribed leaves are not
asserted to be forced to those colours by the unprescribed graph.

Both tables have total count 24, zero singleton coefficients, and zero
cubic coefficient. Their completion parity is balanced (12 extensions
of either parity), and their existential pair shadows are all full.
Yet F excludes 000 and 111 while G permits both. The exact discrepancy
lies in the three pair modes, with coefficients -8 each for F and zero
for G. These are precisely the three weight-two eigendirections of M,
each of eigenvalue -4. The six-state boundary cycle is therefore an
actual colouring interface using spectral channels absent from the
coordinate drawings of the cube.

Pinning the three ports to 000 is the decisive continuation: F has zero
extensions and G has three. Neither underlying graph gives an ordinary
high-chromatic obstruction: F is three-colourable (a triangle with trees
attached and a separate star), and G is a forest. The prescribed external
colours are the load-bearing data. Removing that prescription restores
many more possible colourings.

## 6. Connection contract and exact scope

| Source -> target | Map / retained predicate | Lost data and sidecar |
|---|---|---|
| Eight conic labels -> binary port assignments | chosen coordinate labels in F2^3 | only the cube and Walsh basis transfer, not complex curves or plane colourability |
| Joint relation -> pair feasibility | existential projections | 35 joint relations share the full shadow; retain joint support |
| Counting table -> proper marginals | coordinate sums | one cubic coefficient; retain it or its permitted interval |
| Two fragments -> complete boundary count | product / Fourier convolution | requires exact graph union, disjoint interiors, matching names and fixed prescriptions; extra private constraints require a larger boundary |
| Six-state NAE -> actual six-palette interface | fixed-colour leaf lists, exact unit coordinates | exterior prescription cannot become a universal obstruction |
| K3 relation quotient -> colouring table | no complete transfer | M kills genuine future-observable singleton modes |

The old HN target is unchanged: a genuine six-colour obstruction must
defeat every admissible exterior colouring, not a prescribed one. The
complete boundary relations in the old compiler and THM-4614's
[scoped extension theorem](../../../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md)
remain applicable. This packet supplies an exact interface and a stronger
geometric loss witness, not an escape from their quantifiers or field gates.

Reproduce `python -B .../cube_boundary/check.py` and repeat with `-O`, from
the repository root. The checker exhausts all 256 cube relations; all
65,536 pairs under exact gluing; all 6,561 tables with entries in {0,1,2}
for marginal fibers; and all eight boundary assignments in both exact
unit graphs. It verifies complete distance graphs, leaf prescriptions,
extension counts, Fourier formulas, and mode-specific hostile controls.
