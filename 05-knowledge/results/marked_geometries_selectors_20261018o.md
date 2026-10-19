# Marked selectors and concrete bridges at 5, 27, 21 and 9

**PROVED:** the elementary selector obstruction and conditional list gadget.
**FINITE-EXACT:** explicit incidence maps and group actions on the specified
small geometries. **CITED:** their classical cubic-surface interpretation.
**OPEN:** an authenticated unrestricted six-colour geometric selector for the
four Hamming charts, and six-colour forcing. No theorem of plane colouring
or Collatz convergence is asserted.

The useful new object is a marked chain of incidence maps, not an equality
of cardinalities. Deriving the ten-label Witt design at one point gives a
nine-point affine plane. Marking its origin gives four opposite pairs. Those
pairs produce the two Hamming extensions already needed by the ten-label
matching construction. The orientation of the pairs, and the particular
four-chart partition, are additional information.

## 1. Inheritance and incoming work

Closest proved mechanism: [the ten-label Hamming packet, §§2,6,7](ten_label_matching_hamming_20261018n.md).
Hostile: all single port words were already covered by matching projections;
extra codes improve joint relations, not that union. Corrected near miss:
Hadamard gauge normalization was necessary to retain the linear origin
(MISTAKE-642). Least-used sidecar: an oriented representative of each
opposite pair in a derived affine geometry.

The board shares the ROOT anchor with [the actual clock certificate](root_actual_clock_certificate_20261018o.md):
source weights, integer clocks, joint boundary relations, finite geometries,
orientation gauges, and source-fixing retractions. The geometry is the niche;
the recovered number atlas is the wildcard.

Incoming 190ea97bf4 and 269c5e124b sharpen the problem. The
[boundary-cell atlas](boundary_fibre_atlas_20261018d.md) constructs each of
the four charts separately but proves that gluing all four imposes their
intersection, which is empty. The [cubic attachment](surviving_repair_geometry_20261018d.md)
has a source-fixing retraction, so all old colourings extend even after
arbitrarily enlarging its named base sets. The
[variable-label compiler](variable_pair_fibres_20261018d.md) retains the
actual odd cycles, pin-path constraints and all original source colours.
These results supersede another search within the same retracting layers.

## 2. A 27-point object whose local boundary is a ten-label matching

Over F2 use

    q(x)=x0*x1+x2*x3+x4+x4*x5+x5,
    B(x,y)=q(x+y)+q(x)+q(y).

Its nonzero zeros form a set P of size 27. Its lines are triples
{x,y,x+y} with B(x,y)=0; there are 45. This is the classical elliptic
quadric model of GQ(2,4). Its relation to the 27 cubic-surface lines and
five tritangent planes through a line is described in
[Brouwer--Van Maldeghem, *Strongly Regular Graphs*, §10.10](https://homepages.cwi.nl/~aeb/math/srg/rk3/srgw.pdf).
The [Veldkamp-space paper, §3](https://arxiv.org/html/0903.0715) also records
the elliptic-quadric model and its 36 doily sections. These are classical
identifications, not newly discovered classification counts.

For each **marked** p in P, its ten other collinear points split as

    {x,x+p}, five pairs.

The checker verifies on all 27 choices that the induced graph on these ten
points is exactly five disjoint edges. Thus the map to the ten-label matching
carrier preserves an actual binary relation. It loses the remaining sixteen
points and all incidences through them; retain p and the original six-bit
labels to reconstruct the ambient context. Changing p changes the actual
ten-point boundary, so it is not an operation on one fixed physical port set.

For each h with q(h)=1, the polar section P intersect h-perp has 15 points
and 15 lines. The 36 sections are distinct. A saved incidence bijection sends
one marked section to the six-label doily: points are duads of six labels,
lines are perfect matchings. This reconnects the old Richelot/Pauli carrier
through actual incidence rather than merely the number 15.

An independent model uses the 27 classes

    E_i; H-E_i-E_j; 2H-sum_(j!=i)E_j,

with H^2=1, E_i^2=-1 and other basis products zero. Their intersection-one
triples reproduce all 45 lines, and a saved labelled isomorphism is replayed.
This also corrects a historical error: for distinct line classes L,M,

    (L-M)^2=-2-2(L.M).

Thus a root difference corresponds to **skew**, not meeting, lines. The old
`categorification_23.py` mixed the blowup-count and degree conventions for
del Pezzo surfaces as well. Both displayed claims are repaired, with their
lineage in MISTAKE-643. The 27/45/36 incidence counts survive.

### The 27 -> 9 projection has a concrete central coordinate

Brouwer--Van Maldeghem also give an affine F3^3 model in §10.10. In the
complement graph, distinct (v,z),(w,t), with v,w in F3^2, are adjacent when
v=w or t-z=area(v,w), where area(v,w)=v1*w2-v2*w1. This model is a
Heisenberg group with multiplication

    (v,z)*(w,t)=(v+w,z+t+area(v,w)).                       (1)

Our checker authenticates an incidence map from the binary quadric to these
literal 27 addresses, all associativity cases, and the following projection.
Forgetting z gives the same nine-point affine plane used below. Nine vertical
three-point lines collapse to its points; the other 36 lines project onto
its twelve lines, each with three central lifts. At a point the five pairs
of neighbours are one vertical opposite pair plus four horizontal direction
pairs. Thus 27=3*9 and 5=1+4 describe this one object, rather than two
unrelated counts.

The missing coordinate is an accumulated oriented area. Its commutator is
(0,2*area(v,w)). A matrix A in GL(2,3) acts by
(v,z)->(Av,det(A)z); this is a group and incidence automorphism, checked
on every pair for all 48 matrices. The same determinant will select the
Hamming orientation in §4. Here the restored coordinate is modulo three;
it is not the integer arrival clock of ROOT.

This is also a precise reason to retain the contact predicate. The lift has
vertical edges, so it is not a loopless matching thickening of the quotient
graph; the incoming odd-cover support theorem does not apply. And two
nonadjacent vertices have five common neighbours. Two distinct planar points
have at most two common unit neighbours, by intersecting their unit circles.
Consequently this entire 27-vertex graph cannot be an injective planar
unit-distance graph. It can organize constraints, but its abstract edges
cannot simply be installed as physical contacts.

## 3. The actual 21 -> 9 -> four-pair map

The old [HYP-2572, Mathieu residue](../hypotheses/HYP-2572-unit-distance-n22-mathieu-residue.md)
already built PG(2,4), with 21 points and five points on each line. Its use
as a general unit-distance classifier remains a proposal, not a realization
theorem. We reuse the finite plane itself.

Normalize nonzero triples over F4 by their first nonzero coordinate. On this
21-point plane, take

    U={ [x:y:z] : x^3+y^3+z^3=0 }.

Every nonzero F4 element has cube one, so U consists exactly of points with
two nonzero coordinates. It has 3*3=9 points. The 21 lines meet U in one
point (nine tangents) or three points (twelve secants). A saved isomorphism
identifies these twelve triples with the affine lines of F3^2.

Independently, derive the inherited S(3,4,10) Witt design at literal point 9.
Its twelve remaining three-element blocks are the same affine plane, through
a second saved coordinate map. Mark the preimage of (0,0). For the four
directions choose representatives

    (1,0), (0,1), (1,1), (1,2).

The eight nonzero points are these representatives and their negatives.
The four opposite pairs are **exactly** the pairs obtained by deriving the
Witt design at the two marked points. Adjoining that marked pair gives the
ten-label matching from the preceding session.

This is an incidence-preserving bridge among the old Mathieu residue, the
current ten-label design, and the matching surgery. The affine coordinate
map, origin, ordered directions and signs are retained. No full automorphism
group of the derived affine plane is asserted to extend to the original
ten-label design or to physical port isometries.

Two nine-point geometries must not be confused. This affine plane has twelve
lines and every pair of points is collinear; its collinearity graph is K9.
A nine-point grid GQ(2,1) has six lines and collinearity degree four. Equal
point counts do not identify their predicates. Nor can making every unital
collinearity a planar unit distance work: four pairwise equidistant planar
points are already impossible.

## 4. Orientation chooses the Hamming extension; an atlas needs more marking

Order the eight points by the four displayed representatives followed by
their negatives. Let e_i be the mask of opposite pair i. The space E of even
unions of these four pairs has dimension three. Adjoin one transversal t
choosing one point from each pair. The parity of the within-pair choices
gives two possible codes E+<t>, each H8. In the retained coordinate order,
one is span(170,204,240,255). Adjoin the original marked pair to obtain
the exact ten-coordinate H8+I2 used previously.

**FINITE-EXACT action statement:** all 48 matrices in GL(2,3) act on these
eight literal ports. Determinant 1 preserves each Hamming extension;
determinant -1 exchanges them. The parity of within-pair sign changes is
the determinant character. Each matrix and induced permutation is saved.
This supplies a genuine orientation/determinant connection, not an analogy
between a field of nine elements and nine arbitrary labels.

Preserving the code still does not preserve the particular four-chart atlas.
On eight coordinates the inherited atlas consists of four cosets of
span(255,240). Under SL(2,3), which has 24 elements, this atlas has an orbit
of **four different atlases**, with stabilizer of size six. The induced
permutation group on those four atlases is A4 of order twelve. All four
partitions and all twelve actions are exported and replayed.

There are therefore two distinct fours: four charts in one partition, and
four partitions in this local group orbit. They must not be collapsed.
An action can preserve all 16 Hamming words yet move the subspace defining
which eight-bit words share a chart. A geometric selector needs the original
ports and the chosen partition, not merely a code automorphism.

## 5. A stronger selector obstruction and an authenticated local gate

The four charts can be written using selector bits alpha,beta and free
bits p,q,k:

    (p,p+alpha,p+beta,p+alpha+beta,
     q,q+alpha,q+beta,q+alpha+beta,k,k).                 (7)

All sums here are in F2. Fixing alpha,beta gives one eight-word chart;
varying all five bits gives precisely the 32-word ten-port code C.
This is a reusable truth table with literal ordered ports.

**Selector obstruction (PROVED).** No graph in which every vertex has a
fixed colour list of size at most two can realize C as its projected binary
boundary relation, even with arbitrary hidden vertices and pins.

Choose a Boolean coordinate for each two-element list. An edge forbids
some of the four assignments at its ends and is a conjunction of 2-CNF
clauses. Pins are unary clauses. Every such constraint is preserved by
coordinatewise majority: if a forbidden pair appeared in the majority,
two of the three assignments would supply each of its two entries, and
those two index sets intersect. One original assignment would violate the
constraint. Projection preserves majority by taking the majority of three
full extensions. But 0,170,204 belong to C while their majority, 136, does
not: it is an unsupported weight-two word. Contradiction.

This extends the fixed-ternary-label obstruction to arbitrary fixed
two-element lists, including lists from different palette pairs. It does
not forbid adaptive label choices: a selector must genuinely use a larger
domain somewhere. The incoming odd-cover majority theorem preserves a fixed
label fibre for a different reason; neither majority result creates a selector.

There is a small positive primitive. Give three input ports Boolean colours
0,1. Join each to a different vertex of an internal triangle, whose lists
are {0,1,2}. Its boundary relation is NAE (not all equal). Equal input bits
leave only two internal colours and cannot colour the triangle; mixed bits
permit all three. Exhaustive extension counts authenticate both directions.

This graph has a literal planar unit drawing. Internal coordinates are
(0,0),(1,0),(1/2,sqrt(3)/2); the three ports are
(-1,0),(2,0),(1/2,1+sqrt(3)/2). All fifteen pair distances are checked
exactly in Q(sqrt(3)); there are exactly the three triangle edges and three
spokes. **The lists are an additional assumption.** With all six colours
available, every Boolean input triple extends, including the two rejected
by NAE. Thus this is not yet an unrestricted six-colour gadget.

It is nevertheless a composition-ready logical primitive. Direct truth
tables prove

    NAE4(a,b,c,d) iff exists w: NAE(a,b,w) and NAE(1-w,c,d).

Taking d=0 implements a three-literal OR. Four such clauses exclude the
four odd-parity assignments of three bits, producing an even-parity gate.
Together with complemented ports and a fixed Boolean zero, those gates
express (7). This compilation is abstract list-colouring logic; it does
not assert that the separate drawings can be glued at arbitrary prescribed
locations or that the lists can be forced geometrically.

The incoming odd-degree one-neighbour theorem explains a concrete obstacle
to forcing those lists. A point exterior to the old conjugation-stable field
in an odd-degree extension has at most one old unit neighbour. Direct old
palette contacts can therefore exclude at most one named colour there,
whereas the triangle primitive needs a three-colour list. Indirect coupled
exterior constraints remain possible, but the current six cubic layers have
no mutual contacts and retract. That mechanism, rather than a failed SAT
chart alone, dictates the next geometric search.

## 6. A typed number atlas recovered from earlier sessions

| Number | Recovered scenarios | Mechanism and missing coordinate |
|---|---|---|
| 5 | golden discriminant; anchor -5 in (9n+5)/8; five presentations of the strong four-vertex tournament; five tritangent pairs at a cubic line | The [5/8/9 note](five_eight_nine_transfer_20261005.md) keeps field reduction, legal valuation word, and presentation fibres distinct. The quadric supplies an actual five-pair boundary; it does not identify those five objects with tournament presentations. |
| 27 | three odd multipliers in the [charged terminal family](collatz_charged_terminal_family_20261012a.md); the universal pair-chain relation X-1=27(Y-1) in [THM-4601](../../01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md); 27 cubic lines; the three-channel bound B_n^2<=27^n in [the Fourier packet](cube_unitary_clock_packets_20261018n.md) | Ordered affine composition, a marked rational anchor, a finite incidence configuration, and a symmetric-power norm are different operations. The Fourier cubic and Eisenstein source share a coefficient field, but no authenticated source-preparation map is supplied. |
| 21 | PG(2,4) points; Fano flags and Paley-seven arcs in [the earlier cardinality bridge](cardinality_bridges_20261018m.md); source ell-mass coefficient 41/21 | The first is 4^2+4+1; the second has an explicit order-21 equivariant map with a base flag; the third is a geometric shell sum. Their labels and actions differ. The current unital construction gives a new actual inclusion into the first. |
| 9 | F9 from the golden discriminant modulo 3; affine-plane points above; nine product Richelot kernels after a symplectic split; nine perfect matchings of the cube; two-step multiplier 3^2 | F9's additive plane is F3^2 after a chosen basis, but multiplication by phi has order eight and norm -1. Product kernels need a chosen splitting, and cube flows need factorial weights. These structures cannot be substituted for one another by a nine-element bijection. |

The older [HYP-2892 unital frame](../hypotheses/HYP-2892-lrc14-unital-pair-frame-hodge-refinement.md)
and [HYP-2894 realizability hostile](../hypotheses/HYP-2894-lrc14-unital-pair-slot-realizability-guardrail.md)
are particularly relevant: 28 unital points did not identify a natural
S8-invariant design on K8 pair slots. This session keeps explicit incidence
maps and does not replace the marked arithmetic source by a uniform design.

The [historical number-type reflection](../../07-reflections/numbers-as-type-patterns.md)
supplied provenance for the 27-line tangent. Its verbal assignments such as
“lossy cubed” have no defined invariant that proves a transfer. The corrected
[Fano/snark synthesis](../../07-reflections/fano-primes-snarks.md) already
separates Eisenstein splitting, cubic lines and four-flows. We retain that
correction rather than deriving a weaker or false version anew.

## 7. Reproduction and next useful objects

```sh
python3 04-computation/experiments/marked_geometries_selectors_20261018o.py --out /tmp/geometry.json
PYTHONHASHSEED=31 python3 -O 04-computation/experiments/marked_geometries_selectors_20261018o.py --out /tmp/geometry-optimized.json
cmp /tmp/geometry.json /tmp/geometry-optimized.json
```

The [receipt](marked_geometries_selectors_20261018o.json) retains all point
labels, incidence maps, 48 local matrices, four chart partitions, extension
tables and exact coordinates. There are 91,850 active checks. Incidence maps
are checked edge by edge after construction; the Picard and binary-quadric
models are independent descriptions. Stable integer graph IDs remove an
arbitrary hash-order choice of isomorphism from the receipt.

Next obligations are precise: force and transport the three-colour list of
the authenticated NAE primitive, or construct a different unrestricted
six-colour selector; preserve the chosen chart partition as well as the
Hamming code; and defeat the existing source-fixing repair on an actual
edge before spending another search on the same cubic layers. The ROOT
counterpart is to bound the cost of the actual weighted readout, not just
the norm of an ambient carrier.

Research moves used from [META-PATTERNS](../../00-navigation/META-PATTERNS.md):
“The same representation is not the same carrier” and “Compose witness
fibres before taking their scalar shadow.” The code-preserving A4 atlas
motion and the non-unit-realizable Heisenberg graph provide explicit
counterindications to dropping the partition or metric sidecar.
