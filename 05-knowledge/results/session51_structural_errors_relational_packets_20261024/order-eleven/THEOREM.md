# Order-eleven refinement error: isolate the regular sector

24 October 2026. **PROVED:** the exhaustive nonregular strong stable-cell
classification, exact finite-quotient formulas, SCC transfer, and prospective
trace congruence. **FINITE-EXACT:** the small quotient and crossing calculations,
with standalone literal-presentation and unrestricted automorphism checks.
**INHERITED CERTIFIED INPUT:** completeness of the fifteen regular-nine classes
and the lower-order strong error distributions. **NOT COMPUTED BY THIS PACKET:** the regular-eleven
error, the order-eleven upper trace, and the unrestricted order-eleven class count.

All uses of "unknown" or "remaining" below refer to the input boundary of this
packet, not to a claim that another calculation has not subsequently supplied
the coefficient. The residual theorem remains valid unchanged when that input
is supplied.

## 1. Statement, definitions, and scope

For an ordinary tournament T, directed color refinement starts with one color
and repeatedly retains the current color together with the outneighbor counts
into each current color. Let C(O) be the eventual stable cell containing an
actual automorphism vertex orbit O. Define

    e(T) = sum_O ( |C(O)| - |O| ).

Write E_s,n for the sum of e(T) over ordinary strong n-vertex tournament classes,
E_n for its sum over all ordinary n-vertex classes, and E_regular11 for its sum
over regular eleven-vertex classes. Converse tournaments are not identified.
Regular tournaments are strong, so the last collection is a strong sector.

**Main result.** The positive-error nonregular strong order-eleven classes have
exactly the following actual stable-cell patterns and contributions:

| Actual stable cells | Positive-error classes | Error sum |
| --- | ---: | ---: |
| 9+1+1 | 12, with errors 18, 36, or 72 | 630 |
| 7+3+1 | 2, each error 14 | 28 |
| 7+1+1+1+1 | 22, each error 14 | 308 |
| 5+5+1 | 88 rigid classes, each error 40 | 3520 |
| 3+3+3+1+1 | 260 rigid classes, each error 18 | 4680 |

The patterns are pairwise disjoint because these are actual stable cells, not
merely a supplied equitable partition. Consequently

    E_s,11 = E_regular11 + 9166,
    E_11   = E_regular11 + 19480.                         (1)

The full positive-defect class polynomials are

    strong order 11:
      R_11(t) + 24 t^14 + 263 t^18 + 2 t^36 + 88 t^40 + 7 t^72;

    all order 11:
      R_11(t) + 48 t^14 + 348 t^18 + 8 t^36 + 256 t^40 + 28 t^72,  (2)

where R_11(t) counts positive-error regular-eleven classes by error. The known
parts contain respectively 384 and 688 classes. No numerical value for R_11,
E_regular11, U_11, or a_11 is claimed.

If U_11 is the pointed stable-cell upper trace, then

    a_11 = ( U_11 - E_regular11 - 19480 ) / 11.

If P_10 is the number of almost-regular ten-vertex classes and r_11 the number
of regular eleven-vertex classes, the remaining correction is

    E_regular11 = 11 (P_10 - r_11).                       (3)

Thus any independently obtained complete upper trace must satisfy

    U_11 = 10 mod 11.                                    (4)

This is a prospective necessary check, not a predicted count or a substitute
for the unknown regular-eleven calculation.

Inheritance: the closest mechanism is the Session50 equitable-cell and signed
triangle holonomy theorem; the canonical hostile is the nontransitive regular
seven; the corrected near miss is counting cyclic-majority triangle blocks as
actual cells when refinement merges them; and the useful sidecar is the complete
crossing action, not just crossing degrees. Session51's two-C5 degree-two sector
is the new edge-rigidity predecessor. The conceptual board is actual-cell
recognition, internal vertex orbits, signed holonomy, weighted quotient
symmetry, and SCC transfer. Each sector below records how these interact.

## 2. Restrictions and the complete odd-partition list

Every stable cell induces a regular tournament. If its size is m, the constant
internal outdegree is (m-1)/2 by counting its internal edges, so m is odd.
For distinct cells of sizes a,b, constant crossing row and column sums r,s obey
ar=bs. If gcd(a,b)=1, the crossing is uniform. In particular every singleton
attaches uniformly to each other cell.

The twelve odd partitions of eleven are

    11;
    9+1+1;
    7+3+1; 7+1+1+1+1;
    5+5+1; 5+3+3; 5+3+1+1+1; 5+1+1+1+1+1+1;
    3+3+3+1+1; 3+3+1+1+1+1+1; 3+1+1+1+1+1+1+1+1;
    1+1+1+1+1+1+1+1+1+1+1.

These are possible size patterns, not a claim each is realized. Pattern 11 is
exactly the regular-eleven sector, and is deliberately retained as unknown.
All other patterns not listed in the main table have zero error:

- A regular five-cell is the unique cyclic tournament C5. Outside 5+5+1 its
  other cells have coprime size, so its internal transitive C5 action extends
  independently.
- A three-cell is a directed cycle. Nonuniform crossings between two such
  cells are a permutation matrix or its complement. In oriented F3 coordinates
  the matching has form x -> epsilon x+c, with epsilon in {+1,-1}. A rotation
  vector k preserves it exactly when k_j=epsilon k_i. A forest of these edges
  permits a transitive C3 action on every triangle component. With at most two
  triangles there is no cycle, so there is no error.
- In the only three-triangle partition, positive error requires an unbalanced
  signed support cycle. This forces all three triangle crossings nonuniform.
  Uniform or incomplete-support cases are zero-error and need not be counted.

Every regular tournament is strong: in a tournament with multiple ordered
strong components, all degrees in an earlier component exceed all degrees in
a later component. Hence equal degrees rule out such a decomposition.
This also ensures that uniform substitutions are strong exactly when their
quotients are strong.

## 3. Uniform nine- and seven-cell substitutions

### 3.1 Pattern 9+1+1: contribution 630

All crossings are uniform, and the quotient must be a directed triangle.
The two singleton names are interchangeable, so there is one rooted quotient
class. The three cell degrees are 5, 1, and 9, so they really are the stable
cells. Every internal automorphism of the regular-nine cell extends, and every
full automorphism restricts to one. The error is exactly that of the internal
regular-nine tournament.

The byte-pinned Session50 regular-nine certificate proves fifteen classes, with
error histogram

    error 0: 3; error 18: 3; error 36: 2; error 72: 7.

Their total error is 630. The present checker independently reconstructs one
nine-vertex representative from each inherited component's almost-regular
order-eight parent, verifies all fifteen full automorphism groups and orbit
multisets, and checks pairwise nonisomorphism. Completeness is inherited from
the authenticated predecessor atlas; these new checks do not manufacture
completeness from a closed list.

### 3.2 Pattern 7+3+1: contribution 28

Again crossings are uniform. The three distinct cell sizes intrinsically name
the quotient vertices, so both directed-cycle orientations survive as different
ordinary classes. In one orientation the cell degrees (on cells 7,3,1) are
(6,2,7); in the other they are (4,8,3). Thus the stated stable cells are actual.
The full group is Aut(R_7) times C3, and the error is e(R_7).

The independent regular-seven enumeration gives three classes with full group
orders 3,7,21. Only the order-3 class is defective: its vertex orbits are 3,3,1
inside a seven-cell, giving error (7-3)+(7-3)+(7-1)=14. The two weighted cycles
therefore contribute 2*14=28.

### 3.3 Pattern 7+1+1+1+1: contribution 308

Substitute a regular-seven tournament into the distinguished vertex of an
ordinary five-vertex quotient. If the quotient is strong, the root's quotient
outdegree is 1,2, or3, so seven-cell degrees are 4,5, or6. A singleton degree
is at most3 if it loses to the seven-cell and at least7 if it beats it.
Thus refinement separates the seven-cell at once.

The remaining four singleton vertices either become four actual singleton
cells or contain an actual three-cell plus one singleton. A two- or four-cell
is impossible by parity. A merged triple is a directed cycle, uniformly
attached to the seven-cell and to the remaining singleton. Collapsing it gives
one of the two strong 7+3+1 weighted cycles just counted.

There are 544 labelled strong five-vertex quotients. This small count follows
without a target atlas from the labelled source-SCC recurrence, using strong
counts 1,0,2,24 through order four:

    2^10 - 5*2^6 - binom(5,3)*2*2 - 5*24 = 544.

Exactly 16 collapse: choose the singleton outside the triple in four ways,
choose the triple's directed-cycle orientation in two ways, and choose one of
the two weighted cycles. The other 528 are actual 7+1+1+1+1 presentations.
The root-preserving S4 action on them is free, because a quotient automorphism
would otherwise permute actual singleton colors. There are 528/24=22 quotient
classes. Every substitution has full group exactly Aut(R_7), so the error sum
is 22*14=308. The checker independently verifies all 1024 labelled quotients,
the two collapsed orbits of size8, and all 22 actual orbits of size24.

## 4. Pattern 5+5+1: contribution 3520

Give both five-cells the cyclic coordinates of C5. Its full automorphism group
is translations: the unique source in a vertex's two ordered outneighbors is
its successor, so every automorphism preserves successor.

Let M be the binary crossing matrix A -> B. It has constant row and column
sum d. Swapping the names A,B replaces d by5-d, so a nonuniform crossing has
a unique normalization d=1 or2. This normalization uses crossing degree, not
total vertex degree after decorating by the singleton. Write t_A,t_B in {0,1}
for whether the respective cells beat the singleton. Then the three degrees are

    on A: 2+d+t_A; on B: 7-d+t_B; on the singleton: 10-5(t_A+t_B).

The core is strong. Strongness after adding the singleton is exactly the mixed
condition t_A != t_B.

For d=1, both mixed attachments yield actual 5+5+1 cells: the degrees are
(3,7,5) or(4,6,5). These orientations cannot be identified because the crossing
normalization intrinsically names A and B.

For d=2, the attachment (0,1) has degrees (4,6,5), hence actual 5+5+1 cells.
But attachment (1,0) has degree5 at every vertex: its actual stable partition
is 11, not 5+5+1. ALL those presentations belong to the unknown regular-eleven
sector. Counting their C5-by-C5 presentation orbits as extra ordinary classes
would be invalid, since hidden decompositions can become invisible there.

The translation stabilizer of a nonuniform matrix is either trivial or a
noncoordinate line in F5^2. A coordinate-axis element would force constant
rows or columns and then uniformity. A nontrivial line is transitive on each
cell, giving error0; a trivial stabilizer makes the full tournament rigid,
giving 10*(5-1)=40. Actual-cell recognition makes this the full group, not
just a known subgroup.

There are 120 degree-one matrices. There are 2040 degree-two matrices: their
bipartite 2-regular graph has component row-size partition 5 or3+2, so

    R(5,2) = (5!)^2 [ 1/10 + (1/6)(1/4) ] = 2040.

Every full-slope nonidentity translation fixes binom(5,d) matrices; each
nontrivially stabilized matrix is fixed by four such translations. There are
16 nonidentity full-slope elements. Consequently the rigid-class numbers are

    d=1: (120 - 4*5)/25 = 4;
    d=2: (2040 - 4*10)/25 = 80.

There are two surviving decorations of d=1 and one of d=2. Therefore there
are 2*4+80=88 rigid classes, with error sum 88*40=3520. The accompanying
nonuniform zero-error classes number 2*4+8=16. The uniform crossing sector
has just the weighted directed cycle, with independent C5^2 action and zero
error. The checker treats it as a separate positive-control example.

## 5. Pattern 3+3+3+1+1: contribution 4680

Positive error requires all three triangle crossings to be nonuniform. Their
row sums are1 or2. Orient the majority quotient Q on triangle cells by the
row-sum-two direction. For each singleton h=0,1 let t_hi record whether
triangle i beats it. Also retain the edge between the two singletons.
There are

    8 majority tournaments * 2^6 attachments * 2 singleton edges = 1024

named skeletons. They retain every count needed by refinement and strongness.
Signs and phases of the three matching permutations are a separate sidecar.

### 5.1 Actual cells and the complete collision rule

The singleton degrees are 3 times an integer plus its outdegree within the
two-singleton edge. Their residues modulo3 differ, so the two singletons have
different degree colors and can never coalesce.

A stable cell containing exactly one singleton and one or three whole triangle
blocks has even size and is impossible. If it contains one singleton and two
triangle blocks, it has size7. Its crossing to the remaining size3 triangle
would have to be uniform by coprimality, but each of those two included
triangle blocks has a nonuniform crossing to the remaining triangle.
That is impossible. Thus both singletons remain actual singleton cells.

Two triangle blocks cannot coalesce alone, because their combined size6 is
even. All three coalesce precisely when both singleton attachments are constant
across the triangles and the majority quotient is cyclic. Indeed constant
attachments cannot separate them, and a cyclic majority gives equal core
degree4; a transitive majority gives distinct core degrees3,4,5. Conversely
coalescence forces equal singleton relations and equal core degree.

Therefore the only possible stable patterns in this presentation family are
actual 3+3+3+1+1 and the controlled collapse 9+1+1. The checker validates this
rule over all 1024 skeletons both on literal eleven-vertex tournaments and on
an independent compressed five-block refinement matrix. The signs and phases
cannot affect this calculation, since all within-block and between-block
counts are constant and fixed by the skeleton.

### 5.2 Strongness and the number 98

The nine-vertex core is strong because every triangle is strong and every pair
of triangle cells has crossing arcs both ways. Relative to this core, a singleton
has one of three attachment types: source (one constant string), sink (the other
constant string), or mixed (six strings).

For any fixed Q, the strong skeletons consist of:

- both singletons mixed: 6*6*2=72 choices;
- exactly one mixed: two positions, six mixed strings, two uniform types, and
  one necessary singleton-edge orientation, giving24 choices;
- one source and one sink: two positions, and the singleton edge directed from
  the sink to the source, giving2 choices.

Thus there are98 strong skeletons for each Q, or784 total. The only strong
collapses have cyclic Q, one uniformly attached source, one uniformly attached
sink, and the closing singleton edge. There are 2*2=4 such named skeletons.
They are already members of the 9+1+1 sector and must not be counted again.
The number with actual 3+3+3+1+1 cells is therefore780.

Equivalently the finite breakdown is

| Majority | Strong actual | Strong collapsed | Nonstrong actual | Nonstrong collapsed |
| --- | ---: | ---: | ---: | ---: |
| Transitive, six labelings | 588 | 0 | 180 | 0 |
| Cyclic, two labelings | 192 | 4 | 48 | 12 |

### 5.3 Exact class quotient and full automorphisms

The block-name group is S3 on the triangles times S2 on the singletons. On
actual skeletons this group acts freely: any skeleton automorphism would
preserve refinement and hence every actual cell color individually. There
are780/12=65 skeleton classes, split49 transitive and16 cyclic.

In oriented triangle coordinates each matching is x -> epsilon x+c over F3.
There are eight sign assignments and27 phase assignments. For a fixed skeleton,
internal rotations act on phases by a linear map whose kernel obeys
k_j=epsilon_ij k_i. Four negative-product sign assignments have zero kernel,
so each has exactly one phase orbit and no nontrivial automorphism. The other
four have kernel dimension1, and each has three phase orbits. These twelve
balanced classes are transitive on every triangle cell and have error0.

Because actual cells are canonical, no further block permutation enters the
full group of a fixed representative. Thus each of65 skeletons gives exactly
four rigid classes and twelve zero-error classes. The rigid classes number
260, and each has nine fixed triangle vertices, contributing9*(3-1)=18.
Their total contribution is260*18=4680. The transitive/cyclic split of rigid
classes is196+64.

As an independent check, for each of the65 skeleton representatives the code
constructs all216 sign/phase presentations, partitions them under all27
internal rotations, and tests the full eleven-vertex automorphism group by
unrestricted degree-pruned arc backtracking. It recovers all1040 resulting
classes:260 rigid and780 with group order3. It does not merely assume that
the displayed rotations comprise the entire automorphism group.

## 6. Completion and the SCC distribution transfer

The twelve odd partitions are exhaustive. Coprimality, internal transitivity,
and the signed-forest mechanism dispose of every pattern outside the five
nonregular exceptional sectors and the single unknown regular-eleven sector.
The sector total is

    630 + 28 + 308 + 3520 + 4680 = 9166.

Within sectors all allowed presentations and isomorphisms have been accounted
for. Across sectors, actual stable-cell sizes prohibit overlap. The three
important collision exclusions are:

1. Sixteen strong rooted-five quotient presentations collapse from proposed
   7+1+1+1+1 to actual7+3+1.
2. Every d=2 two-C5 presentation with attachment(1,0) is regular eleven and
   is withheld from the nonregular residual.
3. Four strong three-triangle skeletons collapse to actual9+1+1.

For the passage from strong to arbitrary tournaments, distinct strong components
occupy disjoint degree intervals. Stable cells and automorphism orbits therefore
remain inside individual components and e is additive. Ordered strong-component
decomposition gives

    E(z)=A(z)^2 E_s(z).

The inherited lower coefficients are E_s,7=14, E_s,8=0, E_s,9=716,
E_s,10=3992, and E_s,k=0 for k<=6. From a_0,...,a_4=(1,1,1,2,4),

    [z^0,...,z^4] A(z)^2 = (1,2,3,6,13).

Hence the lower-component padding at order eleven is

    2*3992 + 3*716 + 13*14 = 10314.

Adding9166 proves the total constant19480 in(1). In particular the multiplier
of E_s,9 is3, the degree-two coefficient, not6, the degree-three coefficient.

The smallest defective strong component has order7, so at order11 there can
be at most one defective strong component. The same A^2 factor therefore
applies separately to every positive defect value. The inherited lower strong
polynomials are

    order7: t^14;
    order9: t^14 + 7 t^18 + 2 t^36 + 7 t^72;
    order10: 4 t^14 + 32 t^18 + 84 t^40.

This proves the two polynomials in(2), not merely their first moments.

## 7. The isolated regular-eleven problem and sharper support boundary

Deleting a vertex from a regular-eleven tournament gives an almost-regular
ten-vertex tournament with five vertices of degree4 and five of degree5.
Conversely any such parent has a unique regular extension: degree4 vertices
beat the new root, which beats degree5 vertices. An isomorphism of parents
extends uniquely over the distinguished roots. Almost-regular parent classes
therefore correspond exactly to pointed regular-eleven states.

A regular-eleven class with q vertex orbits contributes q pointed states and
has one stable cell of size11, so its error is11(q-1). Summing proves(3).
This leaves a precisely specified lower-order reroot problem, but this packet
does not enumerate its parents or obtain its trace.

Every tournament automorphism group has odd order: a nontrivial involution
would swap some pair and reverse their directed edge, and an even-order group
would contain an involution. Therefore every vertex orbit has odd size. Since
eleven is odd, q is odd. Positive regular-eleven errors can consequently lie
only in

    {22,44,66,88,110}.

This is an allowed support, not an existence claim for each value. These values
are disjoint from the known defect values14,18,36,40,72 in(2). Thus the listed
688 total-order-eleven classes at those five known error values are exact even
without resolving the regular sector. The unknown polynomial can be written

    R_11(t)=r_1 t^22+r_2 t^44+r_3 t^66+r_4 t^88+r_5 t^110,

with unknown nonnegative integers r_i. In particular E_regular11 is divisible
by22, and E_11 is congruent to19480, hence10, modulo22. The upper-trace
consequence remains U_11=10 mod11 because its exact trace is11*a_11.

## 8. Certificate contract and limits

The source object is the actual stable-cell decomposition, internal oriented
regular cells, and literal crossing matrices or matchings. The target is the
class-weighted refinement-error distribution. Permitted internal actions and
weighted quotient isomorphisms preserve ordinary isomorphism, strongness,
vertex orbit sizes, and e. Crossing degrees alone lose matrix stabilizers;
unsigned triangle skeletons lose matching signs and phases; a chosen equitable
partition loses whether refinement really recognizes its cells. These are the
three required sidecars, each tested independently here.

The checker enumerates small quotients and exceptional families, never all
2^55 labelled order-eleven tournaments. It verifies full automorphism actions
on1232 nonregular strong eleven-vertex representatives, plus the inherited
fifteen regular nines and the three newly regenerated regular sevens. Some
zero-error families are disposed of analytically rather than enumerated.

No unrestricted order-eleven atlas or target count is an input, output, or
calibration. No scalar recurrence for all future errors, complexity improvement,
unmarked metagraph sufficiency, or literature priority is claimed. The genuine
advance is that all nonregular order-eleven error is closed, leaving only an
explicit regular-sector problem and a prospective trace congruence.
