# Independent audit of the order-ten refinement correction

Status: PROVED structural classification, FINITE-EXACT small crossing checks.
The proposed coefficient survives independent derivation. No ordinary order-ten
atlas was read or constructed. No full order-ten count is predicted here.

## Verdict

The exhaustive stable-cell cases give

- two five-cells: 84 rigid classes, each error 40, total 3360;
- one seven-cell and three singletons: four strong substitutions of the unique
  bad regular-seven class, each error 14, total 56;
- three triangle cells and a singleton: 32 rigid classes, each error 18,
  total 576.

Therefore E_s,10 = 3992. The inherited SCC identity and independently inherited
E_s,7 = 14, E_s,8 = 0, E_s,9 = 716 give

    E10 = 3992 + 2*716 + 6*14 = 5508.

The coefficients 2 and 6 are the coefficients of z and z^3 in A(z)^2, with
A(z) = 1+z+z^2+2z^3+... . They are not arbitrary padding multipliers.

## Complete partition argument

The ten odd partitions of ten are

    9+1; 7+3; 7+1+1+1; 5+5; 5+3+1+1; 5+1+1+1+1+1;
    3+3+3+1; 3+3+1+1+1+1; 3+1+1+1+1+1+1+1; 1^10.

These are actual stable cells, not just an equitable block presentation.
Coprime cell sizes have uniform crossings. Thus 9+1 and 7+3 cannot be strong.
A unique five-cell has only uniform crossings to the remaining cells, so its
C5 action extends independently. At most two triangle cells have forest
interaction support and no unbalanced rotation constraint. Therefore all
patterns other than the three listed in the verdict have zero strong error.

This argument does not claim every odd partition occurs, and does not turn
an arbitrary equitable presentation into the stable partition.

## Two five-cells

Use the unique regular-five tournament on each block. If each first-block
vertex sends r arcs to the second block, column sums are also r. The two
vertex degrees are 2+r and 7-r, always different. Thus the blocks are the
actual stable cells and cannot be interchanged by an automorphism.
For r=0 or 5 the tournament is not strong; r=1,2,3,4 are strong. Relabelling
the two blocks sends the matrix M to 1-M^T, pairing r with 5-r. Normalize to
r=1 or 2. This identification is not tournament converse.

All isomorphisms within a normalized sector are internal C5 rotations, so
matrices are quotiented by C5 x C5, not S5 x S5 and not independent reflections.
For a nonconstant crossing, a stabilizer with one zero projection is impossible:
a nonzero translation of five rows/columns would make those rows/columns equal,
forcing a constant column/row after using equitability. Hence a nontrivial
stabilizer is exactly one diagonal subgroup of order five. There are four
possible nonzero slopes. Two distinct slopes together generate the full
translation group and would force a constant matrix, so they never overlap.

For r=1 there are 120 permutation matrices. Exactly 4*5=20 are invariant
under a slope (one offset out of five). They give four size-five orbits;
the other 100 give four size-25 orbits. The latter four classes are rigid.

For r=2 there are 2040 matrices. One independent combinatorial derivation:
2-regular bipartite graphs have either one alternating ten-cycle, giving
5!4!/2=1440 choices, or a six-cycle and four-cycle, giving
C(5,2)^2 * 6 = 600 choices. Exactly 4*C(5,2)=40 matrices have a nontrivial
translation stabilizer. They give eight size-five orbits; the other 2000
give eighty size-25 orbits. Each rigid class has ten fixed vertices in
five-cells and error 10*(5-1)=40.

The independent program exhausts the 120 and 2040 matrix universes, uses
literal translations to form all 96 classes, and independently searches all
arc-preserving vertex maps compatible with stable colors for every class.
It verifies the r versus 5-r block-swap identity on every representative.

## Seven-cell and three singletons

Every crossing is uniform, giving a four-vertex quotient with the seven-cell
marked. The expansion is strong if and only if that quotient is strong.
There is one ordinary strong-four type with trivial automorphism group, hence
four classes after marking its vertex. Independently, all 64 labelled
four-tournaments modulo the six permutations fixing the mark produce twelve
rooted quotient classes, exactly four strong.

Within each strong expansion the seven-cell has degrees 3,4,5 or6, whereas
singleton degrees lie in 0..2 or7..9. Any merged pair of singletons would be an
even stable cell and is impossible; merging all three would require identical
adjacency to the seven-cell, making the quotient nonstrong. Thus its actual
partition is 7+1+1+1. The seven-cell action is exactly the internal action.
Only the inherited nontransitive regular-seven type contributes, with error14.
The program expands that concrete hostile in all four quotient classes and
checks the full automorphisms, not only the quotient action.

## Three triangles and a singleton

Positive error requires an unbalanced signed support cycle. On three triangle
cells this forces all three crossings nonuniform, each a matching or its
complement. Their union is already strong. The singleton must have both
incoming and outgoing block incidences; the six mixed attachment masks are
exactly the strong possibilities.

The specified four blocks are equitable, so the actual stable partition can
only coarsen them. A union of two triangles has even size six and cannot be
stable. Merging the singleton with one triangle gives even size four. Merging
it with two triangles gives seven, but that cell would have a uniform crossing
to the remaining triangle by coprimality, contradicting nonuniform crossings.
Merging all three triangles gives nine, but then equitability to the singleton
would require all three incidences equal, contrary to the mixed mask. Merging
all four blocks gives even size ten. Consequently these are the actual four
stable cells, including when two of the initial degrees coincide.

The majority tournament on the triangle cells is either transitive or cyclic.
Up to cell permutation there are six mixed masks for the transitive majority
and two for the cyclic majority: eight skeletons total. Each skeleton has no
nontrivial cell automorphism (transitive order is rigid; a mixed mask breaks
the cyclic majority's C3 action). For each skeleton, three matching maps
have eight sign choices. Four balanced signs give three phase orbits each;
four unbalanced signs give one phase orbit each. Thus each skeleton has
sixteen ordinary classes, four rigid. This gives 128 classes, 32 rigid.

The program independently forms the 10,368 presentations using all twelve
nonuniform crossings at each of three pairs and all six mixed root masks.
It quotients by all 162 cell permutations and internal rotations. Full
arc-preserving automorphism searches give 96 classes of group order three
and 32 rigid classes. Their majority split is 72+24 balanced and 24+8 rigid,
for transitive and cyclic majority respectively.

The cyclic-majority sector must not be dropped merely because it was regular
at order nine: a mixed singleton resolves its nine-cell into three cells.
Its eight rigid classes contribute 144 of the total576. The program includes
an explicit control whose root deletion has a single stable nine-cell but
whose strong order-ten extension has three triangle cells and error18.
Conversely a source/sink singleton must not be counted as strong; both are
explicit hostile controls. Deleting a support edge kills the signed cycle
and its error, also verified explicitly.

## Independence and reproduction

`audit.py` was written before reading the Session51 producer implementation.
It imports only Python standard-library modules. It uses no parent atlas,
order-ten graph list, canonicalization package, or producer routine.
The inherited proofs are pinned separately; their results are declared inputs.
Normal and optimized executions use explicit exceptions, not removable asserts.

    python3 audit.py > normal.out
    python3 -O audit.py --output optimized.json > optimized.out
    cmp normal.out optimized.out
    cmp receipt.json optimized.json

`receipt.json` stores all 228 representatives of the three audited sectors,
their actual automorphism orders and literal vertex orbits. Zero-error classes
are retained as positive controls, not discarded from the finite universe.
The additional all-equitable-crossings pass tests all 14^3*8=21,952 triangle
crossing/root presentations, including uniform crossings and nonstrong root
attachments. It verifies by full automorphism search on every strong case that
no positive strong error was omitted by the complete-support restriction.

## Information boundary

This classifies only the entire order-ten refinement-error coefficient, not
all order-ten tournament classes and not A000568(10). The correction still
needs the independent pointed upper trace to yield a count. The proof relies
on inherited actual stable-cell restrictions and the seven/nine coefficients;
it is not a recurrence from size alone, group orders alone, or an unmarked
metagraph. The finite crossing computations certify the small residual
universes rather than replacing the exhaustive partition proof.

## Frozen producer comparison and distribution extension

The producer packet was frozen before inspection at manifest SHA256
905dd04456f9b4cf302e8b71c7bdd915db209ba081b4735b0f7a733378866f43.
All fifteen manifest-listed payloads were byte-verified and physically copied
into `producer/`; no producer file was edited. Both producer checkers were
executed normally and with `-O` from that copy, and all four outputs matched
the respective pinned outputs exactly. Its structural proof and source agree
with the independent derivation above. No correction was required.

The new `distribution_audit.py` independently reconstructs all765 root-deletion
destinations from the pinned85 almost-regular-eight parents. It does not import
the producer's canonicalizer or distribution checker. It only compares
order-eight tournaments. The fifteen regular-nine classes have defect
histogram 3 at0, 3 at18, 2 at36, 7 at72, and diagonal trace135. Combining that
histogram with this audit's own strong-ten representatives yields exactly

    12 t^14 + 46 t^18 + 4 t^36 + 84 t^40 + 14 t^72.

There are160 positive-error ordinary classes, first moment5508 and raw second
moment229416. No class below order14 can contain two defective strong
components, since a defective strong component needs at least seven vertices.
Thus the A(z)^2 placement factor applies to each positive defect value here;
there is no missing collision term at order10. The producer's general raw
moment identity M2=A^2 M_s,2+2 A^3 M_s,1^2 follows correctly by applying
(t d/dt)^2 to 1/(1-S(z,t)) and evaluating at t=1.

Six semantic mutation controls were run on physical copies in both normal
and optimized modes. Every altered program failed at the intended mathematical
check, and the manifest verifier independently rejected its changed bytes:
removing the degree-two five-cell sector; removing triangle-cell permutations;
dropping cyclic-majority strong classes; setting actual error to zero; omitting
one SCC side; and halving the order-nine distribution padding. This is twelve
semantic failures and twelve byte-integrity failures, not merely a hash test.
`hostile-results.json` records the exact failure labels. The audit and producer
normal/optimized receipts all match; the final `replay.py` reproduces them
from another manifest-only physical copy.
