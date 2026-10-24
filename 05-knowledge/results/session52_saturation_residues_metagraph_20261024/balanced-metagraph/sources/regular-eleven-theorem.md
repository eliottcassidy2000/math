# Regular-eleven correction from almost-regular parent rerooting

24 October 2026. **PROVED:** the restricted parent normal form, atlas construction,
reroot trace identity, full-group recovery, and converse correspondence.
**FINITE-EXACT:** the attached bounded computation and distributions. Independent
audit is a separate packet. No literature-priority claim is made.

## 1. Result, inheritance, and boundary

The complete, action-equipped nine-vertex atlas certified in
`session51-parent-nine` and `session51-audit-parent-growth` yields:

    almost-regular ten-vertex classes M10 = 13333;
    regular-eleven parent-only reroot trace tr(D) = 13453;
    regular-eleven classes R11 = tr(D)/11 = 1223;
    regular-eleven error Ereg11 = 11(M10-R11) = 133210.

The independent nonregular-eleven classification supplies Es11=Ereg11+9166
and E11=Ereg11+19480. Consequently

    strong order-eleven total correction Es11 = 142376;
    unrestricted order-eleven total correction E11 = 152690;
    a11 = (U11 - 152690)/11.

This closes the correction, NOT U11 or a11. No full ten-vertex atlas is
constructed or imported. No eleven-vertex canonicalization is performed.
Only a small, explicitly selected class of ten-vertex parents is normalized.
The source retains entire adjacency matrices and complete permutation actions;
this is not prediction from unmarked metagraphs or scalar counts alone.

The closest inherited mechanism is Session50's regular-nine parent-only trace.
The hostile obstruction is a nontransitive regular cell, which stable color
size alone cannot resolve. The corrected near miss is the cyclic-cell decoder
that forgets crossing phases or the full crossing matrix. The indispensable
sidecar here is the FULL row matrix under all intrinsic internal rotations.
The working concept board is: degree-filtered extensions; intrinsic stable
cells; complete matrix normal forms; pointing/unpointing; and converse symmetry.

The source-atlas producer manifest SHA-256 is
161449ddf5c2a07b7afefd213438b85c74897850b89050dde66428819b0838d6;
its independent audit manifest is
2878857405e9d1e6f6d999db9e35deabc420afbe48861f5b5f277a92f30245b0.
The full byte-pinned source contract is in `SOURCE_PINS.json`.

## 2. Exact universe and exhaustive input filter

An almost-regular ten-vertex tournament has exactly five vertices of outdegree4
and five of outdegree5. Deleting a specified vertex gives a nine-vertex parent
with all degrees in {3,4,5}. Conversely let P be a nine-vertex tournament and
S the set of old vertices beating a new root. The extension is almost regular
if and only if:

- every degree3 vertex of P belongs to S;
- no degree5 vertex belongs to S;
- degree4 vertices may be included or excluded;
- |S| is4 or5.

The last condition gives root outdegree5 or4. The total degree sum45 then
forces five vertices of each degree. The code additionally tests that exact
multiplicity, not merely the two allowed degrees.

Isomorphic rooted extensions correspond precisely to Aut(P)-orbits of these
eligible subsets. The set is invariant under the complete supplied parent
actions. Enumerating the least subset in each such orbit is therefore exhaustive
and duplicate-free as POINTED almost-regular types. Forgetting this root may
still identify several states; the normal form below performs that quotient.

The source atlas has191536 classes. Exactly10252 pass the score filter, with
135418 eligible labeled masks and133130 eligible subset-orbit states. These
are emitted results, not expected-value tests in the production executable.
The initial cost gate found the135418-mask upper bound before the full run.

## 3. A complete restricted parent normal form

Stable directed refinement starts with outdegree colors and splits by vectors
of outcounts into the currently ordered color classes. Canonical lexicographic
color numbering makes each final cell intrinsically named. Since refinement
only splits, each cell stays in one of the initial size5 degree classes.

Every stable cell induces a regular tournament: every vertex has the same
internal outdegree, and counting internal edges gives degree(m-1)/2. Thus m
is odd. Here m is1,3,or5. The regular tournaments on these sizes are the unique
consecutive cyclic tournaments C1,C3,C5, and each automorphism is a rotation.

For C3 the unique outneighbor of v is its intrinsic successor. For C5, the two
outneighbors form a directed edge, and its source is the intrinsic successor.
To prove uniqueness for regular order5, write v's outneighbors as a->b and
inneighbors as c->d. Degree2 forces b to beat both c,d, c to lose to a, and d
to beat a. Thus the cyclic order is v,a,b,c,d, with each vertex beating its
next two vertices. This also proves that the intrinsic successor forms a
5-cycle. Every automorphism commutes with successor and is therefore a rotation.

List the actual cells in invariant color order, give each its intrinsic cyclic
order, and enumerate all independent rotations of these orders. For EVERY
candidate relabeling retain and compare the entire tuple of outgoing row masks.
The lexicographically least full tuple is the normal form. There is no
crossing-degree, spectral, sign-only, gain-only, or coarse fingerprint quotient.

An isomorphism preserves the intrinsically ordered stable cells, their successor
relations, and hence belongs to this enumerated rotation family. Conversely
equality of full normal-form matrices supplies an actual isomorphism.
Therefore the normal form is an iff invariant for the specified universe.
The minimizers f_i from an input to its normal form give its complete canonical
automorphism group as f_i composed with f_0 inverse. Every such map preserves
all arcs; every automorphism gives one unique minimizer. The code verifies
permutation validity, all arcs, identity, distinctness, and group closure.
These checks supplement, rather than replace, the completeness argument.

Each initial five-cell can split only as5,3+1+1,or five singletons, giving
rotation factors5,3,or1. Hence at most25 candidates are required per parent.
For the lower controls, initial size3 or4 yields at most9. The executable's
normalizer is compile-time restricted to orders6,8,10 and runtime restricted
to the corresponding almost-regular scores and cell structures.

Deduplication obtains13333 almost-regular ten-vertex classes. Its full-group
order histogram is {1:13282,3:38,5:8,9:5}. The total number of pointed source
states over each canonical parent is independently checked to equal the number
of vertex orbits of its complete recovered group.

## 4. Complete labeled mass, independently checked

The deduplicated parent mass is

    sum_P 10!/|Aut(P)| = 48251508480.

On the source side it is also

    sum_Q (9!/|Aut(Q)|) * number_of_eligible_labeled_masks(Q).

Both count labeled almost-regular ten-vertex tournaments with a fixed new
label10: the first through unpointed classes and the second through deleting
that fixed label. The source atlas is already proved complete and distinct;
the complete normal-form iff establishes distinctness in the output. Mass is
an additional coverage check, not an independent substitute for these proofs.

`degree_mass.py` independently counts labeled tournaments realizing a prescribed
outdegree vector. Delete its first labeled vertex of degree d; choose d of the
other labels to lose to it. The other labels lose one unit of residual degree
exactly when they beat it. Sum the recursively counted residual degree vectors.
Equal-degree labels are grouped only to multiply identical subproblems by the
appropriate binomial coefficients; no unlabeled quotient is taken. Sorting a
degree vector is valid because relabeling preserves the number of realizations.
Landau inequalities provide only a safe impossible-sequence rejection filter.

It obtains191474240 for a specified assignment of five4s and five5s, hence
binom(10,5)*191474240=48251508480 over all assignments. Independently, the
all-five degree vector on eleven labels gives the same48251508480. Equality
also follows from the unique regular-extension bijection below. These are
fresh calculations, not imported target values.

## 5. Unique regular extension and the reroot matrix

Every almost-regular even parent P on2k vertices has a unique regular extension
on2k+1 vertices: its degree(k-1) vertices beat the root, which beats its degreek
vertices. Parent isomorphisms extend uniquely over the distinguished root.
Conversely deleting any vertex of a regular tournament yields just such a
parent. Thus almost-regular parent classes are exactly pointed regular-child
types, one for each vertex orbit of an unpointed child.

For each canonical ten-parent i, construct this extension T_i, delete each of
its eleven roots, normalize ONLY the resulting ten-parent, and set

    D[i,j] = number of roots whose deleted parent has normal-form index j.

Every row sums to11. Equality of deleted-parent types is exactly equality of
rooted regular-child types, by unique extension. Therefore the nonzero support
of a row is the complete set of pointed states over one child; all rows in that
component are identical, and their entries are its vertex-orbit sizes.
If a component has q states, its diagonal sum is11 and its stable-cell error
is sum_over_orbits(11-|O|)=11(q-1). Hence

    R11=tr(D)/11;       Ereg11=11(M10-R11).

These are trace/component identities, not target canonicalization. Components
are checked disjoint and exhaustive, all within-component rows identical,
every component's diagonal sum and row sum eleven, and all deletions in the
complete targeted parent atlas.

The parent full group is the stabilizer of the root in the regular child.
Consequently |Aut(T_i)|=D[i,i]*|Aut(P_i)|. This gives the full group ORDER,
not merely a subgroup order, without searching the target's permutations.
It agrees for every state in a component. The resulting histogram is

    order1:1205; order3:6; order5:7; order9:1; order11:3; order55:1.

The full regular labeled mass sum11!/|Aut(T)| again equals48251508480.

## 6. Error distribution and completed correction

The regular-eleven orbit-size multisets and multiplicities are:

- eleven singletons:1205;
- eight singletons and one triple:2;
- five singletons and two triples:5;
- one singleton and two fives:7;
- one orbit of size11:4.

Thus the regular positive-error polynomial is

    7 t^22 + 5 t^66 + 2 t^88 + 1205 t^110.

The possible44 sector is empty. Its first moment is133210. Adding the separately
proved nonregular total polynomial gives the COMPLETE positive-error polynomial
at order eleven:

    48 t^14 + 348 t^18 + 7 t^22 + 8 t^36 + 256 t^40
      + 5 t^66 + 28 t^72 + 2 t^88 + 1205 t^110.

There are1907 defective ordinary classes in total, with first moment152690.
The independent nonregular strong polynomial similarly combines with this
regular polynomial; `combined-correction.json` records both exact distributions.
The regular sector is strong because a nontrivial ordered strong-component
decomposition forces disjoint degree intervals, impossible at equal degree.

Every tournament automorphism group has odd order: an involution would swap
some pair and reverse its edge. All vertex orbits are therefore odd, and at
odd total order their number q is odd. Accordingly every regular-eleven error
is divisible22, an independently predicted support check.

### 6.1 Structural exclusion of the absent error44 sector

In fact the five observed orbit patterns are the only possible patterns for
ANY regular-eleven tournament, independently of the census. Actual automorphism
orbits have odd size and induce regular subtournaments. Between orbits of sizes
a_i,a_j, constant outgoing degrees r_ij satisfy

    a_i*r_ij + a_j*r_ji = a_i*a_j,
    sum_(j != i) r_ij = (11-a_i)/2.

If the sizes are coprime, their crossing is uniform. Enumerate the odd-part
partitions of eleven; the following are impossible:

- 9+1+1: a fixed vertex has degree0,1,9,or10, never5.
- 7+3+1: the fixed vertex has degree0,3,7,or10, never5.
- 5+3+3: the five-orbit has internal degree2 and must beat exactly one triple.
  The other triple beats all five vertices, giving its vertices degree at
  least1+5=6, a contradiction.
- 7+1+1+1+1: every fixed vertex has degree at most3 or at least7.
- 5+3+1+1+1: each fixed vertex has external degree0,3,5,or8, so its degree
  within the three fixed vertices would have to be2 or0. Three such numbers
  cannot sum to3, the number of internal edges.
- 3+3+3+1+1: a fixed vertex has degree a multiple of3 plus0 or1, never5.
- 5+1+1+1+1+1+1: each fixed vertex has degree0 or5 within the six fixed
  vertices. Their sum15 would require three degree5 vertices; at most one
  vertex of a tournament can beat all other vertices.

The remaining partitions are exactly11,5+5+1,3+3+1+1+1+1+1,
3+1+1+1+1+1+1+1+1, and eleven singletons. In particular no regular-eleven
class can have five vertex orbits, so error44 is impossible analytically.
These restrictions concern ACTUAL automorphism orbits. Feasibility of an
arbitrary equitable partition does not assert that it is the orbit partition;
the regular-five all-singleton partition is a simple warning. The computation
establishes the multiplicities, not the preceding necessity statement.

## 7. Converse controls and a general fixed-point corollary

Reverse every parent arc and normalize its full matrix to obtain the converse
involution c on parent states. The computation verifies c(c(i))=i and
D[c(i),c(j)]=D[i,j], so c also induces the exact converse involution on child
components. It never merges ordinary classes during construction or tracing.

There are229 self-converse almost-regular ten-parents and229 self-converse
regular-eleven classes. The regular self-converse error histogram is

    error0:4; error22:3; error66:5; error88:2; error110:215.

The equality229=229 has an all-order explanation, not a numerical coincidence.
For a self-converse odd tournament T, the group of automorphisms and
anti-automorphisms has order2|Aut(T)|, with |Aut(T)| odd. An involution in that
group must be an anti-automorphism a. It fixes exactly one vertex: the number
of fixed vertices is odd, while two fixed vertices would force their edge to
reverse. The involution a permutes Aut(T)'s vertex orbits. Each invariant orbit
has odd size and hence contains an a-fixed vertex. Therefore exactly one orbit
is invariant. Such invariant orbits are exactly the self-converse vertex-pointed
types, since any two anti-automorphisms differ by an automorphism.

Thus each self-converse odd T has exactly one self-converse pointed type;
non-self-converse T has none. For regular odd T, the unique extension identifies
these pointed types with self-converse almost-regular parents. This proves the
equality of their counts at every odd order. For even tournaments no
vertex-pointed self-converse type exists, since an anti-involution has no fixed
vertices. No claim of novelty is attached to this elementary group argument.

## 8. Controls, reproducibility, and limitations

The same engine, with no hardcoded class totals, gives:

    order7:  M6=5,  tr(D)=21,  R7=3,  Ereg7=14;
    order9:  M8=85, tr(D)=135, R9=15, Ereg9=630.

The latter agrees with the archived Session50 certificate, including orbit
multisets. Source rows and actions are pinned byte-for-byte, and the smaller
input text is reconciled with archived certified JSON. Normal/O2 and
O3/NDEBUG production runs are byte-identical across all receipts, parent
atlases, state rows, and components. Normal and Python -O degree-DP runs agree.
Six hostile source corruptions are rejected in both compilation modes. All
mathematical guards use explicit exceptions and survive optimization.

The main order-eleven run makes293126 ten-parent normalizations, evaluating
340206 full matrix candidates. It constructs no unrestricted order-ten atlas,
uses no eleven-canonicalizer, and never imports or calibrates against R11 or
a11. It DOES retain and normalize full almost-regular parent adjacency. The
normalizer is special to this small-cell universe and is not a claimed general
polynomial-time tournament-isomorphism algorithm. No forward scalar recurrence,
unmarked-metagraph sufficiency, U11 evaluation, or a11 prediction is proved.
