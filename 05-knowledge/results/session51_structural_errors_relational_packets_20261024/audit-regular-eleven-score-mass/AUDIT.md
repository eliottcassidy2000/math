# Independent audit: regular-eleven score-mass closure

24 October 2026. Verdict: PASS. The labelled score recurrence and the
orbit-stabilizer subtraction establish the proposed regular-sector closure.
This is an atlas-free score-count/full-group sidecar, not a theorem that bare
metagraphs or a scalar sequence determine the unrestricted order-eleven count.

## Frozen inputs and independence

The audited producer is `session51-regular-eleven-score-mass`, manifest
`42f8d34f2791b513c4ef906f836c6b9bd8a411c7db91c2d3c4e385715f49e6d9`.
Its thirteen manifest payload files and original manifest are copied under
`producer/`. The inherited complete nonrigid classification is
`session51-regular-eleven-nonrigid`, manifest
`3fd0d159b9e356fb5d86906e657e9ef992061408f2391419ceaf5d248813acaa`.
Its manifest and every listed payload are copied under `nonrigid/`.

`independent.py` was written separately before inspecting the producer's
frozen code. It does not import either producer's implementation. It uses no
almost-regular-ten parent list, regular-eleven list, unrestricted-eleven list,
canonicalization engine, or target census length. The previously reported
13,333-parent holdout is not an input or comparison oracle in this packet.
The new code derives the numerical closure without embedding either 1,223
or 13,333 as an expected count. The nonrigid completeness theorem is an
explicit inherited, source-pinned hypothesis, not reconstructed from a total.

## Labelled counting proof

Let N(d) count tournaments on a fixed ordered label set with prescribed
outdegree vector d. Permuting d preserves the count by a vertex-relabeling
bijection. Therefore sorting d is valid as a cache key, but does not identify
any tournaments within the set counted by N(d). This is the key distinction
between memoization of equivalent subproblems and quotienting the objects.

Distinguish one fixed vertex v of degree a. Its actual outneighbor subset X
has exactly a elements. If a remaining vertex of required degree r belongs
to X, it loses to v and keeps residual required degree r. Otherwise it beats
v and its residual degree is r-1. Every tournament has exactly one X and
one residual tournament; every permissible pair extends uniquely. Thus
summing the residual counts over all literal subsets X, with unit weights,
counts every labelled tournament once. There is no factor for choosing v.

The primary audit implementation distinguishes the largest-degree vertex and
enumerates these literal subsets. It uses only necessary degree feasibility
filters: bounds, total edge sum, and the fact that every k vertices contain
k(k-1)/2 internal edges. It does not assume that those filters are sufficient.

A second independently written implementation distinguishes the smallest-degree
vertex and groups equal degrees. For each degree r occurring m times, choosing
x named outneighbors has multiplicity binomial(m,x), and the residual multiset
contains x copies of r and m-x copies of r-1. The product of these binomial
weights exactly partitions the literal subsets into equal-score cases. This
second method deliberately does not use the primary prefix-sum pruning. Both
implementations retain arbitrary-precision integer counts.

The normal and optimized computations agree on:

- orders 1, 3, 5, 7, 9, 11: 1, 2, 24, 2,640, 3,230,080,
  and 48,251,508,480 labelled regular tournaments;
- the fixed assigned vector (4,4,4,4,4,5,5,5,5,5): 191,474,240.

An independent direct edge-orientation enumeration covers every labelled
tournament through order six, totaling 33,867 orientations. Both recurrences
match all 637 nondecreasing candidate score vectors, including impossible
vectors with zero count, and all 3,271 realized fixed label assignments.
The exhaustive comparison uses counts of each exact assigned score vector,
not sums over permutations, so it detects spurious factorial normalization.

## Fixed-root deletion identity

Deleting a particular root from a regular-eleven tournament leaves five
vertices of degree four and five of degree five. Conversely degree-four
vertices must beat the new root, while the root must beat all degree-five
vertices. This unique extension is regular and has the specified root label.

There are binomial(10,5)=252 disjoint choices of which five labels have degree
four. Hence

    L11 = 252 * N(4^5,5^5)
        = 252 * 191,474,240
        = 48,251,508,480.

The factor 252 is a labelled score-assignment multiplicity. It is not an
ordinary-class count and must neither be omitted nor introduced inside each
fixed-score subproblem.

## Coverage and ordinary group mass

The inherited theorem classifies all nonrigid regular-eleven tournaments by
full automorphism orbit partitions. Its elimination of the other odd orbit
partitions, the exceptional full-orbit exclusions, and the hidden transitive
Aut55 presentation are necessary to its completeness. The present audit does
not weaken that premise to a selected subgroup or equitable presentation.

Using only the pinned representatives, the independent code revalidates the
literal tournament edges and regular degrees, computes unrestricted ordinary
isomorphisms by injective adjacency search, and independently recovers all
eighteen full automorphism groups and vertex orbit partitions. Common-neighbor
signatures are used only as necessary rejection filters; every accepted map
is a full adjacency-preserving bijection. All 153 pairs are nonisomorphic.
The hidden q=3 presentation is isomorphic to exactly one transitive record and
is excluded from the final list once, not counted twice.

The recovered spectrum is:

- three classes of group order 11 and one of group order 55, each q=1;
- seven classes of group order 5, each q=3;
- four classes of group order 3 and one of group order 9, each q=7;
- two classes of group order 3, each q=9.

These are ordinary orientation-preserving classes. A positive relabeling
control passes, and four of the eighteen representatives are not isomorphic
to their converses, which the isomorphism procedure correctly keeps separate.

Each ordinary class T has exactly 11!/|Aut(T)| distinct labelled versions.
The complete nonrigid labelled mass is therefore

    M = 3*(11!/11) + 11!/55 + 7*(11!/5) + 11!/9 + 6*(11!/3)
      = 151,764,480.

Every remaining labelled regular tournament belongs to a rigid class, and
all such classes have labelled orbit size 11!=39,916,800. The exact remainder
is 48,099,744,000, so division gives 1,205 rigid ordinary classes. Adding the
eighteen nonrigid classes gives 1,223 regular ordinary classes. Divisibility
is checked but not substituted for completeness: integrality alone could
never rule out a missing nonrigid class.

The inherited pointed-class bijection and regular refinement definition give
78 nonrigid pointed types and error 660. A rigid class contributes 11 pointed
types and error 110. Consequently

    P10 = 11*1,205 + 78 = 13,333,
    E_regular11 = 110*1,205 + 660 = 133,210.

The full regular error polynomial is

    4 + 7t^22 + 5t^66 + 2t^88 + 1,205t^110.

The source-pinned nonregular/SCC constants independently give E_s,11=142,376
and E_11=152,690. These yield only the residual formula

    a11 = (U11 - 152,690)/11.

U11 and unrestricted a11 are not evaluated here. No larger-order completeness
or metagraph sufficiency is inferred from this calculation.

## Integrity, hostile controls, and relocation

`crosscheck.py` matches the separately obtained counts, fixed-profile factor,
nonrigid mass, six closure fields, and complete error histogram to the frozen
producer. It checks that all three producer source copies equal the applicable
files under the complete pinned nonrigid dependency. The nonregular offsets
are read from the pinned inherited proof, rather than a target-eleven atlas.

All checks use explicit exceptions and survive Python optimization. Normal
and `-O` outputs agree byte-for-byte for independent counting, crosschecking,
and fourteen independent semantic hostile controls. The controls reject:
a missing class, a nonliteral isomorphic duplicate, a subgroup mistaken for
the full group, incorrect vertex orbits, incorrect error, an invalid edge,
a valid tournament made nonregular, fractional rigid mass, negative rigid
mass, inconsistent pointed mass, omitted binomial weights, reversed residual
updates, and factorial division or multiplication at sorted cache states.
These are semantic checks, not merely hash failures. The frozen producer's
own semantic/manifest hostile suite is also retained and replayed.

`replay.py` copies only this packet's manifest-listed payloads and manifest to
a new directory. All scripts run with an unrelated current working directory.
It reexecutes both modes of independent computation, crosscheck, semantic
controls, integrity checks, the producer, and the inherited classification.
Both inherited packets' own relocation/hostile replays run as well. Changed
receipt bytes, changed producer code, changed inherited receipt bytes, duplicate
manifest paths, and path traversal are rejected in normal and optimized modes.
The final integrity check runs again after restoring the temporary mutations.

The relocation receipt is written outside the manifest's payload list to
avoid a self-referential digest cycle. No external workspace packet is needed
by a relocated replay. The audit's source and payload digests are in MANIFEST.json.
