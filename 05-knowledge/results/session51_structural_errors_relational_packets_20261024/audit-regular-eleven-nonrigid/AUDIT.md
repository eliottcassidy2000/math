# Independent audit of the nonrigid regular-eleven classification

24 October 2026. PROVED structural completeness; FINITE-EXACT full-group and
ordinary-isomorphism certification. The result is independent of a known
regular-eleven class total and of the 13,333-parent reroot calculation. The
latter is used only as a subsequent holdout. No publication is performed.

## Result and gain in the prediction interface

Exactly eighteen ordinary nonrigid regular-eleven classes occur, with

    error0:4; error22:7; error66:5; error88:2;
    full Aut order3:6; order5:7; order9:1; order11:3; order55:1.

They contribute78 pointed types. Their unpointing deficit is

    Delta = sum_nonrigid (11 - number_of_vertex_orbits) = 120.

Every rigid class contributes11 pointed types and deficit0. Therefore, writing
P10 for the number of almost-regular ten-vertex classes,

    R11 = (P10+120)/11;
    E_regular11 = 10P10-120;
    E_s,11 = 10P10+9046;
    E_11 = 10P10+19360.

The last two use only the independently audited nonregular constants9166 and
19480. No U11 is supplied and no unrestricted a11 is predicted. This replaces
the need to reroot every regular extension by a complete small exceptional
classification plus the single parent count P10. It is not an all-order
scalar recurrence and does not discard the source adjacency/action data.

The eighteen classes contribute total error660. Thus P10=11r+78 and
E_regular11=110r+660 for the number r of rigid regular classes. This also proves
P10 congruent1 modulo11, E_regular11 congruent0 modulo110, and E11 congruent10
modulo110. The same modulo110 statement does not follow for U11.

## Independent finite computation

`independent.py` is newly written and uses none of production's support code.
It generates all2,040 five-by-five matrices with every row and column sum3 by
literal row-subset enumeration. It computes every translation orbit, finding
80 trivial-stabilizer orbits and8 order-five-stabilizer orbits. Production
instead starts from40 affine-line-invariant matrices. The independently found
8 exceptional matrix orbits give seven error22 classes and one hidden
vertex-transitive class.

The independent code generates the five triangle-pair candidates from full
attachment subsets and two affine matching signs. It separately exhausts all
64 labeled order-four tournaments, checks the unique strong-four class by
full isomorphism, and constructs both permitted crossing matches. Finally it
generates ALL32 circulants and uses unrestricted ordinary isomorphism between
them, not a multiplier-orbit theorem or a circulant-isomorphism assumption.

Full isomorphisms and automorphisms use dynamic minimum-domain bijection search.
Exact adjacency and common-outneighbor signatures are rejection filters only;
every surviving full bijection is returned. Candidate deduplication retains
orientation, not converse. All47 raw candidates give18 distinct classes. Normal
and optimized Python outputs match. `crosscheck.py` independently recomputes
full groups on every accepted producer representative and proves a complete
one-to-one class correspondence, not just agreement of the count or histogram.

## Exhaustiveness proof audit

The objects are FULL automorphism orbits, not arbitrary equitable blocks.
Each is odd and internally regular; every crossing is biregular, and coprime
orbit sizes force uniform crossing. The degree obstructions exclude all odd
partitions of eleven except

    11; 5+5+1; 3+3+1^5; 3+1^8; 1^11.

In particular five actual orbits are impossible, hence regular error44 is
analytically impossible rather than merely absent from a computed sample.
The exclusions are recorded and audited in the separate correction packet.

For5+5+1, the singleton intrinsically names A and B by s->A and B->s. Both
five-orbits induce C5. The A->B matrix has degree3. Its stabilizer in C5^2 is
nontrivial and projects to both factors; an axis element would force all rows
or columns identical, impossible for degree3 margins. It is thus a noncoordinate
order-five line. The affine-line construction covers every full-orbit instance.
Our complete2,040-matrix enumeration independently verifies that no other
nontrivial stabilizers were discarded. Full target groups are required because
one displayed5+5+1 presentation actually belongs to a vertex-transitive class.
The producer keeps and correctly removes that hidden-symmetry overcount.

For3+3+1^5, each fixed vertex must beat exactly one triangle and have internal
fixed-set degree2. The fixed set is therefore C5. If A->B has crossing row
degree d and A beats fixed subset S, then |S|=4-d. Swapping the triangle names
uniquely normalizes d to0 or1. The uniform case gives one four-subset orbit.
The matching case gives two three-subset orbits and two affine matching signs,
with phase removed only by actual triangle rotations. Thus one plus four
candidates exhaust this sector. Their unrestricted full groups have exactly
the advertised two triple orbits and five fixed vertices.

For3+1^8, the fixed vertices induce four degree2 vertices L and four degree5
vertices H. Exactly two edges run L->H. Internal L degrees are at most2 and H
degrees at least1. Each side is the unique strong4, except potentially a
triangle->sink on L or source->triangle on H. In either exceptional case the
triangle has zero matching incidence, identical external adjacency, and an
independent C3 rotation, contradicting that those vertices are fixed. This
is the essential full-orbit exclusion. Both sides are therefore strong4 and
rigid. The two crossing edges match two named degree-one L vertices to two
named degree-two H vertices. Exactly two matchings remain, and full groups
confirm the intended nine orbits.

For one orbit, transitivity on eleven vertices implies an automorphism group
order divisible by11. Cauchy's theorem supplies an element of order11, an
11-cycle. In its coordinates the graph is one of the32 circulants. Full
ordinary-isomorphism testing reduces these to4 classes. This proves exhaustive
coverage without using a prime-order circulant CI theorem. The order55 member
is the same hidden-symmetry candidate found in the five-cell sector.

Every other class is rigid and contributes deficit zero. Summing the eighteen
exceptional deficits gives120, with no assumption about how many rigid classes
exist. Hence the compressed formulas above are independent of the full atlas.

## Complete-atlas holdout, after the structural result

`atlas_control.cpp` reuses the earlier audit's incoming-refinement unrestricted
parent normalizer and its distinct full-automorphism search, with the original
source bytes archived. It never calls the renamed unused original main.
No target-eleven canonicalization is performed.

For each of the independently generated18 classes, all11 deleted-root parents
are normalized against the complete13,333-parent atlas. All198 destinations,
their multiplicities, and full group orders agree. The eighteen candidates
hit EXACTLY every nonrigid full-atlas component, with no omissions or extra
components. The computed deficit120 then recovers1,223 regular classes and
error133,210 as holdout results, not inputs to the structural generator.
O2 and O3/NDEBUG receipts agree byte-for-byte.

## Hostiles and portability

Controls verify a genuine relabeling, require that ordinary isomorphism does
not fold a non-self-converse graph into its converse, and reject a missing
nonrigid class, a nonliteral isomorphic duplicate, a one-sided invalid arc, and
a valid tournament made nonregular by reversing one arc. They run normally
and optimized. The frozen producer's own semantic controls also remain included.

Only files listed in MANIFEST.json are copied for `replay.py`; both independent
Python modes, both C++ control modes, cross-bijections, semantic hostiles, the
producer checks, and manifest-corruption rejection run from that fresh copy.
Source manifests and their relevant payload entries are pinned independently.
The replay receipt remains outside its own manifest to avoid a hash cycle.
