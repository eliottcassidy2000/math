# Independent review: all unit ternary clone sources through order nine

Status: **PASS**, no correction. The complete finite catalog was reconstructed independently. Its scope is sources of order at most nine; no larger-order amplitude exclusion follows from it.

## Frozen source binding

Source directory: `../../session197-hybrid-constructions/unit-clocks-through-nine/`.

- `PROOF.md`: `100a774901dc03ec3e684607560d67feb1b31ed2d6dc3e3f4f62a897056491fa`.
- `census.cpp`: `e9c609e952d84fda5f25d4218d10a161b5d99d04fa37da69ed8e9de8bb5ed9d8`.
- `verify.py`: `d135cf1075420385d46e7698562a03ce73f16cf1c6dd84d3a1a25d4fdfce3277`.
- `graphs.py`: `3204eeb1e9d2b05cac9bcc21c98f25697d10f925d6ff743393e945ce4e73615c`.
- `catalog.json`: `87ffadc1351902637e10daf79ed68d4f9a316e925f73ec123914f9644f695cb4`.
- `VALIDATION.json` and `VALIDATION_OPTIMIZED.json`: `97c736a29169173f02be48f071009df13481a4b6a90be12f28078bdf29f0a450`.

The source was copied before replay, including its temporary-build location. Normal and Python `-O` runs reproduce the frozen receipt. All original source hashes remained unchanged.

## 1. Mathematical completeness

The matching/run-degree proof is valid. Distinct internal root-to-root fragments contain vertex-disjoint crossing edges. For the reverse inequality, a maximum matching supplies internal fragments; extending one edge through the remaining vertices preserves the order of its matched endpoints. Trimming before the first O vertex and after the last I vertex leaves legal outer prefix/suffix pieces and a legal internal O-to-I path. The zero-matching case is separately handled by an ordinary Hamilton path.

The largest run index therefore has positive leading exponential coefficient. An eventual positive pure ternary ray forces matching rank two, and independence of bases one, two and three gives exactly a3=A, a2=3A, a1=3A+beta. No larger leading base can cancel, and no m=0 interpretation is needed.

For order nine, both root sides have at least two vertices. The degree splits 2+6, 3+5 and 4+4 cover the possibilities up to converse. Converse preserves every positive clone count. Its unbalanced degree counterparts must be restored explicitly, as the source does; balanced converses are already in the complete balanced physical universe.

Every row/column isomorphism of a crossing pattern corresponds to a root-preserving relabeling, and all orientations inside both sides are enumerated. Thus quotienting only the crossing matrix loses no rooted tournament. Source canonicalization exhausts every permutation in its invariant root-incidence/outdegree classes, so it is exact.

Universal source/sink padding preserves every clone count, giving the complete amplitude list for smaller sources as well. This reasoning does not bound arbitrary larger sources.

## 2. Independent full generation

`check.cpp` enumerates every labelled bipartite matrix of each required shape, tests its matching rank, and constructs its full orbit under explicit row and column permutation group actions. This differs from the source's sorted-row tuples and column-minimization generator.

It independently obtains:

    shape  labelled rank-two matrices  pattern orbits
    2+6    3963                        42
    3+5    3448                        38
    4+4    2912                        37

Every internal orientation is then enumerated, giving all 3,215,360 presentations, with unit-hit counts 0, 40 and 72 in the three cases. Unbalanced converses are added.

The independent rooted canonicalizer uses iterative invariant color refinement, starting with outdegree and a distinguished root color, followed by exhaustive permutations within each final class. It does not use the source's root-incidence/outdegree partition alone. Every resulting full adjacency class is compared against the source catalog in both directions.

The result is exactly sixteen rooted classes, twelve strong and four weak, with profiles:

    (60,-1):4; (100,-1):2; (156,-1):2; (180,-1):4;
    (72,+1):1; (120,+1):2; (200,+1):1.

## 3. An independent coefficient recurrence

The source DP revisits an explicit root state. The independent counter instead removes the root and counts exterior permutation words with optional root insertions in their gaps.

An initial exterior vertex v contributes 1+z when v lies in O, otherwise 1. An exterior transition u,v contributes

    [u->v] + z [u in I and v in O].

A terminal exterior vertex in I contributes the extra factor 1+z. Each selection of gaps gives exactly one valid root-skeleton word, and every such word has exactly this exterior order and gap selection. This gives an independent polynomial Hamilton-path DP on only the eight exterior vertices, with no explicit root state.

The five coefficients are packed into 24-bit lanes of an unsigned 128-bit integer. Each coefficient is bounded by 8! times binom(9,k), at most 5,080,320, below 2^24, so ordinary additions cause no lane carry. Coefficients above degree four do not affect the lower lanes; rank two in the full census guarantees a4=0. The implementation checks that vanishing for every presentation.

As an independent primitive control, the gap counter agrees with literal root-word enumeration for all 5,405 labelled rooted tournaments through order five, padded by universal sources to the fixed nine-vertex format. The empty exterior word at the one-vertex boundary is explicitly counted once.

## 4. Actual graphs and neutral padding

For every one of the sixteen catalog graphs, `controls.cpp` independently counts:

- Its root-deletion Hamilton paths.
- Its actual positive clone graphs at m=1,2,3,4, giving sixty-four physical counts.

These counts use inclusion-exclusion over unrestricted fixed-length directed walks, not either source or independent root-run DP. Every result agrees with A*3^m+beta and the catalog's deletion value.

For each weak graph, literal reachability SCCs give exactly one singleton and one strong eight-vertex component containing the root. The singleton is an actual universal source or sink. Removing it yields one of exactly two rooted strong eight-vertex classes, each occurring once with source padding and once with sink padding. All four have the negative amplitude60 profile. Direct clone counts on those eight-vertex cores verify that the padding is neutral. Thus the weak-source description is checked explicitly, not inferred merely from four false strongness flags.

## 5. Replay and artifacts

The complete independent census and physical controls were each compiled and run under optimized flags and under `-O0 -D_GLIBCXX_ASSERTIONS`; their respective receipts agree byte for byte. The source normal and optimized Python replays also agree with the original frozen receipt.

- `check.cpp`: SHA256 `59b3306a7b5c9506b369ab96b70d9752b70d4336c2539e8e073f3104f16d1be2`.
- `independent.txt`: SHA256 `1ef7b9e82ab2c503abc89fd1c952fca73d2c84cafe7b7222998754732e84ffc3`.
- `controls.cpp`: SHA256 `f2ce1e7c1fe8e1bb6e2f9c7ffb1fc9ddc36f3bb97fa9d4d05f96cf3dc0a3ddb0`.
- `controls.txt`: SHA256 `3dd38e0648f7db087c4ec6ae10f2dd16e6544841644912e8900ad098491254b3`.

The C++ files include the earlier independently checked audit helper `../check_c3_complete.cpp`, SHA256 `f5d974967398108a01b98e298bbc1a266db41e6e1b94b9abd52331ce5ccb2946`, solely for elementary graph/canonicalization support. Their small text inputs `source-catalog.txt` and `source-run-catalog.txt` are direct row/profile transcriptions of the frozen source catalog. No external graph bank, repository mutation, network access or unverified numerical library is used by the independent full census.
