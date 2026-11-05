# Independent global H541 supergraph and rebasing audit

**PASS.** No tournament of Hamilton count 541 can contain the specified minimum-clock seed as an induced subtournament. This is genuinely the unrestricted induced-supergraph statement, not merely the earlier three-port-language exclusion. The separate equal-count rebasing graph has both claimed one-birth children, with the same four exposed ports. No returning grammar is inferred from those two children.

## Frozen source binding

The source is `session197-hybrid-constructions/clock-supergraph-barrier/`:

- `PROOF.md`: `0083018e9156f1105b73dbafc4531eb54769c0ed68c75e90f50ea3e1e43411a6`.
- `census.cpp`: `733e7eb551f3be49f811622cd6b79d30e02d29eb7bbf053c85f449590bcc603e`.
- `verify.py`: `cf80900fc73fee79aa286f283fa02d413f68a93252113be0321522d0d1c9f10d`.
- `VALIDATION.json`: `f502eb547e2303662e82c43475c491766801626ff3a28fac81fb884426bfe017`.

The source core helper has the same frozen hash `3204eeb1e9d2b05cac9bcc21c98f25697d10f925d6ff743393e945ce4e73615c` as the earlier pruning packet. No source file was changed.

## Source proof and representation checks

The seed is the strong eight-vertex graph `[2,12,145,149,99,143,175,19]`, with count 179. Any induced copy lies wholly inside one strong component of a containing tournament. At total count 541, primality forces that component's count to be 541: it cannot have count one because induced monotonicity gives at least 179. Other components are singletons. Removing them preserves the entire embedded seed and the relevant count. Thus arbitrary singleton padding is handled before applying a strong-order bound.

The official Busch Theorem 1 was checked directly in the preceding audits and has exactly the required scope. The integer inequality `5^11<=541^3<5^12` bounds the remaining strong component by order twelve. At most four added vertices therefore suffice for a hypothetical witness. No stronger bound is applied to the full weak graph.

The source enumerates every orientation of each new vertex, including sources and sinks, so it does not accidentally require the embedded seed to have a strong one-vertex extension chain. It fixes all seed labels pointwise and permutes only added vertices when canonicalizing. These are exact seed-fixed isomorphisms, with at most 24 tested permutations. Transporting any future new vertex along such an isomorphism gives one of the unrestricted extension masks of the representative. Deduplication therefore loses neither an embedding type relevant to this search nor a possible continuation.

The cap is safe at every stage: each intermediate graph is induced in the final graph, and canonical vertex insertion injects its Hamilton words into those of any extension. The same reasoning justifies the C++ counters' early rejection when an induced subset already exceeds the target. All branch checks are exact integer counts.

The full source replay reproduces the inventories `1,33,68,105,144` at orders eight through twelve. All 301824 arbitrary extension candidates are checked by both endpoint recurrences. Normal and optimized source runs match the frozen receipt byte for byte.

## A separate seed-preserving strong-extension proof

The independent census uses a different, smaller search justified as follows. Let `S` be any current strong induced graph containing the seed, inside a hypothetical strong target `U` of count at most 541. If every vertex outside `S` is unmixed to `S`, strongness of `U` supplies vertices `a,b` with

    S -> a -> b -> S.

Indeed, the outside vertices split into those dominated by `S` and those dominating `S`; both sides must exist and some arc must run from the first to the second. The induced cyclic two-vertex extension has exact count

    H(S)+P_S(2) >= 3H(S)+2^|S|-2 >= 791.

The exact formula follows by separating the three one-block Hamilton patterns from the two-nonempty-block pattern. The lower bound uses the two whole-set terms of `P_S(2)` and at least one path for every proper nonempty subset. Induced monotonicity contradicts the target cap.

Thus at every proper stage some outside vertex is mixed to `S`, and can be added while preserving strongness and the original induced seed. This proves completeness of a mixed one-vertex chain **under this cap**, without misusing the unrestricted strong deletion lemma for a fixed embedded seed.

This independent search deduplicates by full exact isomorphism rather than by fixing seed labels. That is safe here because all subsequent mixed orientations are allowed and transport under any isomorphism; the existence of an induced seed copy is preserved. The independent checker additionally locates an actual seed copy in every retained graph.

The resulting complete strong inventories are:

| Order | Mixed extensions tested | Accepted before deduplication | Strong isomorphism classes |
|---|---:|---:|---:|
| 8 | initial seed | initial seed | 1 |
| 9 | 254 | 31 | 27 |
| 10 | 13770 | 6 | 3 |
| 11 | 3066 | 0 | 0 |
| 12 | 0 | 0 | 0 |

No target 541 occurs. Every one of the 31 retained graphs has its count recomputed by vertex-cover inclusion-exclusion over directed walks, independent of the source endpoint DPs. The three order-ten counts are 511, 517 and 529. Optimized and unoptimized C++ audit builds, together with normal and optimized Python checks, give identical receipts.

## Physical rebasing checks

The different strong parent `[2,64,11,51,39,135,61,95]` has count 179 and is not isomorphic to the original seed. Using common ports `[0,1,2,5]`:

- The left birth losing only to vertex 2 gives newborn out-mask 251 and count 539.
- The right birth beating vertices 0,1,5 gives newborn out-mask 35 and count 541.

The literal rows, strongness, counts and unchanged old-port list were checked. Parent and both children were counted independently by directed-walk inclusion-exclusion and full vertex-permutation enumeration. The original seed's positive count-539 birth was also checked separately.

This proves two scalar children from one actual replacement graph and one common exposure. It does not imply that either child returns to the same two-child class, nor that the original seed can be changed into the replacement by adding vertices. The global obstruction correctly shows why that latter inference would be false.

## Reproduction

Run `check_strong_extensions.py`; use `STRONG_AUDIT_DEBUG=1 python -O check_strong_extensions.py` for the unoptimized checked C++ build. Source replay copies are confined to this audit directory.

- Independent C++ checker: `149d3333901e58f5deb580b560b5843d2aa4226588e67d96042f616d62d37e29`.
- Independent Python checker: `9cb04f5d3a41348a5f1ee80d6bb4af43fb42c8be493894b7bd3fc87aef80fe37`.
- Independent receipt: `5cdf017365e7119f07cb6d77e48b0f53f8562d4137bc9133423df82e377c91f4`.

The exclusion concerns the exact induced seed and its converse, not every count-179 tournament or an unrestricted gap at 541.
