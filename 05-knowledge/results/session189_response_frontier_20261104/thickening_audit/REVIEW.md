# Independent review: dyadic transitive thickening

**APPROVED for the stated source-relative original-module theorem.**

The proof in `thickening/PROOF.md` has been reviewed independently, including the all-order structural argument rather than only its finite controls. Its frozen SHA-256 is `b2dfa4eb00502ec6b667e8d7f3a3bd305b8fa97a1daa195551121158e26d6ea3`. The final source includes the explicit qualifier that the actual augmented-profile identity holds modulo 2^L_Q. This precision repair changes no conclusion or retained controls; no further repair was required.

## 1. Exact conclusion and dependencies

For every finite nonempty named source T, let L_T be the maximum of 1 and the reviewed Seidel exponents 1+v2(E(Q)) for its original prime quotients. For every positive multiple N of 2^L_T, the source T[TT_N] has a response different from every distinct target obtained by XOR of modules of that inflated source. The module family is fixed in the original inflated source and may nest or cross; newly created sequential modules are not allowed. This is not a claim about the uninflated source, arbitrary tournaments, or pairwise separation of two unrelated family targets.

The Seidel theorem and its independently reviewed substitution closure are genuine inputs. They provide every complete normalized child profile modulo the needed modulus, a nonzero actual child-uniform local query, and an odd ambient multiplier. The reviewed Session188 hereditary-profile and highest-node theorems are separate inherited inputs. The finite controls below do not replace any of those all-order theorems. Exact local dependency hashes are retained in `DEPENDENCIES.json`.

## 2. Structural proof checks

At a prime node, a module crossing an inflated child has quotient support containing at least two names. Its support is a quotient module, so it must meet every child. A vertex omitted from the partially met child would then have to see all other quotient children in one direction, making its quotient name a source or sink. A prime quotient of order at least three is strongly connected, so this is impossible. A proper union of several entire children is likewise excluded by primality. Thus the maximal proper strong children and the prime quotient persist, with exactly the inflated physical child sets.

At a linear node the canonical children are ordered strongly connected components. A nontrivial strongly connected quotient remains strongly connected under arbitrary nonempty substitution, including connectivity within one cell through a directed quotient cycle. Singleton components expand into the N ordered singleton components of TT_N. Recursing gives no new prime quotient types. The proof distinguishes strong modules from strongly connected tournaments correctly and handles singleton and prime-free sources using L_T=1.

The profile closure then applies to each actual prime child A[TT_N], not merely to the individual TT_N cells. For a hypothetical equal response, the actual all-zero port query gives scalar H equality; inherited hereditary equality gives the complete q-profile on every fixed physical subset. A highest linear change is already excluded. At a highest prime change the quotient is the original Q or its converse, corresponding target child profiles agree, and the Seidel separator applies.

Other highest changed nodes are disjoint. The selected ports imitate the selected node at every exterior physical name, so each other node remains a module in that very same actual port extension. Replacing each other node by its equal-profile target preserves the queried count exactly. All changed edges lie in the highest nodes, so those replacements reach the full target. This validates simultaneous original-module XOR and does not assume that selected modules remain modules under an arbitrary sequential ordering.

## 3. Independent exact computations

`audit.py` imports no producer code. It uses a different module criterion: for every pair inside a candidate module, the two outgoing rows must agree on its complement. It independently derives strong modules by the non-overlap definition and maximal proper strong children, distinguishes linear quotients by the absence of a directed triangle, and explicitly verifies that every nonlinear strong quotient is prime.

The audit independently checks all physical inflated arcs and cell names, every old and inflated module, every strong module, and every prime quotient with its exact physical children. There are 155 literal inflations, 635,370 nonempty subset checks, and 90 prime-node/child comparisons, with maximum physical order 18. The complete independent data are retained in `structural_recount.json`.

A separate endpoint subset recurrence and unordered-cover partition recurrence compute complete profiles modulo 16 for all 11 labeled base tournaments of orders one through three inflated by TT4. In all cases q_1,...,q_5 are (1,14,4,8,0), and all q_h with h>=6 vanish by the factor h!. These actual closure controls are retained in `profile_closure_controls.json`. They corroborate substitution closure; neither these N=4 cases nor the producer's N=2,3 structural cases are claimed to satisfy every energy-based exponent or to enumerate a full large module-family orbit.

Normal and optimized (`python -O`) audit runs pass byte-identically. All checks use explicit exceptions rather than optimization-sensitive assertions. The complete retained data, source, proof, and receipt hashes are fixed by the publication allowlist and checksum manifest.

## 4. Reproduction

Place the producer at sibling `thickening/` (the original local name `session189-thickening/` is also supported). Run:

    python3 audit.py > normal.json
    python3 -O audit.py > optimized.json
    cmp normal.json optimized.json

Existing retained independent artifacts are compared exactly, not silently regenerated with changed contents. The proof review is an all-order argument; these finite runs are independent structural and profile controls, with precisely the scope stated above.
