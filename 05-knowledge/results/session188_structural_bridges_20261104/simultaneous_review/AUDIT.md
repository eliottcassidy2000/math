# Independent audit: simultaneous module reversals and pair hulls

**Verdict: PASS within every explicitly stated scope. No unresolved mathematical discrepancy found.**

Reviewed artifacts are `simultaneous_response/PROOF.md`, `ANCHORED.md`, their two checkers and full response banks, and `returns_next/PROOF.md` with its finalized D1–D5 checker. The incoming E1–E4 statements used by these proofs were inspected in the retained parity proof, including the distinction between exact scalar invariance and modulo-four ambient transfer.

## Disjoint and exposed-laminar separation

For a homogeneous module, reversing each internal path while preserving ordered cover slots proves exact q_k invariance. In a lifted query exposing one of several disjoint modules, every other selected module is still homogeneous with respect to both new ports. Its reversal therefore preserves the actual query count exactly. Incoming E4 then supplies the nonzero mod-four contrast for the exposed module. This also proves pairwise distinctness of the full disjoint reversal cube.

The singleton-port averaging lemma is sound. In a Hamilton word with x–A contacts ignored, reversing maximal A-runs preserves every other edge. The numbers of valid singleton attachments are |A|, 1, |A|-1, and 1 for zero, incoming-only, outgoing-only, and both contacts. In the last case the endpoints are distinct before and after reversal. The argument handles singleton runs and arbitrary outside attachments.

The all-order anchor battery is the stated finite-difference consequence of incoming E1: the two optional coefficients are d and n-1+d. At odd order either zero or two optional queries occur; at even order exactly one occurs. The n=2 E1 boundary is not silently used: the endpoint parity is one and provides its own separator.

For general nested pairs and exposed laminar families, the proof correctly lifts the whole-outer-module battery **first**, using only the inherited congruence. It then applies exact averaging in the actual ambient extension. Every selected proper module avoids the anchor, and all its indexed singleton queries share the same outside mask. Laminarity preserves modularity throughout the successive reversals. The proof now explicitly notes that every actual contrast is even, so a sum equal to 2 modulo 4 supplies an individual contrast equal to 2 modulo 4.

This does not cover a laminar family in which all selected maximal modules are fully covered by proper selected members. The stated unresolved boundary is accurate.

## Profile monotonicity and two crossing modules

D2's source ordered sum visits each occupied block at most once in each path slot. Reorder those segments by a fixed Hamilton order in the induced target quotient, reversing internal segments in conversed blocks. Physical block names and retained slot labels make the map injective. A transitive target quotient has the inverse map on every cover. For k=1, the separate quotient-Hamilton-order choices yield the stated multiplier.

For exactly two properly overlapping original modules, the intersection and union are modules. The three parts P,I,R have one uniform transitive order, not arbitrary pairwise incidences. Simultaneous edge XOR converses P and R internally, leaves I unchanged, and turns the three-part quotient into a directed triangle. Nonnegative occupied-run coefficients then preserve all profile gains and give the strict scalar lower bound from the positive first coefficient H(Q).

The uniformly split nested-pair construction is also valid: the fixed extra port lies between its two homogeneous parts before reversal, producing a cyclic three-block quotient afterward. The second port is a global sink and is forced last. This argument genuinely needs its uniform-cut assumption; the more general anchored proof supplies the broader nested result separately.

## Pair-generation and incidence independence

The maximal-proper-module proof of D5 is correct. Induced submodules of a module lift to ambient modules. If maximal proper modules are disjoint, a pair from distinct ones has full hull. If two maximal proper modules overlap, their union is the whole current module. A hypothetical proper module containing a vertex from each extreme part extends to a maximal D; union closure forces both extreme parts into D, and every intersection vertex distinguishes those extremes and must also enter D. This contradiction proves pair generation.

The chosen pair e_M lies in N exactly when M is contained in N. Ordering by inclusion gives a square zeta minor with diagonal one, hence an integer determinant one and independence over every field. The direct maximal-selected-module argument correctly proves that **no nonempty family of distinct nontrivial original modules has identity edge XOR at any order**. It does not prove full-response injectivity for all interacting families.

The sharp module-count characterization is also correct: a directed triangle's three pairs have the same hull, whereas every pair in a transitive tournament generates its distinct interval.

## Independent exact computation

`independent.py` imports no producer response functions. It generates complete Boolean polynomials from literal Hamilton-word attachment patterns, then compares all 270,316 named integer coordinates in every full response bank through order five. It separately enumerates original modules and their edge XOR families.

Results match exactly:

- 1,099 tournaments and 130,942 nonempty source/family presentations;
- 3,830 singles, 824 disjoint multi-families, 12,228 other laminar families, and 114,060 crossing families;
- 10,650 pair hulls and 3,830 modules, including every inclusion-zeta incidence;
- no identity edge XOR and no exact full-response collision in this finite universe.

For all four supplied larger controls, the independent engine reconstructs full polynomials and complete truth tables by literal word patterns, rather than the producer's subset DP. It matches every saved vector and literal query. In particular the six-vertex whole-plus-three-pairs control has 36 quadratic, 132 cubic, and 132 quartic nonzero differences. All 64 pair-uniform probes are zero modulo four, exactly 24 differ as integers, the query (3,12) gives 509 versus 545, and the nonuniform query (1,5) restores residue 2. These are not full-response collisions.

`anchored_independent.py` reconstructs the banks again, uses literal Hamilton permutations for endpoint parity, and validates every saved certificate: 15,818 exposed laminar families, 33,808 anchor certificates, and 3,828 exact averaging checks. It separately retains 6,396 two-vertex boundary certificates. The 1,064 fully covered laminar families are counted as outside the anchor theorem.

`overlap_independent.py` checks the claimed quantitative bounds for all 2,166 proper-overlap pairs and 3,580 uniformly split nested pairs through order five. Scalar and query values come from independently reconstructed full vectors, while block/quotient Hamilton counts come from literal permutations.

The independently constructed profile controls in `../interval_review/independent.py` additionally verify D2 on 320 arbitrary quotient/block/converse cases, with induced-subset Hamilton counts obtained from literal permutations.

## Optimization and wording fixes

All independent normal/optimized receipts match byte for byte. The finalized D1–D5 producer checker was also replayed in both modes with identical output. The supplied search and anchor normal/optimized receipts agree.

Two small presentation fixes were requested and verified: the six-vertex full basis is stated as 484 coordinates of bidegree **at most** (2,2), and the anchor proof explicitly uses even actual contrasts when extracting a single modular separator. The pure quartic basis at that order has 225 coordinates. Neither fix changed a computed result.

Artifacts include the three independent checkers and paired receipts in this directory, plus the D1–D5 normal/optimized replay outputs. No unbounded conclusion rests on the bounded census alone.
