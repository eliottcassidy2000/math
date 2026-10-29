# Toward a signed classification of tournaments and their native metagraph

**Status:** audited elementary proofs, recovered classical/repository mechanisms, and exact finite computations. This is progress on an operation-sensitive classification program. It does **not** complete the unrestricted prime-tournament atlas, classify the entire native metagraph, or settle the open Hamiltonian inequalities.

Source pin: `3b6e3ef82382da222fc59aa02914eea72d5a105b`.

## 1. What a useful classification must contain

The finite-simple-groups analogy becomes productive only after choosing what counts as a subobject and which operations must be preserved. Four different objects occur here:

| layer | atomic or recursive data | what it classifies | essential missing input |
|---|---|---|---|
| ordinary modular tree | oriented prime/linear quotient plus child isotypes | full tournament, recursively | unrestricted prime quotient library |
| signed cycle-index recursion | positive/negative permutation-cycle distributions and compatible child colours | specified attachment and converse orbit counts | does not determine tournament isotype |
| switching/bijoin structure | rank-one cuts, rooted normalized charts and their transporters | switching-stable decomposition and class symmetry | not the native legal-move graph |
| native metagraph | actual directed Hamilton paths, line colours, ordered ports and legal intermediate states | legal complement transitions | unrestricted legal transition rules are not supplied by these strata |

A decomposition formula is not an enumeration of its prime inputs. Existing A000568 and aggregate prime-symmetry formulas remain valid; the missing library here is a structural classification of prime isotypes, not the absence of an enumeration formula. An exact observable is not automatically a classifier. A quotient path is not automatically a legal path of the original system. This packet makes those boundaries constructive, with explicit tests and repairs.

## 2. Main advances and recovered mechanisms

### A. Sign changes have a precise recursive route

The [signed modular cycle-index formula](grammar/RESULT.md) extends the inherited child-isotype stabilizer recursion. For each signed quotient permutation, an l-cycle contributes the child index with sign ε^l and substitutes x_j↦x_(lj).

A tournament anti-automorphism has only even cycles of length2·odd and at most one fixed point. Therefore, in a negative recursive term, only one possible fixed-child branch continues with a negative state. All other child-cycle factors are positive. On an odd tournament the negative chain ends at its unique fixed vertex; on an even tournament it ends at a quotient without a fixed child. This concerns a chosen anti-map, not a canonical global branch of the unmarked object.

The formula gives exact marked one-vertex attachment counts and converse-fixed counts. It also gives a special unpointed descent for odd self-converse tournaments, using the unique automorphism orbit of vertices fixed by anti-maps. These are proofs with explicit scopes, not a claim of historical priority for self-converse enumeration.

Two hostiles make the classifier boundary sharp:

- Two order-six tournaments both have automorphism order3 but have24 versus32 marked attachment orbits.
- A transitive and a strong order-four tournament have identical **both** signed cycle indices. Signed indices are sufficient for the stated orbit counts, not for isomorphism classification.

### B. A complete extremal native-word classification at seven

The prior [tree-gauge packet](../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md), §8 Theorem3, already proves the switching-quotient generator injectivity from order8. We recover its consequence: with labels fixed, two native moves cannot yield a nontrivial vertex switch at order≥8. This obstruction is inherited, not rebranded as new.

The [extremal refinement](metagraph/SIGN_GRAMMAR.md) classifies every order-seven exception: the two Hamilton paths must be edge-disjoint and partition K_(4,3). Up to relabelling there is one undirected two-path skeleton. Its512 admissible within-part fillings split into64 blue-blue and448 black-black cases; mixed colours are impossible.

More precisely, three independent affine parity equations on nine free orientation bits give eight fibres of64. The zero syndrome is blue; the seven nonzero syndromes are black. The actual normalized skeleton symmetries fix these syndrome coordinates. An abstract resemblance to seven Fano points does not authorize a GL(3,2) action on the tournament configurations.

The same endpoint parity-graph method gives necessary bounds for words of arbitrary length, with separate cut and converse-cut cases. It retains every path guard and excludes the empty-cut cancellation case from inappropriate size bounds.

### C. A connected quotient need not connect native states

The [Cayley quotient](metagraph/CAYLEY_QUOTIENT.md) is recovered from the same prior packet and independently retested. Labelled switching classes form an elementary binary group. A directed path can be uniquely normalized inside each class; complementing its nonpath edges induces a fixed translation.

This quotient has two components when n≡1 mod4 and one otherwise, for n≥3. But at n=3 the actual labelled native graph has two disjoint K_(1,3) components while the switching quotient is connected. From n≥8 there are quotient edges that have no legal lift from a specified tournament, even after allowing all possible path witnesses for that edge.

At order 6, an exhaustive labelled census finds two 16,384-vertex native components separated by coordinate parity, while the switching quotient is connected. The all-word parity guard explains why an odd-size cut cannot be realized by native moves when n≡2 mod4.

This is why the switching classification cannot silently solve the native metagraph problem. The preserved sign information and the lost guards must both be stated.

### D. A literal simple group appears, in the correct category

The [signed-cut analysis](switching/SIGNED_CUTS_AND_SIMPLE_SYMMETRY.md) gives elementary rooted proofs and exact group identifications. An ordinary module can be destroyed by switching; rank-one cut structure is preserved. Choosing a root converts signed cuts to modules of a normalized deletion.

The switching-class symmetry group has root stabilizer equal to the ordinary automorphism group of that deletion. Its Sylow2 order divides n; allowing global converse gives a2-part dividing2n. These necessary action constraints sit inside the stronger classical Babai–Cameron switching-group classification, explicitly credited in the report and already used in prior repo work.

For a source adjoined to Paley7, exhaustive8! testing identifies the **actual permutation groups** as PSL(2,7) of order168 and, with converse, PGL(2,7) of order336. An independent conjugacy-class certificate proves simplicity of the former. This is a genuine connection to simple groups, not a resemblance between group orders. It concerns switching-class symmetry, not ordinary tournament symmetry or native metagraph automorphisms.

### E. Pfaffian recursion needs mixed ports, not just scalar energies

The [parity interface theorem](pfaffian/PARITY_INTERFACES.md) makes the sign changes algebraic. A rank-one skew cut factors the closed Pfaffian into ordinary Pfaffians for even/even blocks and rooted cofactor responses for odd/odd blocks. Switching and global reversal have explicit, different sign factors.

A four-vertex hostile has equal separate scalar root energies but joint Pfaffians1 versus3. The exact repair for an even interior A with specified ports U is the pair

    (Pf(A), Uᵀ A⁻¹ U).

An odd tournament interior retains one residual odd mode in addition to the mixed port matrix. For an even invertible interior, scalar closure is safe for isotropic port columns, including proportional columns, but not for arbitrary later attachments. This refines the inherited rooted-response warning; it does not prove H≥disc or an integral Smith-form classification.

## 3. How the pieces fit together

The proposed organizing theory is **a classification relative to the operations**, with three explicit obligations at every node:

1. **Object:** retain the oriented quotient and child types that determine isomorphism.
2. **Interface:** retain the ordered attachment data that future operations can observe.
3. **Symmetry:** retain the actual compatible transporter action, not only its order, rank or marginal statistics.

This is a research framework, not a newly proved universal theorem. The packet supplies several exact instances and a completed local stratum:

- signed modular recursion provides the symmetry rule;
- rooted switching charts identify the right cut decomposition;
- Pfaffian transfer proves a finite interface closure for specified ports;
- the K_(4,3) theorem completely classifies the extremal two-step switching pattern;
- the Cayley/native discrepancy proves why the guard layer cannot be discarded.

The independent audit retains both major hostile controls: the twelve-vertex pair with67 identical cyclic triples but automorphism orders81/27, and the all-strong run-order pair with Hamilton counts2,899,675 and2,906,091. Neither is erased by these new summaries.

## 4. Precisely what remains open

The [source and target map](history/PROVENANCE_AND_TARGETS.md) distinguishes already achieved classifications from open frontiers. In particular, the prior converse-exception family theorem is a closed predicate-specific classification; an unrestricted prime atlas is not.

Useful next mathematical targets are:

- a native transition grammar that updates marked ports and module splits after arbitrary run interleavings;
- mixed-word integral glue, rather than only powers of one operator or first-prime-level ranks;
- minimal read-only attachment carriers, rather than summaries that are sufficient only after adding reconstructive operations;
- realized prime-cycle-index data and deletion-defect strata, without treating a prime census as a free oracle.

These are remaining research questions, not pending verification or claimed breakthroughs. No complete classification or solution of a major open conjecture is announced.

## 5. Reproduce and inspect

Run `python -B verify_all.py` from this packet. It runs every listed mathematical checker normally and optimized, checks identical output, and compares to the saved receipts. The code is standard-library Python. Full producer and independent test universes are stated in their reports; finite cases are never extrapolated into unproved all-order claims.

- [Independent audit](audit/AUDIT.md)
- [Verification record](verification.json)
- [Pinned historical sources](history/SOURCE_BLOBS.json)

The source tree was kept unchanged. Publications use isolated new paths and preserve incoming work. Exact remote bytes, manifest hashes and local relative links are checked after publication. Empty reported CI statuses are not presented as a CI pass.
