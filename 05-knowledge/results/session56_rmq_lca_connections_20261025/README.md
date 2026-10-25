# Session 56: RMQ/LCA as exact interfaces, not universal summaries

25 October 2026. Second verified checkpoint; research is ongoing. The user's request is to connect range-minimum queries and least-common-ancestor algorithms with older niche themes and current arithmetic/tournament work. Classical algorithms are credited, and no external novelty or A000568 speedup is claimed.

## 1. Overlap-safe signature aggregation

[Proof and implementation](range-signatures/THEOREM.md); [independent audit](audit-range-signatures/AUDIT.md).

The overlapping-block RMQ trick transfers to fixed selected collision masks because AND is idempotent. It does NOT automatically transfer to families of alternative choices. Meet convolution of a family with itself can choose two different members from the same factor and fabricate a repair.

The exact CRT hostile uses the ordered coprime factors [2,5,3]. Modulus 30 has 112 square-magic arrays, none with nine distinct entries. The overlapping query duplicates the factor 5, obtains 784 arrays, and falsely permits distinctness. Component identities and shared cell labels are essential. Correct disjoint sparse tables and segment trees are supplied, with operation-order and provenance guards.

Finite families are idempotent under meet convolution exactly when meet-closed; a reduced finite antichain is idempotent exactly when empty or a singleton. Complexity is stated in monoid operations AND interface/payload costs, not as a free constant-time solve of the entire choice problem.

The producer's assertion-based verifier intentionally rejects optimized Python. The independent exception-safe audit passes normally and with -O; no unsupported producer optimized pass is claimed.

## 2. A genuine tournament LCA carrier

[Structural interface](tournament-hierarchy/REPORT.md); [independent audit](audit-tournament-hierarchy/REPORT.md).

In a canonical strong-module hierarchy, the least common module of two vertices determines the quotient where their paths separate. The actual arc is recovered from the two child ports and that quotient tournament. LCA indexes the hierarchy; it does not supply the missing quotient orientation.

All 1099 labelled tournaments through order 5 are independently checked. Three prime order-five tournaments share the same module-star shape but have 9,13,15 Hamiltonian paths. A separate minimal fixed-Hamilton-order hostile has identical ordered degrees and cut flux, yet 11 versus 13 paths. Reversing a directed four-cycle explains the degree-preserving change.

The producer is a linear, assertion-based reference verifier, validated in normal mode. Its optimized output is not treated as assertion validation; the independent audit validates both modes. Its O(log n) indexing interface is proved and independently implemented in the audit; optimal classical LCA indexing is cited separately. SCC ordered-join composition and prime substitution require different counting rules, with quotient automorphisms retained. None of these results solves native all-size metagraph wiring or replaces the data needed for A000568.

## 3. A source-pinned list of niche connections

The [ten-card historical atlas](niche-atlas/ATLAS.md) supplies a source, target, map, retained/lost coordinates, hostile example, cheapest test and status for each connection. Topics include persistent versions and named ports; Rule30 innovation ultrametrics; Fibonacci/Berggren address distortion; Zeckendorf carry cocycles; finite Collatz itineraries; tropical range composition; and irreducible higher-arity interactions. Applications marked OPEN are not promoted merely because their elementary control passes.

[Independent atlas, valuation and fibre audit](audit-valuation-atlas/REPORT.md) checks the current source pins and records two repaired wording errors. Cartesian structure preserves argmin indices, not minimum values unless the values are retained separately.

## 4. Valuation trees and exact quotient fibres

- [Valuation/LCA carrier](valuation-lca/REPORT.md): low-digit-first p-adic tries turn equality depth into LCA; adjacent common-prefix values answer pair queries by RMQ. CRT needs named partitions across primes. A finite marked prefix gives only a capped fuel bound, not global Collatz coverage.
- [Cartesian fibres](cartesian-fibres/PROOF.md): a shape T represents exactly n!/product(subtree sizes) permutations. Shape weights normalize to one; treating shapes uniformly loses source multiplicity. A precise cardinality bound limits tree-only encodings of all tournaments, but is not a lower bound on computing A000568.

## 5. Algebra and the minimum number of hierarchies

- [Partition projectors](partition-projectors/THEOREM.md): zeta coordinates diagonalize meet convolution. Formal signed Möbius projectors exist, but a nonnegative real idempotent is only zero or a deterministic state. Independent resampling and correlated reuse are different operations.
- [Hierarchy covers](hierarchy-covers/REPORT.md): the minimum number of same-label rooted hierarchies containing a given nonempty-cluster family equals its crossing graph's chromatic number. Affine lines on nine points require four trees. A synchronous common-refinement filtration can still fit one tree; it does not retain all independent prime clusters as nodes.
- [Independent algebra/cover audit](audit-projectors-covers/REPORT.md) checks the signed and probabilistic statements, tree-cover constructions, actual Session55 nine-cell CRT example and odd-cycle hostile. The hierarchy producer's assertions are validated in normal mode; its optimized output is not treated as assertion validation. The independent exception-safe audit validates both modes.

## Verification and scope

Exact proofs, finite universes, independent certificates and source pins accompany both lanes. No whole external article is republished. Earlier Session 55 results remain unchanged. Publication uses temporary storage because the original overlay is full; no Git pack data is deleted. The session will have a separately announced end after further connection cards and audits are packaged.
