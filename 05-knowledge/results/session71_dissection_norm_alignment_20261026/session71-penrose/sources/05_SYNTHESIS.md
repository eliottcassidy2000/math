# Operation-relative minimal carriers: what this session closes

Session70, 2026-10-26. This synthesis distinguishes new deductions, inherited mechanisms, cited classical inputs, and finite checks. It is not a claim that all the compared open problems are equivalent or that a matching numeral identifies their proof obligations.

## 1. The organizing question has three genuinely different answers

A carrier is not smallest in isolation. Fix the objects, the named operations, the failure convention, and the observations. Two states can be identified exactly when every allowed continuation gives the same observation, including failure. The coarsest such equivalence is the inherited operation-kernel construction; its concrete classes depend on the experiment.

This session exhibits three sharply different instances:

| Category and experiment | Exact retained state | What can really be discarded |
|---|---|---|
| Tournament with a fixed named Hamilton path; every named permutation can be queried for legality | The complete tournament presentation | No two distinct free-bit masks have the same full guard language |
| Mandatory disjoint nontrivial paths with distinct named endpoints; unrestricted fresh connector paths | The perfect matching of exposed endpoint pairs | Internal shapes and lengths, for this Boolean feasibility observer |
| The inherited factorial graceful family; named edge-overwrite requests, no-ops allowed, only validity observed | Color-to-component permutation | All independent left/right realization signs |

The first result is a reconstruction theorem, the second a boundary minimality theorem, and the third a genuine nontrivial guard quotient. It would be wrong to turn the first into a universal statement that guards recover every graph, or to turn the third into a refutation of tournament reconstruction.

## 2. Why a whole Hamilton-path language can recover unused arcs

For any finite directed graph D on a fixed named vertex set, let L(D) be its spanning directed Hamilton permutations. Define its traceable-arc core C(D) as the union of the arcs occurring consecutively in those permutations. Then

    L(C(D)) = L(D),    C(C(D)) = C(D).

The proof is two inclusions: core arcs are original arcs, and every accepted original path uses core arcs. Thus equality of full languages is equivalent to equality of traceable-arc cores. A tournament arc not appearing in any Hamilton path disproves the naive identity C(T)=T, not this exact equivalence of information.

The derived all-order tournament theorem uses a classical structural input about nontraceable arcs in a strong tournament (Busch–Jacobson–Reid), alongside Camion's Hamilton-cycle theorem. The detailed report and independent audit state and source those inputs precisely. The reconstruction mechanism is:

1. A missing arc x→y has a separating vertex z. In T−z, x is in the first strong component and y in the last, with at least two intervening components.
2. A vertex in an intervening component can be deleted without destroying strong connectivity of T. Consequently every arc incident with that vertex lies on a spanning Hamilton path, by cutting a Hamilton cycle of the vertex-deleted tournament and adjoining the vertex at the appropriate end.
3. Two distinct intervening vertices supply two internally disjoint length-two core paths from x to y. For any deleted vertex other than x and y, at least one detour avoids it. Replacing each missing arc by such a path shows that the traceable-arc core preserves all reachability after any one vertex deletion, within a strong tournament.
4. The separating z therefore exposes the orientation of the missing pair from the core alone: x reaches y in C(T)−z, while y does not reach x. Any other tournament with that same core has the same reachability and must use the same orientation.
5. For a general tournament, its ordered strong components are recovered from its core. Hamilton paths concatenate componentwise, so the strong reconstruction applies within each component and the component order determines every cross arc.

This is a decoder and an all-order proof, not extrapolation from a through-seven census. Its classical input is explicitly CITED; no literature-priority assertion is made for the resulting reconstruction deduction.

On fixed-path presentations X_n, there are exactly 2^M states, M=(n−1)(n−2)/2, for n≥1. Since first-step whole-permutation guards already distinguish every state, any exact operation-stable quotient retaining those guards has singleton cells. M bits are the exact fixed-length information cost for valid presentation states in this gauge. These counts exclude the distinguished failure symbol; an explicit totalization of the partial operations adjoins a separate failure sink. This says nothing by itself about the optimal number of oracle queries, efficient prime generation, or the size of a weighted response representation for a weaker observer.

There is also an exact algebraic distinction. Over any field, the unital POINTWISE algebra generated by all named Boolean guards is the entire function algebra on X_n: products of separating guards or their complements construct each singleton indicator. Its dimension is 2^M. This does not say that the raw guards linearly span that space. Already at n=6 there are at most 6!=720 named guards but 1,024 states, so their raw linear span cannot be the full function space. Pointwise closure, linear response compression, and implicit circuit storage are separate requirements.

One-step elementary relocation guards are a weaker observation. The guard report supplies a first collision at order four. Full legal relocation traces recover the full language using the inherited connectivity theorem, so time/continuation can recover information absent from one local observation. The proof must keep left/right operation names and their transported coordinates; fresh anonymous path sampling is a different experiment.

## 3. The Vitali connection becomes a native correlation theorem

The previous C4 construction already had 2^k distinct alignments at order 7k+2 with identical separate rooted endpoint types and colour. This session sharpens it in two ways.

First, the common strong source is rigid at every k. Loss of a common source-target frame is therefore not explained away by source automorphisms. An equivariant-selector obstruction on an unpointed symmetric object and an alignment ambiguity between independently recorded endpoint models are different mechanisms.

Second, every Boolean conjunction of hidden alignment bits is an actual simultaneous one-step re-marking guard. The guard-by-state matrix is the Boolean zeta matrix, with determinant one over the integers. Its rank is 2^k over every field. Restricting to conjunctions of at most t named blocks leaves rank sum_{j≤t} binom(k,j), with precisely the complementary number of missing linear coordinates for arbitrary signed laws.

Uniform even- and odd-parity laws agree on every proper-block marginal but disagree on a genuine global guard. This supplies the previously missing native realization for the old finite Boolean-moment obstruction. It does not invoke classical Vitali nonmeasurability: the tournament universes here are finite. Complete named moments reconstruct the full law by Möbius inversion, including explicit nonnegativity constraints for actual probability laws.

Do not conflate costs. An individual alignment needs k bits, and singleton guards read them. A linear representation capable of evaluating all guards of arbitrary laws needs dimension 2^k (or 2^k−1 nonconstant affine coordinates). The exact parity gap shrinks with k; no uniform approximation-hardness claim follows.

## 4. Square-sum and graceful work explain different parts of the same loss

The relevant square-sum graph is Q_n on 1,...,n with an undirected edge exactly when the distinct endpoints sum to a square. Its path-cover transfer already retains named endpoint pairing. For unrestricted connector contexts, different perfect matchings of 2r exposed endpoints are distinguished by a connector that closes a cycle in one state but merges components in another. The minimum number of states is (2r−1)!!, giving r log2 r+O(r) fixed-length bits.

The unrestricted qualification is essential. A given square-sum graph may not realize every pairing or every distinguishing connector. The actual arithmetic birth recurrence therefore inherits a sufficient joint interface, not automatically this unrestricted minimum. Singleton endpoint coincidences and multiplicities remain explicit.

The graceful-tree recovery supplies a second caution: the arbitrary frozen-color width-two connectivity question was already refuted by an actual 12-color (13-vertex) fiber with five completions split 3+2. The upper-tail connectivity theorem survives, and neither statement proves every prescribed abstract tree graceful. In the factorial family, the new guard quotient retains the color-to-component permutation but genuinely forgets all left/right signs under its stated no-op-permitted validity-only protocol.

## 5. Petersen, K5, and K3,3 are port tests, not universal obstructions

Contracting Petersen's five spokes gives K5. A spoke cluster has two local vertices, and quotient edges attach to different ports. A one-visit lift would require alternating outer/inner edge types around a five-cycle, impossible by parity. That alone does not rule out all spanning traversals: clusters may be revisited. The full one-/two-run response restores the exact six 2-factors, each two disjoint five-cycles. Global connectivity is the final missing condition.

Deleting one Petersen vertex and suppressing its three degree-two neighbors gives K3,3. A quotient Hamilton cycle is spanning in the expanded graph only if it uses all three marked suppressed edges. Of the six quotient Hamilton cycles, exactly two do, corresponding exactly to the two spanning nine-cycles of the vertex-deleted Petersen graph.

Planarity, Hamiltonicity, and gracefulness separate even on these examples. K5 and K3,3 are both nonplanar and Hamiltonian, but K3,3 admits an explicit graceful labeling while K5 fails a simple odd-difference count. Petersen has Hamilton paths but no Hamilton cycle. None is intrinsically a tournament, and completing missing pairs by arbitrary orientations destroys the original feasibility question.

The common compositional rule is to sum over compatible JOINT boundary assignments, retaining local traversal alternatives and the global acceptance predicate. It is not to multiply marginal counts or identify all graph obstructions.

## 6. What stays open

- Economical generation and classification of unrestricted modular-prime tournament atoms.
- Optimal query complexity and useful algorithmic compression, despite the exact named-state lower bound.
- Minimal carriers for specific restricted arithmetic attachment languages and anonymous fresh-sampling protocols.
- A compact universally colour-, complement-, and re-marking-compatible aggregate algebra. Session69 already closed the ordinary uncoloured all-arity response recursion; that result is preserved.
- The Graceful Tree Conjecture and other prize-level questions remain untouched by equivalence claims.

The concrete closure here is the full named guard reconstruction question and the three precisely typed carrier results above, supported by native witnesses, independent checks, and reproducible artifacts.
