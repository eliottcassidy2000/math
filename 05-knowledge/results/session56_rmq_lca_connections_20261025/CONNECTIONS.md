# Seventeen tested connections: RMQ/LCA, arithmetic, and tournament structure

25 October 2026. This is the requested list of connections, not a list of claimed open-problem solutions. **PROVED** means the stated bridge has a proof in the linked packet; **RECOVERED** credits an earlier result; **OPEN APPLICATION** means the proposed use is not yet realized in the native problem. Classical algorithms and theorems are credited, with no claim of external novelty.

## 1. Range minima and named witnesses

**PROVED / classical bridge.** A stable Cartesian tree preserves the index of every interval minimum; LCA is the corresponding query. It does not preserve values, ties, or arithmetic admissibility unless those remain separate data. This is a precise model for earlier work that had to retain the owner of a terminal event or exposed port rather than only its value. [Cartesian proof](cartesian-fibres/PROOF.md); [historical card 1](niche-atlas/ATLAS.md).

**Test:** [1,1] and [5,5] have the same index-query structure but different minimum values.

## 2. Cartesian trees, hook weights, and support versus multiplicity

**PROVED / classical mechanism rederived.** A binary shape T has exactly n!/product(subtree sizes) permutation preimages. The reciprocal-hook weights sum to one, while unweighted shapes are counted by Catalan numbers. This is a concrete explanation of why an exact quotient for one query does not automatically preserve enumeration weights. [Fibres and capacity](cartesian-fibres/PROOF.md).

**Boundary:** one n-node Cartesian tree cannot encode all tournament classes for all n (already fails at n=7); this is an encoding-capacity bound, not a lower bound on computing A000568.

## 3. RMQ overlap and the prohibition on reusing alternative choices

**PROVED and implemented.** Minimum and fixed-mask AND are idempotent, so overlapping blocks are harmless. Meet convolution of alternative local states is generally not idempotent: a repeated component can secretly choose two different states. The actual CRT query on [2,5,3] has 112 arrays, but naive overlap gives 784 and fabricates distinctness. Disjoint sparse tables and segment trees repair this without pretending the alternatives are deterministic. [Range algebra](range-signatures/THEOREM.md).

**Boundary:** costs include interface size, witness data and integer bit lengths, not just the number of abstract monoid operations.

## 4. Correlated copying, database joins, and shared ports

**PROVED / classical joining specialized.** Copying one chosen state twice is not the same as sampling two independent states. Natural join enforces equality of shared component choices, then evaluates each weight once. It recovers exactly the 112 genuine arrays in the preceding example, including zero or signed weights without unsafe division. [Correlated gluing](correlated-gluing/THEOREM.md).

Full choice IDs are sufficient, not universally minimal: a certified continuation-equivalent quotient may compress them if its correlated weights are retained correctly. A running-intersection join tree gives an all-size sufficient global compatibility criterion; cyclic constraints can defeat pairwise consistency. This is not an unconditional native metagraph realization theorem.

## 5. Möbius projectors and the difference between algebraic and physical states

**PROVED / semilattice algebra.** The partition-lattice zeta transform diagonalizes meet convolution. Signed Möbius columns give exact projectors. But every nonnegative real convolution idempotent is zero or one deterministic state; normalizing a nontrivial choice family into probabilities does not make independent overlap valid. [Projectors](partition-projectors/THEOREM.md).

**Test:** correlated X=Y preserves any distribution, whereas independent X,Y with X meet Y distributed like X force a Dirac law. Group averaging has different multiplication and is not a counterexample.

## 6. p-adic valuations as low-digit LCA queries

**PROVED.** A trie read from least significant digits gives LCA depth v_p(x-y). After sorting by those addresses, pairwise common-prefix depth is the minimum of adjacent prefix depths, an exact RMQ reduction. Repeated records, negative integers and truncated depths need explicit treatment. [Valuation carrier](valuation-lca/REPORT.md).

**Application recovered:** repeated Collatz-block fuel becomes comparison with a marked rational 2-adic path. This computes a stated itinerary's admissible repetitions, not universal Collatz coverage.

## 7. Crossing hierarchies and graph colouring

**PROVED.** Put an edge between two requested clusters when they intersect without nesting. The minimum number of same-label rooted hierarchies representing all clusters equals the crossing graph's chromatic number. Twelve affine lines on nine points require exactly four trees. [Tree-cover theorem](hierarchy-covers/REPORT.md).

**Boundary:** a synchronous intersection of prime filtrations is itself hierarchical and can use one tree, but loses the independent prime blocks as nodes. A five-cycle shows that maximum crossing-clique size is not always the answer.

## 8. Tournament modules and quotient-labelled LCA

**PROVED interface; classical decomposition.** The lowest common strong module tells where two tournament vertices first separate. Their two child ports and the quotient tournament at that module determine the arc. The tree alone does not. [Hierarchy queries](tournament-hierarchy/REPORT.md).

**Test:** three prime five-vertex tournaments share the same module-star shape but have 9,13,15 Hamiltonian paths. Quotient automorphisms also prevent replacing unlabelled prime substitution by naive scalar power-series substitution.

## 9. Electrical divergence, cycle space, and legal tournament moves

**PROVED all-size scoped application.** Ordered scores determine divergence and all net cut fluxes, but not circulation. If any forest of tournament arcs is locked, each fixed labelled-score fibre is connected by directed triangle and quadrilateral reversals. A quadrilateral is indispensable in a two-state five-vertex fibre with path counts 11 and 13. [Cycle-move theorem](cycle-moves/REPORT.md).

**Boundary:** locked arcs and scores survive; Hamiltonian counts modulo 8 need not. This adds a lawful degree-fibre move structure, not a solution of every native wiring guard.

## 10. Score fibres, graphical zonotopes, and unrooted forests

**PROVED / classical zonotope mechanism with exact application.** Realizable labelled outdegrees are the integer points of a graphical zonotope. Its Ehrhart polynomial is the forest enumerator, so the number of score fibres equals T_G(2,1), the number of unrooted spanning forests. For G=K_n minus the locked path, composition compression and an independent subset-tree calculation certify the sequence through n=14 without enumerating all orientations. [Forest counts](score-forest-counts/REPORT.md).

**Test:** a triangle has seven such fibres/forests, while det(I+L)=16 counts rooted forests with a different weight. Neither count is A000568.

## 11. The user's Cayley sequence and mirror bases

**RECOVERED from Session 31, independently rechecked.** The coefficients 1,3/2,8/3,125/24,54/5,16807/720 are a_m=(m+1)^(m-1)/m! for m>=1. With T=z exp(T), the rooted-forest generating function is F=exp(T)=T/z. For n>e, the distinct mirror base is exactly F(log(n)/n). [Recovered proof and source lineage](cayley-mirror/PROOF.md).

**Boundary:** at n=2 the principal branch returns the diagonal value 2, not the mirror 4. This session's connection is to the rooted/unrooted observer distinction in card 10, not a new discovery of the bridge or an MRB-defect formula.

## 12. Prefix minima, Dyck-style feasibility, and affine debt budgets

**PROVED exact summary.** Additive credit paths have summary(total increment, minimum prefix). For affine updates x↦lambda x-delta, the exact summary retains endpoint multiplier, endpoint debt and the maximum normalized prefix debt. Chronological composition is associative and supports disjoint range queries with exact clocks and congruence guards. [Debt range theorem](debt-range/RESULT.md).

**Test:** a loss of 2 followed by a gain of 2 has identity endpoint map but still requires initial credit 2. The native Collatz ledger includes positive charge; deleting it changes the process into a lower bound, not an equal orbit model.

## 13. Zeckendorf carries as update data

**RECOVERED arithmetic; OPEN native data-structure application.** A Fibonacci colour or charge can be updated correctly only with its carry cocycle. A prefix tree or persistent index does not make naive colour XOR a lawful addition rule. [Historical card 7 and checks](niche-atlas/ATLAS.md).

**Test:** two decompositions of the same integer can have different naive colours when the carry is discarded. The full integer/charge observer, rather than the colour alone, is the arithmetic interface.

## 14. Rule30 innovation scales and non-isometric tree equivalence

**RECOVERED scoped result; OPEN application.** A low-bit phase tree can be topologically correct while its ambient state metric uses a different innovation-depth scale. LCA supplies phase depth; a regrading map supplies physical ultrametric depth. [Historical card 5](niche-atlas/ATLAS.md).

**Boundary:** matching trees does not establish a new fractal dimension or Rule30 theorem.

## 15. Fibonacci/Berggren addresses and density distortion

**RECOVERED exact hostile.** The same ray can have a linear arithmetic address in one embedding and an exponential ternary-heap address in another. Ancestry survives while harmonic sums and densities change. [Historical card 6](niche-atlas/ATLAS.md).

**Lesson:** tree depth is not automatically integer size, arrival time, or available credit. Those maps must be retained before transporting asymptotics.

## 16. Tropical kernels versus scalar minima

**PROVED elementary algebra; OPEN application to full native kernels.** Min-plus matrix composition is associative but not generally idempotent. A range index may compose exact finite transition kernels over disjoint blocks, but overlapping them repeats transitions and costs. [Historical card 9](niche-atlas/ATLAS.md).

**Test:** the all-ones 2-by-2 kernel squares to the all-twos kernel. Taking a minimum also loses witnesses and multiplicities unless supplied separately.

## 17. Persistent histories and irreducible higher-order compatibility

**RECOVERED interface/hostile; OPEN engineering application.** A version-tree LCA identifies where editing histories diverged, not whether their mathematical boundary states glue. Pairwise interaction tables can also miss a first saving that appears only at arbitrary higher arity. Named boundary states and continuation semantics remain essential. [Historical cards 4 and 10](niche-atlas/ATLAS.md).

**Boundary:** these abstract counterexamples do not refute a particular native theorem without a map into its actual objects.

## Common mechanism and remaining obligations

These connections support a concrete methodological statement, not a universal theory of mathematics: a summary is exact relative to specified queries and continuation operations. Reusing a summary under a new operation can lose owners, correlations, weights, cycle data, or arithmetic scale.

The live obligations remain native all-guard tournament wiring/A000568 progress, uniform Collatz coverage, and genuine integer-square existence. The tested bridges give usable tools and sharp failure boundaries; none silently discharges those separate problems.
