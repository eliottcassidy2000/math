# Independent audit: finite representatives for rooted four-states

**Verdict: PASS for the elementary theorem and nonsingular numeric inverse corollary.** The finite checks below are implementation controls, separate from the proof. Repository reviewed at 9f9bac37c83c839cb5fd3fe2220991232a87b2b2.

## Proof audit

The producer's conventions D=H(A), E=paths ending at root, U=ordered covers (P,R) with P possibly empty and R starting at root, d=H(A-root) are consistent throughout. A singleton root has (1,1,1,1).

1. Deletion SCCs are totally ordered. For each block B, write H=H(B), Q=sum over subsets S of H(S)H(B-S), and alpha=sum of the same covers whose nonempty second path begins at a root-outneighbor. The quantity V=U-d counts exactly covers whose root path has at least two vertices. In a new block, an already active root path and the other path append their block segments with Q choices. A newly active root path leaves the root for this block, leaving the entire previous deletion graph to the other path, giving alpha*d choices. Thus (d,V) becomes (Hd,QV+alpha*d). Empty first paths cause no exception: they simply start in the new block if its assigned segment is nonempty.
2. Q>=2^|B| follows from Hamilton path existence in every induced subset, including the unique empty path. Alpha is positive exactly when the root beats some vertex: choosing that vertex as the second path singleton already proves positivity, without any stronger endpoint Hamilton-path theorem.
3. Before first activation all blocks beat the root and all later blocks. Thus the inactive prefix is an ordinal factor C with all four coordinates multiplied by H(C). Removing any singleton components of C changes H(C) by factor 1 and preserves all four coordinates. This does not justify deleting arbitrary later singleton components.
4. Camion supplies |B|<=H(B) for nontrivial strong B. Such blocks have order at least 3; their H-values have product d. Their total order is at most d if present. If d=1 there are none.
5. After first active block, V>=1. The final V is at least 2 raised to the total order of subsequent blocks. Hence that suffix has order <=floor(log2 V). The surviving singleton components consist of at most the first active singleton plus this suffix. Add the root and the nontrivial-component budget: n<=d+floor(log2 V)+2 for d>1, and n<=floor(log2 V)+2 for d=1. If V=0 the root is a sink, all states equal d, and deleting every deletion singleton gives n<=d+1, or n=1 when d=1.

The supplied bound is safe even though it double-counts nontrivial suffix blocks. The proof does not need to optimize it.

## Decision and inverse scope

The integer-coordinate checks and fixed-root enumeration give a terminating membership decision: every physical tuple has a bounded representative, and every accepted witness is physically checked. The producer's find_state has been inspected; all 12 distinct states through order 4 return correct witnesses. Five negative boundary tests returned None: (1,1,1,0), (1,2,1,1), (1,0,1,1), (2,0,2,1), (1,0,3,1).

For a specified invertible rational matrix that is a genuine rooted-graft transfer, computing its unique inverse candidate and deciding membership is valid. Do not broaden this conclusion to arbitrary singular rational matrices: their affine preimage can have unbounded d and U, so the present bound does not give a finite enumeration of that preimage. The result is also not a structural graft recognition algorithm for a supplied graph, nor a decreasing recursive inverse, nor scalar universal coverage.

## Prior result and optional strengthening

The existing result at 05-knowledge/results/camion_busch_gaps_polyhedra_collatz_20261002.md section 4 (lines 190-197) already proves scalar H=m decidability using deletion of all singleton SCCs and the Camion sum-product bound, and records a Busch logarithmic refinement. It does not track a root, E, or U; deleting all singleton SCCs is not justified for a rooted four-state. The new scoped ingredient is simultaneous inactive-prefix compression plus the V suffix budget. Session197 already proves exact numerical membership and an O(K^3) inventory under U<=K for its physical sources, explicitly defined in section 1 as rooted tournaments with exactly one root-start Hamilton path. Sections 2–3 use the complete strong start-one grammar and weak transitive tails, and supply an injective three-response signature. That restricted theorem does not subsume arbitrary rooted four-state membership. Its tracked proof was recovered with git show HEAD:session197/source-signature/PROOF.md; initial sparse working-tree absence was not a missing repository source. A portable snapshot is inputs/S197_SOURCE_SIGNATURE.md. No broad priority claim is made.

The optional inequality |B|<=3 log_3 H(B) follows algebraically from the cited Busch minimum formula: for 2a+3b=|B|-1 with a+b>=1, 3 log_3(3^a5^b)-|B|=a+(3 log_3 5-3)b-1>=0. Summing gives the claimed logarithmic nontrivial budget. This audit did not independently obtain the primary Busch paper; retain the producer's cited-input qualification. The elementary theorem avoids that dependency entirely. If implementing the logarithmic bound, use certified integer comparisons rather than floating-point floors.

## Independent finite controls

`audit.py` is standalone, uses literal recursive lists of every directed path, and imports neither producer code nor repository helpers. It independently computes SCCs via transitive closure and recomputes all rooted prefix states from literal paths. It verifies Q>=2^blocksize, alpha positivity, every prefix recurrence, the total nontrivial budget, all four coordinates after compression, sink diagonal, and the representative bound.

Coverage: all labeled tournaments with root 0 through order 5 (1,099 cases, hence every rooted isomorphism class there), 130 deterministic cases each at orders 6 and 7, and an order-12 stress case with eight added inactive source singletons. Both normal and optimized Python executions passed with byte-identical receipts. No assert statements are used. Digest of independently checked cases: 4173dc35932ebadccd04a6a997fe51295f13537e2517200a4af6f61770586128.

The producer's separate exhaustive n<=6 verifier and receipts were inspected, not rerun by this audit. Its first and end path dynamic programs are independent of this literal-path audit's counting method.

## Audited source bindings

The final producer proof snapshot is inputs/PRODUCER_PROOF.md, SHA-256 13dd50a9a7abdc786394ae67b076cfc28d5cd934875421fd7c9af4002b36cb32. All reviewed producer artifacts and prior sources are copied under inputs/ with relative paths and exact SHA-256 bindings in AUDITED_INPUTS.json. These immutable snapshots make the audit portable independently of the sparse checkout.
