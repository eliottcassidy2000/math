# Session117: Ramsey and circuit inheritance

Source snapshot: `session116_source` at `b435f833bc4c5b9d2c5c04f764c2160fabd6a792`. Bounded materialized-file pass; no claim of exhaustive repository or literature coverage. Computation scripts below are named by their sources but `04-computation` is not materialized in this sparse worktree. No scripts were executed in this inheritance pass.

## 1. Strongest exact reusable mechanisms

### Tournament doubling: exact mixed-chain predicate (PROVED)
`01-canon/theorems/THM-483-zigzag-law-refutes-the-sandwich.md` supersedes the sandwich speculation in `THM-455-trans-subtournament-doubling-and-tower.md`.

For skew doubling D, `trans(D(T))=z(T)`, where z is a shuffle of an ascending chain and a descending chain with the explicitly retained cross orientations. A chain in D contains at most one twin pair. Bounds are `trans(T)+1 <= z(T) <= 2 trans(T)`.

Hostile: alternating family A_l on 2l+1 vertices has `trans(A_l)=l+1` and `trans(D(A_l))>=2l+1`. Thus no constant increment bound exists. The source says exact doubled value 2l+1 only verified for l<=5; the remaining uniform upper gap is one vertex. **Do not repeat the source corollary's equality delta=l as an all-l proved result:** its own honesty section only proves delta>=l generally.

Source/target/map: T -> D(T), with transitivity preserved through a mixed-chain certificate; scalar trans(T) loses the shuffle compatibility. Sidecar is z(T), or a certifying compatible shuffle, not only the two chain lengths.

Actual finite next datum in the source: tower `z(T63)` / next-level largest transitive chain. Distinguish full doubled order126 from dephased tower order127 before reusing its shorthand. Existing tower values 3,5,7,11 at orders7,15,31,63 are computational source claims, not asymptotic control.

### Graph doubling: exact split-subgraph predicate (PROVED in note)
`05-knowledge/results/collatz_catalan_ramsey_20260927.md`, Theorem4:

- `omega(D(G))=cs(G)`, maximum order of induced complete split graph K_a join independent_b.
- `alpha(D(G))=max(alpha(G)+1,sep(G))`, maximum separated independent/clique structure, with twin-overlap case separately retained.

Source states exhaustive verification through n6, random n7–8. Its graph +2 increment also fails (reported n9). Scalar `(omega,alpha)` is not sufficient; cross adjacency is precisely what gets lost. Paley doubled graphs give poor Ramsey bounds, not records. Open D18 asks closure under iterated doubling. This is the cleanest algebra-to-combinatorics bridge available here.

### Book Ramsey: exact finite witness target (historical control verified; target OPEN)
`05-knowledge/results/ramsey_book_offdiag_macmini_S171.out` and `ramsey_book_n12_witness_macmini_S171.sets`.

Graph vertices Z_(2n-1) x {0,1}, within-level symmetric U0,U1, arbitrary cross D. For a lower bound `R(B_(n-1),B_n)>=4n-1`, seek graph on4n-2 vertices with every edge common-neighbor count <=n-2 and every nonedge common-nonneighbor count <=n-1. These are directly checkable complete certificates.

Control n12,N46 has reported maxima10,11 and witness sets. Open n100,N398 historically reached nonzero class energy182–184; this is NOT a witness or improved lower bound. Source claims best size seed (94,104), O(m) updates, and bad QR behavior; independently derive correlation identities before treating the QR claim as a proved obstruction. Its prose budget expression is visibly unfinished.

Named scripts: `04-computation/ramsey_book_offdiag_anneal_macmini_S171.c`, `04-computation/ramsey_book_offdiag_verify_macmini_S171.py`. No current running-search status inferred from the old report.

## 2. Circuit inheritance: boundaries rather than a gate breakthrough

`05-knowledge/hypotheses/HYP-3116-lrc14-circuit-lower-bound-missing-input-ledger.md` explicitly describes proof-payload and missing-input discipline, not a Boolean gate lower bound. Exact gap/endpoint-cover functionals, monotone proof DAGs and finite fitted literals cannot transfer automatically to unrestricted fan-in-two circuit size. The ledger itself warns that a literal fitting a finite bank lacks a uniform family theorem.

A literal NP family suitable for a new correctly typed experiment is `TT(n,k)(adjacency bits)=1 iff tournament contains an ordered transitive k-set`: certificate is k distinct vertices and pairwise orientations, verifiable in polynomial time. Similarly `CLIQUE(n,k)` has a standard vertex-set certificate. Naming these families is not a lower bound. Parameters, gate basis (e.g. arbitrary binary Boolean gates), free constants/negations, number of outputs, uniformity, and input encoding must be fixed before comparing bounds. Counting almost-all truth tables, dimension of a representation, communication/state lower bounds, monotone bounds and formula bounds do not prove unrestricted Boolean gate bounds.

`05-knowledge/results/boolean_circuit_info_s20am.out` is unsafe historical speculation: it prints sum_(d=1)^(n-2)(n-1-d)2^d=22 at n5, then still claims formula2^(n-1)-2=14; the correct sum is2^n-2n. It reports 64 inputs mapping to7 outputs and entropy loss3.4874 bits, then claims a perfect/no-information-loss channel. Neither low Walsh energy nor low polynomial degree establishes circuit depth/gate lower bounds. Do not promote this file's informal circuit conclusions.

`05-knowledge/results/ramsey_matroid_23.out` is another unsafe historical numerical analogy: it lists R(5,5)=43 and later says exact value unknown. Use none of its classical bounds without current primary-source verification.

Conditional separate branch: `01-canon/theorems/THM-4559-algebraic-spherical-sets-are-ramsey-conditional-on-the-tensor-classification.md` explicitly depends on unrefereed classification #172. Its separability-idempotent deduction is a precise algebraic bridge but remains conditional, concerns Euclidean Ramsey configurations, and supplies no Boolean gate lower bound.

## 3. Five cheap decisive tests

1. **Alternating-family one-vertex gap.** Classify all maximum transitive subsets of A_l, then ask whether two such subsets can form a legal doubled shuffle of size2l+2. Enumerate l<=6 as controls; seek uniform contradiction or a concrete chain. Stop when exact formula or smallest hostile obtained. Keep twin condition.
2. **Split-state sufficiency hostile.** Enumerate graphs n<=5 and bucket by `(omega,alpha,cs,sep)`. Compute doubled cs/sep. A collision disproves four-scalar recurrence; retain the smallest pair and identify cross-incidence sidecar. If no collision, extend only until a structural proof/counterexample is available. This is algorithmic-state classification, not circuit complexity.
3. **Independent book certificate audit.** Parse n12 witness and count all ordered pair codegrees directly, without cyclic-correlation formulas. Compare all class counts with a separately derived cyclic evaluator. This validates exact incremental update signs and produces a reusable checker before optimizing n100.
4. **Book budget/QR hostile.** For small odd m, exhaustive symmetric U and arbitrary D; derive per-class sum constraints and compare minima with reported slack formula. Separate a necessary average bound from an impossible cap assertion: mean near98.5 alone does not contradict mixed98/99 caps. Stop with exact obstruction conditions or retraction of overbroad QR restriction.
5. **Literal small-circuit certificate.** Fix TT(n,k) or another explicit NP family, exact gate basis and input count; enumerate small truth tables and SAT-synthesize circuits. Upper witnesses require independent truth-table checks; lower bounds require UNSAT proof checked separately. Use finite results only to test a proposed measure's subadditivity under every gate/composition, not to extrapolate all-n lower bounds.

## 4. Synthesis priority

Prefer exact mixed-structure states and proof-carrying pruning over scalar analogy. The central retained datum across the Ramsey doubling and owner-window themes is compatibility across blocks: symmetry reduces equivalent branches only if stabilizers, boundary adjacency and certificate reconstruction are retained. This is a concrete route to better exact algorithms. A transfer to unrestricted NP circuit lower bounds remains an additional unproved theorem, not an interpretation of these results.
