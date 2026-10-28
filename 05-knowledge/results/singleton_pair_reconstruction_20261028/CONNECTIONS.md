# Information, older enumeration, and two changed conjecture landscapes

2026-10-28. This is a targeted continuation, not a claim to have read the full
archive. Proven connections below have an explicit map; proposals remain OPEN.

## Recover before extending

The [previous history audit](../strong_attachment_tail_compiler_20261028/HISTORY.md)
and [incoming Session87 recovery](../session87_four_vertices_enumeration_20261028/HISTORY_AND_TRANSFERS.md)
route March's odd-partition/GMP counters, April's CRT and rounding failures,
June's unit-tail kernel, and October's future-response DP. The new350-prefix
uses the last of these, not a rediscovered partition enumerator. THM-514
[`the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md`](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md)
and THM-516
[`burnside-core-kernel-phi-reframe.md`](../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md)
remain the structural ancestors. Their small outer kernels do not prove a
polynomial construction bound.

October's [XOR disagreement corollary](../converse_tournament_count_channels_20261021/session22-converse-review/COROLLARY_REVIEW.md)
deliberately did not claim a doubled-order structural map. Recovering Eplett's
1979 odd-order representation now supplies that missing classical bridge for
the consecutive-pair interpretation. The even-order version lacks its unique
fixed pole; we have not extended the odd construction to it.

The [new weighted-fiber transfer account](../session88_symmetry_certificates_20261028/operators/TRANSFER.md)
explains why reset and root switching carry actual mark multiplicities.
The [incoming refinement](../root_ports_closed_carriers_20261028/CARRIERS.md)
and [full block transport](../session88_symmetry_certificates_20261028/blockfill/BLOCK_FILLING.md)
are both credited in the lazy implementation. The improvement here is a
smaller stored full-response representation through9; it is not a claim that
we invented individualization/refinement or improved its worst-case bound.

## An exact bridge to codes, and the limit of that bridge

For any fixed deterministic singleton/complement context w, its action on the
original free bits x is affine over F2. Blue is equality between reflected
edge coordinates, so its pullback is a system A_w x=b_w. This connects the
blue observations to binary linear constraints without forcing their
dependence into an ordinary graph cycle matroid.

But the actual blue observer returns only the **zero-test**
1[A_w x=b_w], not the syndrome A_w x-b_w. Equality of all permitted zero-tests
does not mean that a difference lies in every linear kernel. In the seven-vertex
control the syndromes after one left singleton are distinct and nonzero, while both
blue flags are false. `cross_domain.py` records the literal masks. Thus the
Hamming-code/square-class lane's full linear-value observer and the present
conjunction observer are different targets. This is a precise version of the
Vitali warning: which conjunctions can later be selected matters. Here the
restricted geometry cannot expose the hidden distinction; named guards can.

Source: free edge bits. Target: affine blue systems. Map: substitute the affine
context into reflection equations. Preserved predicate: blue. Lost information:
which equation fails. Sidecar for a richer target: the full syndrome. Cheapest
hostile: masks1072/4176 after J_L. Their empty-context syndromes are both5216;
the initial proposed empty-context witness was rejected and repaired before
publication. The distinguishing syndromes after J_L are1820815 and1982543.

The LCA/Frobenius comparison has a different source and map: actual group
frames, transported by products along a tree. It recovers relative phases
only because compatible witnesses were retained. The new pole-normalized
Paley control keeps a group of order3 inside the full21, so even a bijection
of isomorphism classes need not preserve full symmetry groups.

## ETH: separate the announced result from the transfer proved here

**CITED, checked 2026-10-28:** the [October26 primary announcement](https://openai.com/index/building-scientific-understanding-with-academic-communities/)
reports an ETH refutation and a Lean verification, with an external working
group examining the result. This session has not independently audited the
proof or built its formalization. The released
[October25 manuscript, Theorem1.1](https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/main_full_version.tex)
states the stronger exact fixed-k #CNF bound, for v variables and m clauses,

    (v+m)^{O_k(1)} 2^{O_k(v log log v / sqrt(log v))}

with polynomial space. The following consequence is **CONDITIONAL on that
stronger counting theorem**, not a consequence of the logical statement
"ETH is false" alone.

**Orbit-aware CNF transfer — CONDITIONAL, reduction PROVED.** Suppose a
relabeling-invariant tournament property has a fixed-width CNF on its
binom(n,2) arc variables with polynomially many clauses and no auxiliary
assignment multiplicity. For each odd permutation type, orient one variable
per unordered-pair orbit; every other arc is that variable or its negation.
Substitution preserves bounded clause width and is a bijection between these
assignments and fixed tournaments. Even permutation types contribute zero.
Apply exact #CNF to each resulting instance and sum by Burnside with weights
1/z_lambda. Conjugacy invariance justifies one representative per type.

There are exp(O(sqrt(n))) partitions, at most binom(n,2) variables per instance,
and only O(n²) output bits. Thus the conditional running time is2^{o(n²)}
times a polynomial, with polynomial working space. For example, forbidding
each embedding of a fixed finite orientation supplies bounded-width clauses.
If auxiliary variables are used, their number and unique extension property
must be proved; a generic polynomial reduction does not preserve this bound.

For unrestricted A000568 the fixed count is already the closed expression
2^e. The classical odd-partition scalar algorithm runs in
2^{O(sqrt(n))} times polynomial bit work, so invoking #SAT would be weaker.
The opportunity is **additional arithmetic or attachment guards**, not an
automatic speedup for the bare sequence. Explicitly listing all classes still
needs at least a(n)>=2^{binom(n,2)}/n! outputs. Neither ETH's failure nor faster
SAT removes this unconditional output requirement or restores erased alignment.
No conclusion P=NP is drawn.

## The three-dimensional Jacobian counterexample as an owner-loss control

**PROVED in existing canon / exact certificate rechecked:** THM-1300
[`jacobian-counterexample-dixmier-A3-explicit.md`](../../../01-canon/theorems/THM-1300-jacobian-counterexample-dixmier-A3-explicit.md)
contains the polynomial map with determinant -2 and three distinct rational
points mapping to(-1/4,0,0). `cross_domain.py` recomputes both facts symbolically.
This is independent verification of the inherited Claude-discovered example,
not a discovery claim by this session. The repo's current status is false in
dimension3 and above, with dimension2 still open.

[Tao's July21 explanation](https://terrytao.wordpress.com/2026/07/21/a-digestion-of-the-jacobian-conjecture-counterexample/)
relates the construction to multiplying a designated linear factor by a
quadratic factor. Forgetting the factor designation can leave several global
choices despite local invertibility. That is a useful comparison to forgetting
a root or alignment. It is not an isomorphism between the two theories, and
the comparison does not give a planar Jacobian counterexample.

The operational distinction is exact: local certificates, complete global
fibers, and a chosen point in a fiber answer different questions. Our paired
tournament construction retains an alignment; its lazy root-switch carrier
retains an actual transport; its singleton boundary profile drops only data
that all declared future operations provably fail to test.

## Four OPEN targets with decisive next tests

1. **Is K minimal?** At7, K has32272 cells; the horizon7 quotient has24754.
   Classify the remaining gap using the explicit closed K transition rules.
   The exhibited collision is permanently merged; no larger horizon can
   reconstruct every presentation. Also test whether Psi is already determined
   by the two endpoint-rooted types and b beyond7, as happens in the full7-bank.
2. **The weakest menu that restores reconstruction.** Add one selected interior
   port, or one degree-conditioned fresh re-marking, and test the seven-vertex
   witness first. Prove an all-size separation theorem only after specifying
   the selector. The named-path theorem is an upper benchmark, not the answer
   for every intermediate menu.
3. **Even doubles with an honest gauge.** Replace the odd fixed pole by an
   explicitly retained antipodal owner and switching action. Determine the
   smallest sidecar giving both a class bijection and transporter recovery.
   Test n=2,3 against the inherited F_n(2,0) counts before a general claim.
4. **Count constrained atoms without listing them.** Choose a fixed forbidden
   orientation or bounded joint-port guard; compile its fixed-permutation CNF
   parsimoniously and compare small exact counts to the actual banks. This
   is a concrete testbed for the conditional counting transfer; practical
   superiority to the specialized enumeration remains unproved.

The unifying proposal is a program of **operation-specific exact interfaces**:
separate the relation, its coherent transport, and the predicates realizable
by future operations. S2--S4 now give a nontrivial all-size theorem for one
such interface. They do not establish a universal compression theorem for
arithmetic, LRC(14), Collatz, or all joint-port attachments.
