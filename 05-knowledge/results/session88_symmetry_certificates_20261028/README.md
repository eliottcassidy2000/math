# Symmetry through operations, reconstruction, and certificates

Session 88, 2026-10-28. This session revisits old mechanisms across enumeration, number theory, graceful labels, conditional expectation, automata, and marked gluing. Its aim is a coherent theory of **what an exact representation must preserve**, together with algorithms that stop solving a larger problem than the requested output requires.

Classical results are credited, including Fulkerson–Gross reconstruction, the Hunter spanning-tree bound, Smith theory, conditional expectation and two-projection theory. The contributions are the explicit connections, scoped derivations, implementations and hostile tests below. No universal equivalence between the source problems is asserted.

## 1. Compute the equivalence block, not its same row repeatedly

The user's proposed change works. A complete **named** root-switch row identifies its entire ordinary tournament block. Compute it once, fill the block, and transport the verified maps to its other pointed states.

At order eight, the implementation retains exactly the same full Q, ordered tournament bank, parent maps, every root's complete transporter set, and full automorphism permutations:

| Work/output | Full-row reference | Block filling |
|---|---:|---:|
| Rows computed | 54,256 | 6,880 |
| Root classifications | 434,048 | 55,040 |
| All requests, including bank preparation | 493,328 | 114,320 |
| Cached classification searches | 48,326 | 11,074 |
| Matched local median time | 12.255 s | 3.323 s |

The time comparison uses the same output contract and three alternating runs. It is a finite 3.69x observation, not an asymptotic claim. A subgroup-propagation variant reduces root queries to47,391 but did not improve runtime in its initial illustrative pass. Knowing fewer answers to compute is not automatically cheaper than managing the inference.

A bounded [stabilizer-first seed experiment](blockfill/seed_order/STABILIZER_FIRST.md) saves only14 more queries at eight, with one block becoming worse. Its explicit order15 automorphism-group example shows why: recovering the full group can still leave a named root type undiscovered. This is a useful failed implication, not a large extra speedup.

An immutable factored representation stores54,256 weighted entries instead of429,936 explicit sparse Q entries and applies Q to a vector with linear arithmetic work in the pointed-state count. Exporting every explicit matrix entry still incurs its output cost.

The dictionary is essential: two vertex-transitive tournaments on seven vertices have the same anonymous block `[7]`, but Hamilton-path counts175/189, automorphism orders7/21 and different ordered-pair orbits. A numeric row shape is not an object identifier.

Read [the exact block-fill construction and proof](blockfill/BLOCK_FILLING.md) and [reproduction instructions](blockfill/README.md).

Adjacent work reports a191,536-member order-nine bank using compute-once rows, and a separate order-nine root-selector improvement. Those outputs omit parts of this session's full Q/transporter contract; their timings cannot be multiplied or substituted for ours. The [concurrent review](concurrent/SOURCE_STATUS.md) records those limits and independently replays small marker/frame witnesses rather than claiming a second order-nine run.

## 2. A certificate can replace repeated discovery without replacing the source

One explicit deletion-to-parent bijection for **each root** certifies a complete named row. The verifier checks those maps against the original adjacency data instead of performing canonical searches. The row, one known parent-mask stabilizer, and its certified root transversals also reconstruct the full child automorphism group.

The general orbit version is precise: a generator-labelled transport tree proves every listed state is reachable, while closure under every action generator proves no orbit member is missing. Both are required. An incomplete list or a quotient that identifies distinct states proves less.

For integer linear equalities of rank r, failure has a different kind of short certificate:

- rational inconsistency uses a row dependence supported on at most r+1 source equations;
- rationally consistent but integrally impossible data admit a modular character supported on at most r source equations.

The character must annihilate **all** original columns. A selected-minor shortcut wrongly rejects `[2,3]x=1`; `4x=2` shows that a prime-power modulus can be necessary. Source identities and boundary signs make these certificates composable at an LCA without copying or cancelling the wrong constraints.

The arithmetic applications keep their quantifiers: a path obstruction rejects one proposed graceful signed assignment, not all graceful labelings. An authenticated prime-power factorization and a square or adjacent-square witness decide square-σ for a supplied integer; this does not make construction of the factorization free.

Read [proof-carrying symmetry and decision witnesses](certificates/CERTIFICATES.md) and [source-faithful arithmetic dual certificates](arithmetic/DECISION_CERTIFICATES.md).

## 3. Recover a joint pattern from overlap data instead of expanding Boolean space

An old graceful-gap mechanism records intervals by their Gram matrix G of intersection sizes. We replace the displayed Boolean inversion by an exact geometric query.

For a required set of interval names A and forbidden set B, first intersect the required intervals, then clip each forbidden interval to that intersection. Their lengths and pairwise overlaps are computable from G. For intervals,

    size(union) = sum(lengths) − maximum-spanning-tree overlap weight.

Thus one partial membership-pattern count costs O(m²) integer arithmetic/comparison operations. A positive-prefix trie finds all s realized full patterns in O(s m³), and s≤min(N,2m+1), giving O(m⁴) arithmetic work without visiting absent patterns. Bit costs still depend on log N.

The stronger uniqueness mechanism is **classical Fulkerson–Gross, Theorem2.1 (1965)**: if G has an interval realization, every binary-set realization with the same ambient size has the same column multiset. A short list of interval endpoints can therefore certify the promise. Zero-slack set inequalities explain why the missing higher correlations are forced, rather than assumed independent.

This is a useful counterweight to the old unbounded-memory examples: special structure can make pair data determine joint data. It still does not recover named gap addresses or the original column order. Arbitrary sets or circular arcs do not inherit the interval formula.

Read [the MaxST query, sparse reconstruction and classical rigidity connection](certificates/INTERVAL_PATTERNS.md). Independent tests cover106,392 interval families and over three million exact-pattern queries, plus arbitrary-set realizations of interval-realizable Gram matrices.

## 4. Abundancy becomes an angle in an explicit divisor model

This is the most concrete new cross-topic model in the session. Let N=AB with A,B>1 and gcd(A,B)=1. Give each divisor its numerical value as weight. The full divisor set is a product of the two coordinate divisor sets; restricting to proper divisors removes its top corner `(A,B)`.

Record only whether neither coordinate is top, only A is top, or only B is top. The three weights are

    x=s(A)s(B),  y=A s(B),  z=s(A)B.

The two coordinate conditional averages are projections P,R. On the nonconstant two-dimensional flag space, their principal angle satisfies

    cos²(theta)=N/σ(N),       tan²(theta)=s(N)/N.

It is independent of the chosen unitary split. Perfect numbers correspond to45 degrees; abundance and deficiency lie on opposite sides. This is a declared weighted resampling model, **not** an identification with aliquot iteration or a method for finding perfect numbers.

Let E be global averaging and S=PRP. On the three-flag space,

    E = [σ(N) S − N P]/s(N).

On the full divisor-coordinate space, additional zero modes require

    E = [σ(N) S² − N S]/s(N).

N=12 gives an explicit failure of using the flag formula on the full space. The noncoprime split36=6·6 overcounts divisor sources and produces the wrong mass and angle, so the unitary guard is essential.

Two integers with equal abundancy give isometrically equivalent **flag projection pairs**, but this is not a permutation of named divisor states or an integral-lattice equivalence. For example, at N=30 the splits2·15,3·10,5·6 give weights `(9,18,15)`, `(8,24,10)`, `(6,30,6)`. They share the angle, while their minimum averaging denominators are14,21,7. The real spectral quotient has discarded meaningful arithmetic structure.

The constants-removed projections also give a finite two-generator Temperley–Lieb representation with parameter sqrt(σ(N)/N). This is only the stated finite algebra relation, not a Jones tower, subfactor index, or physical model.

Read [the weighted transfer theory and divisor construction](operators/TRANSFER.md).

## 5. Changing the operation language can remove a real obstruction

For two positive weighted conditional expectations, a finite word of the two native operations reaches the common averaging projection exactly iff the two projections commute. The joint contingency table decides commutation; marginal totals do not.

For N=21 the proper-divisor weights are1,3,7. Local averages use denominators4 and8, while the global average needs11. No finite native word can introduce that denominator. Nevertheless the exact signed expression is

    E=(32 PRP−21P)/11.

For N=6 no new denominator prime is needed, but the operators still do not commute. The signed expression `E=2PRP−P` succeeds where native positive words do not. The difference is subtraction and the permitted coefficient ring, not a contradiction.

A characteristic-three hostile goes further: all normalizations can exist, the projections can fail to commute, and yet three operations reach the common projection by isotropic cancellation. The positive-real theorem must not be transferred unaltered to finite fields.

This is the session's precise interpretation of “redefinition is a weapon”: change the stated operation language and prove what becomes possible. Merely renaming an object does not remove its information debt.

## 6. The older connections that make the synthesis more than one example

The [seven-bridge archive map](archive/ARCHIVE_SYMMETRY_RESEARCH_MAP.md) recovers and checks source statements from March–August, including:

1. Burnside cycle types and weighted divisor-incidence algebra.
2. LRC collision bounds as conditional-variance/Dirichlet loss.
3. Marked recursive tilings as a span with a phase torsor, not a quotient graph.
4. Walsh/species/Krylov states for selected future means.
5. Gaussian prime toggles and source-dependent Berggren/LCA transport.
6. Mahler safe/reset states, finite certificates and infinitely distinguishable followers.
7. Eisenstein discriminant/divisor gates reconstructing fixed taxicab solution orbits.

One derived theorem is especially sharp: a positive weighted incidence-Gram response vector recovers the full multiplicity state and preserves compatible symmetry orbits and stabilizers. For odd permutation cores, the entire one-step add-cycle response spectrum determines the core. But `(3,11)` and `(5,9)` have identical current mass, length and edge-orbit exponent; adding a5-cycle separates them immediately. A terminal fingerprint is not necessarily a future state.

The old mechanisms remain subject to their original domains and corrections. Source recovery is distinguished from new independent reproduction. No LRC, Collatz, graceful-tree, or prime-atlas conjecture is solved by the shared language.

The concurrent review also makes the pair-versus-joint distinction observable: two P7 triples agree on every ordered pair-orbit key, but their explicit ten-vertex marker attachments have Hamilton-path counts4,281 and4,315. Relative frames on a *fixed free orbit* can instead be carried by tree gains. Left-gauge and right-gauge conventions require different multiplication orders because the21-element Paley group is nonabelian. Neither fixed-frame compression nor the bounded complement/reset tests establishes a universal closed attachment carrier.

## What is unified, and what remains open

Three interfaces recur:

- **An exact continuation state** must support the declared future operations, including names, source weights and coefficient ring.
- **An exact counting state** may merge histories when every allowed future contribution agrees; it need not reconstruct each history.
- **A decision certificate** is checked against the original input and may be far smaller than either state. Its construction cost is separate from verification.

The proofs and implementations show where these interfaces coincide, where geometry forces hidden correlations, and where a missing phase, source bit, denominator or transporter separates them. This is a reusable framework for symmetry, not a claim that all mathematical structures are isomorphic.

The smallest closed attachment carrier remains operation-relative. An unrestricted prime-class atlas is still open. Current class-count formulas, compact fixed-frame laws and decision certificates do not supply missing representatives or prove all-order closure under richer operations.

## Reproduction and scope

Run the scripts in each linked subdirectory; generated JSON receipts are local outputs, not required committed audit data. New code uses Python's standard library and explicit exception checks, so `python -O` does not remove the tests. The block-fill loader checks the exact inherited Session87 source hash.

Independent audits include interval reconstruction, small full-output tournament banks, all root/group transports in those universes, weighted partition pairs,35,308 small integer systems, positive incidence-state inversion, and108 additional unitary divisor splits. Producer and independent-auditor coverage are stated separately in the notes. Local timings are not universal performance promises; no CI pass is inferred from local verification.
