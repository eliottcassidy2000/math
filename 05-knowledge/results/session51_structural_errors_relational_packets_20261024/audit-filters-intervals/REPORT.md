# Independent dual audit: prime filters and implication intervals

Date: 24 October 2026. Status: **PASS in each lane, with the stated scope retained.** No mathematical correction or publication blocker was found. This audit did not publish or edit either producer's mathematical source.

## 1. Prime-offset lane: PASS

Audited packet: `copies/prime-offsets/`. Source manifest SHA-256:

`224415411c63d2a4045a65560b972d41a310c5e71ce09cbced0035c5f10cd746`

The authoritative actual hash is recorded in `input-integrity.json`; the audit manifest also pins the entire copied packet. The producer manifest covers ten declared payload files. The optional generation log is not a declared producer payload, but it is additionally pinned by this audit. Runtime bytecode is excluded.

### Proof audit

1. The indexing distinction is correct: F starts at n=1, L starts at n=2, the supplied lower prefix beginning with 7 omits L_2=3, and sorting/deleting duplicate F-values changes the indexed object.
2. If a prime endpoint exceeds p_n, a common prime divisor of d and P_n would be a proper endpoint factor. Thus d is coprime to P_n, d >= p_(n+1), and any composite d is at least p_(n+1)^2. This is necessary, not sufficient, for endpoint primality.
3. For n >= 3 the endpoint defining L_n exceeds p_n: Bertrand supplies a prime between p_n and 2p_n, strictly below P_n-1 since P_n >= 6p_n. L_2=3 is the genuine exception, and L_1 is undefined. This validates the finite membership cutoff without either Fortune conjecture.
4. The cutoff 71 is complete globally: only p_n < 71, hence n <= 19, can contribute regular values <= 71. Handling the isolated L_2 term separately closes the negative side. An initial prefix alone would not establish this, but the proved lower bound does.
5. Both endpoints at d=17 are prime and every d=2,...,16 is excluded. At d=289 both endpoints are again prime despite a composite offset. The d=19 one-sided and d=23 two-sided false positives refute the inverse wheel implications. Offset 1 is deliberately excluded.
6. Lucas's criterion is exact here: complete recursively prime factorization of N-1 and the order tests force order N-1 modulo each prime divisor of N. The proof does not depend on probable primality. Strong-Miller-Rabin failures certify compositeness, and passing tests are never used to prove primality.
7. The Ogg comparison keeps distinct predicates. Independent class-number and curve counts agree; the global Ogg list remains CITED, not proved by this finite run. Exact set intersections and complements follow from the proved offset cutoff.
8. Prime supports versus wheel supersets, existential versus conjunctive holes, and branch-conditioned full witnesses remain distinct. No native-routing or supersingular correspondence is inferred from numerical overlap.

### Independent arithmetic and certificate controls

`independent_check.py` does not import either producer checker. It validates all 200 Lucas certificates in increasing numerical order, checks all 1,734 endpoint classifications and first-success minimality, and independently reconstructs the 19 rows.

A separate census enumerates all 29,360 nonsingular short-Weierstrass models over primes 5 through 71, rather than selecting one curve per j. Direct square-multiplicity point counts produce exactly the same supersingular rational j-values. All twists agree on supersingularity. Primitive reduced forms are enumerated with a deliberately loose independent bound and agree with both the order-class-number expression and the published field-class-number expression. Characteristics 2 and 3 remain the separately cited elementary cases.

Thirteen new literal-file certificate mutations are rejected in separate normal and optimized subprocesses, including forged minimality, missing recursive factors, incorrect products, invalid bases, and two composites that pass Fermat tests (341 and 561). All 26 failures occur through active mathematical checks, not import or file errors. The producer's seven mutation controls also pass in both modes.

## 2. Implication-interval lane: PASS

Audited packet: `copies/implication-intervals/`. Final producer manifest SHA-256:

`aee23e3dab48b5bd931142c37d8ee50a65e72ac9ab4879e4e226c79e733ebda4`

Theorem SHA-256: `957767df1e3b6463da50815535614642be1c306317262e1b0f091445a71599ba`.

Receipt SHA-256: `58e115a3e412ed9fb363ef471638a5e71fe363b88c53702b0e4eac97ed64b9c6`.

### Proof audit

1. Every restricted binary relation on two differing coordinates accepts 00 and 11. Its only possible exclusions are the two mixed pairs, exactly the two implication directions. Fixed coordinates, unary relations, and repeated-coordinate relations introduce no unsupported extra restriction. SCCs synchronize, and requirement closure gives precisely ideals of the correctly oriented condensation poset.
2. Full acceptance additionally requires zero aggregate defect in the COMPLETE abelian guard, including exact integer margins. The zero-defect subfamily need not inherit the distributive lattice operations.
3. For accepted I subset J, intermediate ideals are exactly ideals of the induced poset on J minus I. Defect subtraction is valid because I is neutral. Covers are therefore relative neutral atoms, not arbitrary neutral subsets and not necessarily global first-event atoms.
4. Positive physical costs make cover refinement weakly improve both bottleneck coordinates. The simultaneous Pareto frontier, rather than separate minima, is the right budget invariant.
5. A route partitions SCC blocks by their first switch event. Each packet is neutral, and every interpacket dependency points toward an earlier packet, so a directed packet cycle is impossible. Conversely a requirement-first order of an acyclic neutral contraction gives accepted prefixes and no repeat switches. This proves both directions with actual packet costs.
6. Merging SCCs of a packet dependency graph retains neutrality and produces an acyclic condensation. The cost is the sum over the WHOLE merged packet-SCC. The L D(Q) s and L D(Q) w bounds follow only after exact-margin neutrality, finite-group extraction, and the actual packet-SCC size bound L are authenticated.
7. The chain construction telescopes to nonzero proper prefix charges in fixed C2 squared, although each original implication SCC is a singleton. Its accepted fibre has only the two endpoints for arbitrary n. The even-length population condition does not create new accepted states.
8. The paired lift strengthens the obstruction: original SCCs have size two and preserve every pair population exactly, while every physical charge increment is nonzero. Thus bounded margin-neutral synchronization blocks do not rescue the Davenport bound without packet-cycle control. The relation graph remains a tree.
9. Resource reservations, record non-return, exact changed-clock costs, continuing retired-resource history, and atomic-arrival semantics are explicitly scoped. Nothing authenticates native electrical charge labels or additional microscopic transitions. The audit preserves these limits.

### Independent finite universes

The independent checker uses literal relation tables, sets, and packet permutations rather than importing producer compilation or routing routines:

- Every combination of arbitrary binary tables on three coordinates, with every unary product restriction, yields 166 distinct solution families. All 2,352 accepted ordered endpoint pairs agree with their semantic implication ideals, including fixed coordinates and switch gauges.
- Every directed graph through four labelled vertices is checked, including cyclic graphs: counts 1, 4, 64, 4,096. SCC-ideal lifting agrees with literal implication satisfaction. Distinct labelled posets are counted as 1, 3, 19, 219.
- All those labelled posets, all total-zero C2-squared charge assignments, and declared exact integer-margin patterns give 42,681 instances and 87,902 neutral partitions. Relative intervals, cover refinement, packet-order equivalence, cycle compilation, accessibility, and conditional cost bounds pass.
- Physical costs satisfy clock total at least slot total. A separate actual tradeoff has dependencies 2=>1 and 4=>3 and block cost pairs (1,5),(1,7),(3,5),(3,7). Its exact Pareto frontier is {(6,12),(4,14)}. The coordinatewise minima (4,12) are NOT attainable: the tempting third neutral pairing has a cyclic contraction.
- Paired-lift hostiles are checked through logical length 16. Every legal prefix preserves total and all individual pair populations; only the two endpoints are charge-neutral, and no physical increment is zero.

These finite universes audit the proofs; they do not establish all-size results by extrapolation.

## 3. Actual-copy, source and mutation integrity

- Both producer checkers were run from real copied packet directories, normal and `-O`, and their output bytes agree with the respective frozen receipts.
- The independent checker itself is run in normal and optimized modes with byte-identical output. An additional whole-audit temporary-copy replay from an unrelated working directory passed both modes, all 26 independent certificate rejections, and standalone input verification; see PORTABILITY.json. There are 591,433 active independent mathematical checks.
- Every producer manifest entry is checked independently. The interval manifest additionally enforces its complete payload inventory and rejects source edits, inherited-input edits, missing theorem, extra ordinary file, and extra nested MANIFEST.json. The portable checker performs another temporary copy with an unrelated current working directory.
- All inherited source pins match the exact Git blobs, including the interval inputs' Git object SHA-1 values. Source and copy equality is recorded in `input-integrity.json`.
- The bounded inbox command reported missing `.machine-id`; no inbox registration or shared-checkout mutation was attempted. This does not affect mathematical verification. An intermediate shared-checkout inspection showed unrelated navigation edits; they were preserved. The final status inspection was clean on main relative to its cached origin/main.

## 4. Boundaries and release recommendation

Recommend both lanes for publication under their existing labels. Prime arithmetic is FINITE-EXACT with elementary cutoff/wheel proofs; Ogg's classification is CITED; both Fortune conjectures remain OPEN. Interval statements are abstract all-size proofs with finite-exact controls, not native current-bank supply or universal small-move theorems. No LRC closure, new literature priority, or conjecture resolution is supported.

The all-neutral case and identical endpoints carry the usual empty-route convention with zero work. This is already consistent with the inherited endpoint theorem; no nonempty packet need be invented for an empty interval.

No publication was performed by this audit worker.
