# Words, necklaces, and an all-order observer classification

**PROVED, independently audited twice.** For the classical five-letter construction of unique-Hamiltonian-cycle tournaments from the strong four-vertex seed, the native ordered-two-path observer is singular **exactly** when the insertion word contains `be` or `ee`. Every avoiding word gives a unimodular observer. More strongly, on an avoiding word **every ordered-k-path observer** is permutation-equivalent to a unit upper-triangular matrix.

This closes the arbitrary-order spectral criterion left open in the pinned source at `d3132eafbfa9c43c6d9f1110566b1e9c3a3286e1`, `05-knowledge/results/operation_closed_synthesis_20261029/THEORY.md`. It does not solve unrestricted tournament enumeration, nor classify the Smith forms of all singular observers. Independent mathematical reviews are [audit 1](audit/ALL_ORDER_SPECTRAL_PROOF_AUDIT.md) and [audit 2](audit2/REPORT.md). These are human-readable proof audits with exact computational checks, not proof-assistant certificates.

## 1. Definitions and exact grammar

A tournament has one oriented edge on every unordered pair of distinct vertices. M_k has rows indexed by ordered k-tuples of distinct starts and columns indexed by ordered k-tuples of distinct ends. Its entry counts spanning covers by k ordered, nonempty, vertex-disjoint directed paths with those starts and ends. Singleton paths are allowed. All statements concern these native integer matrices.

Start with vertices 0,1,2,3, directed Hamilton cycle 0→1→2→3→0, and chords 0→2, 1→3. At every stage the Hamilton cycle is unique. Type A has one maximum-outdegree vertex a. Type B has two, a0→a1. These maxima have indegree one. Write pred for predecessor on the unique cycle.

Each new vertex x beats every old vertex except its parent p, which beats x; insert x after p in the cycle.

| Type | Letter | Parent | New type |
|---|---|---|---|
| A | a | pred(a) | A |
| A | b | a | B |
| B | c | pred(a0) | A |
| B | d | a0 | B |
| B | e | a1 | B |

The seed is type B. A legal word of length m has order n=m+4. We use the classical bijectivity/completeness of this constructor, not a new classification claim: Michael R. Shelor, *On Tournaments Having A Unique Hamiltonian Circuit* (1977), §§2–5. The actual implementation is [actual_grammar.py](singularity/actual_grammar.py); the pinned finite class data are included verbatim for reproducibility.

## 2. Why avoiding words are unimodular, for every k

For the marked cycle with successor s, build a vertex digraph D_C: each noncycle tournament arc u→v contributes s(u)→v. The [acyclicity proof](triangularity/AVOIDING_WORD_TRIANGULARITY.md) shows that D_C is acyclic for every avoiding word, by tracking the exact outgoing-edge change under insertion.

Choose an integer potential f strictly increasing along D_C. For a spanning path cover, telescoping gives

    Σ_used arcs [f(v)−f(s(u))]
      = Σ_ends f(s(t_i)) − Σ_starts f(s_i) ≥ 0.

Equality holds exactly when every used edge belongs to the Hamilton cycle. Cutting that cycle immediately before the prescribed starts gives the unique equality cover. The induced row-to-column correspondence is a bijection. Ordering by the potential sum makes the matrix unit upper triangular. This works simultaneously for all 1≤k≤n, without a separate determinant calculation.

## 3. Why forbidden words have an explicit kernel

The [bridge–sink star lemma](bridge/BRIDGE_SINK_STAR_LEMMA.md) proves a general sufficient condition. Partition the vertices into R,Q. In R, a has sole in-neighbor s. The only R→Q edge is q→r, with q equal to a or s. If Q−r has a universal sink c and c→r, then

    M_2 ((e_a−e_s) ∧ e_r) = 0.

Here e_i∧e_j means the native ordered-pair vector e_(i,j)−e_(j,i), so this is an integer kernel certificate with no rational-to-integral identification.

At either a `be` or an `ee` birth, the required Q is a directed triangle. Choose the **last** forbidden pair in the word. Every remaining legal parent stays in Q and is different from c. Thus the same sink, bridge, and nonzero star survive the entire suffix. This proves singularity at every order. It does not claim that the whole kernel, or its dimension, is monotone under arbitrary insertions.

There is also an [intrinsic cut characterization](singularity/FORBIDDEN_WORDS_ARE_STRONG_BRIDGE_CUTS.md): forbidden words are exactly those admitting a partition into two strong subtournaments, each of size at least three, with exactly one edge against the prevailing cross direction. The size restriction matters: singleton blocks would make the assertion false.

There is a second intrinsic certificate: avoidance is equivalent to the Hamilton cycle being the **unique directed cycle cover**, or permanent of the adjacency matrix equal to one. Normalize the Hamilton-cycle matching to the diagonal; every other permutation term corresponds to a directed cycle in D_C. Conversely, a forbidden strong bridge partition gives the disjoint cycle cover of its two blocks. This elementary corollary is proved in §6 of the acyclicity note and was separately reviewed in audit 2.

## 4. The recurrence is now an exact observer count

The language transfer and generating function were already proved in the source. The new spectral identification makes them count unimodular M_2 classes at every order. Let u_n count these classes for n≥4. With m=n−4,

    Σ_(m≥0) u_(m+4) x^m = (1+x)/(1−2x−x²+x³),
    u_4=1, u_5=3, u_6=7,
    u_n=2u_(n−1)+u_(n−2)−u_(n−3)  (n≥7).

Thus the counts at orders 4 onward are

    1, 3, 7, 16, 36, 81, 182, 409, 919, 2065, 4640, 10426, 23427, ...

In particular **2065 at order 13** is now a theorem consequence, not an extrapolation from orders 10–12. The source's minimal four-state acceptance automaton and six-state legality-aware output machine now decide this observer predicate at all orders. A short decision state still does not reconstruct a tournament or retain every later observable.

The [word refinements](languages/ACTUAL_GRAMMAR_REFINEMENTS.md) also count forbidden occurrences, including overlaps. If t marks length and z marks occurrences, then

    Z(t,z)=(1+t−tz)/(1−2t−t²+t³−zt(1−t)²).

These counts do not assert that the number of forbidden pairs determines nullity or Smith form.

## 5. Construction words are not tournament necklaces

Rotating the legal insertion word `ecb` gives `cbe`; their two-path determinants are respectively −1 and 0. Even reversing `de` to `ed` changes the isomorphism class. Therefore a necklace quotient of insertion words destroys information.

A different [faithful necklace model](necklaces/FAITHFUL_CHORD_NECKLACES.md) encodes the chord rows along the intrinsic Hamilton cycle. Its constrained cyclic words classify cycle-marked tournaments, and classify ordinary tournaments when the cycle is unique. Converse is a twisted reflection, not bare reversal. The alphabet grows with order, so this is not a finite-alphabet shortcut for A000568.

This supplies a useful separation: construction histories support a small finite-state singularity decision, while faithful cyclic coordinates retain the full object. The same coarse grammar type and determinant can miss the next insertion response. The original order-five d/e pair exhibits precisely this: appending e yields de unimodular but ee singular.

## 6. Verification and provenance

The proof is all-order; finite tests corroborate it. The included checkers cover the pinned 4180 legal words through order 12 for the vertex-DAG criterion; all ordered-k observers on the small avoiding family; 861 explicit last-birth kernel certificates through order 11; independent arbitrary bridge–sink examples, including the two-vertex Q boundary; and an independent fixed-cycle necklace census through order 7. The two audits separately check matching bijectivity, signs, singleton paths, and the last-forbidden suffix guard.

Run `python3 verify_all.py` for the packaged replay; its receipt states precisely which tests ran. Earlier source statements describing the criterion as open are historical. The current status is the audited theorem above. General clone closure, arbitrary singular Smith forms, and unrestricted A000568 finite-state reconstruction are not consequences.
