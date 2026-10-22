# Separate addendum: all-order labelled reconstruction audit

Verdict: PASS. This extends, rather than silently edits, the earlier frozen response audit. The earlier statement that full-response injectivity was open is superseded by this new argument. General two-port injectivity remains unresolved.

## Argument independently checked

Let k≥2, choose distinct core labels u,v, and order the other k−2 core labels arbitrarily. Let R be the alternating path beginning at v, traversing all k−1 base vertices and all other core vertices except u, and ending at the last base vertex z. There are d=2k−3 distinct cross edges in R. Randomize these and the distinct cross edge {u,z}; leave all other arcs fixed. Give each R-edge sign +1 along the forward direction of R, −1 in reverse. Sum Hamiltonian counts against the product of these d signs over all 2^(d+1) assignments. The extra edge has no character sign.

Exchange the finite sums over assignments and vertex words. If a word omits an R edge, flipping that unused edge leaves its validity unchanged and reverses its character, so its total is zero. A word containing all undirected R edges must contain R as a contiguous block, in one of its two directions: internal vertices already have both incident word edges prescribed, and the path spans all vertices except u. This assertion uses only path combinatorics and does not assume a core Hamiltonian path after deletion.

Exactly four candidate words remain:

- uR: valid at the core endpoint exactly when u→v, has all R signs positive, and leaves {u,z} unused. Contribution 2Q_uv.
- reverse(R)u: valid at the core endpoint exactly when v→u, has all d signs negative, and leaves {u,z} unused. Contribution −2(1−Q_uv), because d is odd.
- Ru: one assignment of {u,z} is valid, with positive character. Contribution +1.
- u reverse(R): the opposite assignment of {u,z} is valid, with negative character. Contribution −1.

Thus the signed sum is exactly 2(2Q_uv−1), with no hidden normalization. If expressed as a normalized Fourier coefficient, divide by 4^(k−1). All unspecified attachment arcs and all base-base arcs are irrelevant to these four words, so arbitrary fixed choices cannot affect the identity. In particular the required ordered transitive base is permitted.

## Boundaries and interpretation

- k=2: R consists of v and one base vertex, d=1; the same four-word cancellation works. There are four response queries per ordered pair.
- k=0,1: only one labelled tournament exists, requiring no reconstruction queries.
- A single arc uses 2^(d+1)=4^(k−1) evaluations, each with exactly k−1 added base vertices. Querying one ordering for each unordered core pair yields a sufficient total binom(k,2)4^(k−1), without claiming no overlaps or optimality.
- Equality of response functions on fixed labelled masks implies equality of every labelled core arc. Consistent relabelling gives equivariance. Do not replace the labelled/marked response by a multiset or scalar average and retain this conclusion: the construction uses port identities and ordered masks.
- k−1 is sufficient, not minimal. For k>3 this does not resolve whether two ports alone suffice. The previously verified finite two-port census through six remains a separate result.
- No novelty or efficient bit-complexity claim is supported by this audit. Exact signed sums may involve large counts and exponentially many queries.

## Independent computational evidence

reconstruction_independent.py implements the contrast without importing producer code or any subset DP. Every Hamiltonian count is computed by full vertex permutation enumeration. It additionally accumulates contributions per word and checks that no word outside the four prescribed candidates survives.

- k=2: all labelled cores, all ordered distinct endpoints, all fixed-arc choices: 4 contrast tests.
- k=3: all labelled cores, all ordered distinct endpoints, all fixed-arc choices: 192 contrast tests.
- k=4: 20 deterministic samples, varying cores, ordered endpoints, remaining-core order and all other fixed attachment arcs.

All 216 tests pass. Normal and optimized runs agree. Producer faithfulness.py is separately rerun, including its 428 arc tests and 27,656 direct-DP response evaluations; those controls supplement rather than replace the independent permutation tests.

## Relabelling and converse

For a core-label permutation π, define (πQ)_(πu,πv)=Q_uv and apply π to every mask coordinate. Relabelling the core vertices is a bijection of Hamiltonian words, so F_(πQ)(πw)=F_Q(w). Faithfulness therefore identifies core isomorphism classes exactly with response orbits under this consistent permutation action, not with unordered collections of response values.

Reversal of every tournament arc reverses each Hamiltonian word. To restore the standard increasing base order, also reverse the base labels. Consequently F_(Q^op)(w)=F_Q(complement(reverse(w))). Both complement and word reversal are necessary in this convention. Equality with this transformed response is equivalent, by faithfulness, to labelled Q=Q^op, impossible for k≥2. Equality up to a consistent core relabelling instead characterizes self-converse cores; these two notions must not be conflated.

A further 160 independent permutation checks verify relabelling and converse formulas at k=1–4, recorded in equivariance-converse.json.

## Final source freeze and optional recognition

Reviewed final faithfulness-proof.md SHA256 e64ae1527002562379b8de71448eb355aa27217be71d3c32dba77ad79a5ebfe4, including its relabelling/converse quotient and optional rational-series recognition sections. All five entries of FAITHFULNESS-SHA256SUMS pass. The combined relabelling/converse orbit statement is sound because the two actions commute, and faithfulness applies after alignment. No metagraph-wiring or enumeration-speed claim follows.

The optional recognition algorithm is correct: extraction of ±2 signs gives a unique candidate, but acceptance requires exact equivalence to the candidate's slot series. The direct-sum realization has dimension at most n+D, yielding witness bound n+D−1. Mere passage of the sign tests would not suffice. The k=0,1 direct-comparison branch is necessary and is included. This is a proved decidability consequence, not an audited implementation of a new recognition program.
