# 101, 35, 175: continuation spectra, squarefree clocks, and symmetry

**Status:** independently audited elementary proofs and exact reproducible computations. Collatz convergence, exclusion of all nontrivial integer cycles, and unrestricted tournament reconstruction remain **OPEN**. This packet makes no literature-first claim.

Source snapshot: `e871dd4513c7dcaa75cb29bb9847aafac35a795c`. The triggering example was already proved in [the preceding cycle-transport packet](../collatz_fermat_reframes_20261029/README.md), C3. Here we recover its older ancestry and strengthen the surrounding mechanisms.

## What changed

1. **The same denominator and the same root set still hide different clocks.** The repeated words115 and133 both have slope27/128, denominator101, and exactly the same phase-zero residues modulo101². Nevertheless, the sharp longest realizable all-squarefree repetitions are **4985 and8600** respectively. The first bound is inherited; the second is proved here, including all-prime avoidance and a fixed-word squarefree sieve.
2. **The35/175 example extends to all appended valuations.** For115 the new numerator is185=5·37; the two factors can cancel together. For133 it is221=13·17; each factor can cancel separately, but they can **never cancel together** under a one-letter append, because their required exponent phases disagree modulo4. Both cancellation sequences and reduced-denominator sequences have explicit rational generating functions and linear recurrences.
3. **A single unmarked cycle already exposes the attachment problem.** The rotations133,331,313 all describe the same rational cycle of denominator101. Appending the same letter1 gives denominators **175,175,7**. Forgetting the attachment position cannot commute with this operation.
4. **A compact exact coordinate is not free compression.** Given multiplier, word length and total valuation, the exact carry reconstructs every letter by successive2-adic valuation. The one retained integer encodes the full marked word.
5. **There is an intrinsic tournament boundary.** Pairwise affine swap comparisons form a total preorder, so their strict quotient is always a transitive tournament. The entire weighted elementary comparison matrix is blind to the multiplier; a common suffix reveals it. This identifies an exact useful observer and proves what it cannot represent.
6. **The same175 lives in different integer lattices.** The focal triple101,35,175 is a unimodular basis of its prime-unit group. The older actual27-orbit triple91,175,325 is rationally independent but has Smith invariants(1,1,4), and its full product is a square. Rational independence, binary dependence, and integral saturation answer different questions.

## Read the results

- [Historical map](history/PROVENANCE_MAP.md): six prior101 campaigns, twelve35 campaigns, and twelve175 campaigns, conservatively grouped and with exact source scopes. [Pinned source records](history/SOURCE_BLOBS.json) identify the text actually read. Counts measure this bounded recovery pass, not an exhaustive corpus census.
- [All-clock cancellation and denominator generating functions](spectra/PERIODIC_CANCELLATION_SPECTRA.md).
- [Multiword continuation, prime-adic centres, and carry decoding](continuations/CONTINUATION_STATE_AND_SYMMETRY.md).
- [Squarefree return clocks and the new8600 boundary](squarefree/README.md).
- [Intrinsic swap order and false composition symmetries](symmetry/INTRINSIC_SWAP_ORDER_AND_SYMMETRY_LOSS.md).
- [Integer lattices versus square classes](lattices/VALUATION_LATTICE_AND_SQUARE_CLASSES.md).
- [Direct tournament and addressed-face controls](tournaments/TOURNAMENT_AND_FACE_CONTROLS.md).
- [Independent audit](audit/AUDIT.md).

The older map does more than collect matching numbers. It recovers the September101 prime-square obstruction, an October two-edit word with raw gap175 and reduced denominator35, the rotational tournament with175 Hamilton paths, addressed35-dimensional Kneser layers, genuine35 ROOT-transport failures, and the prime-incidence triangle in27's orbit. Each connection names its operation and what the representation loses.

## The exact anchor

For a marked valuation word u at multiplier3, let k be length, H total valuation and C its ordered carry. Its periodic source is C/(2^H-3^k).

| word | k,H | carry C | periodic source | append1 and reclose |
|---|---|---:|---|---|
|115|3,7|19|19/101|185/175=37/35|
|133|3,7|31|31/101|221/175|
|313, a rotation of133|3,7|49|49/101|275/175=11/7|

“Same clocks” here means the aggregate pair(k,H), not the same ordered valuation word. An actual step on any of these original rational cycles preserves101. Appending and solving a different fixed-point equation is a different operation.

For any appended v>=1, the common gap is Delta_v=2^(7+v)-81. The exact cancellation gcds are

    g_115(v)=1+4·[v=1 mod4]+180·[v=25 mod36],
    g_133(v)=1+12·[v=9 mod12]+16·[v=7 mod8].

The13 and17 clocks are incompatible because9≠7 mod4. This is a concrete instance where two individually possible conditions have an impossible conjunction. No independence assumption is allowed.

## A unifying formulation: composition must preserve responses, not just objects

The general principle is classical in transition systems and automata: given a monoid M acting on the state set (including the identity operation) and an observable f, two states are operationally equivalent exactly when f(mx)=f(my) for every allowed continuation m. This equivalence is stable under continuation and is the coarsest stable equivalence refining the original observation. It is not claimed as a new abstract theorem here.

The new usefulness is in computing its concrete ingredients:

- For a finite menu of arithmetic suffixes, full numerator responses are exactly residues modulo the lcm of their closing gaps. Denominator-only responses are coarser prime-adic distance profiles, and no minimality claim is made for that coarser observer.
- For one-letter suffixes at all clocks, the continuation response reduces to divisor-labelled cyclic subgroup phases. Joint feasibility is a compatible fiber product over the shared exponent clock. For any finite collection of these congruence channels, pairwise phase compatibility already implies joint compatibility by generalized CRT. This is a genuine compression theorem for this clock model, not a general claim that pairwise tournament-port data suffice.
- For squarefree repetitions, the carrier is the affine return together with intermediate zero-root times and the endpoint convention. The unordered roots alone fail in the115/133 example.
- For tournament attachments, marked sites and a common transporter play the corresponding role. The historical81/27 and joint run-order controls remain hostile examples; this packet does not assert an isomorphism from them to Collatz.

This is more specific than saying “retain more information.” It tells which coordinate restores each operation, provides exact finite response partitions for fixed clocks and a finite operation menu (not necessarily stable under arbitrary further composition), and proves that some scalar encodings secretly retain the entire source word.

## What this does and does not offer a convergence argument

The arithmetic machinery separates operations and gives exact tests. It does not produce a well-founded descent for arbitrary changing words. The squarefree bound applies to fixed repeated blocks with every displayed node squarefree. An arbitrary integer orbit can change blocks or contain square factors. The periodic cancellation generating functions vary a final letter of a fixed prefix; they do not enumerate all Collatz trajectories.

A useful next target is therefore not another unqualified denominator invariant. It is a source-preserving transition rule between distinct block types with an exact carry update and a proved decreasing quantity. It must survive the same-cycle175/175/7 example, the incompatible13/17 clock phases, and the nonunit-carry prime-square hostile(1,13). Those are explicit obligations, not discharged claims.

## Reproduction

All mathematical checkers use Python's standard library. From this packet directory:

    python -B verify_all.py

The runner executes all nine scripts normally and with `-O`, compares stdout byte-for-byte, and compares parsed output to the stored receipts. [The verification record](verification.json) lists the tested scripts and hashes. Independent scripts recompute the denominator, squarefree-clock, and lattice claims without importing the producing implementations.

The source snapshot is read-only. No theorem IDs were reserved. Concurrent work is preserved. Publication checks compare every delivered byte with the remote Git tree; absence of CI checks is not a CI pass.
