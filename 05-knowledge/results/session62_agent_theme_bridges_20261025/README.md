# Session 62: exact continuations of other agents' mechanisms

**Status: two audited checkpoints.** This packet starts from live repository boundaries rather than repeating the previous six/seven classification. It distinguishes proofs, bounded exact controls, and open applications; no literature-wide novelty claim is made.

## Wrapped physical incidence: the missing periodic boundary is closed

The inherited cube-plus-bridge theorem stopped before physical line labels wrapped. The [new proof](incidence/cyclic_common_kernel.md) retains the actual line partitions and replaces the boundary equations by an anchored path forest. It gives a constructive common-kernel basis for **every nonempty subset of a prime cycle with `p>3`**, over **every coefficient field**.

For a physical prime `p>=17`, an interval of length `1<=L<=p`, and `delta=1` exactly when the coefficient field has characteristic `p`, the three-direction common kernel satisfies

```
dim J = (L-1) + (L-10)_+ + (L-p+6)_+ + (L-p+3)_+
        + delta * 1_{L>=p-6}.
```

The characteristic-dependent mode first survives at `L=p-6`; full cyclic length has dimension `2p-2+delta`. Physical primes `5,7,11,13` are covered by the general path theorem and explicit tables, not by applying the displayed large-prime formula outside its range.

The [independent audit](audits/incidence/AUDIT.md) reconstructs physical point/line incidence without importing the producer and adds exact Smith-form controls. These finite tests support the all-range proof; they do not replace it. The coupled-kernel identity `dim ker M=2 dim J` holds under the hypothesis that all three weights are **nonzero**. Reusing the physical weights `(-3,-9,1)` in characteristic three violates that hypothesis and gives a recorded counterexample.

## Permuted pairs: the missing information is a local coefficient system

The [unordered-pair transport proof](holonomy/TWISTED_PAIR_LIFTS.md) extends the ordered-pair model by retaining both endpoint-flip gains and pair permutations. For a connected graph of cycle rank `beta`, the parity observer forgets precisely a rank-two local-system cohomology space of dimension

```
2*(beta-1) + f,
```

where `f` is `2`, `1`, or `0` for trivial, transposition, or transitive pair holonomy. Thus the old `2*beta` loss is not valid for arbitrary pair permutations. A three-cycle holonomy repairs every even gain on a single cycle, but does not erase higher-cycle obstructions. Compatible endpoint choices are not the same as trivializing the entire affine transport.

The [independent audit](audits/holonomy-clone/HOLONOMY_AUDIT.md) reconstructs the affine equations and Fano parity halves without importing the producer. The geometric comparison concerns marked label actions, not a plane-colouring construction or an identification of different eight-state automata.

## Clone cancellation: one element and the whole group see different ranks

The [clone theorem](clone/README.md) proves that, for a finite odd-order permutation action on gaps, the subset-orbit polynomial vanishes at `t=-1` to order equal to the **minimum number of cycles of one group element**. The rank of legal invariant toggles instead counts the orbits of the **whole group**. These quantities need not agree: a genuine 24-vertex tournament realizes the regular `C3^2` action on nine gaps, with cancellation order three but only one invariant toggle direction.

The proof includes the exact operation sidecar, equality boundary, real exterior-invariant interpretation and an unbounded family of actual tournament realizations. The [independent audit](audits/holonomy-clone/CLONE_AUDIT.md) separates the abstract group calculation from the intrinsic tournament construction. No faster enumeration or black-path transport follows automatically.

## Palettes: Boolean states, weighted continuations and integral information

The [palette continuation theorem](palettes/RESULT.md) retains joint boundary relations and their relative colour alignments. Over the rationals, with all singleton-list letters and named queries allowed, three linear interfaces have different sharp dimensions: `q` for a nonempty prefix, `q+1` for the whole path-counting series including the empty path, and `q^2` for two-port response. These are not the inherited Boolean automaton's `q+2` states.

The [independent audit and extensions](audits/palettes/AUDIT.md) identify the exceptional coefficient characteristics and the exact integer lattice: the quotient by the integer span of **nonempty-word** matrices is `(Z/(q-1)Z)^q`; including the identity reduces this to `(Z/(q-1)Z)^(q-1)`. The [independent Smith controls](integer-response/integer-response-normal.json) verify `q=2,...,8`; the general statement has an explicit lattice proof, not a numerical extrapolation.

A single homogeneous linear scalar register attached to the deterministic Boolean state loses necessary responses. This does **not** exclude an arbitrary `q+2`-dimensional weighted automaton: the `q+1` construction already fits. Separate palette marginals likewise cannot replace one common exterior assignment. No new plane-colouring bound or realizable equality gadget is claimed.

## Scope and inheritance

The source nonwrapping theorem was correct in its stated range: this is an extension, not a retraction. Kernel arrays are not automatically binary errors, realizable sections, or norm-graph bounds. No LRC(14), plane six/seven, integer-square hourglass, or unrestricted black-connectivity closure follows.

[Concept board](CONCEPT_BOARD.md) · [Source pins](RECOVERY_SOURCE_PINS.json) · [Reproduction and publication gates](PUBLICATION_CHECKS.md) · [Exact file manifest](MANIFEST.json)

Run `python3 verify_packet.py` from any directory to verify the published file set and bytes. Mathematical reproduction commands and exact finite universes are in the lane note and audit. Remote commit verification is recorded separately; local checks do not imply external CI success.
