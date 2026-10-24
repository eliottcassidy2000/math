# Session52: marked regularization and circulation transport

This portable mathematical packet extracts exact structural consequences
from a **known** Type-I regularization, and tests the limits of forgetting
its root. Read `THEOREM.md` for statements, proofs, inheritance, and scope.

## Main outcomes

- The known map T→D(T), of orders n→2n+1, is fully faithful with a root:
  Aut(D(T))_z=Aut(T), with a uniquely recoverable paired copy.
- Every source edge reversal lifts to an edge-disjoint directed four-cycle
  reversal. Source Hamiltonian paths biject with distinguished matched
  alternating target Hamiltonian cycles.
- Source fixed-path complement lifts exactly to block reversal, but differs
  from full target noncycle complement on exactly 5(n−1) target edges.
- If T has neither source nor sink, the root is intrinsically recognized
  by positive common-out-neighbor counts with every other vertex. On this
  subdomain the unmarked map and full automorphism group remain faithful.
- Recentring at a source-copy vertex proves D(1+Q)≅D(Q+1). Ordinary class
  faithfulness first fails at source order four (target nine). Converse-merged
  faithfulness first fails at source order five (target eleven), witnessed
  by 1+1+C3 versus 1+C3+1. An explicit order-17 witness is also included.

The construction is credited to Hatem, Mohamed, Müller and Siniora,
arXiv:2605.04846v1 Lemma 4.3. It is not THM-4557's skew-Hadamard doubling.
No general A000568 formula, native full-target C-line correspondence,
port-clock transport, or first-novelty claim is made.

An elementary count corollary gives N_n=A_n−2A_(n−1)+A_(n−2)
no-source/no-sink ordinary classes and J_n=S_n−S_(n−2) self-converse ones.
These provide ordinary, self-converse, and merged regular-class lower bounds
N_n, J_n, and (N_n+J_n)/2 at order 2n+1. They are weak comparison bounds,
not a fast or asymptotically sharp formula for A000568.

## Reproduction

From this directory, using only the Python standard library:

    python3 verify_manifest.py
    python3 check.py --output replay-normal.json
    python3 -O check.py --output replay-optimized.json
    cmp replay-normal.json normal.json
    cmp replay-optimized.json optimized.json

Expected: status PASS, byte-identical replay JSON, 1,099 labelled source
matrices across orders 1 through 5, 10,650 single-arc lifts, 7,887 source
Hamiltonian paths, and 1,100 source-to-sink recenterings for Q orders 0..5.
Nine hostile mutation controls must be rejected. Normal and optimized
outputs must be identical. Total wall time is normally well under a minute.

## Files and audit scope

- `engine.py`: matrix and independent arc constructors, exact generic
  isomorphism search, independent source permutation checks, path tools.
- `check.py`: exhaustive, explicit witness, and hostile controls.
- `normal.json`, `optimized.json`: complete deterministic evidence, including
  full adjacency matrices and verified isomorphism permutations for witnesses.
- `normal.out`, `optimized.out`: captured standard output of those runs.
- `SOURCE_PINS.json`, `sources/`: source URLs, attribution boundaries, and
  the pinned contrasting repository theorem.
- `MANIFEST.json`, `verify_manifest.py`: local artifact integrity and replay
  portability. The manifest hashes packet files, not the external website.

The generic isomorphism search does not use the structural theorem to prune:
only adjacency, common-neighbor counts, and their invariant multisets are
used, followed by a full-matrix bijection check. No repository, third-party
package, network, or external tournament atlas is needed. The finite first-loss
thresholds are exhaustive; the all-size laws rest on the written proofs.
