# Verification and reproducibility

**PASS — frozen checkpoint1, 2026-11-02.**

## Commands

From the repository root (the replay itself uses an unrelated working directory):

```sh
python3 -B 05-knowledge/results/session140_collatz_hexagonal_geometry_20261102/replay_checkpoint1.py
(cd 05-knowledge/results/session140_collatz_hexagonal_geometry_20261102 && sha256sum -c SHA256SUMS)
```

The replay copies only selected inputs into a fresh temporary directory and
runs all three scientific verifiers under both ordinary Python and `-O`.
Their outputs must byte-match frozen receipts. It verifies source integrity
in both modes, then requires eight failures: a wrong five-step formula,
a wrong mod55 class, same-size proof corruption, and the independent audit's
rejection of that corrupted proof, each in both modes. Mutations affect only
temporary copies and are restored before a final integrity check.

[Replay receipt](REPLAY_CHECKPOINT1.json): **16 executions,
six scientific positive runs, two source-integrity positives, eight expected
rejections**. No assertion is removed by optimization. A rejection of a
mutated certificate is not a counterexample to the actual mathematical theorem.

## Exact finite universes

- Collatz producer: forward/inverse arithmetic for n<=10000; independent
  inverse fibres for n<=500; complete hostile fibres at root7; selected
  macro-trees to depth14 at six roots in both maps (16384 leaves each);
  independent odd predecessor controls below2000; known-cycle hexagon and
  finite square-grid matching contraction.
- HN producer:50685 checks, including all3413 functional maps on at most
  five vertices; all30 rooted quotient triangles (12 exact,18 modular-only);
  all40 mod55 unit classes; exact small-field trace and quartic controls;
  all36 pair distances of the recovered nine-vertex configuration;13 source
  snapshot hashes. Inherited Selmer/modularity proofs are not replayed.
- Independent auditor:1242745 checks, including103024 arithmetic/color
  samples and dyadic boundaries through exponent1024; complete forward-only
  predecessor fibres at200 roots; eight depth12 selected macro-trees;
 1441 labeled finite trees and77305 repeated-attachment triples; all30
  quotient triangles;21504 finite suspension systems. Suspension controls
  do not promote the pending general principles text.

The infinite claims rely on the written proofs and their independent
mathematical audit, not extrapolation from these finite universes.

## Integrity and documentation

- [Selected source manifest](SOURCE_FILES_CHECKPOINT1.json) pins39 verbatim
  producer/auditor files. All13 historical HN snapshot SHA-256 and Git blob
  identities are checked against their original commits; seven Collatz
  inheritance paths retain Git tree pins without republishing their proofs.
- [Source/status ledger](SOURCE_STATUS.md) separates new proofs, inherited
  recovery, classical citation and absent formalization.
- [Source/link verification](SOURCE_INTEGRITY_CHECKPOINT1.json) validates
  active packet links and preserves historical snapshot contexts explicitly.
- [Scoped documentation gate](DOC_VERIFICATION_CHECKPOINT1.json) checks the
  180000-byte startup and450-line frontier budgets. The central mistakes
  ledger is neither hydrated nor read; the full documentation checker is
  therefore **not claimed as run**. No Session113 approval state is changed.
- The upstream order42 frontier and its incoming wording remain unchanged
  below the compacted introduction. That result is not newly independently
  audited by this release.

## Formal and literature limits

Kühn Theorem2 is CITED with its disconnected-graph proviso and primary-text
review, not reproved or Lean-built. The external HN six-bound is CITED / not
locally Lean-built; no transitive dependency/axiom audit is claimed. Long
inherited 3-adic rigidity, descent and modularity proofs are not recertified.
No Collatz convergence, exact HN value, novelty priority or new LRC result
is asserted.
