# Verification and reproducibility — checkpoint2

**PASS — 2026-11-02.** The first-checkpoint science and audit remain frozen.

```sh
python3 -B 05-knowledge/results/session140_collatz_hexagonal_geometry_20261102/replay_checkpoint2.py
(cd 05-knowledge/results/session140_collatz_hexagonal_geometry_20261102 && sha256sum -c SHA256SUMS)
```

[Latest replay](REPLAY_CHECKPOINT2.json) performs **19 top-level
executions** from a fresh relocation and unrelated working directory:

- Eight new scientific runs: principles, curvature, independent checkpoint2
  and independent native-port controls, each in normal and optimized Python.
- Two positive source-integrity runs.
- Eight required rejections in both modes: wrong period-clock sum, wrong
  hyperbolic edge constant, same-size principles-proof corruption, and the
  independent auditor's detection of that corrupt proof.
- One complete first-checkpoint replay, byte-matching its frozen receipt;
  this additionally executes its 16 internal runs (six scientific positives,
  two integrity positives and eight required rejections).

All 85 replay inputs have SHA-256 pins. Explicit exception
checks remain active under `-O`. Mutations affect temporary copies only.
Finite controls are not substitutes for the written universal proofs.

## New finite universes

- Principles: 131,071 distinct odd macro-tree labels through depth 16;
  4,330 map/clock systems; 30,301 lattice sites; 10,000 site-conjugate steps;
  100,000 bitlength-coloring edges; 42,525 exact rational two-variable jet
  checks, plus boundary-marginal, subdivision, affine-carry and torsion hostiles.
- Curvature producer: 255 vertices, 254 edges, 32,385 vertex pairs and 32,131
  edge pairs; sharp distance/equality, geodesic circle and band controls.
- Independent checkpoint2: 101,445 checks ;127 hyperbolic vertices and 8,001
  pairs; 7,875 edge pairs; 64 non-quasi-isometry hostile depths; a 31-vertex
  exact rational proper unit-tree control; 30,301 addresses; 10,000 exact-clock
  steps; 2,625 three-variable jet controls.
- Narrow J1–J4 audit: 30,940 checks; 8,001 source pairs; merged 425-vertex native
  tree; 20 generating-polynomial comparisons; 385 Gram determinants; both
  fixed-port hostiles. **Leaf-Gram rank 2^h requires h>=1; h=0 has rank 0.**

## Preservation, sources and visual artifact

[Integrity report](SOURCE_INTEGRITY_CHECKPOINT2.json) verifies all prior
immutable files, 20 targeted historical principles source Git/SHA identities,
the separate concurrent joint-port source pin, and active publication links.
The original joint-port packet and its manifest are **unchanged**; the
independent addendum explicitly corrects its h=0 rank boundary. This is not
a whole-packet or external-input seal replay.

[Curated file distinction](CURATION_CHECKPOINT2.json) preserves the curvature
producer's original selector/hash list as provenance, while publishing only
text, SVG, generating source and exact receipts. The PNG was visually checked
but is not published. The SVG is XML-validated; sampled drawing coordinates
are illustrative, not mathematical certificates. There are no caches or
external-paper binaries in the payload.
The verbatim Matplotlib SVG has 1,986 trailing-whitespace path-coordinate
lines; its original hash is preserved. The [style receipt](STYLE_VERIFICATION_CHECKPOINT2.json)
records this narrow generated-file exception; every other path passes the
strict Git whitespace gate. No mathematical check is waived.

The current [scoped documentation gate](DOC_VERIFICATION_CHECKPOINT2.json)
uses the published checkpoint1 parent, preserving the fresh concurrent
S61/S138/S139 routes. Startup remains 179,892/180,000 bytes and frontier 450/450
lines. This checkpoint does not edit the frontier. The full checker, which
reads the central ledger, is not claimed as run; that ledger is neither read
nor changed. Pending Session113 approval and the upstream order42 status
remain unchanged.

No local Lean build, external HN proof/axiom audit, or full recertification
of the inherited theorem inventory is claimed. The accepted elementary
proofs, CITED classical inputs and finite controls are separated in the
[source/status ledger](SOURCE_STATUS.md).

---

# First-checkpoint verification (retained history)

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
