# Audit record and exact reproduction boundaries

2026-10-26. **VERIFIED** independent mathematical review and **FINITE-EXACT**
computational controls. This is not a Lean formalization. The all-order
claims follow from the written proofs and explicitly cited structural
inputs, not an extrapolation from the finite receipts.

## Proof review

- The arithmetic reviewer checked the fibre proof, including the invariant
  large-port set, p-subgroup fixation, the sharp 4p construction, uniqueness
  at twelve, and the lifting of ordered path-cover bijections. The primary
  realization paper was independently checked. See
  [FIBRE_AUDIT.md](arithmetic/FIBRE_AUDIT.md).
- The fibre reviewer checked the arithmetic proof, including both directions
  of the Sylow-tower classification, the free action on ordered partitions,
  the carry formula, and the sharp strong-connectivity correction. The
  complete checker was rerun independently. See [fibre/AUDIT.md](fibre/AUDIT.md).
- The metagraph reviewer checked the parent's decorated-tree classifier,
  exact split extension, lift mass and Burnside count, and the random-source
  connection contracts. See [PARENT_AUDIT.md](metagraph/PARENT_AUDIT.md).
- The parent independently checked M1–M7: forced block order proves join
  closure; conjugation transports both the product and the movement relation;
  only the retained hinge can create a two-block cyclic triple; a singleton
  middle block accounts for the extra three-block triple; grouping by rooted
  target classes proves the multiplicity convolution. Colour is explicitly
  excluded from that sufficient interface. The parent also independently
  constructed the order-five pair and enumerated all Hamiltonian words in
  `check_connections.py`, obtaining 11 and 13 without producer imports.

The audits corrected terminology (`Aut(Q)/Stab(0)` is a coset space),
distinguished repository errors from the published realization theorem,
and kept maximal proper strong modules distinct from arbitrary modules.
No uncorrected mathematical objection remains in these scoped results.

## Finite universes and independent routes

| Lane | Universe / positive controls | Hostile control / independent route |
|---|---|---|
| Fibre | All 64 order-four orientations; all 128 specified order-twelve lifts; all 1296 hypergraph symmetries; prime construction controls p=3,5,7,11,13 | Module formula versus direct automorphism backtracking; subset Hamilton counts for all lifts; independently written arithmetic-agent arc code for four representatives |
| Arithmetic | 160 pair-orbit cases (n<=40, p=3,5,7,11); all 8930 fixed objects in the stated small-cube range; all labelled tournaments through five; 152 strong attaining constructions | Pair-orbit generator action versus digit formula; all permutations versus constrained Aut search; same order81 but distinct groups; carry-defect example; Paley7 prime with nonsimple Aut21 |
| Metagraph | 4752 two-factor presentation pairs; 10904 named relocation instances; 116 triple-factor tuples; all 33868 presentations through order seven for factorization; 32 class-pair convolutions | Exact minimum two-factor ambiguity at total order five; H from subset DP and all words; independent parent matrix construction |
| Random recovery | 24 frozen files selected from 12883; five deeper connection probes | Integer checks for Sylvester coprimality/remainder, recurrence seed, C25/C5² subset orbits; unsupported geometric transfer stopped |

There are no imported census filters in the new exact checkers. Larger
ambient tournament universes are not purportedly exhausted. The published
C3 realization theorem supplies completeness of the displayed lift fibre;
the sharp all-order statements have separate proofs. The old A000568 census
was preserved, not rerun or independently re-certified by this packet.

## Provenance and freezing

Initial live base: `fc3768ac5fcd`. `random_recovery.py` reads those exact Git
blobs, preserving the original random source population and hashes even
after repairs or concurrent additions. The correction lineage is the Git
diff of THM-4803, its source notes and mirrored historical reflection, plus
the explicit entry in `01-canon/MISTAKES.md`.

`SHA256SUMS` covers every packet file except itself, including the verifier
and proof notes. The arithmetic independent fibre-audit receipt additionally
pins the reviewed fibre proof's bytes. Normal and optimized receipt pairs
must match. `verify.py --replay` checks the manifest first, then computes into
temporary destinations and compares bytes to the frozen receipts; it never
updates an expectation. Source edits require an intentional new freeze and
review of affected receipts. `agents/check_docs.py` and `git diff --check`
are the repository documentation gates.

**Remaining limits:** no all-prime atlas, no composition theorem for all
coloured metagraph data or arbitrary substitution, and no new LRC(14)
consequence. No historical-priority claim is made for the elementary group
identities or the classical modular decomposition.
