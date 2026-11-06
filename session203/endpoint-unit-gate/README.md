# Fixed rooted unit-correction gates

Status: **PROVED in the stated scope, independently audited.**

A fixed physical context copying exactly one input root, with every other core vertex uniform toward the input's nonroot vertices, cannot realize an even-multiplier recurrence with constant correction on all strong unit-start inputs. The obstruction persists above every source-count cutoff. The two cores may have arbitrary orders and need not preserve one another's induced graphs.

The mechanism is an exact twin contraction: two adjacent free root-incidence bits have mixed difference minus twice an actual tournament's odd Hamilton-path count. The difference of child and any even multiple of parent is therefore 2 modulo 4 on the four-input square.

The edge/triangle positive specialization and the larger family-C replacement instead have the physically valid correction 2 start_root(U)-H(U). The source count is retained, not independently selectable. These are interfaces, not an eventual spectrum theorem.

## Files and reproduction

- PROOF.md: context class, twin lemma, four-source obstruction and positive control.
- TAIL_ADDENDUM.md: arbitrarily large physical strong/unit-start input sources.
- FAMILY_C_AUDIT.md and KNOWN_C_TENSOR.json: independent exact larger-core tensor.
- audit/REVIEW.md: independent audit, including source, scope, and zero-retained-root hostile.
- reviewed-originals/: pre-promotion reviewed proofs. PROMOTION.json records the status-only changes.
- dependencies/: frozen sources corresponding to the historical paths named in the proof. PROVENANCE.json binds their original paths and hashes.

From this directory:

    python3 verify.py
    python3 -O verify.py
    python3 check_known.py
    python3 audit/check.py
    python3 audit/check_tail_companion.py

Both producer and independent normal/optimized receipts agree. Every verifier is standalone and uses exact integer counts. The all-order theorem follows from the proof, not finite enumeration.

No unconditional eventual odd-spectrum threshold is claimed. Target-dependent context choice, multiple sensitive copied vertices, positional incidences, input-edge surgery and other source classes remain outside this theorem.
