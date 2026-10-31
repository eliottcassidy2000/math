# Portable independent minimal-core audit

**PASS:106,801,611 explicit checks in normal and optimized Python.** The two receipts agree byte-for-byte.

The independent auditor imports no producer module or sibling auditor and reads no producer Python. It consumes only the196 new JSON certificates, the existing base manifest and11 chunks, and package hash metadata. It requires no old transition table, Git database, network access or external library.

Reproduce from this package directory:

    python audit_independent.py --package-dir . --base-dir ../square_sum_B_base_interval_20261019 > independent.normal.json
    python -O audit_independent.py --package-dir . --base-dir ../square_sum_B_base_interval_20261019 > independent.optimized.json

Replace the base path explicitly if the package is not installed beside the existing data directory.

The audit checks:

- all2568 base pairs;
- all98 paired and98 terminal universal states;
- all98 paired depth-two expansions;
-196 terminal cycles across every state at depths1 and2;
- eight hostile certificate mutations, all rejected;
- the frozen hashes, symbolic ownership and length contracts, common parity only for reusable paired outputs, and the threshold/prefix identities.

Auditor SHA256:

    3531f81468b7815672dc4eae9466d2ea5f99b5a1d531da145bda2f7e3da967f4

Each retained receipt SHA256:

    49770a428a1e0cc1311e951d78270cb3c010dad91a393b33586ee280088ba74e

The independent mathematical dependency-reduction argument is in `INDEPENDENT_AUDIT.md`. No all-size inference relies on finite expansion counts alone; the universal atom semantics and induction are explicit.
