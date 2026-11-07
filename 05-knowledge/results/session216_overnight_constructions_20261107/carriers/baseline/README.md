# Paired-source coverage investigation

Independently audited proof and exact controls. F(8,137) and eventual odd-spectrum completeness remain open.

- `PROOF.md`: arbitrary common right factors preserve the synchronized two-residue adjusted base exactly when their root deletion has H=1; these factors are precisely O/Z words. This is a scoped maximality result, not a global construction obstruction.
- `probe.py`, `matches.json`: complete symbolic single-nonroot-module test on the six inherited pairs. Its four successful position pairs give exactly old O/Z families.
- `verify.py`, `RESULT.json`, `RESULT_OPTIMIZED.json`: 32 independent literal module profiles and 2,754 exact arithmetic controls, identical with and without Python assertions.
- `amalgamation*.cpp`: asymmetric-rooted-factor searches, with raw text receipts. The broadest search checks 3,467,984 candidates after all left O/Z/J prefixes through length2, finding no two-residue persistent matches involving a nontrivial right factor. This finite failure is not an all-order theorem.
- `PROVENANCE.json`: frozen input provenance; `SHA256SUMS`: payload hashes.

All paths are local packet paths. No Git publishing was performed by this worker.
