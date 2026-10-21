# Independent source63 debt audit

Verdict: **PASS for the frozen packet's scoped proofs and finite-exact source
certificates.** Read [PROOF_AUDIT.md](PROOF_AUDIT.md) for the proof and boundary
review. No producer file, repository file, or remote was modified by the audit.

`check_source_debt_cones.py` independently parses the mathematical note,
reconstructs the 63-exponent residual universe, and discovers each original
source's native trajectory without using its stored word as a driver.
It uses only the Python standard library and explicit active checks.

After publication, run from the repository root:

    python3 05-knowledge/results/source63_separator_debt_20261021/audit/check_source_debt_cones.py --repository . > /tmp/debt-audit.normal.json
    python3 -O 05-knowledge/results/source63_separator_debt_20261021/audit/check_source_debt_cones.py --repository . > /tmp/debt-audit.optimized.json
    cmp /tmp/debt-audit.normal.json /tmp/debt-audit.optimized.json
    cmp /tmp/debt-audit.normal.json 05-knowledge/results/source63_separator_debt_20261021/audit/receipt.normal.json

The script's default packet path is its parent directory, independent of the
working directory. For a separately staged packet, pass `--packet PATH`.
`--repository PATH` checks all upstream hashes against pinned Git objects and
is needed to reproduce the saved receipts' authentication section. Without
that option, the portable mathematical and packet-integrity checks still run,
but the upstream-authentication list is empty and that receipt is different.
No command modifies the repository.

The normal and optimized saved receipts are identical: 1,204,846 active checks,
with independent totals 575,652 odd, 1,152,705 shortcut, and 1,728,357 ordinary
ROOT steps across the 63 original sources. `ADD_ONLY_MANIFEST.json` specifies
the exact new publication paths and hashes. The publisher must refuse an
existing destination rather than overwrite or append to it.
