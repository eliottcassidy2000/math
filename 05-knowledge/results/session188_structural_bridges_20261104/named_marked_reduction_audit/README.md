# Named marked-cover independent audit

Frozen publication allowlist: `REVIEW.md`, `README.md`, `audit.py`, `runs.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`.

Run `python audit.py > normal.json`, `python -O audit.py > optimized.json`, then compare. The audit imports no producer implementation. Its marked tables use unordered set partitions and explicit designated-slot assignments; its physical counts use a separately written predecessor recurrence retaining each first-child run number.

PASS: all78 complete marked-table hashes, explicit repeated-name and singleton-role controls, and all96 physical run vectors across six complete faces. The result remains sufficient finite marked data, not a nonannihilation or reconstruction theorem.
