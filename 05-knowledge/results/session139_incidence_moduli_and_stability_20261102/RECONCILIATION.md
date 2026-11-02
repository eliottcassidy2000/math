# Final checkpoint: concurrency reconciliation

The guarded publication stopped before any tree/ref write when main advanced
to `c16226249f1f27a335cf7c84afdeaa90879b9580`. The final release is rebased on
that exact tree, preserving every concurrent non-owned tree entry. The original
prepared local commit was `143675833c1bffa205eaab0b025ad36730632345`; its scientific
proofs and replay outputs are unchanged. [Machine record](RECONCILIATION.json).

Only four concurrent overlaps needed prose merging: CURRENT-FRONTIER, S111-plus,
S139 and this packet's README. Current-best routing is140, with10..139 open.
The separate276 shared-baseline,314 and466 routes retain their provenance and
credits; their own statements and historical thresholds are not rewritten.
The concurrent modulus/FRACTRAN precision work, S138 routes, result index,
THM-4852 correction, and central ledger remain exactly as in the new parent.
The central ledger is not inspected or edited by this reconciliation.

The earlier full documentation-overlay receipt remains a record for the
pre-reconciliation input. The new scoped check verifies the bounded startup
surface, merged links, explicit overlap provenance and exact preservation of
concurrent tree entries. It does not claim a fresh full-documentation check
of the opaque concurrent ledger. [Scoped receipt](DOC_RECONCILIATION_CHECK.json).
