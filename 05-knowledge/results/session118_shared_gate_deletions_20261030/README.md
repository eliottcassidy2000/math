# Session 118: shared-gate deletion and sequence restrictions

**PROVED local structural and accounting results; no new unrestricted circuit record.**
The user-directed target is a stronger deletion lemma or a global restriction on
worst-case elimination sequences. This checkpoint preserves the exact shared DAG,
not a formula approximation or a list of independently priced gadgets.

- [Four-deletion AP2 lemma](surplus/AP2_FOUR_DELETIONS.md): post-restriction
  nonconstant output forces four distinct deletions, including common partner
  children. Four is sharp even for a minimum circuit on the marked domain.
- [Independent shared-DAG audit and general antichain lemma](audit_ap2/AUDIT.md):
  r incomparable constant seeds force at least 2r distinct deletions, sharply.
- [Complete strengthened weight frontier](frontier/FRONTIER.md): the coefficient/
  penalty region is exactly 0<=c<=4 and b>=max(c,2c-2,4c-8). The unique optimal
  coefficient-four weights are (1,4,2), with dimension charge eight.
- [Global inventory and slack identities](accounting/SEQUENCE_ACCOUNTING.md):
  AP1+AP2<=terminal safe-core dimension; in the small-terminal regime almost all
  dimension reductions are GP and almost all pair removals are P2/PE.
- [Independent algebra audit](audit_accounting/AUDIT.md): exact history and rational
  checks, separate from the producer implementation.

The implication |C|>=4n-8(d+1) requires the full legal-state coverage and termination
assumptions of the pinned external candidate. We prove the local improvement and
its exact algebraic consequences; we do not certify that candidate's explicit
family or claim an established coefficient-four record.

Research is continuing on genuine tight GP/cleanup traces and general constant
propagation. These pending directions are not assumed in the proved checkpoint.
Run each included Python checker normally and with -O; its JSON receipt records
the finite universe. The all-parameter statements rely on the displayed proofs.
