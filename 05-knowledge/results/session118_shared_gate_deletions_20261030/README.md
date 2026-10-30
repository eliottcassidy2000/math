# Session 118: shared-gate deletion and sequence restrictions

**PROVED local structural and accounting results; no new unrestricted circuit record.**
The user-directed target is a stronger deletion lemma or a global restriction on
worst-case elimination sequences. This checkpoint preserves the exact shared DAG,
not a formula approximation or a list of independently priced gadgets.

- [Four-deletion AP2 lemma](surplus/AP2_FOUR_DELETIONS.md): post-restriction
  nonconstant output forces four distinct deletions, including common partner
  children. Four is sharp even for a minimum circuit on the marked domain.
  [Companion sharpness](surplus/ALL_AND_ROWS_SHARP.md) retains the exact three-deletion
  AP1 and two-deletion one-child GP cases, with their compensation caveats.
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

- [Genuine tight-prefix family](worst_cases/TIGHT_CYCLE_HOSTILE.md): arbitrarily many
  exact GP/P2 cycles with shared equations and actual loop eligibility. These circuits
  are nonminimum and have linear dispersal dimension, so the low-dimension/minimum
  regime remains open. [Independent audit](audit_cycles/AUDIT.md).
- [Safe factoring rewrite](audit_cycles/FACTORING_REWRITE.md): shared boundary wires
  are preserved while an exclusive XOR cone is replaced. This saves one ordinary gate
  per block in the hostile family; an explicit legal affine role exchange saves one
  more charged equation per block. The new-coordinate state is minimum, but this
  is not a universal surplus theorem or a free input-coordinate change.

[Research synthesis and remaining obligation](SYNTHESIS.md) connects the positive
lemma, global restriction, counterexample and repaired rewrite. In particular a
complete low-dimension/minimality-sensitive dichotomy for GP/cleanup is still OPEN.

## Reproduction

Run `python3 verify.py` in this folder (or by its repository-relative path).
It verifies file hashes and runs all nine Python checkers normally and under -O,
requiring identical outputs and the recorded exact receipts. All eighteen runs
passed during publication. The finite universes are explicit in each proof/audit.
They corroborate the all-size proofs; they are not Lean certification or an
independent human review. No challenge submission or external proof code was run.
