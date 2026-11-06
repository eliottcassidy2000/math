# Exact append-only signed-section targets

Session 201. **PROVED, independently audited in the stated scope.** This
continuation starts from the independently audited physical signed-section
theorem. It classifies all section targets reachable from a fixed physical
source by appending letters, proves a sharp minimum, and gives an infinite
obstruction to repairing an arbitrary signed representative by append-only
operations. Universal target coverage remains open.

The separately audited SIDECARS.md proves the exact first-growth gap and
an all-order parity obstruction to a proposed formal-boundary collision.
Its source domain is the nontrivial orbit S, excluding the singleton.

Files:

- PROOF.md and verify.py: exact append-only target classification.
- SIDECARS.md and verify_sidecars.py: separately scoped consequences.
- VALIDATION-*.json and SIDECARS-*.json: matching normal/optimized receipts.
- dependencies/SIGNED-SECTION.md: the frozen independently audited dependency.
- audit/: independent proof reviews, positive-coordinate implementations,
  ordinary/optimized receipts, original inputs, and promotion verification.

The target-linear finite decision bound is proved, so each stated bounded
comparison is complete for its declared source/target window.

From this directory, reproduce with:

    python3 verify.py
    python3 -O verify.py
    python3 verify_sidecars.py
    python3 -O verify_sidecars.py
    python3 audit/independent_check.py
    python3 -O audit/independent_check.py
    python3 audit/independent_sidecar_check.py
    python3 -O audit/independent_sidecar_check.py

The audit reports retain their original scratch-directory names as historical
provenance. The publication belongs to session201/append-only-signed-section/.
No mathematical body, implementation, dependency or receipt changed during
promotion; the sidecar gained an explicit statement of its inherited domain.
