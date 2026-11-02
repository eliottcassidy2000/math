# Source status and verification

- The translation lemma and triangle-plus-attachments argument are fully
  proved in `PROOF.md`. No external theorem beyond elementary graph facts
  is required. No novelty claim is made.
- The functional defect formula is reproved briefly in the corollary and
  was proved for all finite graphs in the preceding frozen packet:
  `unit_distance/PROOF.md`, Theorem 1.
- Numerical input only: Boris Alexeev, Dustin G. Mixon, Hans Parshall,
  *The Erdős unit distance problem for small point sets*,
  arXiv:2412.11914v2, dated 2025-02-12, Theorem 1 and Table 1:
  https://arxiv.org/html/2412.11914v2 . This is the same primary source
  already read and pinned in `unit_distance/SOURCE_PINS.json`.
  The computer-assisted global extremality proof is cited, not rerun.
- The table subtracts n and adds one to obtain ε and β. These elementary
  arithmetic operations introduce no new extremal-value assumptions.
- The refinement is proof-only; no standalone checker, floating-point
  experiment, or normal/optimized Python acceptance claim is supplied.
- The original frozen packet is unchanged. Export only `PUBLIC_FILES.json`.
  Verify this supplement with `sha256sum -c SHA256SUMS`; the hash manifest
  covers every allowlisted file except itself, including the selector.
