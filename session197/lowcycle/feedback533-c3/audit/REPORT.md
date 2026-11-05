# Independent review: the H533 feedback-two cycle threshold

**PASS, no correction.** The displayed graph has H=533, exactly three Hamilton cycles and feedback vertex number exactly two. The fixed-port word constructs the literal graph. Combining this with the separately reviewed all-order negative theorem proves that three is the minimum cycle count among strong H533 tournaments with feedback number at most two.

## Frozen binding

Under `session197-lowcycle/`:

- `FEEDBACK533_C3.md`: `db75680e66e686d5f2074e6a7d5bd2e3b22c8ec05c1dc8bce330f99836720484`.
- `verify_feedback533_c3.py`: `01651c832650be28d449c18f6b5c01f54b4a2e0a78e8b39c768dd5fecb933bea`.
- Certificate: `b76a4f50ed6069de2755e4f557417dc0890011c667c6c8c1710b521fea9e35f9`.
- Negative theorem dependency `FEEDBACK533.md`: `24c37bda3f05d45bb98f685f570046313ae1a464573932b145b7d7a58f83a3cb`, independently reviewed in `session197-bridges/feedback533-audit/REPORT.md` (hash `65fd9fc3ba5d41b310063a1152cb8a2d569da97fe2e5da13e5d8eaf0a5903986`).

## Independent proof and controls

The independent checker obtains H=533 by inclusion-exclusion over unrestricted walks, rather than the source endpoint DP or path DFS. It tests all9! permutations of vertices1,...,9 after the fixed vertex0 and obtains exactly the three displayed directed Hamilton cycles. Any one is also a direct strongness certificate.

All1,024 vertex subsets are tested for transitivity by the tournament score criterion. The largest transitive subset has size8; the deletion pair[1,3] is among the minimum feedback sets. This independently proves exact feedback number two, beyond checking one proposed pair or the source's singleton triangle certificates.

The fixed-port word is independently reconstructed into adjacency rows: initial edge1->3, ordered retained ports[1,3], the specified transitive birth order, and masks1,1,1,2,3,0,1,1. All45 resulting arcs match the displayed graph. No port changes or equal-count substitutions occur.

The global minimum claim then follows logically from the named negative theorem: a smaller cycle count in this same strong H533/FVS<=2 class would contradict that theorem. The positive packet does not enumerate all C3 classes, minimize order, or establish unrestricted fixed-port coverage.

## Replays

- `check.py`: `96a7fc97d3ef8ba3a49d21519cd2076acacdaebe9211e664ffb4b9aa14cfe670`.
- Independent normal/optimized receipt: `1dec9cd0566861460a24d24161d953719cee99c16a97aa3e7123d450ee415110`.

The source verifier and its exact negative-proof dependency were copied before execution, since the verifier writes its certificate. Both source modes reproduce the frozen certificate byte-for-byte. Source files remain unchanged.
