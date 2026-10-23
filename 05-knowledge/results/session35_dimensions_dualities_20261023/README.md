# Session35: dimensions, arithmetic and tournament recursions

Audited first checkpoint, 2026-10-23.

- [Dimension and positivity](DIMENSION_POSITIVITY_HORIZON.md): finite dimension positivity horizon and representation-theoretic controls.
- [Lie and Virasoro dimensions](lie/RESEARCH_NOTE.md), with [homogeneous completion](lie/COMPLETION_ADDENDUM.md).
- [Arithmetic](arithmetic/ARITHMETIC.md): 57, Blum/Leyland calculations and explicitly typed fibre correspondences.
- [Tournament recursions](recursions/RESULT.md): exact extension/deletion closure. The marked skew pair is equivalent to the ordinary adjacency characteristic polynomial, not a stronger independent spectral invariant.
- [Tournament zeta](zeta/RESULTS.md): primitive-cycle and spectral statements, with the Paley adjacency/skew normalization distinguished.

Independent audits: [algebra](algebra-audit/AUDIT.md), [Lie/arithmetic](lie-arithmetic-audit/AUDIT.md), [zeta](zeta-audit/AUDIT.md).

Run checkers from this packet directory, for example `python3 lie-arithmetic-audit/independent_checks.py`; all included Python scripts are standalone checks. Python dependencies include SymPy. Run normally for assertion-based checks; `python3 -O` establishes output stability but disables Python assertions. Explicit exception checks remain active. Audit MANIFEST paths are relative to this packet.

This checkpoint excludes evolving collision/Singer addenda and all unapproved attachment-derived work. Finite exact calculations and conditional constructions do not establish a continuum Yang–Mills mass gap, unitary VOA realization, or a new independent invariant.
