# Infinite metagraph families and minimum root defects

**PROVED IN NOTES / FINITE-EXACT controls. Chan–Pak OPEN.**

- [Capped families and operations](FAMILIES.md): for every fixed H cap,
  finitely many strong types plus singleton gaps give exact local birth,
  deletion and arc rules; the capped native part is finite. Retaining fixed
  names gives a finite source-return quotient after free singleton padding.
- [The complete H<=9 slice](FAMILIES.md#4-the-complete-h9-language): six
  SCC-kernel families,15 rooted strong-birth rules,32 directed arc-orbit rules,
  and14 literal native lines. Exactly five neutral components survive at
  arbitrary order. The two-triangle component has sharp upward escape17.
- [All minimum-defect sources, d>=3](MINIMAL_DEFECT.md): the sharp bound
  b>=2^d-2 has a complete equality classification by cyclic transitive collars,
  with the C3 inport alternative at d3. The two-inport ready pentagon is the
  exceptional boundary, not a universal pattern at larger indegree.
- [Sequence export](sequences.tsv): ordinary, rooted, frame, converse,
  attachment-orbit and two/three-named refinements of A000568 through50.
  Closed formulas and the full automorphism cycle-index inventory hold at
  every order. These are not new unrestricted A000568 terms.
- [Independent fixed-defect theorem audit](audit_templates.md): the finite
  transitive-inflation-box proof passes its positivity, threshold and quotient
  checks. Rooted automorphism-group strata and fixed numbers of additional
  names also have eventual polynomial enumerators.

Receipts: [cap dictionary](cap9.json), [independent operation audit](operation_audit.json),
[minimum-defect controls](minimal_defect.json), and
[independent two-inport defect audit](defect_audit.json).

Reproduce from this directory:

    python -B classify_cap9.py --output cap9.json
    python -B audit_operations.py
    python -B minimal_defect.py
    python -B audit_defect.py

The operation audit also reads the pinned concurrent birth-family coefficient
file named in its receipt. No large external census is regenerated. The
all-d minimum-defect proof has additionally passed the arithmetic lane's
independent conceptual review.
