# Session 37: exclusion statistics, parity and arithmetic cuts

This audited research packet distinguishes vertex exclusion, exchange signs, component occupation, and arithmetic completion. Start with the [final synthesis](SYNTHESIS.md).

## Results and independent audits

- [Exclusion and mixed statistics](exclusion/RESULT.md): commuting square-zero path labels, paired exterior realization, and exact polynomial recovery boundaries.
- [Exterior tournament statistics](parity/EXTERIOR_TOURNAMENT_FOCK.md) and [super-PBW and signed curl](parity/SUPER_PBW_AND_SIGNED_CURL.md): distinct component and vertex parity constructions.
- [Arithmetic cut completion](arithmetic/ARITHMETIC_CUT_COMPLETION.md): exact field and cut information, with adjacency/skew distinctions retained.
- [Cycle-gas boundary](root/CYCLE_GAS_REAL_ROOT_BOUNDARY.md) and [one-back-arc family](root/ONE_BACK_ARC_CYCLE_GAS.md).
- [Statistics audit](statistics-audit/AUDIT.md), [arithmetic audit](arithmetic-audit/REPORT.md), and the scoped synthesis reviews in those folders and inheritance-review.
- [Previously published inheritance corrections](../session37_inheritance_audit_20261023/REPORT.md): THM-506 functional recovery through eight vertices, first failure at nine, and filename-qualified bosonic THM-1810 repairs.

## Reproduction and scope

Use the commands in each result and audit from its indicated folder. Run exclusion/export_n8.cpp's compiled executable from the packet root, so its output is exclusion/n8_classes.csv. Run sha256sum -c SHA256SUMS from each file's containing folder; statistics-audit/INPUT_SHA256SUMS is likewise relative to statistics-audit. Source checksum files retain an explicit relationship to their audited inputs; publication portability repairs are logged separately. Normal Python executions validate assertion-bearing checks. Optimized executions establish execution/output stability, not activity of removed assertions. C++ census universes and exact receipts are stated in the linked reports; random witness discovery is not required to verify a fixed witness.

The polynomial pair, skew spectrum, and exterior character are different observables. No physical quantum interpretation, new irrationality or transcendence theorem for a new named constant, or solution of LRC(14) is claimed. External attached-paper critiques and previously withheld material are excluded.
