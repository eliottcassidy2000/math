# Session 144: Jacobian classification, ETH and MRB

**Completed in two independently audited checkpoints, 2026-11-03.** [Read the synthesis](SYNTHESIS.md) for the distinctions between proved classifications, constructed families, external claims and open problems.

- **Jacobian:** [composition tower](jacobian/equivariant_composition_tower.md), [complete fixed-map fibers](jacobian/independent/INDEPENDENT_AUDIT.md), [all-degree inheritance and repairs](jacobian_followup/classification_scope_and_repairs.md), and [universal equivariant boundary law](jacobian_followup/universal_equivariant_boundary_corollary.md). [Arbitrary primitive scale degree is a gauge artifact](jacobian_followup/equivariant_gauge_hostile.md); [a quartic hostile requires an escaping-sheet term](jacobian_followup/independent/QUARTIC_BOUNDARY_AUDIT.md). Full map-equivalence classification and the planar conjecture remain open.
- **ETH:** [source, quantifier and encoding audit](eth/ETH_AUDIT.md), an exact bounded-integer Keller-collision interface, and quantitative reduction budgets. The imported large SAT result is CITED/EXTERNAL, not newly built here. [All eight source byte pins independently match](eth/INDEPENDENT_SOURCE_BYTES.json).
- **MRB:** [complete periodic-filter convergence and finite-part formulas](mrb/RESEARCH.md), stronger entire-perturbation nonuniqueness, and [an exact rational precision evaluator](mrb/audit/README.md). [Packet and commands](mrb/README.md), [independent analytic audit](mrb_independent/MRB_INDEPENDENT_AUDIT.md).

## Verification and immutable scope

The first checkpoint's [Jacobian manifest](jacobian/MANIFEST.json) remains frozen. Its producer and independent programs have 37 and 46 explicit checks; the latter rejects seven corruptions and compares all 1,799 finite-field targets. It can be replayed with `python3 -B jacobian/verify_frozen_packet.py`, normally or with `-O`.

The second checkpoint's manifest and replay guide are linked in [verification](VERIFICATION.md). All new executable acceptance gates remain active under Python optimization. The MRB precision receipts differ only in interpreter-mode/timing metadata before exact normalization; they are not advertised as raw-byte-identical. The original inherited MRB assertion-based scripts are archived as historical evidence, not optimized-mode validators.

The universal statements are supported by written proofs, not extrapolated from finite tests. No new Isabelle/Lean build or full ETH proof replay was run. Documentation checks cover owned bytes, links and unchanged startup objects; the repository-wide checker was not run because it opens the restricted central ledger. No central ledger was read or modified, no empty theorem ID was reserved, and all unrelated prior packets are preserved.
