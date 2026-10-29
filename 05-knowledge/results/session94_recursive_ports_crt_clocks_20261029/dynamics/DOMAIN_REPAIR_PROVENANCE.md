# Propagating an already-proved domain correction

Source inspected at `88c969273f2b6bc4ea2a8c5d29bb36c558eb59d2`.

**Stale live statement:**
`05-knowledge/results/fixed_twists_conductors_reverse_potential_20261025/antihydra/ANTIHYDRA.md`, lines39–46: full countable valuation shift from all odd2-adics and an invalid compactness shortcut. Counterexample `x=1/3` has zero numerator/infinite valuation.

**This is NOT a newly discovered correction.** The more recent proved source already supplies it:
- `05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-reverse-valuation/REPORT.md`, section2, equations3–5: delete complete finite preimage trees of the singular numerators; use inverse-branch convergence and cylinder bases.
- Independent audit `.../session64-audit-dynamics/REPORT.md`, lines38–56, explicitly requires the full countable exceptional set, not just one singular point per cylinder.
- Session64 README lines34–37 routes that result.
- `05-knowledge/results/lca_ports_reverse_20261025/arithmetic/README.md`, section3, lines104–117, independently states the same exclusion.

`PROPAGATE_EXISTING_DOMAIN_REPAIR.patch` makes only the local A3 repair and attributes it to Session64. It does not change finite-cylinder masses, Antihydra's rule, its returned counter, any positive-integer result, or the shared mistakes ledger. It was prepared without mutating the source checkout, then applied after the fresh-read check recorded below. The correction's links use the actual current Session64 paths.

**Bounded correction search:** exact Antihydra source, its package README, current THM-4771, Session64 README, main reverse-valuation report, independent dynamics audit, and LCA/ports/reverse arithmetic README. This was enough to find the prior repair; no claim that the entire repository has no other copy is made.

## Applied checkpoint

The parent applied this exact local repair on main at commit
`77683b7ced760e4fff7a3e2a431a417caf46a925`, after confirming the fresh remote
source equalled the pinned text. Old blob `03a5cc3440a6424937f608b01ec8006305ddbead`
was replaced by `246925c762fb7f7fd4599ca5d066db7435da94e5`. Only A3 changed.
Replay the patch against the initial source pin, not the already-repaired HEAD.
The shared mistakes ledger was not modified.
