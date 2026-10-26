# Concurrent integration and correction lineage

The session began at `fc3768ac5fcd` in an isolated worktree. The original shared
checkout's unrelated changes were left untouched. Startup recovered
`normseven_arithmetic_f1_20261026`, `normseven_descent_moments_20261026`,
`a000568_intransitivity_families_20261026`, Session67 and the HN information
packet before new derivation.

The independent taxonomy work found the false all-order Aut-invariance clause
of THM-4803. Before publication, main advanced through `1886609075` with:

- `d761e4cf0f`: [composition classification](../tournament_composition_classification_20261026/README.md), including the sharp4p threshold, unique order12 hostile fibre, ordered path-cover invariance, Sylow-adapted towers, and marked metagraph composition.
- `c9a8823436`: [cyclic-triple composition](../cyclictriple_composition_20261026/README.md), including complete carousel fibres, weighted flip spectra, and pointed growth.
- `36346d6f19`: orientation/monodromy work with explicit limits on analogy.
- `1886609075`: THM-4804 and HYP-9414 reservations. These are namespaces, not proved dependencies.

The two composition packets independently contain and strengthen the repair
we had prepared. Their current corrected theorem, source notes, mistake entry
and regenerated source hash file are retained unchanged. Our superseded local
patch was preserved outside the deliverable at
`/tmp/norm7-correction-before-integration.patch`; its false pre-repair source
hashes are recorded in [CORRECTION_PREIMAGES.json](CORRECTION_PREIMAGES.json).
That file records lineage, not current source pins. No false statement is
restored to preserve authorship or a duplicate correction.

Our [taxonomy proof and adjacency oracle](taxonomy/TAXONOMY.md) remain an
independent confirmation, also stating the general nonselfconverse-block
version. The wider new arithmetic work does not follow merely from the
shared integer seven: its maps are explicit quadratic restriction, character
transport, and exact evaluation matrices.

The three arithmetic packets line up as follows: the first supplies the
normalization and exceptional-conic equations; the second distinguishes
degree2 and degree4 examples; this session gives the complete finite local
decision rule for every nonzero rational (D,δ) on that specified fibre.
The new rule does not settle global interior points or whole-surface F1.

The final clean rebase also integrated `993de01394`, whose
[THM-4804 source note](../module_shattering_lines_carousel_finite_field_20261026.md)
gives an especially relevant finite-field map. With a chosen normal basis
α,α²,...,α^(2^(m−1)), send a binary decoration w to
`x_w = sum_i w_i α^(2^i)`. Squaring rotates the coordinates, so least period
equals field degree and the surviving rotation stabilizer is
`Gal(F_(2^m)/F_2(x_w))`. Root checked this displayed equivariance and the
resulting orbit/stabilizer interpretation directly. The normal-basis existence
theorem is a classical cited input. This is inherited work, not a new discovery
in our packet.

The incoming commit is labelled a pre-audit checkpoint and its audit section
is unfilled; we do not import its entire note as independently audited. In
particular no F1 surface, multiplication-preserving necklace quotient, or
all-order metagraph connectivity claim is deduced from the normal-basis map.
The chosen basis, the cyclic action and the distinction between a field element
and its conjugacy orbit remain explicit.
