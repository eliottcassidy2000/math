# Session137 — cross-parent spectrum coverage

Source baseline: `bb7cb2d1f61a14632389814469317131d897b640`.

## Main mathematical deliverables

1. **`RESULT.md`** — all-order common-child overlap connectivity, radius n-4; exact two-ear seam reciprocity and lattice; unbounded fixed-core seam divisibility; order-six unit-atlas hostile and its nonstrong-parent repair.
2. **`PHASE2.md`** — complete order-seven two-star and full vertex-replacement atlases; H147 closed-orbit/deletion-deck obstruction; minimum three-arc repair. The independent literal leaf does not depend on a complete class census.
3. **`SYMMETRY_EDIT.md`** — all-order marked-boundary divisor under r-vertex edits; sharper q-arc divisor; completed C3-substitution zero-rigidity census through every strong quotient of order seven; precise formal-network firewall.
4. **`TRIANGLE_CARD_RIGIDITY.md`** — all-order unique homogeneous pair in every deleted-triangle card, canonical family completion, automorphism/card consequences, and reduction of unit isolation to one scalar zero-rigidity equation.
5. **`AUDIT_STATUS.md`** — exact audit scope, independent versus producer-only distinctions, and replay commands.

The full positive-odd H-spectrum except7,21 remains OPEN. So does zero-rigidity for arbitrary strong quotients. None of the finite atlases or arithmetic divisors is extrapolated into a global interval theorem.

## Consequence objects

* `phase2_order*.json`: complete class component membership and literal labelled move certificates.
* `fullstar_order*.json`: deletion decks and common-card spanning forests.
* `deck_barrier.json`: compact H147/H145 hostile certificate.
* `zero_escape6.out`, `zero_escape7.out`: every deletion-orbit row, exact gcd, represented cut count, and original/zero cut masks.
* `zero_closure_controls.out`: actual versus formal zero states with all protected invariants listed.

Primary code is self-contained Python stdlib plus standard C++17. It does not import research engines from the inherited repository. `direct_child_audit.cpp` intentionally uses the producer's reference capacity field, but reconstructs its complete coefficient set using a separate forward-extension DP on literal children. `deck_leaf_audit.py` uses literal permutations and imports no other engine.

`unit_probe.py` and `exchange_probe.py` are exploratory precursors; `audit.py` is the frozen complete first-phase verifier. Compiled executables are local run products, not mathematical source deliverables.
