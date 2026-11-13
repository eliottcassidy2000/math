# Full-count source returns, growing defect, and six-state carriers

**PROVED IN NOTES / FINITE-EXACT; Chan–Pak remains OPEN.** Extended research session, 12–13 November2026. First published checkpoint: `79a49f99d2469ba0ee902ffaf35a6a4107ef37e9`. Proof statuses record the internal independent reviews below, not a literature-priority claim or external referee certification.

[Synthesis](SYNTHESIS.md) · [Inherited and incoming mathematics](INCOMING.md) · [Concept board](BOARD.md)

The main advance is a uniform source return controlling the **full** Hamiltonian-path count. A five-core family with H=45·2^k+7 has a literal two-arc switch from growing defect6·2^k+1 to1. After two births, a larger named replacement restores defect1 at the exact target135·2^k+23. It retains the same physical root and incoming port. The theorem works coefficientwise for arbitrary nonempty substituted tournament modules; the unit offset requires H(module)=1.

The next uniform return, at405·2^k+71, remains open. Separate constructions nevertheless give indefinite neutral returns and complete infinite-family classifications: physical L/T macros with complete-state diamonds, ordinary tournament fibres modeled by Johnson graphs, and a strong rooted six-step block controller. Their scope restrictions are explicit.

| Result | Proof | Independent review |
|---|---|---|
| All fixed-parent indegree-one port fibres; universal five-core repair | [Full-count proof](counts/PROOF.md) | [Root review](AUDIT_ROOT.md), literal and permutations/cuts checks in `counts/audit_counts.py` |
| Uniform successor, retaining the physical root and port | [Module return](compiler/MODULE_RETURN.md) | [Independent return audit](counts/ROOT_RETURN_AUDIT.md) |
| Alternative sparse return by moving the module | [Module migration](returns/MODULE_MIGRATION.md) | [Root review](AUDIT_ROOT.md), `counts/cross_audit.py` |
| Repeatable neutral macros and complete-state diamonds | [L/T language](returns/NEUTRAL_LANGUAGE.md) | [Independent derivation](counts/ROOT_RETURN_AUDIT.md) |
| Exact ordinary-family enumeration and neutral metagraph | [Johnson fibres and canonical peeling](returns/FIBRE_GRAPH.md) | [Fibre audit](returns/FIBRE_AUDIT.md) |
| Restricted birth ancestry and binary/Fibonacci sublanguages | [Ancestry boundary](compiler/ANCESTRY_BOUNDARY.md), [binary language](returns/BINARY_LANGUAGE.md) | Inherited general ancestry credited; exact source-state controls |
| Minimal six-coordinate linear/affine observer; all-periodic E/B obstruction | [Six-state proof](compiler/SIX_STATE_MINIMALITY.md) | [Independent audit](compiler/INDEPENDENT_AUDIT.md) |
| Integral six-cycles, actual strong rooted return, and near-square boundary | [Dynamics proof](dynamics/PROOF.md) | [Independent audit](dynamics/INDEPENDENT_AUDIT.md) |
| Recovered sixes, corrected quadratic orbit, and Collatz clocks | [History and connection ledger](dynamics/HISTORY.md) | Primary sources and exact controls in the dynamics proof |

The general unique-peeling and automorphism mechanism is inherited; the new enumeration refines the selected constructor family by full H and neutral rewrites. For a fixed strong rooted source of order n, the macro family has exactly2^m ordinary classes at order n+8m, with size series x^n/(1−2x^8). Its exact H distribution and neutral-fibre distances are known. For a rooted triangle all outputs are asymmetric, so ordinary and labeled family counts are both explicit. The strict one-/two-root decomposition additionally gives exactly2681 strong ordinary classes with two indegree-one vertices at order9. This contributes structural refinements to A000568; it does not compute new unrestricted terms. A135244 remains an aliquot-inverse sequence, connected through the inherited prime-square extremum, and the near-square branch supplies no new prime-infinitude result.

## Reproduction and scope

From the repository root:

```text
python -B 05-knowledge/results/cycle_opening_six_return_20261112/reproduce_all.py
python -O -B 05-knowledge/results/cycle_opening_six_return_20261112/reproduce_all.py --optimized
python -B 05-knowledge/results/cycle_opening_six_return_20261112/seal_packet.py --check
```

Python with SymPy is required. Exact finite universes, source dependencies, normal/optimized replay receipts, and normalized-text hashes accompany the packet. The replay retains the strongest radius searches once per mode and explicitly labels lower-radius receipts derived by filtering. It does not confuse those projections with independent searches. All-order claims use the proofs, with independent exact controls.

The largest new ordinary-parent inventory is the inherited complete bank through order8; 225 substitution charts are universal,78 with strong parents. The successor search at881 is bounded by its explicit parent-flip radius and defect filters. These are not whole-metagraph exhaustions. No all-odd source selector, unlimited affine return, or Chan–Pak proof is claimed.
