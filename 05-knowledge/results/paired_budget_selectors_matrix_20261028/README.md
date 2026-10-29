# Paired ROOT budgets, five-interior selectors, and aligned tournament pairs

Session: 2026-10-28, branch `codex/session-paired-budget-20261028`.
Inherited live tip: `b38ec1a1eaf7`. No existing research changes were overwritten.

* **PROVED:** [the complete one-step ROOT bound](ROOT_BUDGET.md),
  `D_pair <= min(S/2+3, S/3+4)`, preserving the actual source charges,
  terminal phases, and integer clocks. The coefficient one third and its
  intercept four are sharp. Multilevel closure and universal entry remain open.
* **PROVED:** [at most five private interiors cannot obstruct](FIVE_INTERIORS.md)
  any proper six-colouring constant on four disjoint port pairs in a planar
  unit-distance graph. The next possible selector has at least six interiors.
* **PROVED / FINITE-EXACT:** [signed-core and Gram refinement](SIGNED_CORE_AND_GRAM.md)
  of the concurrently supplied aligned-pair lift. Its complete order-nine
  image is independently checked against all 191,536 ordinary classes.
  The first signed triangle produces determinant 9 or -7 through Eisenstein
  norms. Two strong K4-derived lifts expose the four-probe kernel directly.
* **REFUTED extension:** the one-step budget cannot simply be iterated on
  child charges or repaired using original source mass alone. `depth_two.json`
  retains exact two-generation hostiles, every original source and both clocks.

[Synthesis and next targets](SYNTHESIS.md) connects the finite weighted
controller, groupoid/centralizer sidecar, signed cycle-cover filtration, and
geometric contact obstruction. It records the incoming complete order-ten
generator without attributing that work to this session.

The executable checks and matching JSON receipts live beside the proofs.
The finite checks are distinguished from the event-order and geometric proofs.

Reproduce from the repository root:

```sh
python3 -B 05-knowledge/results/paired_budget_selectors_matrix_20261028/verify.py
python3 -B -O 05-knowledge/results/paired_budget_selectors_matrix_20261028/verify.py
```

`verification.json` pins the source and receipt hashes. See `AUDIT.md` for
scope, proof checks and the same-agent audit limitation.

## Concept board and inheritance

| Lane/object | Closest mechanism | Hostile or lost coordinate | Operation and decisive test |
|---|---|---|---|
| Anchor: paired inverse ROOT rays | Previous packet's single-ray two-thirds bound | Isolated-ray half bound fails; target-seven paired debt is unbounded | Keep the common cutoff; classify its alternating events; test every phase and odd-head state |
| Niche: paired-port unit-distance selectors | Previous four-interior extension theorem | An abstract critical quotient need not have real unit contacts | Classify the nine-vertex critical quotient; recover circle centres |
| Wildcard: square bit matrices and two tournaments | THM-024 anti-involutions; THM-644 folded anti-symmetry | Independent isomorphism classes erase alignment | Normalize each involution pair by its winning endpoint; enumerate common-label orbits |
| K4–K3,3/F64 bridge | Previous exact four-probe and nine-prime bank | Four probes separate nine ports, not all 64 coefficients | Test the same kernel on actual strong tournament constructions |
| Representation and finite control | Shared-label product and controlled-forgetting method cards | Class counts forget stabilizer multiplicity; scalar averages forget phase | Keep the centralizer index; use a weighted finite grammar only after proving its closure |

The least-used sidecar in the anchor is **the order of cutoffs from both
rays**, rather than their separate charge bounds. In the niche it is the
external unit-circle centre of a triple. In the wildcard it is the chosen
anti-involution and its centralizer, not only the doubled isomorphism class.

Method cards actually applied: `00-navigation/META-PATTERNS.md`, “Retain the
product and shared labels before taking scalar totals” and “Controlled
forgetting and unlabeled quotients require a sidecar.” No broad analogy is
being treated as a transfer theorem.

## Follow-up: 57,27,16 and the signed triangle

[The number-bridge packet](../number_bridges_57_27_16_20261028/README.md)
recovers the same norm57 element in the earlier cyclotomic and conductor
work, relates9/−7 to a marked cyclic-matrix elimination, and explicitly
computes the source amplification lost by determinant-only compression.
It also develops the owner's square-sum P4/P3/P2 comparison and a compatible
golden/Eisenstein generator lift with an exact growing binary clock.
These are typed auxiliary constructions; the multilevel ROOT bound above
remains open.
