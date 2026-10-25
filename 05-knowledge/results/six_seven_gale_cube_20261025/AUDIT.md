# Final audit and reproducibility scope

2026-10-25. **PROVED in the explicitly stated scopes, with independent proof
reviews; FINITE-EXACT receipts for the declared universes.** General theorems
are supported by proofs, not extrapolation from the enumerations.

## Independent reviews

- [Canonical budget and Gale review](AUDIT_BUDGET_AND_GALE.md): general-position
  smoothness, complete-intersection degree, canonical budget, exact matrix map,
  root-scaling fields and minus-one twist.
- [Branch and elliptic review](AUDIT_BRANCH_GEOMETRY.md): component connectedness,
  Riemann–Hurwitz, cube multiplicities, the necessary blowup of the base,
  actual divisor relations, component/intersection fields, rational elliptic
  maps and the direction of the Jacobian isogeny.
- [Companion review](COMPANION_AUDIT.md): character transitions and duality,
  D4 resolution, canonical determinant twist, pullback through a common
  resolution, exact free-group criterion, and normalized fiber-product recovery.

The reviewers did not import producer implementations into their independent
checkers. The draft hashes recorded in those reviews preserve the review
lineage. The final [manifest](MANIFEST.json) pins the integrated bytes.

## Integration changes after those reviews

The research folder was renamed to avoid collision with a concurrent frozen
packet. Reproduction paths and status labels were updated. BG1 and the shared
character, cube and Fano results now explicitly credit parallel derivations.
The geometric proofs and their original producer checkers were not changed.

The draft's Kummer target was superseded by incoming Session60 proofs: the
central root K3 is Kummer over C, the K3 deletion classes split 1+6 with
Picard ranks 20/19, and the hourglass has rank58. The incoming quotient
classification is also routed explicitly. These are inherited audited
conclusions, not new consequences of this packet's receipts.

One scope clarification required an extra exact test. The good-to-bad binary
modified complement in the original combinatorics draft is not the geometric
half-turn under that marking. The integrated note distinguishes the two
permutations and supplies an incidence-adapted marking. The separate reflection
checker tests both markings and the half-turn action on all thirty Fano kernels.
The hostile example remains valid for arbitrary binary relabelling; it makes
no claim that the geometric automorphism destroys freeness.

## Exact universes and controls

| Checker | Universe / independent comparison | Receipt |
|---|---|---|
| `verify_budget_gale.py` | 305 gates; exact two-parameter matrix identities and all 720 permutations at u=v=1; fourth-root quartic hostile | `budget_gale.json` |
| `check_character_decomposition.py` | 4,272 gates; canonical dimensions, determinant parity, code restrictions and abstract exterior cubes | `characters.json` |
| `check_geometry.py` | 10,576 gates; all 252 six-old/one-new transitions among the specified nine lines; cube, conic, rational elliptic identities; fourfold and twist hostiles | `geometry.json` |
| `check_free_deck_codes.py` | All 5,040 binary labelings, 1,395 rank-three subgroups of F2^6, 84 six-subsets, 36 seven-subsets and 240 labelled Paley tournaments; good, all-bad and relabelling controls | `free_deck_codes.json` |
| `audit_character.py` | 39 independent gates; three canonical-dimension calculations for n=4,…,16, local D4/sign equations and exterior-cube dimensions | `character_audit.json` |
| `audit_free_quotients.py` | 333 independent gates; point inertia, all thirty Fanos, 120 complementary pairs and 32 free pairs | `free_quotients_audit.json` |
| `check_reflection_markings.py` | 15,334 gates; two explicit markings, actual triple incidence, all thirty kernels and their bad-triple strata | `reflection_markings.json` |

Each of the seven checkers is run normally and under `python -O`: **14
executions**, exact JSON equality and byte-identical output against the frozen
receipts. Exception-based gates remain active under optimization. Run:

```sh
python3 -B 04-computation/six_seven_gale_cube_20261025/verify_packet.py
```

The wrapper validates SHA-256 before executing. It does not rerun the whole
repository or claim to independently reprove Session60's Picard-rank theorem.
Navigation is outside this packet's frozen manifest because it is maintained
concurrently. The manifest excludes itself and the replay receipt to avoid
circular hashes.

## Limits retained

All geometric statements retain characteristic zero and their stated
incidence hypotheses; arithmetic assertions name their field. A complex
isomorphism does not erase a rational twist. A canonical representation with
no invariants does not prove a free action. The raw ramified fiber product is
not identified with its normalization. Abstract exterior powers of a vector
space are not exterior products of differential forms on a surface.

No theorem about an unrestricted tournament recurrence, all-order black-core
connectivity, LRC(14), planar Jacobian, Collatz, RH, P versus NP, F1 descent,
or distinct integer-square existence is promoted by this packet.
