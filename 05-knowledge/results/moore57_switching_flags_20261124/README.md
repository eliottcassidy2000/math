# Moore(57,2): transport obstructions and exact flag repair

**Status:** Moore(57,2) remains **OPEN** as of 2026-11-24. The target is degree 57, diameter 2, girth 5, and 3250 vertices, not group rank 57.

| Object | Takeaway |
|---|---|
| [Root and edge charts](TRANSPORT.md) | Exact matching and two-port equations. The 56-sheet transport must be primitive and nonsolvable; its centered invariant dimensions are all multiples of five. Every regular group lift is excluded, not just the cyclic case. |
| [Gewirtz / PG(2,4)](GEWIRTZ.md) | The collineation-only branch fails an incidence denominator; mixed transport must pass a ternary sign test and a sharp-set composition bound. |
| [Local Smith and lattice tests](LOCAL-SMITH.md) | A radius-two decomposition selects the second of the two published 5-primary critical-group alternatives, conditional on existence. Includes concrete modular rank and discriminant conditions. |
| [Partial flags](FLAGS.md) | Exact midpoint inverse for the generic two-step paths of the linked family; forced edges and a diameter-two thinning to degree 50 at 1085 vertices. Sharp in the no-same-plane repair class. |
| [Replay](flag-check.py) | Run `python3 flag-check.py`; constructs and tests the entire 1085-vertex graph, with no random steps or external dependencies. |

## Inheritance and why the analogy helps

- The [Fermat/Pépin degree obstruction](../s58_fermat_pepin_collatz_20261026/README.md) uses spectral integrality to leave degrees 2,3,7,57. The same idea now works *inside the fiber module*, even without assuming a graph automorphism.
- The [typed 57 atlas](../number_bridges_57_27_16_20261028/ATLAS.md) distinguishes Moore degree from the 57 points of PG(2,7). The tempting 56-sheet projective-coordinate lift was tested as an actual permutation action, and it fails for a specific shared partition, not for numerology.
- The [tournament switching and SL2 work](../chanpak_tournament_sl2_20261124/) supplies the analogy of edge transport modulo changing gauges locally. Transport around loops is retained; a single scalar or an orbit count does not retain individual endpoints. The Moore replacement for a path count is entrywise uniqueness itself. There is no new Chan–Pak consequence here.
- The [linked asymptotic result](https://arxiv.org/abs/2608.03965) achieves ratio approaching one at fixed diameter and growing degree. The smaller experiment preserves diameter exactly while thinning edges, but does not reach the equality target at degree 57.
- For current status, the [Faber–Keegan note](https://arxiv.org/abs/2210.09577) identifies the error in the earlier nonexistence claim; [Less is Moore](https://arxiv.org/abs/2610.04815) provides a recent exact-search framework, not a solution. [Ishida](https://arxiv.org/abs/2606.29183) gives a separate preprint argument excluding involutions of the full graph; our sheet-transport obstruction does not assume it.

**Next exact test.** The catalogue leaves the Gewirtz action as a tractable candidate smaller than A56/S56. The collineation-only case is excluded; first test whether a mixed 56-element sharply transitive permutation set exists in each surviving action. If one exists, test whether 56 such stars can be chosen with inverse reciprocity and whether the two-hop families for their other endpoint pairs are still sharply transitive. Testing separate rows cannot substitute for the coupled equation. The edge chart gives a second independent check and the Smith filters give arithmetic rejection conditions.
