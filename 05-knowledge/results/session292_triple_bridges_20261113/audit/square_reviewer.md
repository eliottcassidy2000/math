# Independent review of the occupancy and finite-observer lanes

Reviewer: `square_reflections` lane. **PASS, scoped to the files and statements below.** This reviewer authored the square-sum lane and does not independently certify that lane here.

Reviewed `triples/PROOF.md` Theorems 1–5 and `manuscript/PROOF.md` sections 1–5 as present on 2026-11-13. The mathematical arguments are sound in their stated rings/carriers. No substantive correction was needed in these files.

## Occupancy lane

- The destination ideal survives every pivot; entering pivot k already supplies x_k, so multiplication by its self-loop kills the extra inverse term. External source/sink states are not eliminated. This restriction is essential.
- Expanding two pivots yields the five listed simple words plus one repeated-pivot term; squarefree multiplication annihilates precisely that extra term. Permuting larger pivot sets does not discard named occupancy information.
- The divided-power identity is integral: every surviving support collection consists of distinct disjoint nonempty subsets, so its orderings contribute exactly k!. Grouping those supports proves the binomial product law, even when different collections share a union.
- The p-adic cutoff is geometric-power precision, not total support precision. The ordinal sum of triangles gives a valid sharpness family since the adjacency determinant is block upper triangular. Nonzero r^(K−1) also certifies failure of the shorter geometric truncation after multiplication by the unit (1−r)^(-1).
- The hostile example in section 3 correctly blocks a claim of bounded-cycle complexity or physical count coverage. No such overclaim survives in the reviewed text.

Independent code `square_reviewer.py` uses frozenset monomials and direct Hamiltonian permutations, not producer imports. It checks every labelled tournament through order four and **all 1,589 pivot orders** across those 75 inputs, including singleton and self-loop boundary controls. Literal counts equal full-occupancy coefficients. It separately expands determinant denominators for one through seven ordinal triangle blocks and checks factorial sharpness and all-subset coefficient one. The producer's own script additionally exercises larger sampled tournaments, all endpoint/subset coordinates, five-term exchanges and modular cutoffs; its replay passes.

## Finite-observer/source lane

- The symmetric escaping quartet example has the claimed roots and compact bound. It is explicitly not asserted to satisfy zeta's arithmetic or growth properties.
- The stabilization argument for V_j is correct: equality at one step gives invariance under every alphabet matrix; starting with nonzero V_0, at most d−1 strict growth steps occur. This proof requires the fixed-dimensional carrier and fixed finite alphabet, as stated. The product-automaton suggestion must retain every legality guard, also explicitly stated.
- The eight-mode Vandermonde certificate and seven-sample hostile are exact. The independent script directly confirms the seven zeros and first difference 3,991,680. The shifted witness and recurrence also pass the producer replay.
- Differentiating F_S(t) gives the displayed derivative and unique continuous maximum. Its domain is t>0 (thus integer applications n>1), not n=1. The note properly keeps the continuous maximum distinct from integer attainment and records the missing absolute-size coordinate.
- These lemmas do not establish the claimed external RH result, a Lean verification, a physical tournament move, or Chan–Pak coverage; the scope labels say so.

Reproduction:

```
python3 05-knowledge/results/session292_triple_bridges_20261113/audit/square_reviewer.py
python3 05-knowledge/results/session292_triple_bridges_20261113/triples/verify.py
python3 05-knowledge/results/session292_triple_bridges_20261113/manuscript/verify.py
```

No external theorem priority or full repository audit was undertaken. The mathematical proofs, not the finite controls alone, establish the all-order statements.
