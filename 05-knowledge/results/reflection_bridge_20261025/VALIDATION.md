# Validation and proof review

2026-10-25. **FINITE-EXACT** finite controls, **VERIFIED** independent proof reviews. All-k conclusions are proved in the linked notes, not extrapolated from a census.

Run from the repository root:

```text
python3 05-knowledge/results/reflection_bridge_20261025/validate.py
python3 agents/check_docs.py
```

The packet validator runs five standard-library checkers normally and with `-O`, rejects nonzero exits, and requires byte-identical standard output and receipts. Checks use explicit exceptions, so optimization does not remove their gates. `validation_receipt.json` records the checker and output SHA256 values. No random sampling is used.

Parent replay passed all ten checker executions. The repository documentation gate passed after retaining the 450-line frontier budget, all local packet Markdown links resolved, and `git diff --check` passed.

| Checker | Explicit universe and decisive controls |
|---|---|
| `check_divisor_ideal_graph.py` | Full exponent boxes k=2..6; 2,160 ordered pairs, 100 action states, 496 arbitrary occupancy masks; all 12 k=2 square states. Checks ideal order/gcd/lcm, actual vertex permutations, rank, and the non-lattice quotient hostile |
| `graph_bridge_check.py` | Six k=1 and fourteen k=3 proper divisors; graph and marked-tournament operations, all relevant source comparisons; all 8! candidate bijections in the inherited eight-state gauge, recovering 32 equivariant and two monotone maps |
| `audit_graph_bridge.py` | Independent producer-free adjacency construction for all 14 k=3 states, 196 ordered source pairs, and 140 tournament triples; direct reachability, directed cycles, closure and blue checks |
| `seven_check.py` | Boolean boundary quotients for 1..8 slots, all seven Fano lines, six coordinate permutations, direct squarefree divisors of 60, cyclic-unit projective normalization m=1..10, 125 exact rational covariance samples and the explicit failed second-layer reflection |
| `geometry_check.py` | Rank-zero actions k=2..12; complete normalized finite-field vectors for k=2..5 over F2,F3 and k=2 over F5; 32,959 checks. Positive action/count controls, the seven-point fixed-count mismatch, and 61-versus-46 deletion hostile |

The graph producer and divisor-ideal reviewer audited each other's proofs. The independent graph replay is saved in [independent_graph_audit.md](independent_graph_audit.md). The ideal proof review confirmed the sharp k+2 coordinate budget, actual square permutations, and the ambient-lattice qualification for gcd/lcm. The parent also checked those formulas directly.

The seven-support reviewer independently audited the projective carrier in [geometric_carrier_audit.md](geometric_carrier_audit.md). Its two requested precision repairs are applied: over Z the deleted objects are coordinate sections, and support specialization direction is explicit. Standard projective chart permutations establish geometric existence; finite point counts alone would not.

The only demonstrated error in a maintained predecessor was the path-word-only formula in the Klein dictionary. Its correction is pinned to the already correct THM-4733 Proposition A, with lineage in `01-canon/MISTAKES.md`. This patch does not re-audit the unrelated detour claims in that predecessor.

Research boundaries: no test establishes an unmarked tournament-class bijection, a global metagraph theorem, a Keller tower theorem, or an F1 model for arbitrary arithmetic equations. The seven-point historical leaf counts are provenance; the explicit Boolean quotient and projective carrier have their own declared domains.
