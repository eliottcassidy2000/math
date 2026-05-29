# TournamentH7 — Submission-Ready Summary

Lean 4 + Mathlib4 formalisation of forbidden H-values in tournaments.

## What is formalised

| Theorem | Statement | Source |
|---|---|---|
| `Tournament.H_ne_seven` | For every tournament T (any n), H(T) ≠ 7 | THM-343 |
| `Tournament.H_ne_twentyone` | For every tournament T (any n), H(T) ≠ 21 | HYP-1753 |
| `Tournament.H_ne_sixtythree` | For every tournament T (any n), H(T) ≠ 63 | HYP-1754 |
| `Tournament.H_pos` | For n ≥ 1, H(T) ≠ 0 (i.e., H(T) ≥ 1) | Rédei 1934 (R1) |
| `Tournament.H_ne_two` | For n ≥ 1, H(T) ≠ 2 (corollary of parity) | Rédei 1934 (R2) |
| `Tournament.H_ne_even` | For n ≥ 1, H(T) is never an even number | Rédei 1934 (R2) |
| `Tournament.H_not_in_forbidden_trio` | H(T) ∉ {7, 21, 63} | bundle |

## Empirical validation of all axioms

Every axiom used in this Lean project was exhaustively verified at n ≤ 6
by `04-computation/lean_axiom_validation_s7.py`. Output saved to
`05-knowledge/results/lean_axiom_validation_s7.out`.

Result: **0 violations across all 34,866 tournaments enumerated**.

Additionally, exhaustive n = 7 enumeration of all 2,097,152 tournaments
(`04-computation/audit_n7_exhaustive_s6.py`,
`05-knowledge/results/audit_n7_exhaustive_s6.out`) confirms:
  - H = 7  has 0 occurrences,
  - H = 21 has 0 occurrences,
  - H = 63 has 0 occurrences.

## Axiom audit

The audit at `TournamentH7/Verify.lean` reveals exactly which axioms each
theorem depends on. The axioms split into three classes:

### Class A: Cited mathematical results (from the literature)

| Axiom | Reference |
|---|---|
| `ocf` / `ocf_extended` | Grinberg–Stanley, arXiv:2412.10572, Corollary 20 (2024) |
| `moonMoser` | Moon, J.W. (1966), Canad. Math. Bull. 9, Corollary 2.1 |
| `moonCamion_oddSize` | Camion (1959), C. R. Acad. Sci. Paris 249 |
| `redei_existence` | Rédei (1934), Acta Litt. Sci. Szeged 7 |
| `redei_parity` | Rédei (1934), ibid. |

### Class B: Elementary combinatorics (could be de-axiomatised in Lean)

| Axiom | Justification |
|---|---|
| `alpha_subset_bound` | Each k-indep set has k distinct vertices |
| `alpha_chain_step` | Each k-indep set has k distinct (k-1)-indep subsets |
| `alpha_binomial_bound` | Distinct k-indep sets are distinct k-subsets |
| `oddCyclesIn_size3` | A 3-vertex SC tournament is a 3-cycle |
| `oddCyclesIn_size4` | The unique 4-vertex SC tournament has score (1,1,2,2) and 2 three-cycles |
| `omegaTriangleLocalises` | Cycles are contained in single SCCs (basic SCC theory) |
| `omegaCliqueLocalises` | Same, for clique on k cycles |
| `oddCyclesIn_upper` | Trivial bound by 2^|S| |

### Class C: Structural axioms with computational citation (project-specific)

| Axiom | Justification |
|---|---|
| `no_alpha_10_0` | Exhaustive n ≤ 7: no tournament has α₁=10 ∧ α₂=0 |
| `no_alpha_8_1` | Exhaustive n ≤ 7: no tournament has α₁=8 ∧ α₂=1 |
| `no_alpha_6_2` | Exhaustive n ≤ 7: no tournament has α₁=6 ∧ α₂=2 |
| `no_alpha_4_3` | Exhaustive n ≤ 7: no tournament has α₁=4 ∧ α₂=3 |
| `H_ne_sixtythree_axiom` | Exhaustive n ≤ 7: no tournament has H=63 |

## Honest status report

This formalisation is **proof-modulo-axioms**. The deep mathematical content
(OCF, Moon's theorems, Rédei's theorems) is axiomatised with citations to the
peer-reviewed literature. The chain-of-subsets reasoning is axiomatised even
though it would be straightforward to de-axiomatise once a computable
representation of `alphaCount` is provided.

The H ≠ 21 and H ≠ 63 results have additional structural axioms (Class C)
whose only justification at present is exhaustive computational verification
at n ≤ 7. For n ≥ 8 these remain conjectures. The H ≠ 21 axioms are
considerably more believable than H ≠ 63 because the H = 21 obstruction
admits a structural reduction (via `omegaTriangleLocalises`) to the H ≠ 7
obstruction in the (4,3) sub-case.

## Build

```
cd 04-computation/lean/TournamentH7
lake exe cache get      # downloads Mathlib oleans, ~6 GB
lake build              # builds all targets
```

After build, `lake env lean TournamentH7/Verify.lean` prints the `#print
axioms` audit confirming the closed set of axioms each theorem depends on.

## Honest assessment for the mathematical community

This is a **structurally interesting** but **not fully rigorous** formalisation:

- **Rigorous part**: the arithmetic of which α-vectors satisfy the OCF + chain
  constraints for h ∈ {7, 21, 63} is fully proved in Lean. The arithmetic
  enumeration for H ≠ 7 reduces the problem to exactly one α-vector; for
  H ≠ 21 it reduces to four.

- **Cited part**: OCF, Moon, Camion, Rédei are accepted as axioms from the
  published literature.

- **Conjectural part**: the structural unrealisability of each surviving
  α-vector for h ∈ {21, 63} at large n is currently only verified
  computationally up to n = 7.

A more rigorous follow-up would prove the Class C axioms structurally,
which appears feasible for the (4,3) case (via reduction to THM-343)
but requires more work for the (10,0), (8,1), (6,2) cases.
