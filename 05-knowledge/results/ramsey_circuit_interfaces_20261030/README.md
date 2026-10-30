# Exact Ramsey capacities, Pascal addresses, and shared-circuit boundaries

2026-10-30. **PROVED scoped / CITED / FINITE-EXACT / CONDITIONAL external
algorithm application.** Same-agent mathematical review with separate
verification algorithms. No claim of literature priority or new unrestricted
Ramsey number, circuit lower bound, or A000568 term.

The main result closes the weighted-template target from the earlier
[open-problem atlas](../tournament_open_problem_routes_20261029/PROOFS.md).
For the hereditary substitution closure of the carousel tournaments
T5,T7,T9, the greatest possible order with transitive number at most k is
EXACTLY

    F(k) = sum_{0<=j<k} 2^popcount(j) = A006046(k).

This is an all-size theorem. A finite bank of1344 joint parity charts
splits every integer attachment configuration into three smaller-budget
configurations with its objective exactly preserved. The direct search
through64 tested18,855,887 candidates; the theorem reduces an exact
single-budget query to O(log(k+1)) arithmetic operations. The formula
also follows for the smaller T5-only palette.

The extremal construction is literal: odd Pascal cells (j,t), j<k, are
vertices; their least significant differing bit-pair determines a cyclic
orientation. Its largest transitive subtournament has exactly k vertices,
and increasing k preserves all old arcs. Its full automorphism groups
and colored-orbit counts satisfy exact binary-budget recurrences. At
budgets1..7 the full group orders are1,3,3,81,27,81,729. These are
proved group laws, with native full-group controls through19 vertices.

## Results and boundaries

| Status | Result | Location |
|---|---|---|
| PROVED; classical exponent credited | Complete12-vertex weighted T5 polytope, sharp exponent and equality supports | [R1](RAMSEY.md#r1-the-complete-five-slot-weighted-polytope) |
| PROVED; finite certificates for T7/T9 | A finite parity compiler for all integer weights | [R2](RAMSEY.md#r2-a-parity-certificate-for-every-integer-weight-configuration) |
| PROVED | Exact capacity, nested Pascal-cell family, logarithmic arithmetic evaluation | [R3](RAMSEY.md#r3-exact-capacity-and-a-logarithmic-arithmetic-algorithm) |
| PROVED | Two-sided alpha/omega bound for the C5 graph-substitution class | [R4](RAMSEY.md#r4-a-two-sided-ordinary-graph-companion) |
| PROVED; FINITE-EXACT controls | Full automorphism and color-orbit recurrences for the extremal family | [R6](RAMSEY.md#r6-full-symmetry-and-decoration-recurrences-on-the-extremal-family) |
| PROVED, standard reduction | Parsimonious B2 shared-circuit to3-CNF compiler with exact size accounting | [B1](CIRCUITS.md#b1-a-parsimonious-compiler-preserving-shared-wires) |
| CONDITIONAL | Primary-bit rate transfer for the reported counting engine | [B2](CIRCUITS.md#b2-preserve-the-primary-bit-clock-in-an-algorithm-to-lower-bound-claim) |
| PROVED obstruction | Exponential future-response rank with a linear B2 circuit | [B3](CIRCUITS.md#b3-exponential-alignment-memory-can-coexist-with-linear-circuits) |

Liu already proved the exponent log_3(2) for the larger U5-free class.
A006046 and its binary maximum recurrence are also classical. The
integer attachment proof and certificate/compiler application are the
refinements pursued here. Paley7 refutes extending our exact capacity
formula to arbitrary tournaments. U5 refutes necessity-free use of our
sufficient parity compiler; that failure is not a capacity counterexample.

The twelve-vertex81/27 pair has the same547 transitive subsets, including
empty, and transitive number6. It therefore agrees on every named
weighted Ramsey query while retaining different full symmetry groups.
This distinguishes a sufficient query observer from full metagraph data.

## Circuit Breaker and the attached paper

The challenge uses unrestricted shared B2 circuits; [CIRCUITS.md](CIRCUITS.md)
audits that model and the attached Hexagon paper. It inherits the repo's
already successful pinned Lean build rather than repeating it. A
subexponential bound measured in auxiliary CNF variables is not
automatically such a bound in the circuit's original inputs. No external
submission or message was sent.

## Reproduction

From the repository root:

```text
python -X utf8 -B 05-knowledge/results/ramsey_circuit_interfaces_20261030/check.py
python -X utf8 -O -B 05-knowledge/results/ramsey_circuit_interfaces_20261030/check.py
python -X utf8 -B 05-knowledge/results/ramsey_circuit_interfaces_20261030/audit.py
python -X utf8 -O -B 05-knowledge/results/ramsey_circuit_interfaces_20261030/audit.py
python -X utf8 -B 05-knowledge/results/ramsey_circuit_interfaces_20261030/seal.py
```

The producer checks15,582 local obligations and seven inherited hostile
obligations; the separate auditor checks35,592 obligations without
importing producer code. [AUDIT.md](AUDIT.md) specifies universes and
limitations; [RECOVERY.md](RECOVERY.md) records inheritance, incoming
work, source pins and the live research portfolio.

Next decisive targets are a uniform parity rule for every odd carousel,
the integer U5-free extension, and a replacement recurrence after adding
Paley7. Actual complement closure and Circuit Breaker's gate-sharing
obligation remain separate open problems.
