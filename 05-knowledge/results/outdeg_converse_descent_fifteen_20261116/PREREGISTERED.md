# Pre-registration ledger — klein-2026-11-16-outdeg

Predictions written in [prediction_outdeg.txt](prediction_outdeg.txt) in three batches, each before the runs it covers
(timestamps in the file). Outcomes:

| # | prediction | outcome |
|---|---|---|
| P-1 | digit-0 characterisation (sink, or sole exit with dead-end port) exact on all tournaments of orders ≤ 8 | **confirmed** (58,650 vertices, 0 mismatches) |
| P-2 | `P_D(2)/h(D)` not always a power of two; an order-5 strong `D` already fails | **confirmed** (3 of 396 are powers of two; order 4 already fails: `(5, 32)`) |
| P-3 | more than half of the strong digit-0 pairs of orders 4–8 have a sink in the base | **confirmed** (279 of 341) |
| P-4 | dead-end descents take ≤ 4 steps at orders ≤ 8; regular tournaments sterile | **refuted** in the first clause (descents of 5–7 steps: 136 cases); second clause true (and now proved) |
| P-5 | square-sum rows of five: internal edge only in the first row; every row pair joined | **refuted** (every row has an internal edge; row pairs (2,5), (3,4), (4,5) have none) |
| P-6 | all 62 strong-`D` counts are near-join values from bases ≤ 8; at least one needs base order 8 | first clause **confirmed**, second **refuted** (all from bases ≤ 7) |
| P-7 | `RT_5` paths: three `Z_5`-orbits of five with patterns (1,1,1,1), (2,2,2,2), (2,1,1,2); double transpositions likewise | **confirmed** |
| P-8 | prime square-sum graph on the 25 primes below 100: no Hamiltonian path, an isolated vertex | **confirmed** (none; 37, 43, 67 isolated) |
| P-9 | square-sum Hamiltonian paths: 0 for `n = 18–22, 24`; 10 for `n = 25` | **confirmed** |
| P-10 | (a) below 125 exactly the seven lack a representation `h(A)h(D) + B_t(A)P_D(2)`; (b) some `v < 125` needs a strong `D` | (a) **confirmed**, (b) **refuted** (the first clause of P-10 as written, "every odd in [23, 125]", was self-contradictory with (a) and is superseded by (a)) |
| P-11 | closure `H(R → near-join) = h(R)·H(near-join)`, 0 failures | **confirmed** (400 cases) |
| P-12 | no strong tournament of orders 4–8 (and 9) has two digit-0 vertices | **confirmed** (0; max 1 at orders 4–9) |
| P-13 | 17, 27 not counts of strong order-9 tournaments; sterile counts in [11, 75] = all odd except 17, 21, 27 | **confirmed** (order-9 strong floor 75) |
| P-14 | the one square-sum edge among the non-supersingular primes below 100 is 61 + 83 | **confirmed** |

Refuted predictions are recorded as refuted; none was revised after the fact.
