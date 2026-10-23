# Additional C5, Paley7, W5 phase examples

Verdict: PASS. `independent_phase_examples.py` independently reconstructs the standard C5, Paley7 and W5 orientations, propagates positive parity equations by relation-graph BFS, cross-checks their ranks by Gaussian elimination, enumerates the remaining phases, and tests every homogeneous subset.

It exactly reproduces `session36-context/REALIZATION_PHASES_EXAMPLES.json`, including clause types: C5 and Paley7 have one positive phase and exactly two realizations, with all negative clauses tautological; W5 has four positive phases, two realizations, four tautologies, one binary XOR and two genuine ternary NAE clauses. All three are indecomposable. Normal and optimized outputs are identical.

The definitions, direct results and verification status are saved in `phase_examples_receipt.json`.
