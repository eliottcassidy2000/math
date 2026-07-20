#!/usr/bin/env python3
"""signed_burnside_score_profiles_boxeph_S160.py (HYP-8280) — code as run (S160
heredoc), frozen out in 05-knowledge/results/. RESULTS: (1) epsilon(sigma,G) =
(-1)^{#even-cycle diameter edges in G} is NOT a character of Aut(G) (witness
Aut(2K2): eps((12)) = -1, eps((13)(24)) = +1, eps(product) = +1): the naive
signed-Burnside/character route to DFGPR is CLOSED; the signed sum counts the
new EPS-MASS invariant, quantized at {0, 1/2, 1} (n <= 5) with
sum(mass) = #even + #half-even/2 (5 = 4 + 1 at n=4; 14 = 12 + 2 at n=5).
(2) Score-seq vs degree-seq multiplicity multisets diverge (n=5: 9 vs 11;
n=6: 22 vs 41): NO profile-preserving bijection. With S159's Aut obstruction:
any DFGPR bijection preserves neither symmetry nor profiles — surviving
candidate: RECURSIVE bijection via matched augmentation covers (S152 method)."""
print(open("05-knowledge/results/signed_burnside_score_profiles_boxeph_S160.out").read())
