# Pre-registered predictions (klein-2026-11-16-monoid)

Written to `prediction_monoid.txt` at 2026-11-16 11:49 MST (P-1–P-6) and before the second directive's runs (P-7–P-9).

| # | prediction | verdict |
|---|---|---|
| P-1 | 33 first at order 6, exactly 1 self-converse class + 2 converse pairs, absent at order 5 | **confirmed** (CONVERSE.md) |
| P-2 | at order 6 some counts only by pairs, some only by self-converse classes | **confirmed** (ten counts only by pairs; 1, 17, 41, 45 only self-converse) |
| P-3 | products of 3, 5, 11, 13, 17, 257 hit every odd residue mod 2^a 3^b 5^c (a ≤ 8, b ≤ 5, c ≤ 4) | **refuted**: with all three primes present the six numbers generate an index-2 subgroup (e.g. 7, 19, 41 missed mod 120); the fifteen generators of THM-5258 (A) do hit everything (replay to 2·10^5) |
| P-4 | the realised numbers 3..199 generate the 2-part Q′ for primes ≤ 63 except 7; realised primitive roots mod p² below 100 exist for every p | **refuted as stated** (numbers below 64 are never units modulo the product of all primes ≤ 63 — the generators must be coprime to the modulus); the corrected check with the counts ≤ 4000 of prime factors ≥ 67 gives full Frattini ranks at every q (THM-5258 (A)); primitive roots mod p² below 100 exist for every p ≤ 63 (`residues.out`) |
| P-5 | indecomposables of odd ∖ {7, 21} below 10^5 = primes ≠ 7, 7p (p ≥ 5), 49, 63, 343 | **confirmed** exactly (11,266 indecomposables) |
| P-6 | no root-moving move r ↦ 2r + 1 at the same charge on the six ready seeds with ≤ 2 added vertices | **confirmed**: 0 hits over 83 seed representatives (PORTMOVE.md) |
| P-7 | pruned accelerated tree: full binary, lengths in {1,2,3,4}, ray 2/4, branch 1–4, no two consecutive sterile forks, labels recovered (depth 18) | **confirmed** (TOWER_TREE.md §2; lengths determined mod 27) |
| P-8 | tower laws: halving length 2(⌊(2m−1)/3⌋ − ⌊(m−1)/3⌋), +2 per +3; three odd-arrow families 4 rows up, 2 longer; first members 1 → 2, 3 → 5, 5 → 8 | **confirmed** (to 10^5) |
| P-9 | P107 on the reconstructed caterpillar/one-parent ready seeds: none | **confirmed**: 0 hits (51 seeds reach no extension of the target count; 32 reach it but never ready at the same charge) |
