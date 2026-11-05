# Pre-registered predictions (opus-2026-11-05-S81), written before the runs

Context: the owner's prompt (merge the attached "Z is Diophantine in Q" paper and the repository's Carmichael themes
into the spectrum programme, prioritising reasoning over verification). Three computations accompany the synthesis.

- Q1 (Carmichael-type multiplicity of strong counts). Let m(h) be the number of strong tournament classes, over all
  orders, with exactly h Hamiltonian paths (finite: H >= Moon's floor). From THM-4744's complete inventories
  (P_d(1) for odd d < 625) the set U = { odd h < 625 : m(h) = 1 } is read off. Prediction: U is contained in [3, 100],
  has at most 12 elements, contains 3, 5, 11, 13 and 41, and every element above 15 is a prime or a prime times 3;
  every element of U is realised by a self-converse class (forced by THM-4744's converse pairing P_d = SC_d mod 2).
  Prediction on the multiplicity function: every value m in {1, ..., 10} occurs as m(h) for some odd h < 625
  (the analogue of Ford's theorem that every multiplicity >= 2 occurs for Euler's function).
- Q2 (Diophantine complexity: minimal rank of a complete sub-grammar). Sub-languages of the cutwidth-2 grammar G(0,2)
  of THM-4878, censused to X = 2047 with the S79 engine (closer bound) and their Hankel ranks computed with the empty
  tournament counting one path: mode 1 = no simultaneous double events (every vertex opens at most one slot and
  terminates at most one); mode 2 = mode 1 with first-in-first-out closing (when two slots are pending only the older
  may terminate: crossing chains); mode 3 = mode 1 with last-in-first-out closing (only the newer may terminate:
  nested chains). Predictions: mode 1 attains every odd number <= 2047 except 7 and 21; mode 2 attains every odd
  number <= 2047 except 7 and 21; mode 3 misses at least 10 odd numbers <= 2047 besides 7 and 21; the Hankel ranks
  satisfy rank(mode 2) < rank(mode 1) <= 17 and rank(mode 3) < rank(mode 1).
- Q3 (orbit-counting congruence). P_n(x) = sum over labelled tournaments of x^H and F_{n,p} the sum over tournaments
  fixed by the rotation by n/p: P_n = F_{n,p} (mod p) for every prime p | n (PROVED by orbit counting; the computation
  for n <= 7 is a check, not a prediction).
