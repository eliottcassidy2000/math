---
id: THM-4527
title: "The Syracuse map in the owner's labels M = (A+1)/2 is F(2N) = 3N, F(4j+1) = 3j+1, F(4n-1) = F(n) (an exact microcosm); its descents K(A) = (A - S(A))/2 satisfy 6K + 1 = (2^v - 3) S(A), so every integer m occurs exactly #{v >= 1 : (2^v - 3) | 6m+1, quotient > 0} times (the 'two copies' are the units 2^1 - 3 = -1 and 2^2 - 3 = 1); and Collatz holds iff 7 Phi_T(n) is an integer for every n >= 1, in which case 7 Phi_T(n) = -2^sigma(n) (mod 7) is a non-residue."
status: >
  PROVED (elementary) + FINITE-EXACT checks (F and K for M, N <= 10^6; the identity for odd A < 4*10^6;
  multiplicities for m <= 10^5; the Paley restatement for n <= 20000). INDEPENDENT AUDIT: OWED.
  Non-consequence: nothing toward Collatz (the restatement is exact, the bridge is the q = 1 shadow and is
  2-adic, see the note section 2.3).
source: opus-2026-10-01-S15 (collatz-functional-uniqueness-20261001), twentieth note, answering the owner's corrected map F(2N) = 3N, F(2N-1) = 2N-1-K_N
depends_on: none (elementary)
related:
  - 05-knowledge/results/collatz_paley_bridge_20261001.md (nineteenth note: the Paley codes, the explained coincidence)
note: 05-knowledge/results/collatz_label_map_paley_snarks_20261001.md
scripts: 04-computation/experiments/collatz_label_map_snarks_20261001.py
output: 04-computation/experiments/collatz_label_map_snarks_20261001.out (ALL CHECKS PASSED)
---

# THM-4527 — the Syracuse map in the owner's labels

**Status: PROVED (elementary); independent audit OWED.**
Full statement and proofs: [`05-knowledge/results/collatz_label_map_paley_snarks_20261001.md`](../../05-knowledge/results/collatz_label_map_paley_snarks_20261001.md), §§1–2.

## Statement

Let `S(A) = oddpart(3A+1)` for odd `A`, and label odd numbers by `M = (A+1)/2`.

1. **Label map.** In labels, `F(M) = (oddpart(3M−1) + 1)/2`. It is determined by three rules:
   - `F(2N) = 3N`;
   - `F(4j+1) = 3j+1`;
   - `F(4n−1) = F(n)`, the microcosm, i.e. `S(4A+1) = S(A)`.

   Writing `F(2N−1) = 2N−1−K_N` gives `K_N = (A − S(A))/2` at `A = 4N−3`. This is the owner's
   `K = 0, 2, 1, 4, 2, 10, 3, 9, 4, 15, 5, …`, and it satisfies `K_{2j+1} = j`, `K_{4m} = 5m−1` and
   `K_{4m−2} = K_m + 6m − 4`.

2. **Difference spectrum.** Let `K(A) = (A − S(A))/2` for every odd `A`, and `v = v_2(3A+1)`. Then
   `6K(A) + 1 = (2^v − 3)·S(A)`. So `m` occurs as a descent exactly `#{v ≥ 2 : (2^v − 3) ∣ 6m+1}` times, and
   every negative `m` occurs once, as an ascent.
   - The mean multiplicity of the owner's `K` is `Σ_{v≥2} 1/(2^v − 3) = 1.34367…`, not 2.
   - The "two copies of each number plus one 0" are the ascents and the shallow descents. These come from the
     units `2^1 − 3 = −1` and `2^2 − 3 = 1` (Gersonides).

3. **Paley restatement.** Let `Φ_T` be the parity code of `T = n/2, 3n+1`. Collatz holds iff `7Φ_T(n) ∈ Z` for
   every `n ≥ 1`. In that case `7Φ_T(n) ≡ −2^σ(n) (mod 7)`, which lies in `NQR_7 = {3, 5, 6}`; here `σ(n)` is
   the number of steps to reach 1.
