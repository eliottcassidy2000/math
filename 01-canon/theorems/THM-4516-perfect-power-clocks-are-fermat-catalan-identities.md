---
id: THM-4516
title: "Perfect-power clocks are Fermat-Catalan identities: for q = p^s an odd prime power, X >= 2, r >= 2, an identity 2^K - q^X = +-m^r is a primitive solution of x^a + y^b = z^c with exponents {K, sX, r} and (THM-4484) makes the mixed shape (K,X) free for y -> y/2, (qy +- m^r)/2; for odd q <= 201, K <= 400 the only such clocks are Catalan 2^3 - 3^2 = -1, the Pythagorean family 2^K + (2^(K-2) - 1)^2 = (2^(K-2) + 1)^2, and the three Fermat-Catalan solutions with a power-of-two term, 2^5 + 7^2 = 3^4, 7^3 + 13^2 = 2^9, 2^7 + 17^3 = 71^2, whose free cycles are those of 3x-49, 9x-49, 7x+169, 13x+343, 71x-4913; 2^3 - 7 zeta_3 = (3 - zeta_3)^2; Beal implies no free shape with all three exponents >= 3"
status: >
  PROVED (elementary, from THM-4484 and the definitions) + FINITE-EXACT (census)
  + CITED (Darmon-Granville finiteness; the ten known Fermat-Catalan solutions; the
  April 2025 survey arXiv:2412.11933v2 of solved signatures). NOT independently audited.
  (1) Dictionary. q = p^s odd, X >= 2, K >= 1, r >= 2. If 2^K - q^X = +-m^r then
  gcd(m, 2p) = 1 and the identity is a primitive solution of the generalised Fermat
  equation with exponent multiset {K, sX, r}; if 1/K + 1/(sX) + 1/r < 1 it is one of
  finitely many (Darmon-Granville). By THM-4484(1) the shape (K,X) is free for
  T_{q, +-m^r}: every necklace of shape (K,X) is an integer cycle of qy +- m^r.
  Conversely every known Fermat-Catalan solution with a pure power of two as a term
  (1 + 2^3 = 3^2, 2^5 + 7^2 = 3^4, 7^3 + 13^2 = 2^9, 2^7 + 17^3 = 71^2) is such a clock.
  (2) Census, odd q <= 201, 1 <= X <= K <= 400, X >= 2: |2^K - q^X| in {1} u {m^r, r >= 2}
  exactly for (q,K,X) = (3,3,2) [-1]; (2^(K-2)+1, K, 2) [-(2^(K-2)-1)^2, K = 4..9 in range,
  the trivial family (a+1)^2 - (a-1)^2 = 4a]; (3,5,4) and (9,5,2) [-7^2]; (7,9,3) [13^2]
  and (13,9,2) [7^3]; (71,7,2) [-17^3]. For q = 3 alone, K <= 3000 adds nothing with X >= 2.
  Hostile extension: all odd q <= 10^4, 2 <= X <= K <= 60: exactly 17 clocks, 13 in the
  Pythagorean family (K = 3..15; Catalan 2^3 + 1 = 3^2 is its K = 3 member) and the four
  Fermat-Catalan readings (3,5,4), (7,9,3), (13,9,2), (71,7,2). Nothing new; likewise
  all odd q <= 10^5, K <= 40 (20 clocks = 16 family members K = 3..18 + the same four).
  (3) Free cycles verified: 3x-49 {65,73,85,103,130}; 9x-49 least 11, 13; 7x+169 nine
  primitive 9-cycles (least 67,71,79,85,93,95,109,121,137) plus 169 x {4,2,1}; 13x+343
  least 15, 17, 29 (and one scaled by 7); 71x-4913 least 73, 75, 79; 5x-9 {7,13,28,14};
  17x-225 least 19 (two scalings).
  (4) 2^9 - 7^3 = (2^3 - 7) Phi_3(8,7) = 1 x 13^2 and 8 - 7 zeta_3 = (3 - zeta_3)^2 in Z[zeta_3]:
  the twisted clocks of a fair 3-split (THM-4515) of a 7x+169 cycle are Eisenstein
  squares above 13, and the untwisted factor is Gersonides' unit 2^3 - 7 = 1 (the
  trivial cycle of 7x+1). Six of the nine primitive 7x+169 cycles admit a fair 3-split.
  (5) Beal shadow: a perfect-power clock with K, sX, r >= 3 would be a Beal
  counterexample, so Beal implies no map py +- m^r (r >= 3) has a free mixed shape with
  K >= 3, sX >= 3, X >= 2. The statement "x^4 + y^3 = z^17 has no solution with
  xyz != 0, gcd(x,y) = 1" is NOT in the solved table of the April 2025 survey (among
  (3,4,n) only n = 4, 5 are solved; smallest open Beal signature (3,5,7)); recorded as
  UNVERIFIED. Its clock shadow is empty: no term can be a pure power of two while the
  other two are powers of one odd base.
source: collatz-necklace-20260929 session (mac-mini), 2026-09-29; owner seed: the (4,3,17) statement
depends_on:
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md (shift criterion: a mixed shape is free iff (2^K - q^X) | d)
related:
  - 01-canon/theorems/THM-4515-fair-consecutive-splits-of-cycle-necklaces-are-circulant.md (the Eisenstein-square reading of 7x+169)
  - 05-knowledge/hypotheses/HYP-3058-lrc14-fermat-catalan-hyperbolic-reciprocal-bound.md (earlier Fermat-Catalan mention, unrelated target)
note: 05-knowledge/results/collatz_necklace_20260929_fair_splits_power_clocks_basins.md
scripts:
  - 04-computation/experiments/collatz_necklace_20260929_power_clocks.py
  - 04-computation/experiments/collatz_necklace_20260929_circulant_check.py (section 4)
outputs:
  - 05-knowledge/results/collatz_necklace_20260929_power_clocks_q201_K400.out
  - 05-knowledge/results/collatz_necklace_20260929_power_clocks_q3_K3000.out
  - 05-knowledge/results/collatz_necklace_20260929_circulant_check_d100.out
  - 05-knowledge/results/collatz_necklace_20260929_power_clocks_wide_q1e5_K40.out
  - 05-knowledge/results/collatz_necklace_20260929_power_clocks_wide_q1e4_K60.out (script 04-computation/experiments/collatz_necklace_20260929_power_clocks_wide.py)
extra_sha256: 757328ca911a56e5d31c29fb22edcc0726e87ffeadfff2beb3239fc85c3a2537 (power_clocks_wide.py), 81b10782f3802fda678888d9b2c0902eb2e936a342de5cc7556c5afc927e8b1c (wide_q1e4_K60.out)
script_sha256: 7dec4393a6248b3be63dd4bd65b28e9988da66c882ae571b7b60f54dc3c95798 (power_clocks), 8e8272b3b2c5027127185721dc99d17ded7f25531eee1aba0542912666108fe0 (circulant_check)
output_sha256: 8288e05f57f1a566191d6a30732645b6ca3d40ce885cf4e234d2f4387f94a1f2 (q201_K400), 485dc20dbead1ceaa1e114dfea28444fb34cf62edaf9c25c4ab75260629182a0 (q3_K3000), b8885a9ad4819611a59f23affbf954d23adc2450353ba684781e6d374bca1d01 (circulant_check_d100)
hash_basis: raw LF bytes
audit: NOT independently audited; perfect-power tests by gmpy2.is_power/iroot; every listed cycle recomputed by iteration.
---

# THM-4516 -- perfect-power clocks are Fermat–Catalan identities

**PROVED + FINITE-EXACT + CITED; not independently audited.** Full note: [collatz_necklace_20260929_fair_splits_power_clocks_basins](../../05-knowledge/results/collatz_necklace_20260929_fair_splits_power_clocks_basins.md), section 2.

## 1. The dictionary

THM-4484: for `T_{q,d}(y) = y/2, (qy+d)/2` a mixed shape `(K,X)` is free (every necklace of that shape is an integer cycle) iff the clock `2^K - q^X` divides `d`. Take `d = +-|2^K - q^X|` itself. If the clock is a perfect power, `2^K - q^X = +- m^r` with `r >= 2`, and `q = p^s` is a prime power with `X >= 2`, then `2^K = q^X +- m^r` is a three-term identity in perfect powers with pairwise coprime bases (`m` odd since `q^X` is odd; `p` does not divide `m` since `p` does not divide `2^K`): a primitive solution of the generalised Fermat equation `x^a + y^b = z^c` with exponents `{K, sX, r}`. Darmon–Granville: finitely many for each signature with `1/K + 1/(sX) + 1/r < 1`. The Fermat–Catalan conjecture (ten known solutions) predicts the complete list.

## 2. Census

Odd `q <= 201`, `1 <= X <= K <= 400`, `X >= 2`:

| identity | `q` | `(K,X)` | clock | free cycles |
|---|---|---|---|---|
| `3^2 - 2^3 = 1` (Catalan) | 3 | `(3,2)` | `-1` | `{-5,-7,-10}` of `3x+1` (THM-4484) |
| `2^K + (2^{K-2}-1)^2 = (2^{K-2}+1)^2` | `2^{K-2}+1` | `(K,2)` | `-(2^{K-2}-1)^2` | `q = 5, 9, 17, 33, 65, 129` in range; trivial family `(a+1)^2 - (a-1)^2 = 4a` |
| `2^5 + 7^2 = 3^4` | 3 / 9 | `(5,4)` / `(5,2)` | `-49` | `3x-49`: `{65,73,85,103,130}`; `9x-49`: least `11, 13` |
| `7^3 + 13^2 = 2^9` | 7 / 13 | `(9,3)` / `(9,2)` | `169` / `343` | `7x+169`: nine primitive 9-cycles; `13x+343`: least `15, 17, 29` |
| `2^7 + 17^3 = 71^2` | 71 | `(7,2)` | `-4913` | `71x-4913`: least `73, 75, 79` |

Nothing else occurs; `q = 3` to `K <= 3000` adds nothing with `X >= 2`, and the hostile extension to all odd `q <= 10^4`, `K <= 60` finds exactly the 13 family members `K = 3..15` (Catalan is the `K = 3` member, `2^3 + 1^2 = 3^2`) and the four Fermat–Catalan readings. Exactly the four known Fermat–Catalan solutions with a pure power of two appear (each prime-power base gives one reading, so `3^4 = 9^2` and `7^3`, `13^2` each appear twice), and `2^5 + 7^2 = 3^4` is the `K = 5` member of the Pythagorean family, hyperbolic only because `9` is a square. The `X = 1` rows (`2^K = q + m^r`; e.g. `2^7 = 3 + 5^3`, the cycle `{1,64,32,16,8,4,2}` of `3x+125`) are one-odd-step shapes and are not counted.

## 3. The Eisenstein square

`gcd(9,3) = 3`, so by THM-4515 the clock of `7x+169` factors along a fair 3-split: `2^9 - 7^3 = (2^3 - 7) Phi_3(2^3, 7) = 1 x 169`, and the twisted clocks are `2^3 - 7 zeta_3^{+-1} = (3 - zeta_3^{+-1})^2`, squares of the Eisenstein primes above `13` (`N(3 - zeta_3) = 13`; verified `8 - 7z == (3-z)^2 mod Phi_3(z)`). The untwisted factor `2^3 - 7 = 1` is Gersonides' unit for `q = 7` (the trivial cycle `{1,4,2}` of `7x+1`). Six of the nine primitive `7x+169` cycles admit a fair 3-split (the two-thirds law of THM-4515), and for them the DFT identities read `N_k (8 zeta^k - 7) = 169 C_k` with `8 zeta^{+-1} - 7` a unit times an Eisenstein square.

## 4. Beal's shadow and the `(4,3,17)` statement

A perfect-power clock with `K, sX, r >= 3` would violate Beal's conjecture; hence Beal implies that no map `py +- m^r` with `r >= 3` has a free mixed shape `(K, X)` with `K >= 3`, `sX >= 3`, `X >= 2`. The examples `13x+343` and `71x-4913` have a square term and are Fermat–Catalan, not Beal. The owner's statement that `x^4 + y^3 = z^17` has no solution with `xyz != 0`, `gcd(x,y) = 1` is not among the solved signatures of the April 2025 survey (arXiv:2412.11933v2, Table 1.1; among `(3,4,n)` only `n = 4, 5`; smallest open Beal signature `(3,5,7)`), and a 2026-09-29 search found no later paper: **UNVERIFIED**. Its clock shadow is empty (no term of `x^4 + y^3 = z^17` can be a pure power of two with the other two powers of one odd base), so nothing in this repo depends on it.

## 5. Boundary

Composite `q` with two or more prime factors gives Fermat–Catalan identities only after regrouping and is outside the census; a census over all odd `q <= 10^4`, `K <= 60` is the cheap hostile probe of the Fermat–Catalan conjecture in this window. No statement about `3x+1` cycles is made beyond THM-4484.
