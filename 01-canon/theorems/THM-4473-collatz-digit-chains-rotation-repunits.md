---
id: THM-4473
title: "Collatz's digit chains are exact Markov laws (last digit: 3 -> 5 certain, all other entries 2^j/15), unlike the vanishing consecutive-prime digit bias; repunit primes are the multidigit prime fixed points of digit rotation; the parity-vector map has the odd 2-cycle {-1/5, 5/7}"
status: >
  PROVED + INDEPENDENTLY AUDITED (items 1-4); FINITE-EXACT (the Q-search
  bound in item 4). Under Haar measure, consecutive odd Syracuse terms
  U(n) = (3n+1)/2^(v_2(3n+1)) form an exact Markov chain modulo 10 and 9,
  and they are i.i.d. modulo 3.
  (1) The last digit 3 always goes to 5. Every other row of the last-digit
  matrix is a permutation of (1,2,4,8)/15; for example 9 -> 9 has
  probability 8/15.
  (2) Mod 3 every row is (0, 1/3, 2/3). The mod-9 stationary law on
  (1,2,4,5,7,8) is (8,16,11,4,2,22)/63.
  The bias is exact and does not decay with size, because v_2(3n+1) is
  geometrically distributed at every scale. By contrast, the Lemke
  Oliver-Soundararajan consecutive-prime last-digit bias (reproduced
  exactly, and extended to 1e11) is a vanishing second-order effect.
  (3) Rotating the k base-b digits of n is multiplication by b modulo
  b^k - 1. Its fixed points are exactly the repdigits c R_k, so the repunit
  primes are exactly the prime fixed points when k>=2, and every other
  prime with k>=2 digits lies on a free orbit of size k (e.g. {337, 373, 733}).
  For k=1 every digit is fixed; the primes below the base are exceptions.
  (4) In Collatz, T^k(2^k - 1) = 3^k - 1 and T^(k+1)(2^k - 1) = (3^k - 1)/2.
  The Bernstein-Lagarias parity-vector map Q has odd fixed points -1 and
  1/3, the 2-adic repunits of bases 2 and -2. Q(-1/5) = 5/7 and
  Q(5/7) = -1/5, a second odd 2-cycle besides {1, -1/3}; -1/5 is the
  base-6 repunit limit. Among odd-numerator rationals p/q with |p| <= 201
  and odd q < 100 these are all odd fixed points and 2-cycles.
source: collatz-procgen-20260922 session, digits/repunits lane (2026-09-24), from the owner's prompt on the non-uniform transitions of primes' last digits {1,3,7,9}, circular primes {337,373,733} and repunit primes; audited and promoted by the session orchestrator 2026-09-24
depends_on: []
related:
  - 05-knowledge/results/collatz_procgen_20260924_inverse_tree_mod192.md (sibling ladder, -1/3 = E(0))
  - 01-canon/theorems/THM-4471-kawasaki-fixed-point-collatz-proof-refuted.md (fixed points)
script: 04-computation/experiments/procgen_repunit_20260924_collatz.py
script_dictionary: 04-computation/experiments/procgen_repunit_20260924_dictionary.py
script_audit: 04-computation/experiments/procgen_repunit_20260924_orchestrator_check.py
output: 05-knowledge/results/procgen_repunit_20260924.out
output_audit: 05-knowledge/results/procgen_repunit_20260924_orchestrator_check.out
script_sha256: 6736412e2edb562414e0126ac936cca7642b9a39e2c576fd036d3fa201a1f3a8
script_dictionary_sha256: 7d9c87bfcea4348358438600b5ced3d02bf7ca53fece2ce5e25209f58898e633
script_audit_sha256: c5ca33151fd8b370aa6893e2682ed575d22b140c77b9c3ca446e83332e279aa3
output_sha256: 3e1ae74599ad0354e107b7d44bff96229bcccbbe6e671ef1295e4fac140686b7
output_audit_sha256: 28a2d8d167c0e901d652c2be82fe3c74210e248fa9c1fd38bf74544050403dce
hash_basis: raw LF bytes
audit: >
  The orchestrator derived by hand:
  * 3 -> 5 (3n+1 = 0 mod 5 and U odd);
  * 9 -> 9 iff v_2(3n+1) = 1 mod 4, with probability 8/15;
  * the mod-9 law, from the 6-periodic inverse powers of 2 mod 9;
  * the Q 2-cycle: the parity vector of -1/5 is 1(100)^inf, which reads
    as 5/7, and that of 5/7 is (1100)^inf, which reads as -1/5.
  Independent code (procgen_repunit_20260924_orchestrator_check.py,
  written without reading the lane's scripts) checks:
  * the exact last-digit law, by enumerating n mod 10*2^18;
  * the mod-9 and mod-3 laws with exact fractions;
  * the rotation fixed points for base 10 and k <= 6;
  * the Q fixed points and odd 2-cycles over odd p with |p| <= 121 and
    odd q < 60;
  * T^k(2^k - 1) for k <= 200.
  The lane's full pipeline (primes to 1e11 and circular primes to length
  16) was re-run: all four parts exit 0, and the output agrees with the
  lane's except for timing lines.
---

# THM-4473 -- Collatz's exact digit chains, rotation fixed points, and repunits

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_repunit_20260924_digits_rotation_repunits](../../05-knowledge/results/procgen_repunit_20260924_digits_rotation_repunits.md).

## 1. The last-digit chain of Syracuse terms

Let `n` be odd and Haar-random, and let `v = v_2(3n+1)`, so that `P(v = k) = 2^-k`, independently of `n mod 5`. Then `U(n) = (3n+1)/2^v`.

| from | to 1 | to 3 | to 5 | to 7 | to 9 |
|---|---|---|---|---|---|
| 1 | 4/15 | 2/15 | 0 | 8/15 | 1/15 |
| 3 | 0 | 0 | **1** | 0 | 0 |
| 5 | 1/15 | 8/15 | 0 | 2/15 | 4/15 |
| 7 | 8/15 | 4/15 | 0 | 1/15 | 2/15 |
| 9 | 2/15 | 1/15 | 0 | 4/15 | **8/15** |

*Proof.*
* `U(n) mod 5 = (3n+1) 2^(-v) mod 5`, and `2^(-v) mod 5` has period 4 in `v`, so each residue class of `v mod 4` has probability `2^(-j)/(1 - 2^-4) = 2^(4-j)/15`.
* `U` is odd, so its last digit is determined by `U mod 5`.
* If `n = 3 (mod 10)`, then `3n+1 = 0 (mod 5)`, so `U = 0 (mod 5)` and `U` ends in 5. ∎

**Mod 3 and mod 9.**
* `U(n) mod 3` is never 0. Its next class is `(3n+1) 2^(-v) mod 3`, which depends only on `v`: 1 with probability 1/3 and 2 with probability 2/3. So the classes are i.i.d.
* Mod 9, the next class is `(3r+1) 2^(-v) mod 9`, which depends on `r mod 3` and on `v mod 6`. The stationary law on `(1,2,4,5,7,8)` is `(8,16,11,4,2,22)/63`.

**Contrast with primes.**
* Lemke Oliver–Soundararajan (PNAS 113 (2016) E4446–E4454) found that consecutive primes avoid repeating their last digit. Their tables are reproduced exactly here and extended to `10^11`: `P(1 -> 1) = 0.193` and `P(9 -> 1) = 0.310`.
* That bias decays like `log log x / log x`, with their second-order term explaining 86% of it at `10^11`.
* Collatz's bias is **exact and first-order**: the halving count has the same law at every size, while prime gaps grow like `log x`.
* The "1/15" of the owner's first snippet is the unit of this matrix.

## 2. Rotation fixed points

* Write `n = d_(k-1) ... d_0` in base `b`. Rotating by one place gives `b n - d_(k-1)(b^k - 1)`, i.e. `b n mod (b^k - 1)`.
* The fixed points satisfy `(b - 1) n = 0 mod (b^k - 1)`, i.e. `n = c R_k` with `R_k = (b^k - 1)/(b - 1)` and `0 <= c <= b-1`: the repdigits.
* For `k >= 2`, a repdigit with `c >= 2` is composite because it is `c R_k` with both factors greater than one. Thus the multidigit prime fixed points are exactly repunit primes `R_k`. For `k=1`, every digit is fixed, including the single-digit primes (2,3,5,7 in base10). In base10 the repunit-prime lengths through1100 are exactly `2,19,23,317,1031`. This boundary correction was audited on2026-09-25 in [the digit lane](../../05-knowledge/results/ternary_digits_20260925.md); it leaves the original multidigit calculations unchanged.
* Every other prime with `k` digits has a free orbit of exact size `k`: rotation acts through `Z/k`, and a non-trivial stabilizer would make `n` a repeated block, which a prime cannot be.
* So "repunits versus `{337, 373, 733}`" is exactly "fixed point versus free orbit": the owner's Brouwer reading is REAL.
* The claim that the digit bias *generates* circular primes is NOT SUPPORTED.

## 3. The repunits inside Collatz

* **Mersenne numbers.** The odd branch is the dilation `x + 1 -> (3/2)(x + 1)` about `-1`. Hence `T^j(2^k - 1) = 3^j 2^(k-j) - 1` for `j <= k`, `T^k(2^k - 1) = 3^k - 1`, and `T^(k+1)(2^k - 1) = (3^k - 1)/2`: base-2 repunits run up to twice base-3 repunits.
* **2-adic limits.** The repunits of base `b` converge 2-adically to `xi_b = 1/(1-b)` for even `b`:
  * `xi_2 = -1` is fixed by `T`;
  * `xi_4 = -1/3 -> 0` (not fixed);
  * `xi_(-2) = 1/3 -> 1`.
  * For even `b`, `T` acts on these limits by the Möbius map `b -> -(b+2)/(b-4)`.
* **The parity-vector map** `Q(x) = sum_j (T^j x mod 2) 2^j` is the inverse of the Bernstein–Lagarias conjugacy.
  * It fixes `0`, `-1` and `1/3`.
  * It has the 2-cycle `{1, -1/3}`.
  * It also has the **odd 2-cycle `{-1/5, 5/7}`**. The orbits are `-1/5 -> 1/5 -> 4/5 -> 2/5 -> 1/5 -> ...`, with parity word `1(100)^inf = 5/7`, and `5/7 -> 11/7 -> 20/7 -> 10/7 -> 5/7`, with parity word `(1100)^inf = -1/5`.
  * Whether this cycle appears in the literature after Bernstein–Lagarias 1996 (for example Hotzel 2003) is UNVERIFIED.
* **The lengths of repunit primes.** Treating them (Mersenne, Wagstaff, base-3 and base-10 exponents) as hidden Collatz structure is NUMEROLOGY: 24 tests, none significant. What is real is that the length of the binary all-ones tail of `n` equals the length of its rising run, `v_2(n+1)`.
