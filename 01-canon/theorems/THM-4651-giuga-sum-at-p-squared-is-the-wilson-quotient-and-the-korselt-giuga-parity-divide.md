---
id: THM-4651
title: "Giuga's sum lifted to p^2 is the Wilson quotient, and the Korselt/Giuga duality splits by parity: (1) S(p) = sum_{k<p} k^(p-1) == p + (p-1)! == p - 1 + p w_p (mod p^2) (Glaisher, Lerch), so Wilson primes are exactly the primes with S(p) == p - 1 (mod p^2) and Giuga's congruence S(p) == -1 holds modulo p^2 exactly when w_p == -1, i.e. at the near-Wilson primes of defect 1 (only 1277 below 2*10^5); (2) an even squarefree number satisfies Korselt's divisibility (p - 1) | (n - 1) only at p = 2, so Carmichael numbers are odd, all known Giuga numbers are even, and a counterexample to Giuga's conjecture is exactly an odd Giuga number that is Carmichael; an odd Giuga number with index k has omega(n) == k + 1 (mod 2), hence at least 10 prime factors when k = 1; (3) every known Giuga number has n' = n + 1, but g - 1 and g + 1 are twin primes only for 30, 858, 1722; (4) among the 255 Carmichael numbers below 10^8 at most half of the prime factors also satisfy Giuga's condition (188461, 2455921, 18900973 reach 2 of 4)"
status: >
  KNOWN: statement 1's congruence (Glaisher 1900; Lerch 1905), Giuga's criterion (Giuga 1950), the characterization of Giuga
  numbers by n' = a n + 1 (Grau and Oller-Marcen 2011), the joint criterion n == p (mod p^2 (p - 1)), and the OEIS sequences
  A079853 and A250406/A250407. PROVED (elementary): the consequences in statement 1, statement 2. FINITE-EXACT: the Wilson-quotient
  scan to 2*10^5 (PARI/gp), the direct check of statement 1 for all primes below 3,000, statement 3 (BPSW probable primality
  above 2^64), and the census of statement 4 (255 Carmichael numbers below 10^8, the known count). COINCIDENCE: the twin-prime
  flanks of 30, 858, 1722. Giuga's conjecture remains OPEN. Not yet independently audited.
session: mac-mini-2026-10-14-wilson
source: 05-knowledge/results/d19_parity_shift_synthesis_20261014.md
scripts:
  - 04-computation/experiments/giuga_korselt_20261014/korselt_giuga_duality.py (+ .out for N = 10^7, _1e8.out for N = 10^8)
  - 04-computation/experiments/giuga_korselt_20261014/wilson_quotient_scan.gp (+ .out)
  - 04-computation/experiments/giuga_korselt_20261014/giuga_twins.py (+ .out)
related:
  - THM-4641 (the Wilson / {2,4,8} / Fermat-multiplier dictionary; this theorem adds the Giuga row)
  - 05-knowledge/results/giuga_korselt_defects_20261014f.md and giuga_block_lattice_20261014f.md (codex: signed Giuga defects, the primary-pseudoperfect completions N(N +- 1) whose factors N +- 1 are twin primes, Sondow-Tsukerman)
  - 05-knowledge/results/collatz_mod6_20260921_grand_circuit_typing.md and bernoulli_boundary_20260925.md (the joint criterion; 30 Giuga not Carmichael, 561 Carmichael not Giuga)
  - THM-4644 (opus S33: Fermat and Wieferich primes as circulant-tournament rigidity and orbit explosion)
---

# THM-4651 — Giuga's sum at p², and the parity divide between Korselt and Giuga

## Setting

* `S(n) = Σ_(k=1)^(n−1) k^(n−1)`. Giuga's conjecture: `S(n) ≡ −1 (mod n)` iff `n` is prime.
* For squarefree composite `n` and a prime `p | n`:
  * **Korselt at `p`:** `(p − 1) | (n − 1)`;
  * **Giuga at `p`:** `p | (n/p − 1)`.
* The **Korselt score** `K(n)` and the **Giuga score** `G(n)` count the primes where each holds.
* Carmichael means `K(n) = ω(n)`; Giuga means `G(n) = ω(n)`. Giuga (1950): `n` is a counterexample iff it is both.
* `w_p = ((p−1)! + 1)/p mod p` is the Wilson quotient.

## Statements

1. **The p²-lift of Giuga's sum is the Wilson quotient.**
   * For every prime `p`: `S(p) ≡ p + (p−1)! ≡ p − 1 + p·w_p (mod p²)`. This is KNOWN (Glaisher, Lerch).
   * Hence (PROVED):
     * `p` is a Wilson prime iff `S(p) ≡ p − 1 (mod p²)`;
     * Giuga's congruence holds modulo `p²`, `S(p) ≡ −1 (mod p²)`, iff `w_p ≡ −1`, iff `(p−1)! ≡ −1 − p (mod p²)`. These are the near-Wilson primes of defect 1 (OEIS A250406/A250407);
     * `w_p ≡ 1` iff `(p−2)! ≡ 1 (mod p²)`, iff `S(p) ≡ 2p − 1` (A079853).
   * FINITE-EXACT, PARI/gp over all primes `3 ≤ p < 2·10^5`:
     * `w_p ≡ 0` exactly at 5, 13, 563;
     * `w_p ≡ −1` only at **1277**;
     * `w_p ≡ 1` at 3, 11, 107, 4931, matching A079853.
   * The congruence itself was checked directly, with both sides computed, for every prime below 3,000.
2. **The parity divide (PROVED).**
   * (a) If `n` is even and squarefree, then `K(n) = 1`: Korselt holds at `p = 2` and fails at every odd `p`, since `p − 1` is even and `n − 1` is odd.
     * Hence Carmichael numbers are odd (Korselt). All known Giuga numbers are even.
     * A counterexample to Giuga's conjecture is exactly an **odd** Giuga number that is Carmichael.
     * Checked on all 405,285 even squarefree composites up to `2·10^6`.
   * (b) An odd Giuga number with index `k` (meaning `Σ_(p|n) 1/p − 1/n = k`) has `ω(n) ≡ k + 1 (mod 2)`.
   * (c) Since the reciprocals of the first eight odd primes sum to 0.99896, an odd Giuga number has `ω(n) ≥ 9`, and `ω(n) ≥ 10` if `k = 1`. Much stronger lower bounds for counterexamples are known in the literature (Borwein, Borwein, Borwein and Girgensohn 1996 and later). This one is only the elementary parity version.
3. **Giuga numbers: derivative and neighbours.**
   * All 13 known Giuga numbers (the 12 of A007850 and the 10-factor one recorded there) have index 1, i.e. arithmetic derivative `g' = Σ g/p = g + 1`. KNOWN in general: Giuga iff `n' = a·n + 1` with `a ≥ 1`.
   * `g − 1` and `g + 1` are both prime for `g = 30, 858, 1722` and for none of the other ten. The largest has `g + 1` prime and `g − 1` composite.
   * So "Giuga numbers are twin-prime centres" is a small-number COINCIDENCE.
   * The structural twin-prime mechanism is different. It is the completion `N(N ± 1)` of a primary pseudoperfect `N` (where `N' = N − 1`) by an adjacent prime `N ± 1`: `30 = 6·5`, `1722 = 42·41`, `2214408306 = 47058·47057`. This is Sondow–Tsukerman's construction, recorded in codex's note. In the derivative language, `N' = N − 1` times a prime neighbour gives `g' = g + 1`.
4. **Giuga scores of Carmichael numbers (FINITE-EXACT, `n < 10^8`).**
   * Among the 255 Carmichael numbers below `10^8`, the distribution of `(G, ω)` is `(0,3): 77, (0,4): 121, (0,5): 18, (1,3): 7, (1,4): 20, (1,5): 8, (2,4): 3, (2,5): 1`.
   * No Carmichael number below `10^8` has more than half its primes satisfying Giuga's condition. The maximum share `2/4` is reached by:
     * `188461 = 7·13·19·109` (at 7, 19);
     * `2455921 = 13·19·61·163` (at 13, 61);
     * `18900973 = 7·13·229·907` (at 7, 13).
   * When Giuga holds at a prime of a Carmichael number, it is mostly a small one (`561` at 3, `1105` at 5, `52633` at 7). This matches the heuristic chance `≈ 1/p`, since Korselt already gives `n/p ≡ 1 (mod p − 1)`.
   * 60 of the 255 have `c − 2` or `c + 2` prime. Below `10^6`, `c + 2` is prime for `c = 561, 2465, 656601`, and `563 = 561 + 2` is the third Wilson prime (THM-4641 row 8, COINCIDENCE).

## Proofs

**1.**
* `k^(p−1) = 1 + p·q_p(k)`, so `S(p) = (p − 1) + p Σ_k q_p(k)`.
* Lerch's formula gives `Σ_k q_p(k) ≡ w_p (mod p)`. Equivalently, `(p − 1)! = −1 + p W_p` gives `S(p) ≡ p + (p−1)! (mod p²)`.
* The three consequences are rearrangements. For the last one: `(p−2)! = (p−1)!/(p−1) ≡ (−1 + p w_p)(−1 − p) ≡ 1 + p − p w_p (mod p²)`.
* The scan multiplies `1, …, p−1` modulo `p²` in PARI. The direct check sums `k^(p−1)` modulo `p²`. ∎

**2.**
* (a) is in the statement.
* (b): `n' = Σ n/p` is a sum of `ω(n)` odd numbers, while `n' = k·n + 1 ≡ k + 1 (mod 2)`.
* (c): `Σ 1/p ≥ 1 + 1/n > 1`, and the first eight odd primes give less than 1. ∎

**3–4.** Exact arithmetic (`giuga_twins.py`, `korselt_giuga_duality.py`). Carmichael numbers are found by a vectorised Korselt test over a smallest-prime-factor sieve, and each is re-checked by factoring. ∎

## Reading

* **"A single value which is both."**
  * Korselt says `n ≡ 1 (mod p − 1)`: in the exponent, `n` acts like 1. Giuga says `n/p ≡ 1 (mod p)`: the cofactor acts like 1.
  * Primes satisfy both trivially. The conjecture says no composite does, and the joint condition is `n ≡ p (mod p²(p − 1))` for every `p | n`.
  * Parity separates the two worlds. The Korselt side lives on odd `n`; all known Giuga numbers are even, and their Korselt score is forced to be 1. A counterexample would have to cross the parity divide as an odd Giuga number.
  * In the census the two conditions never overlap on more than half the primes.
* **Linked lists.** Under the arithmetic derivative `n ↦ n'`:
  * primes map to 1;
  * primary pseudoperfect numbers map to `n − 1`;
  * Giuga numbers map to `n + 1`;
  * `p^p` is fixed.
  * Multiplying a primary pseudoperfect number by a prime neighbour moves it from the predecessor list to the successor list. That is the exact sense in which these numbers "embody one less or one more".
* **Wilson at the next level.** Modulo `p` all primes satisfy Giuga's congruence. Modulo `p²` the defect is the Wilson quotient. So the Wilson primes (5, 13, 563) and the defect-1 near-Wilson primes (1277) are the primes where Giuga's sum is exceptional at `p²`. This ties THM-4641's dictionary to Giuga's conjecture.
