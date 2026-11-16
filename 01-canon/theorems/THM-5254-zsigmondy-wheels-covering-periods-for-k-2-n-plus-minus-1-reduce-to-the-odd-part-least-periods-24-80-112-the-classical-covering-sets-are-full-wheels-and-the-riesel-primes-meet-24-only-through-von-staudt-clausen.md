---
id: THM-5254
title: "Zsigmondy wheels. Let W(N) be the set of primes dividing 2^N - 1, each carrying one residue class modulo ord_p(2); a covering of Z/N by these classes is a covering set of period N for k 2^n +- 1 (Riesel/Sierpinski numbers), and every covering set of period N lies in W(N). (A) PROVED: every dyadic level 2^i has a prime of that order (the prime factors of the Fermat number 2^(2^(i-1)) + 1); nu(2^i) = 1 for i <= 5 and nu(64) = 2. (B) PROVED: for N = 2^a m, m odd, a cover forces sum over non-dyadic p of 1/oddpart(ord_p(2)) >= 1 whenever the dyadic classes do not cover alone, and a cover exists whenever Z/m is covered by one class modulo oddpart(ord_p(2)) per non-dyadic p (nested-dyadic criterion). (C) PROVED: for N = 2^a q, q an odd prime, a <= 5, W(N) covers Z/N iff at least q primes have order with odd part q; the least covering periods are 24 for q = 3 (period 12 fails exactly because 2^6 - 1 = 63 has no primitive prime divisor, Zsigmondy's exception), 80 for q = 5 and 112 for q = 7, and every odd prime is the odd part of some covering period. (D) PROVED: a prime divides both k 2^n - 1 and k 2^m + 1 for one k iff ord_p(2) is even. (E) FINITE-EXACT: base 2, N <= 150, the covering periods are exactly 24, 36, 48, 60, 64, 72, 80, 84, 96, 108, 112, 120, 128, 140, 144 and no other, 132 being excluded by a two-arrangement count, the nested criterion agreeing with every decided period; base 3, N <= 120: 48, 72, 90, 96, 108, 120, base 3 lacking the modulus 2 as base 2 lacks 6; no self-dual (both-sign) wheel among the decided periods below 150. (F) The classical covering sets are the full wheels: Riesel's {3,5,7,13,17,241} = W(24), Selfridge's set for 78557 = W(36) minus {109}, Sierpinski's 1960 set = W(64) with density exactly 1. (G) The primes {3,5,7,13} shared by the Riesel wheel, the prime Coxeter numbers of the pure-A Niemeier lattices and the Leech shell count 196560 = 2^4 3^3 5 7 13 are exactly the odd primes with (p - 1) | 24 (von Staudt-Clausen primes of B_24), forced into every period-24 wheel by (p - 1) | 24 => ord_p(2) | 24; no Leech mechanism enters the covering structure."
status: >
  PROVED (klein-2026-11-16-two3; elementary, not independently audited): (A)-(D) and (F) in the packet's WHEELS.md section 1
  with full proofs; (G) is the von Staudt-Clausen criterion plus the order divisibility. FINITE-EXACT (E): exact-cover search
  on bitmasks (wheels2.py, wheels3.py, brier_sheets.py, exact_one.py; sympy factorisations of 2^d - 1, d <= 150), with the
  five periods the bounded search left undecided settled by the a = 1 identity (90, 126, 150), the density bound (100) and a
  larger search (132). Pre-registered P-W1 confirmed, P-W3 refuted by the base-3 period 108 (packet PREREGISTERED.md).
  CITED: Zsigmondy/Bang (primitive prime divisors; the exception 2^6 - 1), the classical covering sets (Sierpinski 1960,
  Selfridge, Riesel), CLOCKS.md of logarithmic_repairs_covers_20261116 (the two patterns at period 24, the Brier remark),
  THM-484 (involution modulus 24), THM-4838 (von Staudt-Clausen, the "twelve" reading typed NUMEROLOGY), THM-5253 (E).
session: klein-2026-11-16-two3
source: 05-knowledge/results/two_three_wheels_synthesis_20261116/README.md
related:
  - 05-knowledge/hypotheses/HYP-9823-zsigmondy-wheel-conjectures-nested-dyadic-is-without-loss-least-periods-grow-like-sqrt-q-no-self-dual-wheel-and-the-base-b-analogue.md
  - 05-knowledge/results/logarithmic_repairs_covers_20261116/CLOCKS.md
  - 01-canon/theorems/THM-5253-sharp-span-3-growth-lemma-audited-sound-age-structured-caps-every-profile-coordinate-at-most-30-h-history-refined-cones-refuted-and-the-zsigmondy-six-prong-wheel.md
  - 01-canon/theorems/THM-484-involution-modulus-24-source-of-gleason-orders-8-and-24.md
---
# THM-5254 — Zsigmondy wheels

Notation: `W(N) = {p prime : ord_p(2) | N}`, `ν(d) = #{p : ord_p(2) = d}` (the primitive prime divisors of `2^d − 1`;
`ν(d) ≥ 1` for `d ≥ 2`, `d ≠ 6`); a cover of period `N` chooses `r_p (mod ord_p(2))` for `p ∈ W(N)` with
`∪_p (r_p + ord_p(2)·Z) = Z/N`; by the Chinese remainder theorem it yields `k` with `p | k·2^n − 1` whenever `n ≡ r_p`
(and `−k` for the plus sign), and every covering set of period `N` is contained in `W(N)`.

**(A)** `p | 2^(2^(i−1)) + 1 ⇒ 2^(2^(i−1)) ≡ −1 ⇒ ord_p(2) = 2^i`. The Fermat primes 3, 5, 17, 257, 65537 give
`ν(2^i) = 1` for `i ≤ 5`; `F_5 = 641·6700417` gives `ν(64) = 2`. ∎

**(B)** `N = 2^a m`. The dyadic classes (modulo `2^i`, `i ≤ a`) are unions of cosets of `2^a·Z/N`; let `L` be what they
leave. A class modulo `d = 2^i m′` (`m′ | m`, `m′ > 1`) meets each coset in 0 or exactly `|coset|/m′` elements, so it
covers at most `|L|/m′` of `L`; a cover needs `Σ 1/m′ ≥ 1` when `L ≠ ∅` (necessity). Choosing the dyadic classes nested
(`0 mod 2`, `1 mod 4`, …, `2^(a−1) − 1 mod 2^a`) leaves the single coset `C ≅ Z/m`, and by the Chinese remainder theorem
each non-dyadic prime's class can be placed to restrict to any prescribed class modulo `m′` on `C` (sufficiency). ∎

**(C)** With `m = q` prime each non-dyadic class covers exactly one residue of `Z/q` in the nested picture; (B) gives
"iff `c(a, q) = Σ_{i ≤ a} ν(2^i q) ≥ q`" for `a ≤ 5` (where the dyadic classes cannot cover alone). `c(a, 3) = 1, 1, 2, 3`
(`ν(6) = 0`) gives 24; `c(a, 5) = a + 1` gives 80; `c(a, 7) = 1, 2, 4, 5, 7` (`ν(28) = ν(112) = 2`) gives 112;
`c(a, q) ≥ a + 1 − [q = 3]` bounds the least period by `2^(q−1) q` (`2^q q` for `q = 3`). For `a ≥ 6`, `W(64) ⊆ W(N)`
covers. ∎

**(D)** `k·2^n ≡ 1`, `k·2^m ≡ −1 (mod p)` give `2^(m−n) ≡ −1`, so `−1 ∈ ⟨2⟩` and the order is even; conversely
`2^(ord/2) ≡ −1`. ∎

**(E)–(G)** Packet WHEELS.md §2 and §4 (tables, the classical sets, the von Staudt–Clausen identity).

**What is new.** The reduction (B), the exact criterion (C) with the Zsigmondy explanation of the period 24, the
classical-set identification (F), the tables (E) and the typed relation (G). Known and used: Zsigmondy–Bang, the
classical covering sets, CLOCKS.md, THM-484, THM-4838.
