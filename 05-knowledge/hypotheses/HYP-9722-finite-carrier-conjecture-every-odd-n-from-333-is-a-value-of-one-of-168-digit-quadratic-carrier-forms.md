---
id: HYP-9722
title: "Finite Carrier Conjecture: every odd N >= 333 is a value F_C(X,Y) = qXY + tX + sY + a(p(X)+p(Y)) + n + f(X mod 2) + e(Y mod 2), X, Y >= 1, of one of the 168 distinct carrier forms of strong marked cores of orders 3-5 (p the binary inversion surcharge); with the census below 333 this implies Chan-Pak"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-15-carriers, wave 27 of chanpak_waves_20261113.
  FINITE-EXACT (THM-5192 (D)): every odd N in [333, 10^9] is a carrier value. The 32 odd misses above 21 are all at
  most 331; each is the count of a smaller tournament (census). The order-6 cores alone cover [615, 10^8].
  NUMERICAL: the multiplicity is (4.4-5.0)·ln N, matching the first-moment value 4.78·ln N; the window minimum rises
  4, 4, 8, 14, 18 at N = 10^4..10^8. The orders 3-4 union (log-density 1.667) misses a set of density about
  N^(-0.58), unbounded.
  Pre-registered P89 (to 10^9), P92, P93, P96 CONFIRMED, P94 REFUTED (minimum 18 < 20 near 10^8).
session: mac-mini-2026-11-15-carriers
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md#wave-27-carriers-crossings-and-crowding--a-finite-list-of-digit-quadratic-forms-covers-333-to-108
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/unioncar.c (coverage, multiplicities), greedy_cores.c
    (basis), closedform.c (the form), verify_union.c (direct DP)
related:
  - THM-5192 (A), (D); HYP-9702 (two-block completeness, which this sharpens); THM-5182 (pair expansion);
    chanpak_absolute_steering_20261115 and chanpak_viability_shape_20261115 (the fixed carrier, its O(N^(3/8))
    brackets and first moment)
---

# HYP-9722 - the Finite Carrier Conjecture

**The forms.** Let p(0) = p(1) = 0, p(2u+1) = p(u) and p(2u) = p(u) + u + (u mod 2). For a marked block with
cut-transfer data (a, q, s, e, t, f, n), the carrier form is

    F(X, Y) = qXY + tX + sY + a(p(X) + p(Y)) + n + f(X mod 2) + e(Y mod 2),   X, Y ≥ 1.

By THM-5192 (A), F(X, Y) is the Hamiltonian path count of an explicit strong tournament: guard, X-tape, block,
Y-tape, guard, and one apex. The strong marked cores of orders 3, 4 and 5 give 4, 16 and 148 distinct forms,
168 in all.

**Conjecture.** Every odd N ≥ 333 is a value of one of the 168 forms.

Below 333 the census supplies every odd count except 7 and 21. So the conjecture is a finite, explicit,
two-variable arithmetic statement that implies Chan–Pak.

**Evidence.**
- *Coverage.* Exhaustive to 10^9 (4.9·10^10 configurations). A greedy basis of 48 forms, chosen on [333, 10^7],
  covers [333, 10^9] as well. The 3-cycle form
  12XY + 6X + 3Y + 3(p(X) + p(Y)) + 5 + 2(X mod 2) + (Y mod 2) covers 84% of it alone.
- *Order 6 alone.* Its 1,950 forms also cover [615, 10^8]. So the phenomenon is not tied to one lucky list.
- *First moment.* Each form contributes about (2/q)·ln N representations, since #{XY ≤ N/q} ≈ (N/q) ln N. The
  union's measured mean, 4.89·ln N near 10^8, matches Σ 2/q = 4.78.
- *No primality effect.* Primes have the multiplicity of their residue class. The only local bias is mod 3.

**Why the threshold is not a Poisson one.** For the order-3–4 forms, Σ 2/q = 1.667 and every class mod 15 has
mean multiplicity 1.3–2.4·ln N. A Poisson model would then give finitely many misses, yet the zero frequency
falls only like N^(−0.58) (from 3.1·10^(−2) per odd number near 10^3 to 1.5·10^(−4) near 10^7). Representations
come in correlated clumps. The extra order-5 forms (counts 9, 11, 13, 15 with q = 84–108) break those
correlations: no miss above 331 to 10^8.

**What a proof needs.** Each form is a binary problem: a hyperbola qXY + tX + sY perturbed by a 2-regular digital
term a·p.
- *Without the surcharge* this is a divisor problem in residue classes mod q. It fails on every N with q·N + c
  prime.
- *With the surcharge* the multiplicative structure is broken. The natural tools are those for digital functions
  along multiplicative structure (Mauduit–Rivat: carry and Fourier properties). The "freeze the carry" moves of
  the two-tape steering are exactly carry-property arguments.
- *Circle method.* With power-saving Fourier bounds for p twisted by linear phases (a quasi-Riemann-type input),
  the circle method gives density one for each form. Binary problems do not give every N this way.
- *The finite union* would need a lower-tail bound: the probability that all 168 forms miss N must be
  O(N^(−1−ε)). Equivalently, the clumping must be quantified.

**Predictions.**
- *P95.* Every odd N in [333, 10^10] is a carrier value of the 168 forms (P89, registered before its run, is the 10^9 case).
- *P96 (CONFIRMED).* The 48-form greedy basis computed on [333, 10^7] still covers [333, 10^8] (no miss); it in
  fact covers [333, 10^9].
- *P97.* The minimum multiplicity on [10^9, 10^9 + 2·10^5] is at least 20.
