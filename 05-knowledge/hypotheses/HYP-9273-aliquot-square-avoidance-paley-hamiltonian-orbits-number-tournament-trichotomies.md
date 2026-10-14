---
id: HYP-9273
title: "(A) Every aliquot cycle avoids squares and twice-squares, i.e. is parity-homogeneous (contains the opposite-parity amicable problem); (B) the Hamiltonian paths of every Paley tournament P_q form a number of affine-group orbits congruent to 1 mod 8"
status: >
  OPEN (both).
  PROVED around (A): the aliquot sum s(n) = sigma(n) - n changes parity exactly at n in the squares or twice the squares.
  Hence every aliquot cycle has an even number of such members, and an amicable pair has one iff it has opposite parity
  (iff both members are; classical, Gioia-Vaidya 1967, rediscovered by D'Abramo 2005).
  Every cycle satisfies sum sigma(n_i) = 2 sum n_i and prod sigma(n_i)/n_i >= 2^k, with equality iff it is a perfect number.
  FINITE-EXACT for (A): all 248 aliquot cycles with every member <= 10^8 are parity-homogeneous and contain no square
  or twice-square. They are 5 perfect, 231 amicable, ten 4-cycles, one 5-cycle and one 28-cycle; there is no 3-cycle.
  Among odd a <= 10^9, 4476 have s(a^2) a square or twice a square. 4458 of these are primes with s(a^2) = a + 1, and none
  closes an amicable pair, so there is no opposite-parity amicable pair with odd member <= 10^18.
  CITED: Chernykh's complete enumeration of amicable pairs with smaller member < 10^20 (all of equal parity; OEIS A066873),
  which already contains our 10^18 range; Hagis 1970, coprime opposite-parity pairs have both members > 10^60;
  Iannucci-Luca 2009, s(n^2) is a square for at most x/(log x)^(3/2+o(1)) of the n <= x (the same for twice a square,
  Pollack 2011). Our double-flip counts (4476 to 10^9, of which 18 composite) sit far below that bound.
  FINITE-EXACT for (B): q = 3, 7, 11, 19, 23 (orbit counts 1, 9, 1729, 6857869865, 62293308207033, all 1 mod 8) and
  q = 27 (orbit count under the affine group x -> a x + b, a a nonzero square, is 1 mod 8; under the full automorphism
  group with Frobenius it is 3 mod 8).
  CONTROLS (refuting broader versions of B): among all 512 circulant tournaments on Z_19 the orbit counts take all
  residues 1, 3, 5, 7 mod 8, and the doubly regular tournament T_4 of order 15 (skew-Hadamard tower, |Aut| = 21) has
  orbit count 5 mod 8. So B is not a property of circulants or of doubly regular tournaments.
  Falsifiable at q = 31 (out of reach of the direct dynamic programme here, about 33 GB).
session: opus-2026-10-14-S33
source: 05-knowledge/results/number_tournament_bridges_20261014.md
scripts:
  - 04-computation/experiments/aliquot_parity_squares_20261014.py (+ .out, .json)
  - 04-computation/experiments/aliquot_dickson_frustration_20261014.py (+ .out)
  - 04-computation/experiments/number_tournament_bridges_20261014.py (+ .out)
related: [THM-4644, THM-871, THM-307, THM-4641, HYP-9258]
---

# HYP-9273

## (A) Aliquot square avoidance

**Statement.** No aliquot cycle, of any length, contains a perfect square or twice a perfect square. Equivalently, every aliquot cycle is parity-homogeneous. The length-2 case is the open problem of opposite-parity amicable pairs.

**Why it is the right formulation.**
- `σ(n)` is odd iff the odd part of `n` is a square, i.e. iff `n ∈ □ ∪ 2□`. So `s(n) ≡ n + [n ∈ □ ∪ 2□] (mod 2)`.
- Parity is invariant along aliquot sequences except at squares and twice-squares, which have density about `(1 + 1/√2)/√x`.
- **Heuristic.** A cycle member of size about `x` hits such a set with probability about `x^(−1/2)`. Two hits are needed, so mixed cycles should be very rare. This is not a proof.

**Status.** OPEN. Evidence: the exhaustive census to `10^8`, the double-flip census to `10^18`, and the cited enumerations.

**Related mechanism.** Mod 2, the aliquot map is the identity off `□ ∪ 2□`. Modulo an odd `q` dividing every `σ(n_i)` it is the sign flip `n ↦ −n` (Dickson's frustration, note section 3b). (A) says cycles never use the parity defect; Dickson says odd cycles cannot use the sign flip. Both are statements about cycles avoiding atypical `σ`-residues.

## (B) Paley Hamiltonian orbits

**Statement.** For every prime power `q ≡ 3 (mod 4)`, the Paley tournament `P_q` (arcs `i → j` iff `j − i` is a nonzero square in `F_q`) has `H(P_q) / (q(q−1)/2) ≡ 1 (mod 8)`. Here `q(q−1)/2` is the order of the affine group `{x ↦ ax + b : a ∈ (F_q^*)^2}`, which acts freely on Hamiltonian paths.

**Status.** OPEN; six cases (`q = 3, 7, 11, 19, 23, 27`). The broader statements for all circulant tournaments with `p ≡ 3 (mod 4)` and for all doubly regular tournaments are refuted (controls above). The first test is `q = 31`.

**Remarks.**
- **The parity part is free.** By Rédei, `H` is odd, and the affine group has odd order. So the orbit count is odd: the hypothesis concerns the next two binary digits.
- **The 2-adic part.** The mod-8 statement needs finer 2-adic information, in the sense of the repository's 2-adic digits of `H` (OCF: `H = 1 + 2α_1 + 4α_2 + ...`).
- **NUMEROLOGY.** The orbit count `1729 = H(P_11)/55` is the Hardy–Ramanujan taxicab number.
