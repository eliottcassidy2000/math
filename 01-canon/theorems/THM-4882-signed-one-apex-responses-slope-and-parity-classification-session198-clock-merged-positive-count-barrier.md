---
id: THM-4882
title: "Signed one-apex responses aU - bD. (1) Session 198's O/Z source calculus is the one-apex grammar G(1,0) in profile coordinates (D, E, U, d) = (H, 1 - b_m, 1 + bin(b), 1), its g(n) is the one-apex count of the word of value n - 1, and its clock is a left eigenvector of the period matrix for the eigenvalue 1 with a Jordan drift (PROVED, elementary). (2) Lemmas: g odd, g(n) <= n, g(n) >= (n+3)/2 for odd n >= 3; the branch formula F(2m-1) = (2a-b)m - a - b g(m) + b(m mod 2) (PROVED). (3) Classification of S(a,b) = {an - b g(n): n >= 3 odd} for coprime a, b with 1/2 < a/b < 1: all odd integers iff a even, b odd (Session 198); all even integers iff a, b both odd (new: the clock with the cyclic correction K_cyc = K_lin - b); for b even never complete: positive elements in one class modulo 4 when a/b <= 3/4, elements below -a - ab/(4a-3b) in one class when a/b > 3/4 (new); outside the interval one sign is missing (Session 198). (4) Barrier: in an apex grammar H(pW^(h+1)q) >= H(pW^h q), so no fixed block lowers a count and the signed recipe has no positive-count analogue (PROVED)"
status: >
  PROVED: (1) the dictionary and the eigen-functional reading; (2) the lemmas; (3) parts (ii) both odd and (iii) b even,
  and the parity statements; (4) the barrier. KNOWN / CITED: (3)(i) and the sharpness of the interval (Session 198
  rational-response-coverage, 2026-11-05, with its dependencies STRONG_SOURCE_INTERFACE and JOINT_MODULAR_GADGET, used
  as a black box for the modular seeds in (ii)); the one-apex recurrence (THM-4876 (6), HYP-9468, Session 197
  palindrome-mersenne-clocks); the insertion lemma (THM-1370); the one-apex automaton (THM-4878 (1)-(2)).
  FINITE-EXACT: the probe to n < 2^24 and the exact recipe runs (note, section 5); pre-registered P1-P6
  (P4 refuted, the rest confirmed or proved).
  NOT CLAIMED: any positive-count coverage; Chan-Pak remains OPEN.
  AUDIT: PENDING.
session: opus-2026-11-05-S80
source: 05-knowledge/results/signed_responses_slope_classification_clock_barrier_20261105.md
related:
  - THM-4878
  - THM-4876
  - THM-1370
  - HYP-9468
  - HYP-9475
  - session198/rational-response-coverage/PROOF.md
  - session198/rational-response-coverage/dependencies/STRONG_SOURCE_INTERFACE.md
  - session197/modular-gadget-universality/PROOF.md
  - session194/two-apex/PROOF.md
---

# THM-4882 - signed one-apex responses: classification by parity and slope; the clock merged into G(1,0); the barrier

## Setting

g(1) = 1, g(2n) = g(n), g(2n - 1) = g(n) + n - (n mod 2). For coprime positive integers a, b put F(n) = a n - b g(n)
and S(a, b) = { F(n) : n >= 3 odd }. Session 198's sources are rooted tournaments with one root-start Hamiltonian
path and profile (D, E, U, d) = (H, paths ending at the root, ordered two-path covers whose second path starts at the
root, H(root deleted)); its operations O(D,E,U,d) = (D,0,2U,d) and Z(D,E,U,d) = (D-E+U,d,2U-d,d) generate from the seed
(1,0,2,1) the orbit (g(n), n mod 2, n, 1), whose strong points are the odd n. Full proofs: the source note, sections
2-4 and 6.

## Statements

1. **Dictionary (PROVED).** The orbit is the one-apex grammar G(1,0) rooted at the apex: for the chain word b with
   binary value v, (D, E, U, d) = (H(T(b)), 1 - b_m, 1 + v, 1); a bit 1 is O, a bit 0 is Z; g(n) = H(T(b)) with
   v = n - 1 (exact check to n = 4096). The period matrix M_W of Session 198's clock word has, on the strong slice,
   the eigenvalue 2^k on U and the eigenvalue 1 on D and d with the drift (aU - bD) M_W = (aU - bD) - K d; the
   interval (1/2, 1) is the one-apex cone D <= U <= 2(D - E).

2. **Lemmas (PROVED).** g(n) is odd; g(n) <= n; g(n) >= (n + 3)/2 for odd n >= 3. For odd n = 2m - 1 (m >= 2),
   F(n) = (2a - b) m - a - b g(m) + b (m mod 2), and for even m = 2m', F(n) = 2(2a - b) m' - a - b g(m').

3. **Classification.** Let 1/2 < a/b < 1, gcd(a, b) = 1.
   (i) [Session 198] a even, b odd: S(a, b) = all odd integers.
   (ii) [new, PROVED] a and b odd: S(a, b) = all even integers. The period word of the bits of the odd integer
   A = a(2^k - 1)/b (2^k = 1 mod b) starts with Z; a strong input loses K_cyc = K_lin - b = aR - bC_cyc per period
   (C_cyc = Session 198's C plus the cyclic term eps_(k-1) eps_0 = 1), K_cyc > 0 even, and the input right after the
   O loses K_lin; the recipe seed -> O -> W^h with residues 2aU_0 - bD_0 = r + b (mod K_cyc) hits every even r.
   (iii) [new, PROVED] a odd, b even: S(a, b) is never the set of all odd integers. If 4a <= 3b, every positive
   element is -(a + b) mod 4 (the odd-m branch is negative by the lemma, the even-m branch is 4a - 2b = 0 mod 4 plus
   -a - b g, and b g = b mod 4). If 4a > 3b, every element below -a - ab/(4a - 3b) comes from the odd-m branch and is
   2 - a mod 4 when 4 | b, -a mod 4 when b = 2 mod 4.
   Outside the interval: a/b <= 1/2 gives no positive values, a/b >= 1 no negative values [Session 198].

4. **Barrier (PROVED).** For words p, W, q of an apex grammar G(k, 0), H(T(p W^(h+1) q)) >= H(T(p W^h q)) (deleting
   one period leaves an induced subtournament; insertion lemma). Hence no repeated block lowers a count: the signed
   recipe's drop -K per period has no analogue for a Hamiltonian-path count, and positive coverage cannot be
   obtained by repetition of a fixed block from a modular seed; it needs words that vary with the target (S79).

## Verification

`signed_slope_probe_20261105.cpp` (dictionary check, small-witness coverage of 36 slopes to n < 2^24, residue classes
of the positive and negative values), `signed_clock_recipe_20261105.py` (exact drops on 200 random strong points per
slope, 40 random targets per slope hit exactly for 11 slopes of classes (i) and (ii), residue confinement for 15
slopes of class (iii)); pre-registration `signed_slope_prereg_20261105.md`.
