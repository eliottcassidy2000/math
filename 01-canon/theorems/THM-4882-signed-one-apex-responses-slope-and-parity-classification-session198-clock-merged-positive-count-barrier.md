---
id: THM-4882
title: "Signed one-apex responses aU - bD. (1) Session 198's O/Z source calculus is the one-apex grammar G(1,0) in profile coordinates (D, E, U, d) = (H, 1 - b_m, 1 + bin(b), 1), its g(n) is the one-apex count of the word of value n - 1 (KNOWN: Session 197 ONE_CORE_BINARY_INVERSE, Session 198 DIGIT_BLOCK, HYP-9468; restated here in THM-4878's bit order with an exact cross-check), and its clock is a generalised left eigenvector of the period matrix for the eigenvalue 1 with a Jordan drift on the strong slice. (2) Lemmas: g odd, g(n) <= n, g(n) >= (n+3)/2 for odd n >= 3; the branch formula F(2m-1) = (2a-b)m - a - b g(m) + b(m mod 2) (PROVED). (3) Classification of S(a,b) = {an - b g(n): n >= 3 odd} for coprime a, b with 1/2 < a/b < 1: all odd integers iff a even, b odd (Session 198); all even integers iff a, b both odd (new: the clock with the cyclic correction K_cyc = K_lin - b); for b even never complete: positive elements in one class modulo 4 when a/b <= 3/4, elements below -a - ab/(4a-3b) in one class when a/b > 3/4 (new); outside the interval one sign is missing (Session 198). (4) Barrier: in an apex grammar H(pW^(h+1)q) >= H(pW^h q), so no fixed block lowers a count and the signed recipe has no positive-count analogue (PROVED)"
status: >
  PROVED (new): (2) the lemmas; (3) parts (ii) both odd and (iii) b even with the sharp threshold -3(b - a), and the
  parity statements; (4) the barrier (an elementary corollary of KNOWN monotonicity). KNOWN / CITED: (1) the
  dictionary O/Z = one-apex and g = the one-apex count (Session 197 path-polynomials/ONE_CORE_BINARY_INVERSE.md,
  Session 198 source-digit-blocks/DIGIT_BLOCK.md including the first-letter term -E eps_0 of the block transfer,
  Session 199 hidden-structure (3), HYP-9468/session179; audit major 1), restated here with the automaton
  cross-check to n = 4096; (3)(i) and the sharpness of the interval (Session 198 rational-response-coverage,
  2026-11-05, with its dependencies STRONG_SOURCE_INTERFACE and JOINT_MODULAR_GADGET, used as a black box for the
  modular seeds in (ii)); the one-apex recurrence (THM-4876 (6), HYP-9468, Session 197 palindrome-mersenne-clocks);
  the insertion lemma (THM-1370); the one-apex automaton (THM-4878 (1)-(2)); the Jordan reading of additive drifts
  (ribbon lane, Session 194).
  FINITE-EXACT: the probe to n < 2^24 and the exact recipe runs (note, section 5); pre-registered P1-P6: P1/P2
  small-witness clauses REFUTED for the larger denominators (audit major 2), P2's exact coverage PROVED, P3
  CONFIRMED and PROVED, P4 REFUTED, P5 and P6 CONFIRMED.
  CONCURRENT: kind-pasteur-2026-11-05-S12 (THM-4883, HYP-9476, reserved 14:57 after this session's pre-registration
  at 14:22) on the exact 2-adic structure of the even-denominator case.
  NOT CLAIMED: any positive-count coverage; Chan-Pak remains OPEN.
  AUDIT: PASS WITH CORRECTIONS (independent, adversarial, own code: 2 majors, 12 minors, 4 upgrades, all applied;
  MISTAKE-722; 04-computation/experiments/signed_slope_audit_20261105/AUDIT.md). The audit re-implemented the recipe
  (13 both-odd slopes, 23 targets each), materialised 77 clock outputs as tournaments (77/77), checked (iii) on 47
  even-denominator slopes and the dictionary on 254 words.
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

1. **Dictionary (KNOWN; restated).** The orbit is the leading-bit-1 sub-language of the one-apex grammar G(1,0)
   rooted at the apex: for the chain word b with binary value v, (D, E, U, d) = (H(T(b)), 1 - b_m, 1 + v, 1); a bit 1
   is O, a bit 0 is Z; g(n) = H(T(b)) with v = n - 1 (exact check to n = 4096; the identification is Session 197's
   and Session 198's, see the status). The period matrix M_W of Session 198's clock word has, on the strong slice,
   the eigenvalue 2^k on U and the eigenvalue 1 on D and d with the generalised-eigenvector drift
   (aU - bD) M_W = (aU - bD) - K d (K = K_cyc; after an O the drop is K_lin); the interval (1/2, 1) is the one-apex
   cone D <= U <= 2(D - E).

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
   -a - b g, and b g = b mod 4). If 4a > 3b, every element below -3(b - a) comes from the odd-m branch and is
   2 - a mod 4 when 4 | b, -a mod 4 when b = 2 mod 4, and -3(b - a) = F(3) is attained (sharp; audit upgrade 2);
   at (a, b) = (3, 4) both sides are confined (positive values 1 mod 4, values below -3 are 3 mod 4).
   Outside the interval: a/b <= 1/2 gives no positive values, a/b >= 1 no negative values [Session 198].

4. **Barrier (PROVED, an elementary corollary of KNOWN monotonicity).** For words p, W, q of an apex grammar
   G(k, 0), H(T(p W^(h+1) q)) >= H(T(p W^h q)) (deleting one period leaves an induced subtournament; insertion lemma,
   THM-1370; the same monotonicity is used in session195 full-binary-return, session197 minimum-positive-clock and
   Session 194's Theorem 1). Hence no repeated block lowers a count: the signed recipe's drop -K per period has no
   analogue for a Hamiltonian-path count, and positive coverage cannot be obtained by repetition of a fixed block
   from a modular seed; it needs words that vary with the target (S79).

## Verification

`signed_slope_probe_20261105.cpp` (dictionary check, small-witness coverage of 37 slopes to n < 2^24, residue classes
of the positive values), `signed_clock_recipe_20261105.py` (exact drops on 200 random strong points per slope, 40
random targets per slope hit exactly for 11 slopes of classes (i) and (ii), residue confinement of the positive and
negative values for 15 slopes of class (iii)); pre-registration `signed_slope_prereg_20261105.md`. The independent
audit reproduced all of it with its own code, verified (ii) on five further slopes and on 77 materialised
tournaments, and checked (iii) on 47 even-denominator slopes without violation.
