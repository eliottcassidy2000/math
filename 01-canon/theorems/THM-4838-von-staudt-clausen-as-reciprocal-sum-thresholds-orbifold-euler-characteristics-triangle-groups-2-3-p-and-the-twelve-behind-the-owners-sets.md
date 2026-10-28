---
id: THM-4838
title: "von Staudt-Clausen as reciprocal-sum thresholds: the sphere with cone points vSC(2k) has orbifold Euler characteristic -B_2k mod 1, exactly for 2k in {4,6,8,10}; chi(2,3,p) = -B_(p-1) exactly iff p in {5,7,11}, so 2g(X(p)) - 2 = |PSL_2(F_p)| B_(p-1) (g = 0, 3, 26); the owner's {4,2,5,3,7,13} = {d+1 : d | 12} (primes = vSC(12) = genus-zero primes = p^3 | |M|), {2,3,11} = vSC(10) with (p+1) | 12, {2,5,7} with (p-1) | 12, the two prime triples with sigma = 12^2; the vSC sinc product is a Borwein integral with excess -B_2k"
status: PROVED (elementary, on cited von Staudt-Clausen, Riemann-Hurwitz, genus of X_0(p)) + CITED (ATLAS, Adams, Mason/Dummit-Kisilevsky-McKay, Session43, Session83, codex-incidence) + FINITE-EXACT; one NUMEROLOGY remark flagged; independently audited (see the note's section 8)
session: opus-2026-10-28-S67
date: 2026-10-28
note: 05-knowledge/results/borwein_ladder_owner_numbers_vsc_twelve_20261028.md
script: 04-computation/experiments/borwein_ladder_vsc_twelve_20261028.py (+ .out)
related: [THM-4837, THM-4825, THM-4566, THM-4805]
---

# THM-4838 - the twelve behind the owner's sets

Let vSC(2k) = {p prime : (p-1) | 2k} and N_k = B_2k + sum_{p in vSC(2k)} 1/p, an integer by von Staudt-Clausen.
N_k = 1 for 2k <= 12.

**(1) Orbifold reading.** chi_k = 2 - |vSC(2k)| + sum_{vSC(2k)} 1/p (the sphere with those cone orders) satisfies
chi_k = -B_2k (mod 1) for all k. Equality is exact iff N_k = |vSC(2k)| - 2 iff 2k in {4, 6, 8, 10}.
- The corresponding vSC sets are (2,3,5), (2,3,7), (2,3,5), (2,3,11).
- It fails at 2k = 2, 12, 14, 16, 18.
- For 2k >= 20, |B_2k| > 2 sqrt(2k) + 2 >= the possible range.

**(2) Triangle groups.** For primes p >= 3, chi(2,3,p) = 1/p - 1/6.
- It equals -B_(p-1) iff p in {5, 7, 11}: exactness needs |B_(p-1)| < 1/6, so p - 1 <= 10, and p = 3 fails by sign.
- It is congruent to -B_(p-1) mod 1 iff vSC(p - 1) = {2,3,p}: a non-empty sum of reciprocals of distinct primes is
  never an integer.

**(3) Modular curves.** By Riemann-Hurwitz for X(p) -> X(1) branched (2,3,p), 2g(X(p)) - 2 = |PSL_2(F_p)| B_(p-1) for
p = 5, 7, 11 (g = 0, 3, 26).
- These are Galois's exceptional primes. The two lists are both "p <= 11" conditions, so the coincidence is recorded as
  a DICTIONARY only.
- The earlier repository observation "B_6 = 1/42 = Hurwitz {2,3,7}" (`05-knowledge/results/bernoulli_42_s116n.out`)
  is the p = 7 case of (2). The Hurwitz involved is the automorphism theorem, not the quaternion order.

**(4) The owner's sets.**
- {4,2,5,3,7,13} = {d + 1 : d | 12}. As an ordered cycle it is Session43's divisor-height cycle on
  N = 2^5 3^7 5^3 7^4 11 13^2 (|M|/|B| core), recovered by codex-incidence. Its primes {2,3,5,7,13} are all of:
  - vSC(12), the denominator 2730 of B_12;
  - {p : (p - 1) | 12} = the genus-zero primes of X_0(p) (classical);
  - {p : p^3 | |M|} (the extraspecial p^{1+24/(p-1)}, ATLAS);
  - the exponent->=2 primes of N (an observation).
- {2,3,11} = vSC(10): B_10 = 5/66, and 264 = |im J in pi_19^s|. Here p + 1 | 12 for every member.
- {2,5,7} has p - 1 | 12 for every member and is the unique subset of vSC(12) with prod(p+1) = 144.
- {2,3,11} and {2,5,7} are the only three-prime sets with prod(p + 1) = 12^2 (Session83's sigma(66) = sigma(70) = 144).
- The parallel list {p : (p + 1) | 24} = {2,3,5,7,11,23} gives the M_24 eta products.

**(5) Borwein meets von Staudt-Clausen.** For 2k <= 12, the integral
int_0^oo sinc(x) prod_{p in vSC(2k)} sinc(x/p) dx equals pi/2 - pi P(sum U_p > 1), with excess -B_2k:
- exactly pi/2 for 2k = 2, 6, 10;
- pi/2 - pi/43200 for 2k = 4, 8;
- pi/2 - pi (691^5 - 271^5)/(3840 2730^4) for 2k = 12.

NUMEROLOGY: the appearance of 420 = 2 * 2730/13 there is arithmetic, not evidence.

**Novelty.** (1), (2), (3) as stated and (5) are new to the repository, and no published statement was found in a
quick search; they are probably folklore-level. (4) assembles cited facts and credits the concurrent sigma-square-class
lanes (session84, codex-incidence, unifying-carriers).
