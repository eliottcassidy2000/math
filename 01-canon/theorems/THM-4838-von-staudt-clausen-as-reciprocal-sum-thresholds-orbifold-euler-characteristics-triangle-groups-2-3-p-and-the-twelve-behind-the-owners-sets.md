---
id: THM-4838
title: "von Staudt-Clausen as reciprocal-sum thresholds: the sphere with cone points vSC(2k) has orbifold Euler characteristic -B_2k mod 1, exactly for 2k in {4,6,8,10}; for primes p >= 5, chi(2,3,p) = -B_(p-1) exactly iff p in {5,7,11} (and mod 1 iff denom B_(p-1) = 6p, OEIS A092307), so 2g(X(p)) - 2 = |PSL_2(F_p)| B_(p-1) (g = 0, 3, 26); the sinc product over vSC(2k) is a Borwein integral with excess -B_2k (2k <= 12); plus an INTERPRETIVE dictionary for the owner's {4,2,5,3,7,13}, {2,3,11}, {2,5,7} (the 'twelve' reading is numerology)"
status: PROVED (elementary, on cited von Staudt-Clausen, Riemann-Hurwitz, genus of X_0(p)) + CITED (OEIS A092307; ATLAS; Conway-Norton 1979; Adams; Mason/Dummit-Kisilevsky-McKay; BBB 2008 Example 2; Session43, Session83, codex-incidence) + FINITE-EXACT; the unifying 'twelve' reading of the owner's sets is NUMEROLOGY/INTERPRETATION; independently audited SOUND WITH CORRECTIONS (MISTAKE-692)
session: opus-2026-10-28-S67
date: 2026-10-28
note: 05-knowledge/results/borwein_ladder_owner_numbers_vsc_twelve_20261028.md
script: 04-computation/experiments/borwein_ladder_vsc_twelve_20261028.py (+ .out)
audit: 04-computation/experiments/borwein_ladder_vsc_twelve_20261028_audit_S67/AUDIT.md
related: [THM-4837, THM-4825, THM-4566, THM-4805, THM-4644]
---

# THM-4838 - von Staudt-Clausen as reciprocal-sum thresholds

Let vSC(2k) = {p prime : (p-1) | 2k} and N_k = B_2k + sum_{p in vSC(2k)} 1/p, an integer by von Staudt-Clausen.
N_k = 1 for 2k <= 12.

**(1) Orbifold reading.** chi_k = 2 - |vSC(2k)| + sum_{vSC(2k)} 1/p, the orbifold Euler characteristic of the sphere
with those cone orders, satisfies chi_k = -B_2k (mod 1) for all k. Equality is exact iff N_k = |vSC(2k)| - 2 iff
2k in {4, 6, 8, 10}, where the vSC sets are (2,3,5), (2,3,7), (2,3,5) and (2,3,11).

*Proof of "only".*
- Exactness forces B_2k = |vSC| - 2 - sum 1/p. So |B_2k| <= d(2k) + 2 <= 2 sqrt(2k) + 2.
- But |B_20| = 529.1, and for 2k >= 20 the ratio |B_(2k+2)/B_2k| = (2k+1)(2k+2) zeta(2k+2)/(4 pi^2 zeta(2k)) is
  >= 11.7, while 2 sqrt(2k) + 2 grows by a factor < 1.05 per step.
- The cases 2k <= 18 are checked directly: exactness fails at 2, 12, 14, 16, 18.

**(2) Triangle groups.** For primes p >= 5, chi(2,3,p) = 1/p - 1/6.
- It **equals -B_(p-1) iff p in {5, 7, 11}**. Exactness needs |B_(p-1)| < 1/6. The same ratio bound (>= 4.6 for 2k >= 12)
  gives |B_2k| >= |B_12| > 1/6 for 2k >= 12, so p - 1 <= 10.
- It is **congruent to -B_(p-1) mod 1 iff vSC(p - 1) = {2,3,p}**, i.e. iff denom B_(p-1) = 6p. These primes are OEIS
  A092307: 5, 7, 11, 23, 47, 59, 83, ...
  - Proof: a non-empty sum of reciprocals of distinct primes is never an integer.
- At p = 3 both statements fail: chi(2,3,3) + B_2 = 1/3, because 3 is counted twice.

**(3) Modular curves.** By Riemann-Hurwitz for X(p) -> X(1), branched (2,3,p),
2g(X(p)) - 2 = |PSL_2(F_p)| B_(p-1) for p = 5, 7, 11 (g = 0, 3, 26).
- These are also Galois's exceptional primes. Both lists are "p <= 11" conditions, so this is a DICTIONARY only.
- The relevant Hurwitz theorem is the 84(g - 1) bound. The (2,3,7) group is, modulo its centre, the norm-one group of
  the Hurwitz quaternion order over Z[2cos(2 pi/7)]. It is not related to the rational Hurwitz quaternions.

**(4) The owner's sets: facts (PROVED/CITED) and one reading (NUMEROLOGY).**
- Facts about the first set:
  - {4,2,5,3,7,13} = {d + 1 : d | 12}. As an ordered cycle it enters through Session43's prime-to-exponent arrows of
    the |M|/|B| core 2^5 3^7 5^3 7^4 11 13^2 (recovered by codex-incidence). There, 4 = v_7 and 2 = v_13, so the
    equality of sets has no mechanism.
  - Its primes {2,3,5,7,13} are vSC(12) (denominator 2730 of B_12), the genus-zero primes of X_0(p) (Hauptmodul
    (eta(tau)/eta(p tau))^{24/(p-1)}; Conway-Norton 1979 for the pB moonshine functions), and the primes with
    p^3 | |M|.
  - For p = 3, 5, 7, 13 the Monster has normalisers N(pB) = p^{1+24/(p-1)}.(...), and C(2B) = 2^{1+24}.Co_1 (ATLAS).
- Facts about the two triples:
  - {2,3,11} = vSC(10): B_10 = 5/66, and 264 = |im J in pi_19^s|.
  - {2,3,11} and {2,5,7} are the only three-prime sets with prod(p+1) = 144 = sigma(66) = sigma(70) (Session83's
    fibre). This is their only exact characterisation.
  - {p : (p+1) | 12} = {2,3,5,11} and {p : (p-1) | 12} = {2,3,5,7,13} are strictly larger than the owner's triples.
- **NUMEROLOGY:** reading the three sets as "the twelve" is an interpretation, not a theorem.

**(5) Borwein meets von Staudt-Clausen.** For 2k <= 12, the integral
int_0^oo sinc(x) prod_{p in vSC(2k)} sinc(x/p) dx = pi/2 - pi P(sum U_p > 1), with excess -B_2k:
- exactly pi/2 for 2k = 2, 6, 10;
- pi/2 - pi/43200 for 2k = 4, 8;
- pi/2 - pi (691^5 - 271^5)/(3840 2730^4) for 2k = 12.

Context: BBB 2008 Example 2 treats sinc products over prime reciprocals; the vSC version is new to the repository.
NUMEROLOGY: the appearance of 420 = 2 * 2730/13 is arithmetic.

**Novelty.**
- The ingredients were already in the repository:
  - `bernoulli_tournament_primes_s18h.out` (denom B_6 = 42 = Hurwitz number; denominators 30, 42, 66);
  - S33 `number_tournament_bridges_20261014.out` (a chi(2,3,k) table with 1/30 and -1/42);
  - `understand_2357_s90bl.out` (B_10 primes {2,3,11});
  - `bernoulli_42_s116n.out` (the {d+1 : d | n} device and "Hurwitz primes").
- The exact statements (1), the exact half of (2), (3) and (5) are new to the repository. No published statement was
  found by either search; they are probably folklore-level.
- The congruence half of (2) is A092307.
