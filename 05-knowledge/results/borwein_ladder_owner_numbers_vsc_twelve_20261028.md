# The owner's Borwein ladder (15, 113, 80499; 7, 56, 418; 231, 233, 420) and von Staudt-Clausen behind {4,2,5,3,7,13}, {2,3,11}, {2,5,7}

opus-2026-10-28-S67, 2026-10-28.

**Status.**
- THM-4837 is mostly CITED:
  - Borwein-Borwein 2001; Baillie-Borwein-Borwein 2008;
  - Almkvist-Gustavsson 2014; Basel-Baillie 2016;
  - Crandall via Baillie 2011; Schmid 2014;
  - OEIS A092315.
- PROVED: the window criterion, and THM-4838's elementary statements.
- FINITE-EXACT: interval arithmetic and exact rationals, re-certified by the audit.
- HYP-9436 (the crossing formula for every n) is OPEN in general and PROVED for n <= 3000.
- INTERPRETATION: the decode of the owner's 231 and 233. SPECULATIVE: the decode of 420. NUMEROLOGY: the "twelve"
  reading of the three sets.
- Independently audited, **SOUND WITH CORRECTIONS**, all applied (section 8; MISTAKE-692).

## 0. Directive and the concurrent lanes

Owner, 2026-10-28: "merge in prior work you can find and new connections you can reason regarding 15, 113, 80499 and
7, 56, 418 and also {231, 233, 420} and also consider prior work on {4,2,5,3,7,13} and how it relates with {2,3,11} and
{2,5,7}." The same message continues with the attractor-synthesis and abc requests; those are opus-2026-10-28-S68.

Three lanes took the same directive the same morning. Their results are not repeated here.
- **session84** (`session84_guarded_unification_20261028/`):
  - {2,3,11} is the odd-square-bracket doubling set, with multiplier interval (25/13, 25/11].
  - 66 is reached by a Berggren route.
  - The F=U+S frontier.
- **codex-incidence** (`incidence_attractors_20261028/`):
  - Recovers the owner's ordered labels as **Session43's divisor-height five-cycle 4 -> 2 -> 5 -> 3 -> 7 -> 4 with
    13 -> 2**, on the small-prime core 2^5 3^7 5^3 7^4 11 13^2 of |Monster|/|Baby Monster|
    (`session43_observers_phase_guards_20261023/cycle-geometry/PROOF.md`).
  - The sigma square-class chart on {4,2,5,3,7,13,11}: {2,5,7} and {4,7,13} are circuits; {2,3,11} is a parallel pair
    plus a loop.
  - The 364 -> 420 absorber path.
- **unifying-carriers:** arithmetic matroids.

This note supplies the analytic (Borwein) reading of the numbers and the Bernoulli/orbifold reading of the sets.

## 1. Results first

| Owner's numbers | What they are | Typing |
|---|---|---|
| 7, 56, 418 | m(2), m(3), m(4): the least m with 1 + 1/3 + ... + 1/(2m+1) > 2, 3, 4 (OEIS A092315; crossings tabulated by Basel-Baillie 2016) | CITED + FINITE-EXACT |
| 15, 113 | the Borwein breaks: sinc(x/15) in the classical integral, sinc(x/113) in the 2cos(x) variant (Schmid 2014) | CITED |
| 80499 | the Baillie-Borwein-Borwein break: the sinc sum equals the sinc integral for N <= 40248 and first fails at N = 40249, denominator 80499, where 1 + 1/3 + ... + 1/80499 first exceeds 2 pi (BBB 2008, Example 1(b)) | CITED + FINITE-EXACT (re-certified) |
| 231, 233 | match, to three digits, the published first-failure deficits pi/2 - **2.31**00572725e-11 (at 1/15) and pi/2 - **2.33**24296128e-138 (at 1/113) | INTERPRETATION |
| 420 | matches only to two digits the one-sided error at 80499, **4.21**01759371e-226577; "4.20" needs a naive double sum truncated, not rounded | SPECULATIVE |
| {4,2,5,3,7,13} | Session43's exponent cycle of the \|M\|/\|B\| core; as a set it is also {d+1 : d \| 12}, whose primes are vSC(12) = the genus-zero primes of X_0(p) = the p with p^3 \| \|M\| | facts PROVED/CITED; the link between the readings is NUMEROLOGY |
| {2,3,11} | vSC(10) (B_10 = 5/66); the triangle group (2,3,11) has orbifold Euler characteristic exactly -B_10 | PROVED + CITED |
| {2,5,7} | with {2,3,11}, the only three-prime sets with prod(p+1) = 144 = sigma(66) = sigma(70) (Session83) | CITED |

The two halves of the directive share one elementary principle: a **reciprocal sum compared with an integer**. This
is a framing, not a theorem.
- The Borwein integrals break when a sum of reciprocals passes an integer (or 2 pi, via Poisson summation).
- von Staudt-Clausen says the reciprocal sum of the primes p with (p - 1) | 2k misses an integer by exactly B_2k.
- A sphere with cone points of those orders has that sum in its orbifold Euler characteristic.

## 2. THM-4837: the Borwein ladder, with its sources

Let P_m(x) = prod_{j=0}^{m} sinc(x/(2j+1)) and O_m = sum_{j=0}^{m} 1/(2j+1). Let U_j be uniform on
[-1/(2j+1), 1/(2j+1)], independent, and S'_m = U_1 + ... + U_m.

- **(1) Vanishing.** For m >= 1, F_m(xi) := int_R P_m(x) cos(xi x) dx = pi P(xi - 1 <= S'_m <= xi + 1). So for xi >= 1,
  F_m(xi) = 0 iff O_m <= xi.
  - This is the box-spline/frequency-shift principle (Schmid 2014; Majumdar-Trizac 2019; THM-4392 uses the same duality
    in the repository).
  - The edge m = 0, xi = 1 gives pi/2 and is excluded.
- **(2) Ladder.** I_m(t) = int_0^oo (sin t x/x) prod_{j=1}^{m} sinc(x/(2j+1)) dx = pi/2 - pi P(S'_m > t)
  = pi/2 - sum_{k>=0} F_m(t + 1 + 2k). So I_m(t) = pi/2 iff O_m <= t + 1 (Borwein-Borwein 2001).
  - t = 1 is the classical integral and t = 2 the 2cos x variant (Schmid).
  - t = 4 and t = 6 are Almkvist-Gustavsson's cosine families.
  - The telescoping identity is checked exactly for m <= 8.
- **(3) First failure.** With m = m(xi) and eps = O_m - xi, F_m(xi) = pi eps^m (2m+1)!!/(2m)!!. This is
  **Borwein-Borwein's Corollary 1** (2001), also BBB 2008 Theorem 2 "First Bite". Its hypotheses state that m is the
  first failure; for the decreasing sequence 1/(2j+1) they hold at every first failure, as Basel-Baillie (2016) use at
  every integer threshold.
- **(4) Values.**

| threshold | m | 2m+1 | first-failure value | source |
|---|---|---|---|---|
| 2 | **7** | **15** | 2.3100572725e-11 | BB 2001; Wikipedia |
| 3 | **56** | **113** | 2.3324296128e-138 | Almkvist-Gustavsson 2014; Wikipedia |
| 4 | **418** | 837 | 4.5195214277e-1320 | trivial instance of Cor. 1 (Basel-Baillie tabulate this crossing) |
| 2 pi | 40249 | **80499** | one-sided 4.2101759371e-226577 | BBB 2008 Ex. 1(b); Crandall via Baillie 2011 |

- The thresholds 5, 6, 7 give m = 3091, 22845, 168803. The value at 15 reproduces the published rational
  467807924713440738696537864469/935615849440640907310521750000 pi.
- **The 2 pi row is the k = 1 Poisson term.** sum_{n in Z} P_m(n) - int_R P_m = 2 sum_{k>=1} F_m(2 pi k), so the
  cosine transform at one full turn, F_m(2 pi), is BBB's one-sided sum-minus-integral error.
- Crandall's value as reported by Baillie 2011, "about 8.42e-226577", is the two-sided error. **New here:** the one-sided
  value to 11 digits, and that factor-2 remark.
- **(5) Crossing index.** With M = m + 1, O_m = (1/2) ln(4M) + gamma/2 + r(M), and 0 < r(M) < 1/(48 M^2) for all M >= 1
  (GKP (9.89) / DLMF 5.11(ii) envelope). Hence **m(xi) = floor(X) whenever frac(X) >= X/(24 floor(X)^2)**, where
  X = exp(2 xi - gamma)/4.
  - This is a fully proved version of G. Kirchner's 2020 OEIS test. Kirchner used the remainder -127t/(1920 m^4) with
    0 < t < 1, supported numerically, and certified n <= 1000.
  - Here the certified range extends to 2 <= n <= 3000 and to xi = k pi, k <= 500. The audit re-certified this in
    interval arithmetic.
- **(6) An OEIS observation.** The A092315 formula line floor(exp(2n - gamma)/4 + 1/8) ("conjectured", Hasler 2017)
  contradicts the entry's own DATA term a(17) = 81897532160124.
  - a(17) is certified by the envelope and independently by digamma: O_{a-1} - 17 = -5.71e-15 and
    O_a - 17 = +3.97e-16.
  - The +1/8 form differs from the certified value exactly when frac(X_n) >= 7/8: 375 of the n <= 3000
    (17, 31, 45, 57, 59, 63, 86, 118, ...).
  - An OEIS edit is recommended.

## 3. HYP-9436: the crossing formula for every n (Comtet-type)

The conjecture: m(n) = floor(exp(2n - gamma)/4) for every n >= 2. It is the odd, half-shifted analogue of Comtet's
problem for H_k, and it is open in general: a failure needs exp(2n - gamma)/4 within about 1/(24X) above an integer.

Pre-registered (a14c71ec13) before the confirmatory run, P1-P3 were all CONFIRMED:
- P1: no window exception for n <= 3000;
- P2: Kolmogorov-Smirnov D = 0.0192, against the 1% critical value 0.0298;
- P3: the k pi thresholds for k <= 500, with m(pi) = 75 and m(3 pi) = 21553437.

The audit notes that P1-P3 were near-certain predictions; the heuristic hit expectation for n >= 61 is 3.6e-54. Their
confirmation adds little beyond the certificate itself.

## 4. THM-4838: von Staudt-Clausen as reciprocal-sum thresholds

Write vSC(2k) = {p prime : (p - 1) | 2k}. N_k := B_2k + sum_{p in vSC(2k)} 1/p is an integer, equal to 1 for 2k <= 12.

- **(1) Orbifold reading.** The sphere with cone points vSC(2k) has chi_k = 2 - |vSC(2k)| + sum 1/p = -B_2k (mod 1).
  Equality is exact **iff 2k in {4, 6, 8, 10}**, where the vSC sets are the triangle groups (2,3,5), (2,3,7), (2,3,5)
  and (2,3,11).
  - Proof of "only": exactness forces |B_2k| <= d(2k) + 2 <= 2 sqrt(2k) + 2. But |B_20| = 529.1 and
    |B_(2k+2)/B_2k| >= 11.7 for 2k >= 20.
  - The cases 2k <= 18 are checked directly.
- **(2) Triangle groups, p >= 5.** chi(2,3,p) = 1/p - 1/6.
  - **Equality with -B_(p-1) holds iff p in {5, 7, 11}**: |B_2k| >= |B_12| > 1/6 for 2k >= 12, by the ratio bound >= 4.6.
  - The congruence mod 1 holds iff vSC(p - 1) = {2,3,p}, i.e. iff denom B_(p-1) = 6p. These are the primes of OEIS
    **A092307**: 5, 7, 11, 23, 47, 59, 83, 107, ...
  - At p = 3 both fail, because 3 is counted twice.
- **(3) Modular curves.** 2g(X(p)) - 2 = |PSL_2(F_p)| B_(p-1) for p = 5, 7, 11: 60 (-1/30) = -2, 168 (1/42) = 4 and
  660 (5/66) = 50, giving g = 0, 3, 26.
  - The list {5, 7, 11} coincides with Galois's exceptional primes. Both are "p <= 11" conditions, so this is a
    DICTIONARY only.
  - The relevant Hurwitz theorem is the 84(g - 1) bound. The (2,3,7) group is, modulo its centre, the norm-one group of
    the Hurwitz quaternion order over Z[2cos(2 pi/7)], not of the rational Hurwitz quaternions.
- **(4) The owner's three sets: facts, then one reading.**
  - **{4,2,5,3,7,13}** enters through Session43 as the prime-to-exponent arrows of the |M|/|B| core: 4 = v_7 and
    2 = v_13. Its equality with {d + 1 : d | 12} is a set coincidence with no mechanism. The primes of the latter,
    {2,3,5,7,13}, are all of the following (CITED):
    - vSC(12), the denominator 2730 of B_12;
    - the genus-zero primes of X_0(p), whose Hauptmoduln are (eta(tau)/eta(p tau))^{24/(p-1)} (Conway-Norton 1979 for
      the pB moonshine functions);
    - the primes with p^3 | |M| (normalisers N(pB) = p^{1+24/(p-1)}.(...) for p = 3, 5, 7, 13, and
      C(2B) = 2^{1+24}.Co_1; ATLAS).
  - **{2,3,11} = vSC(10).** B_10 = 5/66, E_10 = 1 - 264 sum sigma_9(n) q^n, and |im J| in pi_19^s is 264.
  - **{2,3,11} and {2,5,7}** are the only three-prime sets with prod(p + 1) = 144 (Session83's sigma^{-1}(144) fibre).
    This is their only exact characterisation. {p : (p + 1) | 12} = {2,3,5,11} and {p : (p - 1) | 12} = {2,3,5,7,13}
    are strictly larger.
  - **NUMEROLOGY:** reading the three sets as "the twelve" is an interpretation, not a result.
  - The M_24 list {p : (p + 1) | 24} = {2,3,5,7,11,23} (Mason; Dummit-Kisilevsky-McKay) is context only.
- **(5) Borwein meets von Staudt-Clausen.** For 2k <= 12, sum_{vSC(2k)} 1/p - 1 = -B_2k. So
  int_0^oo sinc(x) prod_{p in vSC(2k)} sinc(x/p) dx is a Borwein integral with excess -B_2k:
  - exactly pi/2 for 2k = 2, 6, 10;
  - short by pi/43200 for 2k = 4, 8;
  - short by pi (691^5 - 271^5)/(3840 * 2730^4) for 2k = 12.

  BBB 2008 Example 2 already treats sinc products over prime reciprocals; the vSC version is new to the repository.
  **NUMEROLOGY:** the appearance of 420 = 2 * 2730/13 is arithmetic.

**Repository prior work behind THM-4838:**
- `bernoulli_tournament_primes_s18h.out`: denom B_6 = 42 = Hurwitz number; denominators 30, 42, 66.
- S33 `number_tournament_bridges_20261014.out`: a chi(2,3,k) table with 1/30 and -1/42.
- `understand_2357_s90bl.out`: vSC(10) = {2,3,11}.
- `bernoulli_42_s116n.out`: the {d + 1 : d | n} device and "Hurwitz primes".

The exact statements (1), the exact half of (2), (3) and (5) are new to the repository. Neither the session's search
nor the audit's found a published statement; they are probably folklore-level.

## 5. {231, 233, 420} in prior repository work

Two readings sit alongside the Borwein matches.
- **Collatz.** 231 -> 233 is a recorded 17-step Syracuse word on the 27/17 clock
  (`collatz_procgen_20260922_order_laws.md`, `crossroads233_20260926_graph.md` section 5,
  `collatz_join_shields_and_lifts_20261004.md`). m_17(231 + 2^28 t) = 233 + 2 * 3^17 t is re-checked here.
- **Aliquot.**
  - 231 = s(345) with sigma(345) = 24^2, and 420 = s(364) with sigma(364) = 28^2. Both are square-sum aliquot targets in
    Session83's census.
  - 364 is session84's least mixed balanced square-sigma source, and codex-incidence's absorber path
    19236--364--420--19180 passes through it.
  - No n with square sigma has s(n) = 233.

**Typing.** The Borwein reading *matches* 231 and 233 to three digits, using the owner's own 15 and 113. It matches 420
only to two digits, and that match is speculative. The Collatz and aliquot readings are genuine prior-work occurrences,
but they do not explain why the three numbers go together.

## 6. What is new, cited and numerology

- **New to the repository:**
  - the one-sided 11-digit value at 80499 and the factor-2 remark about Crandall/Baillie;
  - the fully proved window criterion and the certified range n <= 3000, k pi <= 500 pi;
  - the explicit OEIS observation (a(17) against the +1/8 line);
  - the interpretive decode of the owner's numbers;
  - THM-4838 (1), the exact half of (2), (3) and (5).
- **Cited and published:**
  - the Borwein ladder and its single-corner first failures (BB 2001; BBB 2008; Schmid 2014);
  - the 1/113 value and the higher cosine families (Almkvist-Gustavsson);
  - the crossing table (Basel-Baillie);
  - the 80499 break (BBB 2008; Crandall; Baillie);
  - Kirchner's test; A092307.
- **Numerology:** the "twelve" reading of the owner's sets, and 420 = 2 * 2730/13 in (5).
- **Speculative:** 420 as the third mantissa.

## 7. Reproduction

```text
python3 04-computation/experiments/borwein_ladder_vsc_twelve_20261028.py 3000 500
```

- About 21 s; every claim is an assertion; the output is in `.out`.
- The audit's independent scripts are in `04-computation/experiments/borwein_ladder_vsc_twelve_20261028_audit_S67/`.

## 8. Audit record

Independent audit `04-computation/experiments/borwein_ladder_vsc_twelve_20261028_audit_S67/AUDIT.md`:
**SOUND WITH CORRECTIONS**, all applied here and in THM-4837, THM-4838 and HYP-9436.

- **What reproduced.** All 39 numerical claims, under independent code: exact-rational crossings with rational bounds
  on pi; the window test re-certified in interval arithmetic; a(17) by the DLMF envelope and digamma.
- **MAJOR.**
  - (1) The 80499 = 2 pi break and its size are published (BBB 2008 Ex. 1(b); Crandall; Baillie 2011). The
    "cosine transform at one full turn" is the k = 1 Poisson term.
  - (2) The single-corner first failure (BB 2001 Cor. 1; BBB 2008 Thm 2) and the 1/113 value (Almkvist-Gustavsson 2014)
    are published, and the ladder is tabulated by Basel-Baillie 2016.
  - Both were typed new. BBB 2008 was on the session's own citation list (MISTAKE-692).
- **MINOR.**
  - (3) Kirchner's 2020 test was under-described.
  - (4) The +1/8 "refutation" is visible in the entry's own DATA.
  - (5) THM-4838 (2) scope p >= 3 was false at p = 3.
  - (6) The |B_2k| growth step was missing.
  - (7) The congruence list is A092307.
  - (8) Repository prior work was under-credited (s18h, S33, s90bl).
  - (9) The 231/233/420 decode was typed too strongly.
  - (10) "The twelve" is numerology.
- **NOTES:** the m = 0 edge of (1); premature status lines; the heuristic 3.6e-54; interval re-certification; Monster
  normaliser wording and Conway-Norton; the Hurwitz-order wording; Majumdar-Trizac and Schmid credits; BBB Example 2.
- The pre-registration order of HYP-9436 was confirmed (a14c71ec13 before 7a5cd669d0).
