# The owner's Borwein ladder (15, 113, 80499; 7, 56, 418; 231, 233, 420) and the twelve behind {4,2,5,3,7,13}, {2,3,11}, {2,5,7}

opus-2026-10-28-S67, 2026-10-28.

**Status:**
- **PROVED:** THM-4837 and THM-4838. Both are elementary. They rest on the cited Borwein-Borwein formula, von Staudt-Clausen, Riemann-Hurwitz and the genus formula for X_0(p).
- **CITED:** Borwein-Borwein 2001; Schmid 2014; OEIS A092315 (Sloane, Hasler, Kirchner, Song); Ogg; the ATLAS; Adams.
- **FINITE-EXACT:** interval arithmetic and exact rationals.
- **HYP-9436** (the crossing formula for every n): OPEN in general, PROVED for n <= 3000; pre-registered P1-P3 CONFIRMED.
- **Two NUMEROLOGY items**, flagged where they occur.
- **Independently audited:** see section 8.

## 0. Directive and the concurrent lanes

Owner, 2026-10-28: "merge in prior work you can find and new connections you can reason regarding 15, 113, 80499 and
7, 56, 418 and also {231, 233, 420} and also consider prior work on {4,2,5,3,7,13} and how it relates with {2,3,11} and
{2,5,7}." (The same message continues with the attractor-synthesis and abc requests; those are opus-2026-10-28-S68.)

Three lanes took the same directive the same morning. Their results are not repeated here.
- **session84** (`session84_guarded_unification_20261028/`):
  - {2,3,11} is the odd-square-bracket doubling set, with exact multiplier interval (25/13, 25/11].
  - 66 reached by a Berggren route.
  - The F=U+S balance frontier.
- **codex-incidence** (`incidence_attractors_20261028/`):
  - Recovers the owner's ordered six labels from **Session43's divisor-height five-cycle 4 -> 2 -> 5 -> 3 -> 7 -> 4 with
    13 -> 2**. This is the small-prime core N = 2^5 3^7 5^3 7^4 11 13^2 of |Monster|/|Baby Monster|
    (`session43_observers_phase_guards_20261023/cycle-geometry/PROOF.md`).
  - The sigma square-class chart on the atoms {4,2,5,3,7,13,11}: {2,5,7} and {4,7,13} are circuits, and {2,3,11} is a
    parallel pair plus a loop.
  - A Fano/Hamming extension.
  - The absorber path 19236--364--420--19180.
- **unifying-carriers:** arithmetic matroids and divisor towers.

This note takes the two angles they did not: the **analytic** reading of the first six numbers (Borwein integrals) and
the **Bernoulli / genus-zero** reading of the three sets.

## 1. Results first

| Owner's numbers | What they are | Typing |
|---|---|---|
| 7, 56, 418 | m(2), m(3), m(4): the least m with 1 + 1/3 + ... + 1/(2m+1) > 2, 3, 4 (OEIS A092315) | CITED + FINITE-EXACT |
| 15, 113 | the denominators 2m+1 at those thresholds, i.e. the classical Borwein break (sinc(x/15)) and the 2cos(x) break (sinc(x/113)) | CITED |
| 80499 | the denominator at the threshold **2 pi**: 1 + 1/3 + ... + 1/80497 < 2 pi < 1 + 1/3 + ... + 1/80499 | FINITE-EXACT (interval arithmetic) |
| 231, 233 | mantissas of the first-failure deficits: pi/2 - **2.31**00572725 x 10^-11 at 1/15 and pi/2 - **2.33**24296128 x 10^-138 at 1/113 | PROVED closed form + FINITE-EXACT |
| 420 | the third deficit, at 80499, is **4.21**01759371 x 10^-226577: agrees to two digits; naive double-precision summation gives 4.208, which truncates to 4.20 | interpretation (two-digit match) |
| {4,2,5,3,7,13} | Session43's divisor-height cycle (provenance above); as a set it is {d+1 : d \| 12}, whose primes are vSC(12) = the genus-zero primes of X_0(p) = the p with p^3 \| \|M\| | CITED + PROVED |
| {2,3,11} | vSC(10) (B_10 = 5/66); the triangle group (2,3,11) has orbifold Euler characteristic exactly -B_10; every p+1 divides 12 | PROVED + CITED |
| {2,5,7} | every p-1 divides 12 (so it lies in vSC(12)); the unique subset of vSC(12) with prod(p+1) = 144 = sigma(70) = sigma(66) | PROVED + CITED |

The two halves of the directive meet in one principle: **a reciprocal sum compared with an integer.**
- The Borwein integrals break when a sum of reciprocals passes an integer (or 2 pi).
- von Staudt-Clausen says the reciprocal sum of the primes p with (p - 1) | 2k misses an integer by exactly B_2k.
- A sphere with cone points of those orders has that sum as its orbifold Euler characteristic.

Section 4 makes this exact: the Borwein excess of the von Staudt-Clausen sinc product is -B_2k.

## 2. THM-4837: the Borwein ladder as first non-vanishing points

Let P_m(x) = prod_{j=0}^{m} sinc(x/(2j+1)) and O_m = sum_{j=0}^{m} 1/(2j+1). Let U_j be independent and uniform on
[-1/(2j+1), 1/(2j+1)], and set S'_m = U_1 + ... + U_m. Then sinc(a x) = E e^{i a U x}, which gives:

- **(1)** For xi >= 0, F_m(xi) := int_R P_m(x) cos(xi x) dx = pi P(xi - 1 <= S'_m <= xi + 1). In particular, for xi >= 1,
  **F_m(xi) = 0 iff O_m <= xi.**
- **(2)** For t >= 0, I_m(t) := int_0^oo (sin(t x)/x) prod_{j=1}^{m} sinc(x/(2j+1)) dx = pi/2 - pi P(S'_m > t)
  = pi/2 - sum_{k>=0} F_m(t + 1 + 2k). So **I_m(t) = pi/2 iff O_m <= t + 1** (Borwein-Borwein 2001).
  - t = 1 is the classical integral.
  - t = 2 is the 2cos(x) variant, because 2 cos x sin x = sin 2x (Schmid 2014).
  - The telescoping identity is checked exactly for m <= 8.
- **(3) First failure in closed form.** For xi > 1 let m = m(xi) be the least m with O_m > xi, and put eps = O_m - xi.
  Then 0 < eps <= 1/(2m+1), and
  **F_m(xi) = pi eps^m (2m+1)!! / (2m)!!**.
  - Proof: minimality gives eps <= 1/(2m+1) <= 2 min_j 1/(2j+1). So only the corner of the box [-a_j, a_j] at
    (a_1, ..., a_m) reaches past xi - 1. That corner simplex has volume eps^m/m!, and P = eps^m/(m! prod 2a_j).
  - This is the first-failure case of Borwein-Borwein's formula. The point here is that the **first** failure is always
    single-corner.
- **(4) The owner's numbers.**

| threshold xi | m(xi) | 2m+1 | first non-zero value F_m(xi) = pi/2 - I_m(xi - 1) |
|---|---|---|---|
| 2 | **7** | **15** | **2.3100572725** x 10^-11 |
| 3 | **56** | **113** | **2.3324296128** x 10^-138 |
| 4 | **418** | 837 | 4.5195214277 x 10^-1320 |
| 2 pi | 40249 | **80499** | **4.2101759371** x 10^-226577 |

  - The thresholds 5, 6, 7 give m = 3091, 22845, 168803.
  - The value at 15 reproduces the published Borwein rational, 467807924713440738696537864469/935615849440640907310521750000 pi.
  - A direct quadrature on [0, 697] gives 2.31006e-11, with tail <= 4e-15.
  - The 2 pi row is certified with interval arithmetic: O_40248 - 2 pi <= -1.0075e-5 and O_40249 - 2 pi = +2.3472673396e-6.
  - The 2 pi row has a natural integral: **the cosine transform of the odd sinc product at one full turn,
    int_R cos(2 pi x) prod_{j=0}^{m} sinc(x/(2j+1)) dx, is exactly 0 for every denominator up to 80497 and first
    non-zero at 80499.**
- **(5) The crossing index.** Put X(xi) = exp(2 xi - gamma)/4 and M = m + 1. Then
  O_m = (1/2) ln(4M) + gamma/2 + r(M), with **0 < r(M) < 1/(48 M^2) for all M >= 1**.
  - The bound comes from Graham-Knuth-Patashnik's envelope H_N = ln N + gamma + 1/(2N) - 1/(12N^2) + 1/(120N^4) - eps_N,
    0 < eps_N < 1/(252 N^6). This gives r(M) = 1/(48M^2) - 7/(1920M^4) - eps_{2M} + eps_M/2.
  - Hence **m(xi) = floor(X(xi)) whenever frac(X) >= X/(24 floor(X)^2)**. The upper bound for M >= X and the lower bound
    at M = floor(X) bracket the crossing.

**Prior art credited.**
- Borwein and Borwein, "Some remarkable properties of sinc and related integrals", Ramanujan J. 5 (2001): the
  equality criterion and the first-failure value.
- H. Schmid, "Two curious integrals and a graphic proof", Elem. Math. 69 (2014): the 2cos(x) variant to 1/111.
- Baillie-Borwein-Borwein, Amer. Math. Monthly 115 (2008).
- J. Baez's 2018 post and the 3Blue1Brown video (2022) popularised 15 and 113.
- OEIS A092315:
  - J. Song (2022) already notes that a(2) = 7 and a(3) = 56 are the Borwein breaks.
  - G. Kirchner (2020) gives a(n) = floor(exp(2n - gamma)/4), with a precision check for n <= 1000.
  - M. F. Hasler (2017) conjectures a(n) = floor(exp(2n - gamma)/4 + 1/8).

**New here:**
- the 2 pi break at 80499 and its reading as the vanishing range of the cosine transform at one turn;
- the closed-form first-failure values at 113 and 80499, and the observation that the first failure is always
  single-corner;
- the decode of 231, 233 and 420;
- the rigorous window criterion (5).

**Hasler's +1/8 form is false.** The two OEIS formulas differ exactly when frac(X_n) >= 7/8. There criterion (5) proves
m(n) = floor(X_n), because frac(X_n) is far outside the window. The first counterexample is
**n = 17: a(17) = 81897532160124**, while floor(X + 1/8) = 81897532160125.
- An independent check uses the digamma form O_m = (psi(m + 3/2) - psi(1/2))/2:
  O_{a-1} - 17 = -5.71 x 10^-15 and O_a - 17 = +3.97 x 10^-16.
- The +1/8 form fails at 375 of the n <= 3000 (17, 31, 45, 57, 59, 63, 86, 118, ...).

## 3. HYP-9436: the crossing formula for every n (Comtet-type)

The conjecture is m(n) = floor(exp(2n - gamma)/4) for every n >= 2. It is the odd, half-shifted analogue of Comtet's
problem for H_k. It is open in general, because a failure needs exp(2n - gamma)/4 within about 1/(24X) above an integer.

Pre-registered (commit a14c71ec13, before the confirmatory run):
- **P1** (no window exception for 61 <= n <= 3000): CONFIRMED. The formula is now PROVED for 2 <= n <= 3000. The
  tightest case is n = 2, with frac/window = 101.8. n = 1 holds by direct check.
- **P2** (Kolmogorov-Smirnov uniformity of frac(X_n), n <= 3000): CONFIRMED, with D = 0.0192 against the 1% critical
  value 0.0298.
- **P3** (thresholds k pi, k <= 500): CONFIRMED. Values: m(pi) = 75 (151), m(2 pi) = 40249 (80499),
  m(3 pi) = 21553437 (43106875).

The heuristic expected number of exceptions over all n >= 2 is 0.0063.

## 4. THM-4838: the twelve, as reciprocal-sum thresholds

Write vSC(2k) = {p prime : (p - 1) | 2k}. von Staudt-Clausen gives N_k := B_2k + sum_{p in vSC(2k)} 1/p, an integer.
N_k = 1 for every 2k <= 12.

- **(1) Orbifold reading.** The sphere with cone points of orders vSC(2k) has orbifold Euler characteristic
  chi_k = 2 - |vSC(2k)| + sum 1/p, which is **= -B_2k (mod 1)** for every k. Equality is exact **iff 2k in {4, 6, 8, 10}**.
  - Those vSC sets are the triangle groups (2,3,5), (2,3,7), (2,3,5) and (2,3,11).
  - Exactness needs N_k = |vSC(2k)| - 2, which holds there and fails at 2k = 2, 12, 14, 16, 18.
  - For 2k >= 20, |B_2k| > 2 sqrt(2k) + 2 rules it out.
  - Checked for 2k <= 120.
- **(2) Triangle groups.** For a prime p >= 5, chi(2,3,p) = 1/p - 1/6. It **equals -B_(p-1) exactly iff p in {5, 7, 11}**,
  and **is = -B_(p-1) (mod 1) iff vSC(p - 1) = {2,3,p}** (p = 5, 7, 11, 23, 47, 59, 83, 107, ...).
  - Exactness forces |B_(p-1)| < 1/6, hence p - 1 <= 10. The case p = 3 fails by sign.
  - The congruence holds because a non-empty sum of reciprocals of distinct primes is never an integer.
- **(3) Modular curves.** Riemann-Hurwitz for X(p) -> X(1), branched (2,3,p), gives
  **2g(X(p)) - 2 = |PSL_2(F_p)| B_(p-1) for p = 5, 7, 11**: 60 (-1/30) = -2, 168 (1/42) = 4 and 660 (5/66) = 50.
  - So g = 0, 3, 26: the icosahedron, the Klein quartic and X(11).
  - These p are also Galois's exceptional primes, where PSL_2(F_p) acts transitively on p points (Galois 1832; Kostant,
    Notices AMS 1995). Both lists are "p <= 11" conditions, so the coincidence of {5, 7, 11} is recorded as a
    **DICTIONARY, not a theorem linking them**.
  - The repository's earlier observation "B_6 = 1/42 and the Hurwitz primes {2,3,7}" (`bernoulli_42_s116n.out`) is the
    p = 7 case. The Hurwitz in question is the automorphism theorem (84(g - 1) = 2(g - 1)/|chi(2,3,7)|), not the
    Hurwitz quaternion order.
- **(4) The three sets.**
  - **{4,2,5,3,7,13} = {d + 1 : d | 12}.** Its primes are vSC(12) = {2,3,5,7,13}, the denominator 2730 of
    B_12 = -691/2730. The same five primes are three other things:
    - the **genus-zero primes of X_0(p)** (classical; the eta quotient (eta(tau)/eta(p tau))^{24/(p-1)} is a
      Hauptmodul exactly when (p - 1) | 24);
    - the primes with **p^3 | |M|**, where the Monster has the extraspecial centralizers p^{1+24/(p-1)} (ATLAS);
    - the primes of exponent >= 2 in Session43's N = 2^5 3^7 5^3 7^4 11 13^2. This last item is an observation about one
      quotient, not a theorem.
  - **{2,3,11} = vSC(10).** B_10 = 5/66, E_10 = 1 - 264 sum sigma_9(n) q^n, and |im J| in pi_19^s is 264 = 2^3 3 11
    (Adams). Every p + 1 divides 12: {p : (p + 1) | 12} = {2,3,5,11}.
  - **{2,5,7}.** Every p - 1 divides 12. It is the unique subset of vSC(12) with prod(p + 1) = 144.
  - **{2,3,11} and {2,5,7} are the only three-prime sets with prod(p + 1) = 12^2.** That is, sigma(66) = sigma(70) = 144;
    sigma^{-1}(144) = {66, 70, 94, 115, 119} is Session83's bridge. The twelve thus appears three ways: (p + 1) | 12 for
    {2,3,11}, (p - 1) | 12 for {2,5,7}, and a product of 12^2 for both.
  - The parallel list {p : (p + 1) | 24} = {2,3,5,7,11,23} gives the M_24 eta products eta(tau)^{24/(p+1)} eta(p tau)^{24/(p+1)};
    the p = 11 member is the newform of X_0(11) (CITED: Mason; Dummit-Kisilevsky-McKay).
- **(5) Borwein meets von Staudt-Clausen.** For 2k <= 12, sum_{vSC(2k)} 1/p - 1 = -B_2k exactly. So
  int_0^oo sinc(x) prod_{p in vSC(2k)} sinc(x/p) dx is a Borwein integral with excess -B_2k:
  - it equals pi/2 exactly for 2k = 2, 6, 10 (B_2k > 0);
  - it falls short by **pi/43200** for 2k = 4, 8 (excess 1/30);
  - it falls short by **pi (691^5 - 271^5)/(3840 * 2730^4)** for 2k = 12 (excess 691/2730). The second corner enters
    because 2/13 = 420/2730 < 691/2730.
  - **NUMEROLOGY:** that this produces the owner's 420 (271 = 691 - 420, 420 = 2 * 2730/13) is arithmetic, not
    evidence.

## 5. {231, 233, 420} in prior repository work

Two readings sit alongside the Borwein mantissas.
- **Collatz.** 231 -> 233 is a recorded 17-step Syracuse word,
  231, 347, 521, 391, 587, 881, 661, 31, 47, 71, 107, 161, 121, 91, 137, 103, 155, 233.
  - It runs on the 27/17 clock (`collatz_procgen_20260922_order_laws.md`, `crossroads233_20260926_graph.md` section 5,
    `collatz_join_shields_and_lifts_20261004.md`).
  - The cylinder law m_17(231 + 2^28 t) = 233 + 2 * 3^17 t is re-checked here for t = 1, 2, 3.
- **Aliquot.** 231 = s(345) with sigma(345) = 24^2, and 420 = s(364) with sigma(364) = 28^2. Both are square-sum aliquot
  targets in Session83's census.
  - 364 is session84's least mixed balanced square-sigma source.
  - codex-incidence's absorber path 19236--364--420--19180 passes through it.
  - No n with square sigma has s(n) = 233 (exhaustive below 233^2).
  - So in the aliquot reading the owner's set splits into two square-sigma targets and one non-target.

**Typing.** The Borwein reading explains 231 and 233 to three digits from the owner's own 15 and 113. It explains 420
to two digits from the owner's 80499. The Collatz and aliquot readings are genuine prior-work occurrences, but they do
not explain why these three numbers go together.

## 6. What is new, cited and numerology

- **New to the repository:**
  - all of section 2 except the Borwein-Borwein and Schmid criteria;
  - the 80499 = 2 pi identification, the first-failure closed form at the owner's breaks, and the 231/233/420 decode;
  - the window criterion;
  - the refutation of the OEIS +1/8 conjecture at n = 17;
  - THM-4838 (1)-(3) and (5).

  A quick literature search found no statement of chi(2,3,p) = -B_(p-1) for p = 5, 7, 11, but it is probably folklore
  (the audit searches further).
- **Cited:** Borwein-Borwein, Schmid, GKP, von Staudt-Clausen, Riemann-Hurwitz, the genus of X_0(p), the ATLAS, Adams,
  the M_24 eta products, OEIS A092315, Session43, session84, codex-incidence, Session83.
- **Numerology:**
  - 420 = 2 * 2730/13 in (5);
  - any reading of the Borwein exponents (-11, -138) as the owner's labels.
- **Interpretation, not proof:** 420 as the third mantissa (two digits).

## 7. Reproduction

```text
python3 04-computation/experiments/borwein_ladder_vsc_twelve_20261028.py 3000 500
```

- Standard library plus mpmath and sympy; about 21 s.
- Every claim above is an assertion; the output is in `borwein_ladder_vsc_twelve_20261028.out`.
- The first argument is the n-range of HYP-9436 and the second the k-range of the k pi thresholds.

## 8. Audit record

(pending)
