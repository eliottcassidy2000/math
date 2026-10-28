# Referee audit: opus-2026-10-28-S67 (THM-4837, THM-4838, HYP-9436)

Independent adversarial referee, 2026-10-28. Worktree commit 7a5cd669d0 (origin/main has since added only the
S68 lane, which cites but does not change these files).

**VERDICT: SOUND WITH CORRECTIONS.** Every numerical claim reproduces with independent code, and the proofs are
correct apart from two small gaps (MINOR-5, MINOR-6). The **novelty typing of THM-4837 is wrong, though, and must be
fixed before canon.** The following are all published:
- the 80499 = 2 pi break and its error;
- the single-corner first-failure formula;
- the 1/113 deficit;
- the integer-threshold ladder;
- an explicit-remainder certification of the crossing formula.

Sources: Baillie-Borwein-Borwein 2008; Crandall 2008 via Baillie 2011; Borwein-Borwein 2001 Cor. 1;
Almkvist-Gustavsson 2014; Basel-Baillie 2016; Kirchner 2020 on OEIS.

Two of these sources (BB 2001 and BBB 2008) are cited by the session itself, so this is the recurrent "cited paper
not read" failure.

Files audited:
- `05-knowledge/results/borwein_ladder_owner_numbers_vsc_twelve_20261028.md` (the note);
- `01-canon/theorems/THM-4837-...md`;
- `01-canon/theorems/THM-4838-...md`;
- `05-knowledge/hypotheses/HYP-9436-...md`;
- `04-computation/experiments/borwein_ladder_vsc_twelve_20261028.py` and `.out`.

In this report, "note L" means a line of the note; THM-4837, THM-4838 and HYP-9436 line numbers are those of the
files at 7a5cd669d0.

---

## 1. Findings

### MAJOR-1. The 80499 = 2 pi break, its "natural integral" and its size are published (BBB 2008, Crandall, Baillie)

**Where:**
- note L43, L88-91, L110-111, L209-211;
- THM-4837 title (L3), L41-42, L62-63;
- HYP-9436 L22.

**Evidence.** Baillie, D. Borwein and J. M. Borwein, "Surprising sinc sums and integrals", Amer. Math. Monthly 115
(2008) 888-901 (preprint dated 6 Aug 2007, carmamaths.org/jon/sinc-sums.pdf).
- Their identity (6) is 1/2 + sum_{n>=1} prod_{k=0}^N sinc(n/(2k+1)) = int_0^oo prod_{k=0}^N sinc(x/(2k+1)) dx.
- They write that "(6) holds for every N <= 40248 and fails for all larger integers".
- Example 1(b): "If 7 <= N <= 40248, then sum a_k < 2 pi ... For N > 40248, the equality in the above formula
  fails." It adds that Crandall "rigorously estimates that the error for N = 40249 is minuscule: less than
  10^-226576".

R. Baillie, "Fun with very large numbers", arXiv:1105.3943 (2011, rev. 2015), Sec. 2:
- "holds for every N=0,1,2,..., 40248, but fails for N>=40249";
- "Crandall has shown [9, p. 24] that, at N=40249, the left side minus the right side of (1) is positive, but is only
  about 8.42*10^-226577";
- Sec. 11 prints 40248, 80497 and 2 pi - O_40248 = 1.00752...e-5.

The denominator at N = 40249 is 2*40249 + 1 = 80499.

**Mechanism.** By Poisson summation, sum_{n in Z} P_m(n) - int_R P_m = 2 sum_{k>=1} F_m(2 pi k). At N = 40249 only
k = 1 survives. The session's "cosine transform at one full turn" F_m(2 pi) is therefore exactly BBB's one-sided
sum-minus-integral error.

**Audit numbers** (`audit_borwein_numbers.out` section 2):
- one-sided error F_m(2 pi) = 4.2101759371e-226577 (rigorous enclosure; agrees with the session);
- two-sided error = 8.4203518743e-226577.

Baillie attributes "about 8.42e-226577" to the one-sided identity (1), so that report is off by a factor 2: it is the
two-sided value. The normalisation was cross-checked against BBB's exact Example 4. For sinc^7,
s_7 + 1/2 - i_7 = 2 pi (7 - 2 pi)^6/(6! 2^7) = F(2 pi), the single cosine term, and this is reproduced symbolically.

**Correction.**
- Note L43, typing column: replace "FINITE-EXACT (interval arithmetic)" with "CITED (BBB 2008 Ex. 1(b); Crandall
  2008 via Baillie 2011) + FINITE-EXACT (re-certified)".
- Note L89-91: append "This is the k = 1 Poisson term. It is the error in BBB's sinc-sum identity
  1/2 + sum_{n>=1} prod sinc(n/(2k+1)) = int_0^oo prod sinc(x/(2k+1)) dx, which holds for N <= 40248 and first fails
  at N = 40249, denominator 80499 (BBB 2008, Example 1(b)). Crandall's value, as reported by Baillie 2011
  (8.42e-226577), is the two-sided error sum_Z - int_R. The one-sided error is F_m(2 pi) = 4.2101759371e-226577."
- Note L111 and THM-4837 L62: delete "the 2 pi break at 80499 and its reading ..." / "the 2 pi/80499 identification"
  from NEW. Replace it with: "KNOWN (BBB 2008; Crandall; Baillie). New here only: the 11-digit closed-form value of the
  one-sided error, and the remark that the literature's 8.42e-226577 is the two-sided normalisation."

### MAJOR-2. The single-corner first-failure formula, the 1/113 value and the integer ladder are published

**Where:**
- note L69-75 ("The point here is that the **first** failure is always single-corner"), L85, L110-113, L209-211;
- THM-4837 title, L32-35, L62-63.

**Evidence.**
- **Borwein-Borwein, Ramanujan J. 5 (2001), Corollary 1 (p. 79).**
  - Hypotheses: "If 2a_k >= a_n > 0 for k = 0, 1, ..., n-1 and sum_{k=1}^n a_k > a_0 >= sum_{k=1}^{n-1} a_k".
  - Conclusion: int_0^oo prod_{k=0}^n sin(a_k x)/x dx = (pi/2) prod a_k - pi (a_1 + ... + a_n - a_0)^n/(2^n n!).
  - The hypotheses say precisely that n is the first failure. For the decreasing sequence 1/(2j+1), with any
    a_0 >= a_n/2, they hold at every first failure. That is the session's "observation".
- **BBB 2008, Theorem 2 ("First Bite").** Same statement; Example 1(a) applies it at 1/15.
- **Wikipedia "Borwein integral", General formula.** "...which means that n is the first value when the partial sum
  of the first n elements of the sequence exceed a_0", with 0 < a_n < 2a_k and the same C_n.
- **Almkvist-Gustavsson, "More remarkable sinc integrals and sums", arXiv:1405.1265 (2014).**
  - Item 2 computes the 1/113 failure exactly, giving "1 - 1.484870809*10^-138". Multiplied by pi/2 this is
    2.33243e-138.
  - Item 3 gives the thresholds 5 and 7 (value 1 for n <= 3090 and for n <= 168802) for
    (2cos(pi t) + 2cos(3 pi t) [+ 2cos(5 pi t)]) prod sinc(pi t/(2k+1)). These are exactly the session's I_m(4) and
    I_m(6), since (2cos x + 2cos 3x) sin x = sin 4x.
- **Wikipedia.** "≈ π/2 − 2.3324 × 10^−138" for 1/113.
- **Basel-Baillie, "Sinc integrals and tiny numbers", Elem. Math. 71 (2016) 7-20, arXiv:1510.03200.**
  - They study J_n(a_0) = a_0 int_0^oo prod_{k=0}^n sinc(a_k x) dx with a_k = 1/(2k-1), for a_0 = 1..10.
  - Their crossing indices for a_0 = 2..10 are n = 8, 57, 419, 3092, 22846, 168804, 1247298, 9216354, 68100151
    (for a_0 = 1, n = 2 and t = 1/12). These are the session's m(T) + 1 for T = a_0. The ar5iv extraction doubled
    some digits ("88", "5757"); the values were restored from n = m + 1 and the stated t.
  - They give single-corner deficits t_n = (s(n) - a_0)^n (2n)!/(2^(2n-1) (n!)^2) to 40 digits, at every threshold.

**What survives as new** (all trivial instances of BB Cor. 1):
- the value 4.5195214277e-1320 at 837 in this exact family;
- the 11-digit one-sided value at 80499.

**Correction.**
- Note L74-75: replace with "This is exactly Borwein-Borwein's Corollary 1 (2001), also BBB 2008 Theorem 2 'First
  Bite'. Its hypotheses state that m is the first failure. For the decreasing sequence 1/(2j+1) they hold at every
  first failure, as Basel-Baillie (2016) use at every integer threshold."
- Note L112-113 and THM-4837 L62-63: remove "the closed-form first-failure values at 113 and 80499, and the observation
  that the first failure is always single-corner" from NEW.
- Prior-art list (note L99-108; THM-4837 external): add Almkvist-Gustavsson 2014, Basel-Baillie 2016,
  Crandall 2008 and Baillie 2011.
- THM-4837 title: keep the formula, but drop its implied novelty.

### MINOR-3. The "rigorous window criterion" was anticipated by Kirchner 2020 (OEIS), whose test is under-described

**Where:** note L107, L115; THM-4837 L63; HYP-9436 L72.

**Evidence.** OEIS A092315, comment by Gerhard Kirchner (May 21 2020), quoted exactly:
- "Another formula, see link Asymptotic formulas, formula 1, is OH(m) = (log(4*m)+gamma)/2+1/(2*m)-11/(48*m^2)+1/(8*m^3)-127*t/(1920*m^4), 0<t<1."
- "The test can be carried out with t=0." / "Additionally, the precision can be tested by checking if t=1 makes a
  difference."
- "For n<=1000, the case "Precision too low" does not occur."

His Maxima program checks floor(OH(k,0)) = floor(OH(k,1)) = n and floor(OH(k-1,0)) = n-1 at k = floor(exp(2n-gamma)/4).
That is a two-sided bound certificate of a(n) = floor(X) for every n <= 1000. In his attachment (a092315.pdf; text
extracted with pdftotext) the bound 0 < t < 1 is supported numerically ("t(10) = 0.955, t(100) = 0.995,
t(1000) = 0.9995"), not proved. Basel-Baillie 2016 and Baillie 2011 likewise locate such thresholds with explicit
Euler-Maclaurin error bounds.

**What is new:**
- a fully proved remainder (0 < r(M) < 1/(48M^2), from the GKP/DLMF envelope);
- the clean sufficient condition frac(X) >= X/(24 floor(X)^2);
- the extension from n <= 1000 to n <= 3000 and to k pi, k <= 500.

**Correction.**
- Note L107 and HYP-9436 L72: "G. Kirchner (2020) gives a(n) = floor(exp(2n - gamma)/4) with an explicit-remainder
  test (formula 1, remainder -127t/(1920 m^4), 0 < t < 1, supported numerically) that certifies it for n <= 1000."
- Note L115 and THM-4837 L63: "the window criterion: a fully proved version of Kirchner's 2020 test, extending the
  certified range to n <= 3000 and to k pi, k <= 500."

### MINOR-4. The "+1/8 refutation" is already visible in the OEIS entry's own DATA

**Where:** note L117-122, L213; THM-4837 title, L54-58, L64; HYP-9436 L73-75.

**Evidence.**
- The DATA line of A092315 lists a(17) = 81897532160124. The entry records extensions by Hasler (Jan 24 2017) and
  Kirchner (Jul 23 2020), and Kirchner's b-file covers n = 1..1000. Which extension added a(17) was not checked.
- floor(X_17 + 1/8) = 81897532160125. Hasler's formula line ("a(n) = floor(exp(2*n-gamma)/4+1/8) for all n >= 1
  (conjectured; see also comments in A002387)") therefore contradicts the entry's own data.
- Audit: all 22 DATA terms equal floor(X_n). Among them the +1/8 form fails only at n = 17.
- a(17) is certified here two ways:
  - by the DLMF 5.11(ii) envelope in interval arithmetic: O_a - 17 = +3.9681161e-16 and O_(a-1) - 17 =
    -5.7083785e-15;
  - by digamma.
- The count of 375 failures for n <= 3000 is confirmed.

**Correction.** Retype as "explicit observation: the entry's DATA (a(17)) contradicts the '+1/8' formula line; certified
independently; 375 failures for n <= 3000". It is not a refutation of an open conjecture. Recommend submitting an
OEIS edit to A092315.

### MINOR-5. THM-4838 (2): with scope p >= 3, the mod-1 "iff" is false at p = 3

**Where:** THM-4838 L23-26. The note (L150) correctly says p >= 5.

**Evidence.**
- vSC(2) = {2,3} = {2,3,p} for p = 3.
- But chi(2,3,3) = 1/6 and -B_2 = -1/6, whose difference is 1/3, not an integer. The prime 3 is double-counted.
- See `audit_vsc_monster.out` section 2.

**Correction.** THM-4838 L23: "For primes p >= 5, chi(2,3,p) = 1/p - 1/6. (For p = 3 both the equality, by sign, and the
congruence, because 3 is counted twice, fail.)"

### MINOR-6. The |B_2k| growth step behind THM-4838 (1) and (2) is asserted, not proved

**Where:**
- THM-4838 L21 ("For 2k >= 20, |B_2k| > 2 sqrt(2k) + 2 >= the possible range") and L24 ("exactness needs
  |B_(p-1)| < 1/6, so p - 1 <= 10");
- note L147-149;
- the script checks only 2k <= 120.

**Evidence.** The argument is correct but omits two steps:
- (a) exactness forces B_2k = |vSC| - 2 - sum 1/p, which lies in [|vSC|/2 - 2, |vSC| - 2). Hence
  |B_2k| <= max(|vSC| - 2, 2) <= d(2k) + 2 <= 2 sqrt(2k) + 2.
- (b) monotone growth of |B_2k|.

Audit: |B_2k| > 2 sqrt(2k) + 2 for 20 <= 2k <= 200, and min |B_(2k+2)/B_2k| = 11.70 for 2k >= 20.

**Correction.** Add to THM-4838 (1): "Exactness forces B_2k = |vSC| - 2 - sum 1/p, so |B_2k| <= d(2k) + 2 <= 2 sqrt(2k) + 2.
But |B_20| = 529.1 and |B_(2k+2)/B_2k| = (2k+1)(2k+2) zeta(2k+2)/(4 pi^2 zeta(2k)) >= 11.7 for 2k >= 20, while
2 sqrt(2k) + 2 grows by a factor < 1.05 per step. The same ratio bound (>= 4.6 for 2k >= 12) gives |B_2k| >= |B_12|
> 1/6 for 2k >= 12 in (2)."

### MINOR-7. The congruence list in THM-4838 (2) is OEIS A092307 (denominator of B_(p-1) = 6p)

**Where:** note L150-153; THM-4838 L25-26.

**Evidence.**
- A092307: "Primes p such that there are no primes q, 3 < q < p, such that (q-1) divides (p-1)." Its terms are
  5, 7, 11, 23, 47, 59, 83, 107, 167, 179, 227, 239, 263, 347, 359, 383, ..., identical to the session's list.
- Comments: "Primes p such that 6p is the denominator of some Bernoulli number" (T. D. Noe, Sep 26 2006) and
  "Primes p such that A027642(p-1) = 6p" (Eldar and Ordowski, May 26 2021).
- Since chi(2,3,p) = (6-p)/(6p), the congruence is exactly denom(B_(p-1)) = 6p.

**Correction.** Cite A092307 at note L151 and THM-4838 L25. Type the congruence half as "CITED (A092307) + one-line
reformulation". The exact half (p in {5,7,11}) stays new to the repository. As in the session's search, no
published statement was found (see section 4).

### MINOR-8. Repository prior work behind THM-4838 (1)-(3) is under-credited

**Where:** note L160-162, L214-217; THM-4838 L32-33, L55-56 (credits only `bernoulli_42_s116n.out`).

**Evidence:**
- `05-knowledge/results/bernoulli_tournament_primes_s18h.out` L180-182: "denom(B_6) = 42 = 2*3*7 = HURWITZ NUMBER"
  (the B_6/Hurwitz coincidence). It also tabulates the denominators 30, 42, 66 with primes {2,3,5}, {2,3,7},
  {2,3,11} (L42-46).
- `05-knowledge/results/number_tournament_bridges_20261014.out` L67-75 (S33): chi(2,3,k) for k = 3..8, including
  1/30 for (2,3,5) and -1/42 for (2,3,7) with "Hurwitz: |Aut| <= (2g-2)/|chi| = 84(g-1)".
- `05-knowledge/results/understand_2357_s90bl.out` L177: "For B_10: denom = 66 = 2*3*11. Primes = {2,3,11}."
- `bernoulli_42_s116n.out` L86-96 already uses the {d+1 : d | n} device (for n = 6), which THM-4838 (4) uses for 12.

**Correction.** THM-4838 Novelty: "Ingredients were already in the repository: s18h (denom B_6 = 42 = Hurwitz), S33
`number_tournament_bridges_20261014.out` (chi(2,3,k) table), s90bl (vSC(10) = {2,3,11}), s116n ({d+1 : d | 6}). The
exact statements (1)-(3) and (5) are new to the repository."

### MINOR-9. The interpretive decode of 231, 233 and 420 is typed too strongly

**Where:** note L44 (typing "PROVED closed form + FINITE-EXACT"), L45, L203-204 ("explains"); THM-4837 title (L3),
L43-45.

**Evidence.** The values are proved, but the reading "owner's 231 and 233" is an interpretation: three-significant-digit
matches with published numbers (2.31e-11 in BB 2001 and Wikipedia; 2.3324e-138 in Almkvist-Gustavsson and Wikipedia).

The 420 match is weak (`audit_borwein_numbers.out` section 3):

| summation | mantissa | truncated to 3 digits |
|---|---|---|
| true value | 4.2101759 | 4.21 |
| naive left-to-right double | 4.20837 (session's claim reproduced) | 4.20 |
| naive right-to-left double | 4.20991 | 4.20 |
| math.fsum | 4.21016 | 4.21 |
| numpy pairwise | 4.21023 | 4.21 |

So "420" needs both naive summation and truncation, a post-hoc mechanism. The only published number for this break
is Crandall/Baillie's 8.42e-226577.

**Correction.**
- Note L44 typing: "INTERPRETATION (3-digit match to published deficits)".
- Note L45: "SPECULATIVE (2-digit match; 4.20 needs a naive double sum truncated, not rounded)".
- L203-204: "explains" becomes "matches".
- THM-4837 title and L43-45: drop "the owner's 231, 233 and, to two digits, 420", or mark it
  "(interpretation; 420 speculative)".

### MINOR-10. "The twelve behind the owner's sets" mixes proved facts with a numerological reading

**Where:** THM-4838 title (L3) and L35-45; note L46-48, L163-178.

**Evidence.** Every individual statement is true and re-verified, but none of them explains the owner's sets:
- {4,2,5,3,7,13} enters (Session43) as the prime-to-exponent arrows of |M|/|B|. The label 4 is v_7 and 2 is v_13.
  So "= {d+1 : d | 12}" is a set coincidence with no mechanism.
- {p : (p+1) | 12} = {2,3,5,11} is strictly larger than {2,3,11}.
- {p : (p-1) | 12} = {2,3,5,7,13} is strictly larger than {2,5,7}.
- The only exact characterisation of the two triples is prod(p+1) = 144 (Session83's complete fibre).

**Correction.** In the note table (L46-48) and the THM-4838 title, type the unifying reading "the twelve behind
{4,2,5,3,7,13}, {2,3,11}, {2,5,7}" as NUMEROLOGY/INTERPRETATION. The component facts remain PROVED/CITED.

### NOTES

- **NOTE-11. THM-4837 (1) fails at the edge m = 0, xi = 1.**
  - Where: THM-4837 L25; note L62-63.
  - F_0(1) = int sinc(x) cos(x) dx = pi/2 (conditionally convergent), not 0 as "iff O_m <= xi" claims (O_0 = 1).
    The formula pi P(...) would give pi.
  - Fix: state (1) for m >= 1, or for xi > 1.
- **NOTE-12. Premature status.** The statuses say "independently audited (see the note's section 8)" (note L11;
  THM-4837 L4; THM-4838 L4), but section 8 still reads "(pending)" (note L237). Update after this audit.
- **NOTE-13. HYP-9436 L48:** "about 10^-52" should be "about 4 x 10^-54" (audit: 3.56e-54). P1-P3 were near-certain
  predictions, so CONFIRMED adds little evidence beyond the certificate itself. The pre-registration order is correct
  (section 6).
- **NOTE-14. The session's own certificate was not interval arithmetic.** It used floating point with a heuristic
  10^5 guard (script L241-256). It was re-certified here with mpmath interval arithmetic for 2 <= n <= 3000 and
  k <= 500: no exception, tightest case n = 2 with frac/window = 101.84. "PROVED for 2 <= n <= 3000" stands.
- **NOTE-15. Monster wording.** Note L168 and THM-4838 L40 say "extraspecial centralizers p^{1+24/(p-1)}".
  - In the ATLAS/Wikipedia list, 2^{1+24}.Co1 is C(2B). For p = 3, 5, 7, 13 the listed subgroups are the normalisers
    N(pB) = p^{1+24/(p-1)}.(...). 13^{1+2} is a Sylow 13-subgroup.
  - Cite Conway-Norton 1979 for the pB Hauptmoduln (eta(tau)/eta(p tau))^{24/(p-1)} + 24/(p-1). The genus-zero /
    extraspecial link is classical moonshine.
- **NOTE-16. The Hurwitz quaternion "correction" goes too far** (note L161-162; THM-4838 L33).
  - The (2,3,7) group is, modulo the centre, the norm-one group of the *Hurwitz quaternion order* over Z[2cos(2pi/7)]
    (Wikipedia "Hurwitz quaternion order"; Katz-Schaps-Vishne).
  - Rephrase: "the relevant theorem is the 84(g-1) bound; the s116n phrase points, at best, to the Hurwitz quaternion
    order over Q(cos 2pi/7), whose unit group gives (2,3,7), and not to the rational Hurwitz quaternions."
- **NOTE-17. Credit for (1)-(2).** Also cite:
  - Majumdar-Trizac, PRL 123 (2019) 020201 (the random-walk reading);
  - Schmid 2014 eqs. (12) and (19) (cosine factor = frequency shift of the box-spline transform, general cos(a_0 t)).

  The repository's THM-4392 (LRC box-spline/Poisson duality) uses the same principle.
- **NOTE-18. THM-4838 (5) context.** BBB 2008 Example 2 already treats sinc products over prime reciprocals (odd
  primes 3..29; all primes with a_0 = 1/2, value pi(1 - 1/240), re-checked here). Cite it as context. The vSC version
  is new to the repository.
- **NOTE-19. OEIS status of HYP-9436.**
  - Kirchner's FORMULA line "a(n) = floor(exp(2*n-gamma)/4)." carries no "conjectured" tag. No proof for all n was
    found anywhere, so HYP-9436 "OPEN in general" is correct.
  - Jianing Song's comment writes the sinc arguments inverted, sin((2k+1)x)/((2k+1)x). The note does not copy this.
- **NOTE-20. Process.** BBB 2008 is in the session's prior-art list (note L102; THM-4837 L13), and it contains
  40248/40249 explicitly. Recommend a MISTAKES entry of the "cited paper not read fully" type (cf. MISTAKE-690).

---

## 2. Numeric claims: session value vs independent value

Scripts:
- A = `audit_borwein_numbers.py` (exact rationals by binary splitting; 100-digit rational pi bounds; mpmath interval
  arithmetic; BB 2001 sign-sum formula);
- V = `audit_vsc_monster.py`;
- C = `audit_collatz_aliquot.py`.

| # | claim (location) | session | audit value | method | status |
|---|---|---|---|---|---|
| 1 | m(2) = 7, 2m+1 = 15; O_7 - 2 (note L80, .out L4) | 0.021800422 | 0.0218004218004 | A exact | OK |
| 2 | m(3) = 56, 113; O - 3 | 0.0032870588 | 0.00328705881207 | A exact | OK |
| 3 | m(4) = 418, 837; O - 4 | 0.00069059164 | 0.000690591638868 | A exact | OK |
| 4 | m(5) = 3091 (6183) | 4.1720043e-5 | 4.17200434812e-5 | A exact | OK |
| 5 | m(6) = 22845 (45691) | 2.0676198e-5 | 2.06761975732e-5 | A exact | OK |
| 6 | m(7) = 168803 (337607) | 1.7917843e-6 | 1.79178429102e-6 | A exact | OK |
| 7 | m(2 pi) = 40249 (80499); O_40248 - 2pi, O_40249 - 2pi (note L88) | <= -1.0075e-5; +2.3472673396e-6 | -1.00752472258e-5; +2.34726733956e-6 | A exact + rational pi | OK |
| 8 | m(pi) = 75 (151) (HYP L69) | 75 | 75 (O_75 - pi = +5.53e-3, O_74 - pi = -1.09e-3) | A exact | OK |
| 9 | m(3 pi) = 21553437 (43106875) | 21553437 | 21553437 (O_m - 3pi = +1.1175e-9, O_(m-1) - 3pi = -2.208e-8) | A fsum with error bound + DLMF intervals | OK |
| 10 | first failure at 15 (note L80) | 2.3100572725e-11 | 2.3100572725481e-11 | A exact eps | OK (published) |
| 11 | first failure at 113 | 2.3324296128e-138 | 2.3324296128422e-138 | A exact eps | OK (published: A-G 2014, Wikipedia) |
| 12 | first failure at 837 | 4.5195214277e-1320 | 4.5195214276882e-1320 | A exact eps; exact and loggamma routes agree to 4e-78 | OK |
| 13 | first failure at 80499 | 4.2101759371e-226577 | enclosure [4.2101759371451, 4.2101759371451]e-226577; two-sided 8.4203518743e-226577 | A | OK (cf. Crandall/Baillie 8.42e-226577, two-sided) |
| 14 | BB rational at 1/15 (note L86) | 467807924713440738696537864469/935615849440640907310521750000 pi | equal; also BBB's 6879714958723010531/935615849440640907310521750000 deficit; also via the 2^7-term BB Thm 2 sign sum | A | OK |
| 15 | naive double sum gives mantissa 4.208 (note L45) | 4.208 | 4.20837 forward; 4.20991 reverse; 4.21016 fsum; 4.21023 numpy | A | reproduced, order-dependent (MINOR-9) |
| 16 | direct quadrature 2.31006e-11, tail <= 4e-15 (note L87) | yes | not re-run; anchored by the published rational (row 14). Normalisation of (1) checked by an independent quadrature, m = 4 (A section 6) | A | OK |
| 17 | 0 < r(M) < 1/(48M^2), r = 1/(48M^2) - 7/(1920M^4) - eps_2M + eps_M/2 (note L93-95) | yes | algebra by sympy; direct 50-digit check M = 1..3000; GKP/DLMF bound checked N = 1..199 | A | OK |
| 18 | window test: no exception for 2 <= n <= 3000, tightest n = 2, ratio 101.8 (HYP L65) | yes | yes, by interval arithmetic; ratio 101.84 | A | OK |
| 19 | KS D = 0.0192 vs 0.0298 (HYP L67) | 0.0192 | 0.0192 | A | OK |
| 20 | expected exceptions 0.0063; "about 10^-52" for n >= 61 (HYP L41, L48) | 0.0063; 1e-52 | 0.006288; 3.56e-54 | A | first OK; second NOTE-13 |
| 21 | k pi, k <= 500: no window exception | yes | yes (interval) | A | OK |
| 22 | +1/8 fails at 375 of n <= 3000; first 17, 31, 45, 57, 59, 63, 86, 118 | yes | 375; same list | A (interval) | OK |
| 23 | a(17) = 81897532160124, floor(X+1/8) = ...125, O_(a-1) - 17 = -5.71e-15, O_a - 17 = +3.97e-16 | yes | -5.7083785e-15, +3.9681161e-16 (DLMF envelope + intervals; digamma agrees) | A | OK |
| 24 | OEIS DATA a(1..22) = floor(X_n) | n/a | all 22 equal; the +1/8 form disagrees with DATA at n = 17 | A | MINOR-4 |
| 25 | chi_k = -B_2k (mod 1) for all k; exact iff 2k in {4,6,8,10}; N_k = 1 for 2k <= 12 | 2k <= 120 | 2k <= 200; N_14, N_16, N_18, N_20 = 2, -6, 56, -528 | V | OK |
| 26 | for 2k >= 20, |B_2k| > 2 sqrt(2k) + 2 | asserted | true for 20 <= 2k <= 200; ratio >= 11.70 | V | OK (proof gap MINOR-6) |
| 27 | chi(2,3,p) = -B_(p-1) iff p in {5,7,11}; mod-1 iff vSC(p-1) = {2,3,p} | yes | p < 2000: yes for p >= 5; p = 3 counterexample to THM-4838's scope | V | MINOR-5 |
| 28 | 2g(X(p)) - 2 = |PSL2(F_p)| B_(p-1): 60(-1/30) = -2, 168(1/42) = 4, 660(5/66) = 50; g = 0, 3, 26 | yes | yes (Riemann-Hurwitz and genus formula) | V | OK |
| 29 | genus-zero X_0(p) = {2,3,5,7,13} = {p:(p-1)|24} = {p:(p-1)|12} = vSC(12) | yes | yes for p < 5000 | V | OK |
| 30 | |M|, |B| factorisations; p^3 | |M| iff p in {2,3,5,7,13} | yes | factored from the decimal orders (Wikipedia): equal | V | OK |
| 31 | small-prime core of |M|/|B| = 2^5 3^7 5^3 7^4 11 13^2 (Session43 N) | yes | |M|/|B| = 2^5 3^7 5^3 7^4 11 13^2 * 29*41*59*71; core = 39,046,339,332,000 | V | OK ("core" drops 29, 41, 59, 71) |
| 32 | extraspecial p^{1+24/(p-1)} (ATLAS) | yes | 2^{1+24}.Co1, 3^{1+12}.2Suz.2, 5^{1+6}:2J2:4, 7^{1+4}:(3 x 2S7), 13^{1+2}:(3 x 4S4) (Wikipedia maximal subgroups) | web | OK (wording NOTE-15) |
| 33 | |im J| = 264 (pi_19), 65520 (pi_23); E_10 = 1 - 264 sum sigma_9 q^n | yes | denominators of B_2k/(4k): 24, 240, 504, 480, 264, 65520 | V | OK |
| 34 | M24 eta list {p : (p+1) | 24} = {2,3,5,7,11,23} | yes | yes | V | OK (CITED) |
| 35 | vSC sinc deficits: 0 (2k = 2, 6, 10), pi/43200 (4, 8), pi(691^5 - 271^5)/(3840*2730^4) (12) | yes | identical, via the BB 2001 Thm 2 sign-sum formula (independent of the session's code) | V | OK |
| 36 | sigma^{-1}(144) = {66,70,94,115,119}; only triples {2,3,11}, {2,5,7}; unique subset of vSC(12) | yes | yes | V | OK |
| 37 | {d+1 : d | 12} = {2,3,4,5,7,13}; {p:(p+1)|12} = {2,3,5,11} | yes | yes | V | OK |
| 38 | Syracuse word 231 -> 233 (17 steps), 27/17 clock; m_17(231 + 2^28 t) = 233 + 2*3^17 t | t = 1..3 | t = 0..2000 plus two large t; 27 halvings; the 2^28 modulus is needed (231 + 2^27 changes the word) | C | OK |
| 39 | s(345) = 231 with sigma = 24^2; s(364) = 420 with sigma = 28^2; no square-sigma preimage of 233 | yes | yes. 233 has 7 aliquot preimages, none with square sigma. Complete square-sigma preimages: 231 <- {345, 6169}; 420 <- {364} | C | OK |

The session's own script was also re-run. Its `.out` reproduces exactly, apart from timing lines (runtime 20.3 s).

---

## 3. Proof checks (line by line)

**THM-4837**
- **(1) Fourier formula.** Correct for m >= 1.
  - sinc(a x) = E e^{i a U x} with U uniform on [-1,1].
  - The density of U_0 + S' is (1/2) P(|y - S'| <= 1).
  - Fourier inversion gives F_m(y) = 2 pi g(y) = pi P(y-1 <= S' <= y+1).
  - The vanishing criterion follows because S' has positive density on (-A, A).
  - Independent check: quadrature at m = 4, xi = 3/2 and 7/4 agrees within the tail bound 3.1e-7.
  - Edge case m = 0, xi = 1: NOTE-11.
- **(2)** I_m(t) = (pi/2) P(|S'| < t) = pi/2 - pi P(S' > t) (Dirichlet integral; S' continuous and symmetric for m >= 1).
  The telescoping is additivity of probability. Correct; it is BB 2001 Thm 1 / Schmid.
- **(3) First failure.** Minimality gives eps <= 1/(2m+1) = a_m <= 2 min_j a_j. So the region
  {sum V_j <= eps} lies in the box, and P = eps^m/(m! prod 2a_j). Correct. It is BB 2001 Corollary 1 (MAJOR-2).
- **(5) Envelope.**
  - Writing O_m = H_{2M} - H_M/2 and inserting the GKP envelope gives exactly
    r(M) = 1/(48M^2) - 7/(1920M^4) - eps_2M + eps_M/2 (sympy-verified).
  - The upper bound holds because eps_M/2 < 1/(504M^6) < 7/(1920M^4) for M >= 1.
  - The lower bound holds because r > 1/(48M^2) - 7/(1920M^4) - 1/(16128M^6) > 0 for M >= 1. At M = 1,
    r = 0.01825 < 1/48.
  - The GKP bound 0 < eps_N < 1/(252N^6) follows from DLMF 5.11.2 with the enveloping property of Sec. 5.11(ii)
    (real z > 0). I could not consult GKP itself, so the equation number (9.89) is unverified, but the bound is
    confirmed numerically for N = 1..199.
- **Window criterion, both halves.**
  - M = floor(X) + 1 > X gives O > xi, so m <= floor(X).
  - frac >= X/(24F^2) gives ln(X/F) >= frac/X >= 1/(24F^2), so O_(F-1) < xi.
  - Hence m = floor(X). The "otherwise floor(X) - 1 or floor(X)" clause is also correct.
- **(6)** Correct; see MINOR-4.

**THM-4838**
- **(1)** chi_k + B_2k = 2 - |vSC| + N_k is an integer. Exactness holds iff N_k = |vSC| - 2. The cases 2..18 are
  correct. The growth step is missing (MINOR-6).
- **(2)** chi(2,3,p) + B_(p-1) = N - 1 - sum_R 1/q with R = vSC(p-1) \ {2,3,p}. The step "a non-empty sum of
  reciprocals of distinct primes is not an integer" is correct. The scope is wrong at p = 3 (MINOR-5), and the growth
  step for |B| >= 1/6 is again missing (MINOR-6).
- **(3)** Riemann-Hurwitz for the Galois cover X(p) -> X(1) with group PSL2(F_p), branched (2,3,p), gives
  2g - 2 = |G|(1/6 - 1/p) = |G| B_(p-1) exactly when (2) holds. Correct. The Galois-primes item is correctly labelled
  a dictionary.
- **(4)** Every statement is true (V sections 3-5). The typing of the unifying reading needs correction (MINOR-10).
- **(5)** N_k = 1 for 2k <= 12 gives excess -B_2k. At 2k = 12 two corners enter, since 2/13 < 691/2730 < 2/7, and the
  formula is confirmed. Correct.

**HYP-9436.** The statement and the "what would refute" section are consistent with THM-4837 (5). The statuses are
correct after re-certification (NOTE-14).

---

## 4. Prior-art search log (2026-10-28)

### 4.1 Searches and fetches

**OEIS**
- `oeis.org/search?q=id:A092315&fmt=text` and `oeis.org/A092315`. Quoted:
  - FORMULA lines: "a(n) = floor(exp(2*n-gamma)/4+1/8) for all n >= 1 (conjectured; see also comments in A002387)"
    (M. F. Hasler, Jan 22 2017); "a(n) = floor(exp(2*n-gamma)/4)." (Gerhard Kirchner, Jul 23 2020).
  - Kirchner comment (May 21 2020): quoted in MINOR-3.
  - Song comment (Dec 10 2022): "a(2) = 7 and a(3) = 56 are related to the Borwein integrals." The integrals are
    written with sin((2*k+1)*x)/((2*k+1)*x), i.e. inverted arguments.
  - DATA: a(1..22) as listed in `audit_borwein_numbers.py`.
- `oeis.org/A092315/a092315.pdf` (Kirchner, "Asymptotic formulas", 2020-05-21): formula 1, with t fitted numerically.
- `oeis.org/A002387`. For the H_n analogue, floor(exp(n-gamma) + 1/2) is "conjectured" (Cloitre 2002), with a
  heuristic by Hickerson (2003). Comtet and Boas-Wrench appear only in the references; this audit did not verify
  those references.
- `oeis.org/A092307`: MINOR-7.

**Wikipedia, "Borwein integral."** It gives "π/2 − 2.31 × 10^−11", the 1/15 rational, "≈ π/2 − 2.3324 × 10^−138"
(the 2cos variant at 113), and the General formula with the first-failure statement (MAJOR-2). It does not mention
80499 or 2 pi.

**Web searches leading to the sinc-sum literature**
- "80499 sinc sum integral Borwein": BBB sinc-sums.pdf; arXiv:1405.1265; others.
- "Baillie Borwein Borwein Surprising sinc sums ... first fails": the "(6) holds for every N <= 40248" sentence.
- "226577 sinc OR Borwein OR Crandall": Baillie arXiv:1105.3943 ("about 8.42*10^-226577") and Basel-Baillie
  arXiv:1510.03200.
- "Borwein integral 2cos x 1/113 2.3324 10^-138 Schmid": Almkvist-Gustavsson arXiv:1405.1265, Wikipedia, and a
  golem.ph.utexas.edu post.
- "Crandall Theory of ROOF walks": the preprint was not retrievable (cited via BBB and Baillie).

**Papers read** (PDF text extracted to scratch only; no PDF stored in the repository). pdftotext drops some math
symbols, so in the quotes above "<=", "pi" and some superscripts were restored by hand from context (cf. MISTAKE-684).
- BBB 2008 preprint: Theorems 1-2, Examples 1, 2 and 4.
- BB 2001 (gwern copy): Theorem 1, Corollary 1, Theorem 3, Corollary 2, Remarks 2 / eq. (10).
- Schmid 2014 (gwern copy): eqs. (7)-(12), (19).
- Baillie 2011 (ar5iv): Secs. 1-3, 11.
- Basel-Baillie 2016 (arXiv abstract and ar5iv): Table 1, t_n formula.
- Almkvist-Gustavsson 2014 (ar5iv): Theorem 1, eqs. (10)-(12), items 1-4.
- Majumdar-Trizac 2019 (arXiv:1906.04545 abstract).
- Baez, "Patterns that eventually fail" (Azimuth, 2018): has 2.31e-11 and the 1/111 and 1/113 Schmid/Egan
  discussion. It does not have 80499, 10^-138 or Crandall.

**THM-4838 searches** (no prior statement of chi(2,3,p) = -B_(p-1), or of chi(S^2; vSC(2k)) = -B_2k mod 1, was found)
- "Baez 42 Bernoulli number B_6 1/42 Hurwitz automorphism 84 von Staudt";
- "mathoverflow triangle group Euler characteristic Bernoulli number coincidence (2,3,7) -1/42 ...";
- "Hurwitz B_6 1/42 coincidence Bernoulli Klein quartic 168 Euler characteristic";
- "site:mathoverflow.net Hurwitz 84 Bernoulli 42 von Staudt-Clausen triangle group";
- "chi(2,3,7) ... -1/42 Bernoulli B_6";
- "orbifold Euler characteristic von Staudt Bernoulli cone points sphere". The closest related item is Patchkoria
  2024 (orbifold Euler characteristics of mapping class groups, Harer-Zagier, von Staudt congruences), which is not
  this statement.
- "Galois exceptional primes 5 7 11 ... Bernoulli": found A092307.
- "genus zero X_0(p) ... 2730 / denominator of B_12": the genus-zero list is classical. No source states the
  vSC(12) coincidence, which is elementary.

**THM-4838 reference checks**
- "Monster group extraspecial 1+24/(p-1) ...": Conway-Norton 1979 context; Duncan-Swisher, "Modular functions and the
  monstrous exponents" (relevant to "p^3 | |M|").
- Wikipedia "Monster group" and "Baby monster group": orders and the list of maximal subgroups.
- Wikipedia "Hurwitz quaternion order": NOTE-16.
- DLMF 5.11 (eq. 5.11.2 and Sec. 5.11(ii)): supports the GKP envelope.

### 4.2 Novelty typing required

| item | session typing | required typing |
|---|---|---|
| 80499 = 2 pi break and its integral | NEW | CITED (BBB 2008 Ex. 1(b); Baillie 2011) |
| value at 80499 | NEW | CITED as ~8.42e-226577, two-sided (Crandall via Baillie). The 11-digit one-sided value and the factor-2 remark are new |
| value at 113 | NEW | CITED (Almkvist-Gustavsson 2014; Wikipedia) |
| single-corner first failure | NEW | CITED (BB 2001 Cor. 1; BBB 2008 Thm 2; Wikipedia; Basel-Baillie 2016) |
| thresholds 4-7 / crossing indices | (implicitly new) | CITED (Almkvist-Gustavsson 2014 for 5, 7; Basel-Baillie 2016 for 2-10 as n = m+1; OEIS DATA) |
| window criterion | NEW | rigorous version of Kirchner 2020; new range n <= 3000, k pi |
| +1/8 refutation | NEW | explicit observation of a contradiction already in the OEIS DATA |
| THM-4838 (1), (2) exact, (3), (5) | new to repo, folklore-level | stays, plus the A092307 citation and the repo-ingredient credits (MINOR-7, MINOR-8) |
| decode of 231/233/420 | NEW | INTERPRETATION; 420 SPECULATIVE |

---

## 5. Repository-context checks

| claim | file and line | result |
|---|---|---|
| Session43 defines the divisor-height cycle 4 -> 2 -> 5 -> 3 -> 7 -> 4 with 13 -> 2, and N = 2^5 3^7 5^3 7^4 11 13^2 as the small-prime core of \|M\|/\|B\| | `05-knowledge/results/session43_observers_phase_guards_20261023/cycle-geometry/PROOF.md` L9-15, L35-41 | CONFIRMED |
| codex-incidence recovers the labels from Session43. In its chart {2,5,7} and {4,7,13} are circuits and {2,3,11} is a parallel pair plus a loop. It has a Fano/Hamming extension and the absorber path 19236--364--420--19180 | `05-knowledge/results/incidence_attractors_20261028/README.md` L8-18, L31-34, L50-52 | CONFIRMED |
| session84 has the bracket set {2,3,11} with interval (25/13, 25/11], 66 by a Berggren route, the F = U + S frontier, and 364 as the least mixed balanced square-sigma source | `05-knowledge/results/session84_guarded_unification_20261028/BRACKETS_AND_66.md` L1-24; `FUS_DIVISOR_FRONTIER.md` L12, L203; `README.md` L52 | CONFIRMED |
| `bernoulli_42_s116n.out` contains the "Hurwitz primes" observation | L67-72, L86-96 | CONFIRMED. Its "max order quaternion algebra" phrase is ambiguous (NOTE-16) |
| Session83's census lists (231,345) and (420,364); sigma^{-1}(144) bridge | `05-knowledge/results/session83_aliquot_divisor_grammar_20261028/CENSUS_OUTPUT.md` L110 (class III pairs, target/source); `CARRIER_AND_144_BRIDGE.md` L60; `README.md` L75 | CONFIRMED |
| Collatz 231 -> 233 on the 27/17 clock | `collatz_procgen_20260922_order_laws.md` L104; `crossroads233_20260926_graph.md` Sec. 5 (L224-262); `collatz_join_shields_and_lifts_20261004.md` L42, L109 | CONFIRMED |
| no earlier Borwein-integral or A092315 work in the repository | grep of 00-navigation, 01-canon, 04-computation, 05-knowledge for A092315, Borwein integral, Borwein-Borwein, x/113, 1/111, 80499, odd harmonic | CONFIRMED. Only incidental matches: Erdos-Borwein constant, Borwein-Choi, Bailey-Borwein-Crandall-Pomerance, LRC sinc products (THM-879, THM-4392). For THM-4838, earlier ingredients exist (MINOR-8) |
| S68 lane (same owner message) | commits e458b68aa2, 593105aedd, 4e32454073 (THM-4839) | exists; cites THM-4837/4838 without overlap |

---

## 6. Pre-registration check (HYP-9436)

- **Commit order:**
  - reservation 464ed0c713 at 10:02:09;
  - pre-registration a14c71ec13 at 10:09:21, containing the HYP file and the exploratory script, with no `.out`;
  - checkpoint 7a5cd669d0 at 10:19:06, adding the `.out` and the results.
- `git diff a14c71ec13 7a5cd669d0` on the HYP changes only the status line and appends the Results section. The text
  of P1-P3 is unchanged.
- The script diff only adds sections C3 (Hasler/digamma) and E2 (genus/Monster). This is disclosed at HYP L77-78.
- The pre-registration is valid. Its predictions were low-risk (NOTE-13).

---

## 7. How to rerun

From `04-computation/experiments/borwein_ladder_vsc_twelve_20261028_audit_S67/`, with Python 3, mpmath, sympy and numpy:

```text
python audit_borwein_numbers.py 3000 500    # ~10 s; sections 1-6, ends "FAILURES: none"
python audit_vsc_monster.py                 # ~5 s
python audit_collatz_aliquot.py             # ~5 s
```

The outputs are committed alongside the scripts as `audit_borwein_numbers.out`, `audit_vsc_monster.out` and
`audit_collatz_aliquot.out`. All three report "FAILURES: none".

The scripts take no data from the session code. The only external inputs are hard-coded:
- 100 decimals of pi, checked against mpmath;
- the 22 OEIS DATA terms;
- the decimal orders of M and B from Wikipedia.
