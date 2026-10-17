# S39 independent audit (opus-2026-10-17-S39): THM-4663, THM-4664, HYP-9286, the dictionary note

Independent referee (subagent), 2026-10-17. Every claim re-derived with independent code, with a different method from
the session's wherever possible. Scripts and outputs: this folder (`a1` ... `a10`). (The referee could not write this
file itself; the coordinator transcribed its final report.) THM-4665 (EI17) was written after this audit and is audited
separately.

## Verdicts

| item | verdict | evidence |
|---|---|---|
| THM-4663 (1) spectrum | CONFIRMED | exact closed-walk counts (new field model, mod 3 primes + CRT) reproduce the session's moment file at all 428 primes 5 <= p <= 2999, k <= 8; 40-digit direct sums agree for the 76 primes < 400 (a1) |
| (2) W_2, W_3, W_4, codegree lemma | CONFIRMED | primes 5-59 and q = 9, 25, 27, 49, 121, 125, characteristic 3 included (a3, a10) |
| (3) M_5, M_6, M_8 closed forms | CONFIRMED WITH FIXES | forms rebuilt independently (a2; level-15 a_2, a_3, a_5 = 1, -3, 5; weight-6 newform by a U_2 route, a_2..a_11 = 4, -9, -66, 176, -60); all three closed forms at all 428 primes; M_5 holds at p = 5 too, so "(p >= 7)" is unnecessary; p^2 | S_7 + 1 also at p = 5, 7 |
| (4) pentagons | CONFIRMED | C_5(G_5) = 165 matches, so p >= 5; recounted by meet-in-the-middle at 15 primes <= 59 and networkx at p = 5, 7 |
| (5) hexagons | CONFIRMED WITH FIXES | C_6 and D(p) recounted at primes <= 59; every component of the proof checked by classifying all closed 6-walks at p = 5, 7, 11, 13, 17, 23; one exposition gap (fix 5) |
| (6) q = 49, 121, 169 | CONFIRMED (now exact) | the session's check was float-rounded; exact walk counts in F_(p^4) agree (a5) |
| (7) Bessel moments, p_5'(0) | CONFIRMED (>= 23 digits) | root number +1 for both forms certified by split-point independence to 1e-41; Bessel moments to 1e-40; p_5'(0) = p_4(1) computed without Bessel functions matches 6 L(f15,1)/pi^2 to 7.6e-24; attribution fixes |
| (8) theta(2,2,3) and 30a1 | CONFIRMED WITH FIXES | theta identities exact for p < 200; j, Delta and quartic invariants agree; the plane cubic is singular only mod 2, 3, 5; iota has a fixed point at p = 5, so p >= 7 is needed there |
| (9) R_3(c) = #E'_c, Gamma_1(6) | CONFIRMED WITH FIXES | exact for p < 200; fibre types of both models confirmed |
| (10) moments and twisted sums | identities CONFIRMED; K3 reading WRONG sign | all exact at every prime 5 <= p < 700 (the session's script stopped at 400) |
| THM-4664 (1) decomposition | CONFIRMED | verified exactly in Z[zeta_p] at every alpha != 0 frequency for p = 5, 7, 11; the alpha = 0 cases hold; the fibre parametrization is a bijection for every b (exhaustive, p <= 13) |
| THM-4664 (2) sizes and cone | CONFIRMED | floating point, margins >= 0.256 p^2; max |E| = 2.463p, max |S_Gamma| = 1.447p |
| THM-4664 (3) duality | CONFIRMED | discriminant of the fibre cubic = -c^2 c_3/4, c_3 = THM-4621's / THM-1310's L; top-degree part c(27a^2c + b^3); dual curve alpha^2 gamma = 4 beta^3 |
| THM-4664 "Open" paragraph | WRONG | Gao's G has 20, 24, 32 = 2(p - 1) coefficients with |S| > p^2/2, exactly the lines (0,0,gamma) and (alpha,0,0), not 340, 600, 384 (a8b, evaluation checked against exact substitution) |
| Collatz row | CONFIRMED WITH FIXES | trace formula for every prime 5 <= p < 62, k <= 20 (dynamic programme); Perron-Frobenius bound, the k = 20 fractions and tau_p agree |
| Keller row | CONFIRMED for d = 2, 3 | the actual THM-1300 and Gao maps on F_p^3 reproduce the session's histograms exactly; the d = 4, 5 rows are not Keller maps over Q |
| Legendre row | CONFIRMED | every prime p < 1000; already in the repo as HYP-9274 (Ahlgren-Ono Lemma 2.2 + Theorem 6) |
| HYP-9286 | status overclaims; (b) now well posed | its test case, the third moment, resolves positively |

## Addenda

* Twisted sums exact at every prime 5 <= p < 700; the derivation from W_5 checks algebraically.
* K3 reading: Euler number 2 + 12 + 4 + 6 = 24 (rank-19 trivial lattice); the trace of Frobenius on the 3-dimensional
  piece of H^2, computed directly (split/non-split signs from the Weierstrass nodes), is (p/3)p + a_p(f15) at all 121
  primes 11 <= p < 700; with -a_p it fails at 59. The transcendental trace is +a_p(f15).
* Third moment identified: sum_c a_c^3 = -(1 + (p/3))(4p + 2) - tr(T_p | S_5(Gamma_0(6), (-3/.))), the 2-dimensional
  space S_5(Gamma_1(6)); exact at all 121 primes, and a theorem by Eichler-Shimura-Deligne (a4c, a4d).
* Twisted second moment: lives on a 6-dimensional weight-3 space (observed size up to 2.92 p^(3/2)); "at most two forms"
  cannot settle it; OPEN.
* p_3 identity: provable in two lines (convolve the two-step arcsine density with one unit step; substitute
  t = r^2 + 1 - s^2; the radicand becomes Q_c(t) on exactly the bounded positive interval); matches the BSWZ closed form
  to 22 digits at seven radii (0.2-2.9). The committed script showed only three radii at ~6e-16.

## Required fixes (applied by the coordinator; see MISTAKE-624)

1-13 (THM-4663): status block lists statement 10, the p_3 identity and the mod-6 divisibility, and cites the exact
checks a4b, a5; "quartic nebentypus of conductor 105" -> "quartic nebentypus, level 525, coefficients in
Q(sqrt-1, sqrt6, sqrt14)"; drop the "p >= 7" provisos for M_5 and pentagons; credit the two level-6 Bessel identities as
Broadhurst's conjectures (eqs. 121-122, checked there to 1000 digits) and note that int x J_0^5 and p_5'(0) do not appear
there (p_5'(0) = p_4(1); NUMERICAL here); proof 5's count of admissible second triangles; proof 8's twist statement
(squarefree D with |D| | 30, either sign; witnesses <= 31 for all 15); mark the theta "no fit" as a session run; replace
"same surface with the cusps permuted" by the 2-isogeny L_c ~ E'_c (Phi_2 vanishes identically) and the isomorphism
L_c = E'_(sigma(c)), sigma(c) = (c - 9)/(c - 1) over Qbar(c) with no twist; 6 | R_3(c) PROVED by the torsion points
(0, 4(c-1)^2) (order 3) and (-(c-1)^2, 0) (order 2) on Y^2 = X^3 + (c+3)^2 X^2 + 8(c-1)^2(c+3) X + 16(c-1)^4; p_3 PROVED;
proof 10's "sums to 2p - 1" -> "sums to 2p; replacing the t = 2 term by R_2(0) = p + 1 gives 3p"; K3 sign +a_p;
third-moment identity replaces the non-fit sentence; Schoen -> Schoen's construction, modularity via Eichler-Shimura for
X_1(6) or Gouvea-Yui.
14-15 (THM-4664): Gao counts 2(p - 1); statement 2 typed FINITE (float64, margin >= 0.25 p^2), decomposition exact in
Z[zeta_p] for p = 5, 7, 11.
16-18 (HYP-9286): status scope; Legendre row CITED via HYP-9274; (b) excludes W_7 from "verified"; third moment
identified; untwisted four-bar moments modular by Eichler-Shimura; only twisted moments open.
19-20 (dictionary note): equidistribution only along k with 2^k in <3> mod p; (Boehm-Sontacchi; Lagarias); "finite-place
shadow" is ANALOGY; the characteristic polynomial is irreducible at p = 13, 23 and z^2 + z - 1 is not cyclotomic;
(b)<=>(d) by finiteness, (c)<=>(d) by Parseval; the value-set illustration is FINITE, d = 4, 5 are only the one-variable
family; the script tolerance was vacuous (gaps < 0.012).

## MISTAKE candidates

1. K3 sign (-a_p stated, +a_p true). 2. THM-4664 Open numbers irreproducible (true 20, 24, 32). 3. d = 4, 5 rows are not
Keller maps (Gao's side conditions fail: resultants 141/8000 and 13491239/25920000); THM-4625 (5)'s S_5, S_6 instances
inherit the problem. 4. R_3(1) arithmetic slip (2p, not 2p - 1). 5. Prior-art lapse (Legendre moment = HYP-9274).
6. Citation details (conductor 105; Broadhurst's level-6 conjectures). 7. Collatz characteristic polynomial description.

## Prior art and upgrades

No earlier repo work on Kloosterman moments, four-bar linkages, Gamma_1(6), Bessel moments or walk densities. Literature:
four-bar configuration curves are classically elliptic; the three-Bessel sunset family is the universal curve over X_1(6)
(Bloch-Vanhove), so statement 9 is its finite-field shadow. Upgrades: 6 and 9 and the theta identities now exact;
6 | R_3(c) and p_3 PROVED; the second moment follows from Eichler-Shimura for X_1(6) (S_4(Gamma_1(6)) is spanned by the
hexagon form), and run backwards through W_6 = sum R_3^2 this makes the M_6 closed form and the hexagon count hold for all
p >= 5 without citing Hulek-Spandaw-van Geemen-van Straten.
