# S91 (2026-11-10/11): the local densities of F1, the numbers of an Euler-Sombor abstract, qualifying path counts and source lines

Session `opus-2026-11-10-S91`. Canon: THM-5128, THM-5129. Hypotheses: HYP-9660, HYP-9661. Scripts and outputs:
`04-computation/experiments/s91_f1_local_densities_20261110/`. The Chan-Pak conjecture (every odd number except 7
and 21 is a Hamiltonian path count) is NOT proved here and remains OPEN. This session determines the local
(residue class) densities of the best-understood join family, F1, at each fixed odd prime, each squarefree odd
modulus and every power of 2, with a power saving that is uniform in the region but not in the modulus;
odd prime-power moduli are numerical only, and uniformity in the modulus is OPEN (it needs at least HYP-9660 (B)).
That no congruence excludes any odd number, 7 and 21 included, was KNOWN (the one-apex family alone attains every
residue modulo every odd modulus and every odd residue modulo every even modulus: Session 25; Session193 modular_sync Theorem 1; Session206 exact-modular-image); what is added is
the size of the local densities of F1.

Independent audit 2026-11-11 (seven auditors; record `audit/AUDIT_REPORT.md` in this directory, their scripts and
outputs in `audit/<auditor>/`): PASS WITH CORRECTIONS, all applied; MISTAKE-761. No stated theorem was found false.
Receipts quoted below as `audit/...` are relative to this directory; script and output paths (`bridges/...`,
`rooted_spectrum/...`, `rates_exact/...`, `blind_proof/...`) are relative to the scripts directory above; notes of
other sessions are under `05-knowledge/results/` (Session206: `session206/` at the repository root).

Owner directive: "work a session on the local densities of F1 and see how prior work on the numbers 2, 3, 5, 7, 13,
19, 27 and 1039 relate deeply with symmetry structure we have been studying regarding tournaments and primes and
constructions like collatz which relate deeply to prime geometry. think of qualifying path counts in k-vertex
extensions of tournaments and how chan-pak can be proven by modifying that to Hamiltonian path counts in
tournaments themselves", with four attached papers (section 3).

## 0. Inheritance and portfolio

*Anchor:* the by-value local densities of F1 (obligation 2 of the S90 note). *Niche:* k-vertex extensions and rooted
("qualifying") counts, in the restricted sense of section 4.1 (the words were defined first, and more generally, in
odd_c_seed_symmetry_20261110/EXTENSIONS.md). *Wildcard:* the four attached papers and the eight numbers.

Inherited, and credited throughout (the startup sweep found the first three, which S90 had missed: MISTAKE-758;
the rows on Session206's weighted capacity, the S82 prediction, the congruence statement, the 2-adic balance and
the fans were added by the independent audit: MISTAKE-761):

| What | Where | Status there |
| --- | --- | --- |
| Two favoured classes `0` and `-1`, confinement at `d = U = 0 (mod p)`, limits `1 + (3p-4)/(4p^2)`, `1 - 1/p^2 + 1/(4p)` | THM-4916 (kind-pasteur S15, 2026-11-07), letter-count model, one root change | PROVED; by-value limit "NOT CLAIMED" |
| Two root changes: (A) the census factors tend to the letter-count limits (conjecture); (B) "Fermat boost" at census scale; (C) weighting words by value closes the gap between the letter-count automaton and the census (conjecture) | HYP-9506 (A), (B), (C) | OPEN; pre-registered R2 FAILED |
| 3-adic law (the factor depends on the valuation); "a joint real/residue transfer theorem ... is still required" | Session223 | PROVED at 3; the transfer theorem missing |
| Digit chain on `(v mod 2p, Q(v) mod p)`, uniform stationary law; by-value equidistribution of the one-apex count | Session206 | PROVED for fixed `p` |
| By-value two-class law of a fixed `p`-divisible seed, with the weights `c_0`, `c_1` of the two parities (confinement to the classes `h`, `h - e`) | Session206 weighted-parameter-capacity; THM-4888 (C); HYP-9480 (A) | PROVED for a fixed seed |
| Bank family F2 at `p = 3`: the local factor comes entirely from the triples with 3 dividing `y_h` and `y_R`, split by the parity of `v`; a-priori prediction against the census, with the weights `2 w^1`, `2 w^0` of THM-5128 (D2) | S82 note section 5; `04-computation/experiments/chanpak_mass_exploratory_20261106/local_factor_mod3_prediction.py` | NUMERICAL (prediction against census) |
| Mixing of `H mod M`; joint limit `(UV, R_M)` | Session193 modular_sync Theorems 3, 4; THM-4886 (C) | PROVED |
| No congruence excludes an odd number: the one-apex family attains every residue modulo every odd modulus and every odd residue modulo every even modulus | Session 25; Session193 modular_sync Theorem 1; Session206 exact-modular-image | PROVED |
| Twisted transfer matrices of the count; rates `0.9320` at `1/7` (extremal), next worst `11/93, 1/31, 1/15, 1/63`, `sqrt(3)/2` at `1/3` | S82 note section 7 | VERIFIED, `q <= 120` |
| Readout mass `(X/a) ln(1 + a/c)`, mean `(pi^2/12) ln N` for unordered pairs | THM-4885 | PROVED |
| 2-adic balance of F1; equidistribution of the one-apex count modulo `2^j`; `H mod 4` as an `F_2` polynomial | S78 note, S82 note section 5; Session192, Session193 modular_sync Theorem 3, THM-4886 (C); THM-466, THM-478 | NUMERICAL (observation); PROVED; PROVED |
| `x = y` exactly at `h = 0` (the transitive tournament rooted at its sink) and at the fans, the only unit-offset one-apex roots; the count `2^m + 1` of one back arc; the boost of the primes dividing `2^i + 1` through the multiplier tails | THM-5126 (A) (kind-pasteur S21); Moon 1972, THM-4748 (5), THM-4861; THM-4916 (1c), HYP-9506 (B) | PROVED; KNOWN; PROVED (the identity (1c)), NUMERICAL (the share of the boost carried by these tails, (B)) |
| By-value class means at 13 primes; presence / absence of the excess | HYP-9657 (A) and its audit (S90) | NUMERICAL |
| Mod-4 quadratic form of a vertex cube | THM-5127 (C) (mac-mini, 2026-11-10) | PROVED |

## 1. The local theory of F1 in one page (THM-5128)

Rooted state of the one-apex tournament with word `h`: `x = 2 psi(h) + 1` (its path count), `y = 2h + 1`. F1 is
`N = x1 x2 + 2^s y1 y2 - [s = 0]`: a bilinear form in two rooted states, on two sheets (`s = 0`, `s >= 1`).

1. **One biased class, carried by divisible states.** For any two independent states modulo `p` and a nondegenerate
   bilinear gluing of rank `r`, the class law is exactly
   `1/p + (pi_1 + pi_2 - p^(-r)) ([c = 0] - 1/p) + (second order)`, `pi_i` the mass of states that vanish modulo `p`
   (the identity is a standard character-sum computation; what is used is its reading for joins and the bound on
   the remainder). For F1 and an odd prime `p` (`r = 2`) a state vanishes modulo `p` when `p | gcd(x, y)`; then `N`
   is forced into the class `0` (`s >= 1`) or `-1` (`s = 0`), whatever the other factor is. Everything else is flat
   to first order.
2. **By value the law is first order up to a power saving.** Counting triples with `N <= X`:
   `share(c) = 1/p + ([c = c_s] - 1/p) F(X) + O(X^(-delta))`, `c_s` the distinguished class (`0` for `s >= 1`, `-1`
   for `s = 0`), `F(X)` the share of triples whose smaller state is p-divisible. The proof cuts the larger word
   into cells (top digits fix the size, bottom digits fix `x/y`, the middle digits mix the residues at the rate of
   the twisted transfer operator). The exponent `delta` and the constant depend on `p`: the saving is uniform in
   the region, not in the modulus (uniformity in the modulus is HYP-9660 (B), OPEN). KNOWN before: the two-class
   by-value law of a fixed p-divisible seed (Session206 weighted-parameter-capacity; THM-4888 (C); HYP-9480 (A))
   and, for the bank family F2 at `p = 3`, the first-order prediction with the same weights (S82 note section 5).
   New: the summation over the smaller state with a power-saving error, for F1 and every odd prime, and the
   constants.
3. **The limit is generic, the approach logarithmic.** `F(X) -> 1/p^2`, so the by-value class factors modulo an
   odd prime are `p mu_c = 1 - 1/p^2 + [c = c_s]/p`, and for squarefree odd `m` the limiting density of a class
   modulo `m` is the product of these (PROVED). The product over the primes dividing `N + [s = 0]`,
   `(8/pi^2) prod_{p | N + [s = 0], p odd} (1 + p/(p^2 - 1))`, is *the squarefree-level product (a class average)*,
   not the singular series: inside the distinguished class modulo `p` the classes modulo `p^2` are not equally
   likely.
   - *Prime powers (exact for uniform states; by value NUMERICAL).* For independent equidistributed states modulo
     `p^k` the value of the form has the factor `(1 - p^-2) sigma(p^v)/p^v = (1 + 1/p)(1 - p^(-v-1))` at a class of
     valuation `v = v_p(c - c_s) < k` (an exact elementary count; enumerated by the audit at 9, 25, 27 in
     `audit/proof-5128-CF/a15_padic_factor_exact.out` and at 9, 25, 27, 49 in
     `audit/citations/check_singular_series_prime_powers.out`); averaged over the distinguished class modulo `p`
     (all `v >= 1`) it gives the class factor `1 - 1/p^2 + 1/p`. The by-value limit modulo `p^k` is NUMERICAL:
     audit census at `X = 2^26`, modulo 9, sheet 1, `9 * share = 0.848, 1.229, 1.456` for `v = 0`, `1`,
     `2 or more`, against the limits `0.889, 1.185, 1.296`; the first-order law at 9 predicts
     `0.852, 1.216, 1.457` (`audit/proof-5128-CF/a4_moduli_26.out`). It would follow from the argument of (C), (D2)
     carried out with (B4) for the modulus `p^k`, which is not written out. The singular series proper would then
     be `(8/pi^2) sigma(n)/n`, `n` the odd part of `N + [s = 0]` (mean 1, minimum `8/pi^2`): `8/pi^2` times the
     abundancy index of the odd part. Over both sheets the density is `w_0 S_0 + w_1 S_1` (`w_0`, `w_1` the shares
     of the two sheets), not an Euler product.
   - *The approach.* A single p-divisible state `h` keeps a share of order `1/ln X` of the mass. So at every
     feasible scale the profile is decided by *how early and how densely the p-divisible states occur*.
   - `x = y` exactly at `h = 0` and `h = 2^j`, where `x = y = 2^(j+1) + 1` (KNOWN: these are the transitive
     tournament rooted at its sink and the fans, the only unit-offset one-apex roots, THM-5126 (A) with
     `delta' = 2 q(h)`; the count `2^m + 1` of one back arc is Moon's, THM-4748 (5), THM-4861; the boost of the
     primes dividing `2^i + 1` through the multiplier tails is THM-4916 (1c), HYP-9506 (B). New: the reading of the
     fans as p-divisible states of the first-order law, the density `1/p^2` of p-divisible states and the table of
     first divisible states). Every prime factor of `2^m + 1` therefore has the divisible state `h = 2^(m-1)`:
     `3` (h = 1), `5` (2), `17` (8), `11` (16), `13` (32), `43` (64), `257` (128), `19` (at 256, and already at
     218), `41` (at 512, and already at 266), `65537` (2^15). That state is below the generic scale `p^2` only for
     small `m`, and having it does not imply a boost: what matters is the whole mass of divisible states. Of these
     primes `3, 5, 11, 17, 43, 257, 65537` show an excess far above the generic value (`41` a mild one:
     `+0.032 / +0.031` against `0.024`); `13` and `19` do not. At
     `13` the balanced state `h = 32` is followed by sparse ones (32, then 721, 877); at `19` the divisible states
     start late (218, then 256, 294, 313, 351, 959: the first at `0.6 p^2`). In both cases the weighted mass of
     divisible states stays below `1/p^2` (`A_13 = -0.0038, -0.0049`, `A_19 = -0.00014, -0.00032`:
     `audit/recompute-5128/a5_Ap.out`; the states: `audit/recompute-5128/a4_divisible.out`) and the excess is at
     or below the generic value (`+0.0441 / +0.0446` against `0.077`; `+0.0510 / +0.0495` against `0.053`).
   - Primes whose first divisible state is late (`29`: 1435, `37`: 1350, `31`: 8664, `61`, `73`, `127`: 109283) show
     none, or a deficit: without divisible states the distinguished class is *under*-represented by `1/p`
     (two non-zero vectors of `(Z/p)^2` are orthogonal with probability `1/(p + 1)`, not `1/p`). Even order of 2
     does not make a prime early: `29, 37, 61` divide `2^14 + 1`, `2^18 + 1`, `2^30 + 1`, whose states `2^(m-1)`
     are far beyond `p^2`.
4. **The prime 2 is exact.** Modulo `2^j`, for odd characters, the twisted transfer operator is `sqrt 2` times a
   unitary matrix (`T_k T_k^* = 2I`; an even character is a character of `Z/2^(j-1)`): square-root cancellation
   with no loss, `|sum_{h < 2^n} e(k q(h)/2^j)| <= 2^((n + j + 1)/2)`, and no 2-adic factor: by value every odd
   class modulo `2^j` receives `2^(1-j) #T + O_j(X^(7/8))` of the `#T` triples with `N <= X` of each sheet
   (statement and proof supplied by the independent audit; by value the states themselves are not equidistributed
   modulo `2^j`: the bottom digits of `h` fix both `y mod 2^j` and the size of `x/y`; what is uniform is `x` given
   `y`). The mod-4 law of the one-apex family is a band quadratic form whose symbol over `F_2` is
   `(1 + z)^2 (1 + z + z^2)`; its corank has period 6 in the word length, and the signed sums `S_n`, `D_n` have
   period 12 up to the factor `-64` (`S_(n+12) = -64 S_n`, `D_(n+12) = -64 D_n`, from `M^12 = -64 I` for the `4 x 4`
   signed transfer matrix of `q mod 2`; added by the audit). KNOWN before: the 2-adic balance of F1 as an
   observation (S78 note; S82 note section 5); equidistribution of the one-apex count modulo `2^j` (Session192;
   Session193 modular_sync Theorem 3; THM-4886 (C)); `H mod 4` as an `F_2` polynomial (THM-466, THM-478). New: the
   unitarity with the bound `2^((n+j+1)/2)`, the symbol and the corank.
5. **Tested before looking.** Predictions for the unseen window `[2^24, 2^27)` were computed from the list of
   divisible states and the uniform law of `x/y` alone (the method of the S82 prediction for F2 at `p = 3`), and
   pre-registered (commit 1f97e27805): 8 tests, 279 intervals, all CONFIRMED (re-evaluated by the audit, also with
   the measured standard errors at 4, 3 and 2 se: all inside; `audit/recompute-5128/a10_prereg_full.out`). Three
   primes never used before were predicted to show an excess. The Fermat prime 65537 shows `+0.14 / +0.17`
   (predicted `+0.18 / +0.21`; standard error `0.022 / 0.021` from the measured dispersion, where the registered
   tolerance used `0.029`); its neighbours 65521 and 65539 show nothing. The Wagstaff prime 2731 shows
   `+0.063 / +0.075` (predicted `+0.061 / +0.069`). The Wagstaff prime 43691 shows `+0.018 / +0.046` (predicted
   `+0.061 / +0.069`, standard error `0.018 / 0.017`): its registered interval contains 0, so this prime is not
   resolved (`241` likewise). The constant `0.75` in the registered `se = sqrt(0.75 p / #N)` was a conservative
   bound on the relative variance (measured `0.45 / 0.39` on the window), so the registered "4 se" bands were
   about 5.2 to 5.5 measured standard errors wide.

| p | first divisible state | rate per digit | excess, sheet 0: measured / predicted | sheet 1: measured / predicted | generic limit |
| --- | --- | --- | --- | --- | --- |
| 3 | 1 | 0.8660 | +0.539 / +0.515 | +0.544 / +0.532 | 0.375 |
| 5 | 2 | 0.7418 | +0.3130 / +0.3120 | +0.3280 / +0.3284 | 0.208 |
| 7 | 38 | 0.9320 | +0.233 / +0.174 | +0.230 / +0.175 | 0.146 |
| 11 | 16 | 0.8587 | +0.1690 / +0.1705 | +0.1707 / +0.1749 | 0.092 |
| 13 | 32 | 0.8272 | +0.0441 / +0.0459 | +0.0446 / +0.0448 | 0.077 |
| 17 | 8 | 0.8289 | +0.2107 / +0.2081 | +0.2296 / +0.2288 | 0.059 |
| 19 | 218 | 0.8053 | +0.0510 / +0.0511 | +0.0495 / +0.0497 | 0.053 |
| 29 | 1435 | 0.7894 | +0.0051 / +0.0032 | +0.0020 / +0.0007 | 0.035 |
| 31 | 8664 | 0.9002 | +0.0041 / -0.0061 | -0.0038 / -0.0089 | 0.032 |
| 43 | 64 | 0.8342 | +0.0885 / +0.0877 | +0.0930 / +0.0941 | 0.023 |
| 257 | 128 | 0.8551 | +0.1859 / +0.1815 | +0.2023 / +0.2044 | 0.004 |
| 65537 | 32768 | - | +0.144 / +0.183 | +0.168 / +0.206 | 0.00002 |

(Window `[2^24, 2^27)`; excess = distinguished class over the mean of the other classes; "predicted" used no window
data; generic limit `p/(p^2 - 1)`. At 3, 7, 31 the difference is the second-order term, which the theorem bounds by
a power of `X`. Rates: `rates_exact/rates_by_modulus.txt`, 262 moduli: every odd `m <= 521`, and 1023, 1143; 1039
is in `bridges/b3_f1_data.out`; the rate at 65537 was not computed. `13` and `19` are at or below their generic
limits although they have balanced divisible states: item 3.)

## 2. What this settles, and what it corrects

- *HYP-9657 (A)* ("the classes 0 and -1 are favoured at some primes and not at others; the size is not `1/p`"):
  its first-order part is explained: presence, absence and size are the mass of p-divisible states. The S90 file
  overstated its novelty. The two classes and the confinement mechanism were already PROVED by value for a fixed
  seed in Session206 (weighted-parameter-capacity; THM-4888 (C); HYP-9480 (A)), predicted for the bank family F2
  at `p = 3` and compared with the census in the S82 note, section 5, and PROVED in the letter-count model in
  THM-4916: corrected there (MISTAKE-758 for THM-4916, Session223 and HYP-9506; the by-value credits were added by
  the S91 audit, MISTAKE-761).
- *The census question of THM-4916 / HYP-9506 for one root change.* By value the limit has the same shape as the
  letter-count law but a different split of the boost `1/p` between the classes `0` and `-1`:
  `(0.80359, 0.19641)` for word pairs against `(3/4, 1/4)`. The constants are dilogarithm values
  (`w_0 = (pi^2/12 + 2 Li_2(-1/2) - Li_2(-1/4)) / (pi^2/12)`), because `x/y` is uniform on `[1/2, 1]` and the sum over
  the padding telescopes. So for one root change the by-value limit is not the letter-count limit (the `K = 1`
  analogue of HYP-9506 (A) fails in the split), while a value-weighted law does reproduce the census (the `K = 1`
  analogue of HYP-9506 (C)); the excess of the census over both limits is the `1/ln X` term.
- *The Kronecker question of S90* ("is the mass a class-number-type quantity?"): not at first order. In the limit
  the factor of a class modulo `p` depends only on whether it is the distinguished class (`p | N`, or `p | N + 1`
  on the sheet `s = 0`); the p-adic factor of an integer grows with `v_p(N)` (or `v_p(N + 1)`): measured at 3 on
  `[2^22, 2^24)`, sheet 1, `1.22` for `3 || N` and `1.45` for `9 | N` (NUMERICAL;
  `audit/prior-art-wording/a1_prime_power_factor.out`). At finite size the first-order term depends on the 2-adic
  order structure of `p` through the divisible states; the second-order term (visible at 3, 7, 31) is not
  analysed, and at 7 it cannot be told apart from a quadratic character, the squares modulo 7 being the powers of
  2 (S90).
- *The "Fermat boost"* (HYP-9506 (B), observed on the census of two root changes): for one root change it is the
  first-order term of THM-5128 (C), (D2), quantified, and decaying like `1/ln X`; for two root changes the
  first-order term is not derived (HYP-9660 (A)).

## 3. The numbers 2, 3, 5, 7, 13, 19, 27, 1039

Two sessions treated the same eight numbers first: odd_c_seed_symmetry_20261110/PRIMES_AND_PAPERS.md (commit
4098631a54, 2026-11-10: the table of `ord_p(2)` for `3, 5, 7, 13, 19, 1039`, the even / odd split with `<2>` = the
squares at 7 and 1039, `F_27` against `Z/27Z`, the Paley bridge from 7 to 27 (`189/7 = 27`), a reading of the
multiset-dimension paper) and kind-pasteur S22 (kps_s22_sign_flip_towers_20261110/README.md section 8: a typed
table of the eight numbers; `1040 = 2^4 * 5 * 13` typed there as a coincidence). New here: the identification of
the list with numerals of the Euler-Sombor abstract, the row 27, and the reading of the split through the
divisible states of THM-5128.

**Where they come from (EXACT arithmetic; the statements about the paper are CITED).** All eight are printed in
the abstract of J. Allikvere, "Minimal Trees for the Euler Sombor Index Revisited: A Correction and the First
Thirteen Minima", MATCH Commun. Math. Comput. Chem. 97 (2027) 1281-1293, doi:10.46793/match.97-3.28526 (the
article's own header misprints the year as 2026), one of the four attached papers; the abstract also prints
`0, 4, 10, 11, 18` and `204,990`. More precisely the list is the set of integers occurring in
`delta = sqrt 13 + 2 sqrt 3 - sqrt 7 - sqrt 19`, `eps = 2 sqrt 19 - 5 sqrt 3` and `sqrt 27`, together with the
census size `1,039` (`bridges/b1_origin_list.out`, check C2). The Euler-Sombor index of a graph is the sum over its
edges of `sqrt(i^2 + ij + j^2) = |i - j omega|`, `i, j` the end degrees: the length of an Eisenstein integer.
`7, 13, 19, 27` are the norms at the degree pairs `(1,2), (1,3), (2,3), (3,3)`; `3` is the squarefree part of the
norm `12` at `(2,2)` (bulk weight `sqrt 12 = 2 sqrt 3`), of `27` and of `75` (it is also the norm at `(1,1)`, a pair
that occurs in no graph of the paper); `2` is the degree of the bulk (thread) vertices and a coefficient of `eps`;
`5` is the other coefficient in `eps = 2 sqrt 19 - 5 sqrt 3 = sqrt 76 - sqrt 75`, one of the two gaps
(`0 < eps < delta`) that order the first thirteen minima; `1,039 = 2 + 5 + 13 + 33 + 89 + 240 + 657` is the number
of unicyclic graphs with 4 to 10 vertices used in the verification.

**The same numbers on the 2-clock (EXACT arithmetic; MECHANISM for 3, 5 and 7 only).**
`a^2 + ab + b^2 = (a^3 - b^3)/(a - b)`, so `7 = 2^3 - 1`, `13 = (3^3 - 1)/2`, `19 = 3^3 - 2^3`, `27 = 3^3`. In
base 2:

    3 = 2 + 1,      5 * 13 = 2^6 + 1,      27 * 19 = 2^9 + 1,      7 = 2^3 - 1,
    ord_p(2) = 2, 4, 3, 12, 18, 18, 519   for   p = 3, 5, 7, 13, 19, 27, 1039.

The list splits by the position of `-1` on the 2-clock (the same split is in PRIMES_AND_PAPERS.md). *Fermat type*
(`-1` is a power of 2): `3, 5, 13, 19, 27`. *Mersenne type* (odd order): `7` and `1039`. On this list the split
coincides with the presence of a balanced divisible state. It is not the dichotomy of THM-5128, which is early
against late (and dense against sparse): `3` and `5` are far above the generic value (`3` also has a visible
second-order term), `13` and `19` are below or at it, `27` is a prime power outside THM-5128 (C) (measured excess
`+0.49 / +0.49`, no prediction), `7` is the slowest mixer, and `1039`, of odd order 519, mixes fast (`0.764`) and
shows no bias (a control: predicted `-0.00025 / -0.00031`, measured `-0.0012 / +0.0019` with standard error
`0.0036`; `bridges/b3_f1_data.out`).
- `2^m + 1` is the Hamiltonian path count of the transitive tournament with one reversed arc of span `m + 1`
  (KNOWN: Moon 1972; THM-4748 (5); THM-4861): the one-apex word `h = 2^(m-1)`, a fan. The fans and `h = 0` are the
  only words with `x = y` (KNOWN: THM-5126 (A)). So a prime factor of `2^m + 1` has a balanced p-divisible state:
  `3` at `h = 1`; `5` and `13` together at `h = 32` (`x = y = 65`); `19` and `27` together at `h = 256`
  (`x = y = 513`). These are the first balanced p-divisible states of `3, 13, 19, 27` (`5` has `h = 2` before, and
  `19` the unbalanced `h = 218`); at 3 and 5 the excess is far above the generic value, at 13 and 19 it stays at
  or below it.
- A prime of odd order never divides `2^m + 1`, so it has no balanced divisible state (`7`: first divisible state
  `h = 38`; `1039`: `h = 405,729`). Odd order does not decide the mixing rate. `7` is the slowest modulus found
  (rate `0.932` per digit, shared by all its multiples), followed, among the moduli prime to 7 and up to multiples
  with the same rate, by `93` (`0.914`), `31` and `381` (`0.900`), `15` (`0.892`), `219` and `171` (`0.879`,
  `0.878`), `127` (`0.875`), `51` (`0.869`), `3` (`0.866`); `73` has `0.852`, and the characters prime to `63` have
  `0.889`. On these moduli the powers of 2 cycle quickly (`ord_m(2) = 3, 10, 5, 14, 4, 18, 18, 7, 8, 2`, and `9` at
  73), while `1039`, of order 519, has `0.764`; a small order of 2 is not by itself a criterion: `5`, of order 4,
  is the fastest modulus of the scan (`0.742`; HYP-9660 (B)). KNOWN: the scan for `q <= 120`, the extremality of
  `1/7` and the next-worst list `11/93, 1/31, 1/15, 1/63` (S82 note section 7). New here: the extension to every
  odd `m <= 521` and to `1023, 1039, 1143` (`rates_exact/rates_by_modulus.txt`; 1039 in
  `bridges/b3_f1_data.out`; the value for the characters prime to 63 is in
  `audit/numbers-vs-outputs/check_rates_npz.out`). The second-order term of THM-5128 is visible at `3, 7, 31`.

The doubling map `eta -> 2 eta + b` on `Z/m` under the digit recursion is the same base dynamics as in the Collatz
threads (THM-4553, THM-4629: "characters hear near-returns"). `ord_27(2) = ord_19(2) = 18` with `2^9 = -1` modulo
both is the same statement as `2^9 + 1 = 27 * 19`, and every `2^m + 1` is a fan count: that the repository's
mod-18 / mod-19 pairing appears in the tournament with `H = 513` is a COINCIDENCE (exact arithmetic, no mechanism;
THM-4566 types the 18 / 19 link NUMEROLOGY). The two missing values are `7 = 2^3 - 1` and `21 = (2^6 - 1)/3`, the
repunits `111` in bases 2 and 4 (recorded: seven_twentyone_mersenne_openai_math_20261006.md); `(2^m - 1)/3` for
even `m` are the numbers sent to a power of two by `3x + 1`, `(2^m + 1)/3` for odd `m` by `3x - 1`, and the latter
are one third of a balanced value `2^m + 1`.

**Quadratic forms in the join calculus (EXACT).** A one-arc join is a bilinear form in rooted states,
`H = x_h y_h - x_e y_e + x_R y_R` (signature `(2,1)`); on the diagonal (a one-apex tournament joined to itself with
padding `s`) it is `x^2 + 2^s y^2 - [s = 0]`, the norm form of `Z[sqrt(-2^s)]` at the rooted state: Gaussian for
`s = 0`, `Z[sqrt(-2)]` for `s = 1`. With `x, y` odd these are `1, 3, 5, 1, 1, ... (mod 8)`: a self-join is never
`7 (mod 8)`, and `7` and `1039` are `7 (mod 8)`. `27 = 3^2 + 2 * 3^2` is a self-join (`h = 1`, `s = 1`). The orders
that occur are `Z[sqrt(-2^s)]`, of discriminant `-2^(s+2)`: the orders of conductor `2^floor(s/2)` in the two
fields of discriminant `-4` and `-8`, the 2-power family. For gluings through a core `Z` of interface vertices
(`bridges/b7_selfjoin_eisenstein.py`: every core with `|Z| <= 4` enumerated; the core polynomial confirmed on
explicit tournaments for every core with `|Z| <= 2`, for random cores with `|Z| = 3, 4`, and for one core per
Eisenstein-type form with `|Z| <= 3`; FINITE-EXACT) the literal Euler-Sombor form `x^2 + xy + y^2` never occurs.
With one interface vertex two of the eight cores have the diagonal `x^2 + xy + 2y^2 - x - y - 1`, whose quadratic
part is the norm form of `Q(sqrt(-7))` (the definite forms at `|Z| = 1` have the fields `Q(sqrt(-2))` and
`Q(sqrt(-7))`); a form equivalent to `x^2 + 3y^2` first occurs with two; a form properly equivalent to the
Eisenstein form first occurs with three, `3x^2 + 9xy + 7y^2 = X^2 + XY + Y^2` at `(X, Y) = (x + y, x + 2y)`, at
`|Z| = 3` only for cores whose `Z` is a cyclic triangle (24 cores), and always with lower-order terms. So `7`
enters before `3`, and the lattice `Z[omega]` enters only with a cyclic triangle in the interface (for the cores
enumerated, `|Z| <= 4`: at `|Z| = 4` the only form that is a multiple of a form equivalent to the Eisenstein form
is `5 (x^2 + 5xy + 7y^2)`, on 48 cores, each with a cyclic triangle inside `Z`; added when the corrections were
verified, `audit/verify-corrections/README/v2_cores_eisenstein.out`). Where the polynomial `x^2 + x + 1` of the
cube roots of unity also enters is over `F_2`: in THM-5125 (4b) (nullity of `A^2 + A + I`), and, new, as a factor
of the symbol `(1 + z)^2 (1 + z + z^2)` of the mod-4 form of the one-apex family (THM-5128 (E)), which gives the
period 6. So on the one-apex and F1 families `3` acts through `F_4 = F_2[omega]`, not through the lattice
`Z[omega]`.

**As values and as moduli of F1.**

| n | H-value | one-apex words | `r_F1` | as a modulus: first divisible state, rate |
| --- | --- | --- | --- | --- |
| 3 | yes | h = 1 | 1 | 1, `0.866 = sqrt(3)/2` |
| 5 | yes | h = 2, 3 | 3 | 2, `0.742` |
| 7 | no (THM-1370) | none (`psi != 3`) | 0 | 38, `0.932` (slowest) |
| 13 | yes (order 5) | none (`psi != 6`) | 0 | 32, `0.827` |
| 19 | yes | h = 10, 13 | 4 | 218, `0.805` |
| 27 | yes | h = 17, 23 | 5 | (prime power; rate `0.866`) |
| 1039 | yes (least order 9) | none (`psi != 519`) | 0 | 405,729, `0.764` (58 divisible states below `2^26`, `62.2` expected) |

`1039 + 1 = 2^4 * 5 * 13` lies in the distinguished class of the sheet `s = 0` at `5` and `13` (squarefree-level
product `S_0 = 1.055`), but `1039` is prime, so `S_1 = 0.811` and over both sheets `w_0 S_0 + w_1 S_1 = 0.917 < 1`
(`bridges/b3_f1_data.out`, part F; here `(w_0, w_1) = (0.432, 0.568)` are the shares of the sheets among the
triples, THM-5128 (D1); with the word-pair shares `(0.196, 0.804)` of section 2 the value is `0.859`,
`audit/verify-corrections/README/v1_small_cases.out`): its total local density is below average. F1 misses
`1039`, as it misses 65 of the 512 odd numbers of `[1025, 2049)`: local factors are averages.

**The other three papers** (typed; checks in `bridges/`. The digests of the papers are not committed: the
statements about the papers are CITED from the supplied PDFs, which the audit's citation check re-read).
- *Multiset dimension* (J. Allikvere, arXiv:2607.28813; its results: eight graphs of order 11 with `dim_m = n`,
  refuting the conjecture `dim_m <= n - 1` of Simanjuntak, Siagian and Vetrik; `dim_m = 4` on the `n x n` king
  grid for `n >= 5`; `dim_m(P_3 [x] P_n) = n` for `n >= 6`; a linear two-sided bound at every fixed height).
  Background quoted in its introduction, not a result of the paper: the finite values of `dim_m` are all positive
  integers except `2` (paths have value 1 and every other graph at least 3: Simanjuntak, Siagian, Vetrik,
  arXiv:1711.00225; every `k >= 3` occurs: Khemmani, Isariyapalakul, Int. J. Math. Math. Sci. 2018, Art. 8978193).
  This is a statement of Chan-Pak shape whose exclusion is a one-line symmetry (the two landmarks of a 2-set
  receive the same multiset `{0, d}`), where 7 and 21 need a strong floor and deletion. LOOSE. The strip proof of
  the paper (three-column blindness, a transfer automaton on column words, tight transitions as a regular
  language) has the architecture of the cut-state automata of THM-4878 and of the cells of THM-5128. ANALOGY;
  lost: min-plus optimisation instead of counting. One exact arithmetic fact was added by this session
  (`bridges/b4_strip_automaton.py`, FINITE-EXACT, three independent confirmations of the spectrum; re-derived by
  the audit from the paper's definitions, `audit/citations/check_strip_counts.out`): the recurrent core of that
  automaton (19 states) has characteristic polynomial `x^14 (x - 3)(x - 1)(x + 1)(x^2 + x + 2)`, and the paper's
  counts (`188, 556, 1660, 5004, 44,972` at `n = 5, 6, 7, 8, 10`; the paper states no recurrence) satisfy
  `a(n) = 2a(n-1) + a(n-2) + 6a(n-3)` for `n >= 4` (the roots `1, -1, 0` do not contribute: checked, not implied
  by the polynomial) and `42 a(n) = 32 * 3^n - 9 V_n - 21 U_n`, where
  `((-1 + sqrt(-7))/2)^n = (V_n + U_n sqrt(-7))/2`, so that `V_n^2 + 7 U_n^2 = 2^(n+2)`: the second factor is the
  minimal polynomial of `(-1 + sqrt(-7))/2`, of norm 2, and `|U_n| = 1` exactly at `n = 1, 2, 3, 5, 13`, for every
  `n` by Nagell's theorem on `x^2 + 7 = 2^m` (`m = 3, 4, 5, 7, 15`; CITED). So `Q(sqrt(-7))`, the field of two of
  the one-vertex cores above and of the repository's `{1, 2, 4} mod 7` thread, is also the field of this strip
  count: a COINCIDENCE (exact arithmetic, no mechanism).
- *Edge multiset dimension of hypercubes* (arXiv:2608.09983; already THM-4525, THM-4534). The one-apex count
  depends on its word only through the multiset of gaps between a 1 and a later 0 (THM-4876 (6)): a
  distance-multiset invariant on a path. EXACT restatement, nothing new.
- *Malik's leaf-constrained minimisers*: a `2 x 2` exchange inequality with a classification of equality, as in
  THM-5124, for an additive index. ANALOGY; lost: `H` is not additive over local pieces.
- *Allikvere's correction* itself: the index is the bulk weight times the size plus a finite sum of thread
  penalties (a paraphrase of its Lemma 1; the paper's words are threads, skeleton and penalties), and the error
  came from giving an adjacent-defect link the bulk weight. F1 has the same anatomy (`2^s` for the padding, the
  correction `-[s = 0]` when the two defects touch), and yesterday's over-general sentence about it was repaired
  the same way. ANALOGY; lost: additive versus bilinear.

## 4. Qualifying path counts, sources, and what the domination-defect inequality asks for

Second directive (2026-11-11): "keep going as you were but also draw in new ideas as they become more relevant and
creatively ponder out of the box ways the universal domination-defect inequality allows us to reframe the
constructive task of a chan-pak proof through uniform supply of reusable source interfaces, including their
initial cases". Two other sessions worked the same directive and pushed first: Session262 (Codex; exchange
interfaces, protected supply, selectors, fuelled queues) and kind-pasteur-2026-11-10-S22 (THM-5130 and
SOURCES_AND_DOMINATION.md: domination skeleton plus defect, excess types, sources of the 9-grid at levels 7-11
(one-vertex at 7, 8, 9, 11; a two-vertex selected set at 10), the abundancy index). What follows keeps to what
they do not contain; their results are cited, not re-derived. Canon: THM-5129. Hypothesis: HYP-9661.

### 4.1 Dictionary (interpretation; the identities are EXACT)

The words *qualifying* and *k-vertex extension* were defined first in odd_c_seed_symmetry_20261110/EXTENSIONS.md
(commit 4098631a54, 2026-11-10): qualifying permutations of the `n + k` vertices, with `k` new vertices attached
arbitrarily; an exact compiler `F_{T,k}(x) = H(T_x)` of degree at most `2k`; the aggregate count recovers the
extension mean of THM-5126 (C); refined in session257 (INCOMING_EXTENSION_CONTENT.md). kind-pasteur S22
(kps_s22_sign_flip_towers_20261110/README.md section 4) read the count that qualifies in a one-vertex extension as
the cover count (slope), first. The dictionary below is this session's restriction to the one-vertex interface
with a transitive tail: an interpretation of the owner's words, not the general object.

| the words of the directive | the object |
| --- | --- |
| k-vertex extension of a tournament `B` | `B => TT_(k-1) => z` with `z -> y`: a tail of `k - 1` sinks and an apex over one vertex `y` |
| qualifying path | a directed path `R` of `B` starting at `y`, weighted by `H(B - V(R))`; total weight `U = U_B(y)` |
| the Hamiltonian path count itself | `D + 2 U_even` for the bare apex (`k = 1`), `D + 2^(k-1) U` for `k >= 2`, `D = H(B)` |
| reusable source interface | the pair `(D, U)`: the same base serves every tail length; in `Z_2` the values `D + 2^t U` converge to `D` |
| initial case | the head `k = 1`, which differs from the tail formula by `E = U_odd - U_even`; and the offsets with few slopes |
| the inequality of Session261 in the invoice form of Session262 | `H(T) = H(T - z) + (d(G_z) - 1) + defect`; here `d(G_z) - 1 = 2^(outdeg(y) + t)`, so `U >= 2^outdeg(y)` (kind-pasteur, first) |

A bare apex closes a qualifying path into an odd cycle only if the path has an even number of vertices; one sink
of a tail removes the parity condition (THM-5129 (E)). That is the whole difference between the head and the tails.
The rooted counts themselves have no gaps below a top band: every integer in `[0, 2791]` is a number of
Hamiltonian paths from a fixed start and every integer in `[1, 40162]` is a cover count `U` at order `<= 10`
(helper census, `rooted_spectrum/census10.out`, `q1.out`; FINITE-EXACT, not audited; this extends klein's
observation for fixed-start counts at orders `<= 8`, spectrum_analogies_nearjoin_20261102/README.md). The
exceptions 7 and 21 appear only when `D` and `U` are read together: the difficulty is the coupling, as Session262
puts it ("the base offset is never freely selectable").

### 4.2 What was added (THM-5129)

- **Finite supply (PROVED).** For a fixed offset `D` the slopes, over bases of every order, are finitely many odd
  numbers times powers of two, computable from the strong tournaments with at most `D` paths. KNOWN: the product
  law over strong components (THM-5130 (F2)) and the finiteness of the strong inventory at a fixed path count,
  valid at every order (THM-1370 (1)-(3); Session 24, PROOF section 1, complete for `H <= 73`; THM-4744 (3) for
  every count below 625; THM-4884 (2)). New: the rooted form (the slope set of an offset) and the slope lists.
- **Complete lists at every order for `D <= 73` and for `D = 129, 257` (FINITE-EXACT).** Offsets `3, 5, 9, 17, 33,
  65, 129, 257` have `1, 5, 8, 11, 34, 66, 177, 359` odd parts. For the 9-grid: by `U >= D` a source
  `(2^s + 1, 9 * 2^j)` has `j >= s - 3` (KNOWN: kind-pasteur, SOURCES_AND_DOMINATION.md, Proposition 4), and the
  known sources attain this bound at levels `3, 4, 5, 7, 8` (and `j = 0` at level 2; Session261 at levels 2-6,
  `(129, 144)` and `(257, 288)` at levels 7, 8). The least exponents at levels 2-8 are `0, 0, 1, 2, 5, 4, 5`, and
  the only level where the complete lists say more is 6: least exponent 5 where the bound allows 3 (`(65, 72)`
  and `(65, 144)` are not sources at any order). The family `(2^s + 1, 9 * 2^(s-3))` is present at
  `s = 3, 4, 5, 7, 8, 9` and `11` (11: kind-pasteur), absent at `s = 6` at every order and at `s = 10` through
  order 9: it is not uniform. What the lists for 129 and 257 (THM-5129 (B'); recomputed by the audit, audit-5129,
  after seeing the committed output, so not blind) do give: 177 and 359 odd parts, all with `j <= s`, the absence
  of a dyadic slope, and the window shares below. The strong tournaments with 129 and 257 paths (orders 7, 8, 9
  and 8, 9, 10 only; none of order 11, where the least strong count is 225) are in the inventory of THM-4744 (3);
  they were enumerated again because `U` and `Pi_2` are needed.
- **Permanent omissions (PROVED).** `7, 11, 13, 21, 47, 67` are not `D + 2^t U` (`t >= 1`) at any order: the
  exceptional set of HYP-9453 (1) (klein) is permanent. The four attained ones are heads at a non-sink root.
- **No full window at the offsets computed; lines well supplied in the tested range.** No window `[h, 3h/2]` is
  full at the offsets computed (attained `h <= 73` except `1` and `5`; `129`; `257`); the covered share rises from
  about `0.50` (mean over `37 <= h <= 73`) to `0.75` and `0.72` at `h = 129, 257`; whether windows fill up for
  large `h` is OPEN. Every line `D + 2U = N` with `129 <= N < 16384` carries at least `N/64` sources (rooted order
  `<= 10`), the share of lattice points rises from `0.30` to `0.62`, and in the decile nearest the diagonal
  `U = D` it is `0.94` (odd `U`) and `0.84` (even `U`). Four predictions were registered and pushed before the
  order-10 run; all four held (three of them are one-sided, since the count can only grow with the order, and
  the fourth already held in sample at order 9: THM-5129 (F) 3). These are counts on a finite range (THM-5129 (F)).

### 4.3 Where the local densities meet the sources

F1 is a source family: a one-apex source `(x_1, y_1)` closed by a one-apex module `(x_2, y_2)` in place of the
apex vertex, which multiplies the offset and the slope: offset `x_1 x_2`, slope `y_1 y_2`, head correction 1. Its
two sheets are the head (`s = 0`, value `D + U - 1`) and the tails (`s >= 1`, value `D + 2^s U`). For the F1
source family (nothing is proved here for a general source with the linear readout `D + 2^t U`; HYP-9661 (A) asks
for it) THM-5128 then reads:

- supply is uniform at each fixed odd prime and squarefree odd modulus (and every power of 2), with a power
  saving that is uniform in the region but not in the modulus: each residue class receives its share, and the
  only first-order bias is carried by the triples in which the source or the module is `p`-divisible
  (`p | x_i` and `p | y_i` for one `i`; then `p | D` and `p | U`), where the value is forced into one class;
- that class is `0` on the tails and `-E` on the head. The initial case is visible in the residue statistics;
- the squarefree-level product is at least `8/pi^2` of its mean for every odd `N`, including 7 and 21 (and so
  would be the singular series proper, `(8/pi^2) sigma(n)/n`, which is numerical only).

That no congruence excludes any odd number, 7 and 21 included, was KNOWN: the one-apex family alone attains every
odd residue modulo every modulus (Session 25; Session193 modular_sync Theorem 1; Session206 exact-modular-image).
What THM-5128 adds is the size of the local densities of F1 at fixed odd primes, squarefree odd moduli and powers
of 2 (HYP-9661 (D)). So the residue-class side of "uniform supply" is determined for this family at each fixed odd
prime and squarefree odd modulus (and every power of 2), not uniformly in the modulus: uniformity is HYP-9660 (B),
OPEN, and odd prime-power moduli are numerical only. Beyond that, what a density-one statement for this family lacks
is on the minor-arc (collision) side; completeness is not within reach of F1 alone, which has misses at every
scale tested.

### 4.4 Reframings (typed)

1. **Count lines, not windows (FINITE-EXACT in range; conjecture beyond: HYP-9661 (A)).** The induction
   `N = D + 2U`, `D <= N/3` needs one source among about `N/6` lattice points, and about half of them are sources
   in the tested range. klein has this line over all rooted pairs (HYP-9454 (1) for the sole-exit sweep on the
   lines `3H + 2 sigma`, then HYP-9454 (6) and THM-4865 (e), (f): `h + 2B = v`; digit sets of about 857 members
   for `h` in `[300, 1000)`). What is added is the count `r(N)`, its share and its saturation near the diagonal,
   and the exact slope sets at small offsets (no full window through `h = 73` except at `h = 1, 5`, nor at `129, 257`; whether windows
   fill up for large `h` is OPEN). Residue classes (exploratory; `rooted_spectrum/line_density_classes.out`): on
   `[2^11, 2^13)` the share of lattice points that are sources is `0.578` for `3 | N` and `0.557` otherwise (ratio
   `1.038`; the F1 ratio is `1.375` in the limit); at `5, 7, 11, 13` the ratios are within `0.004` of 1. This is a
   saturating 0/1 statistic at order `<= 10`; whether the multiplicity-weighted supply is equidistributed is not
   tested. The statement to prove is a bivariate local law for `(H(B), U_B(y))`, the rooted version of HYP-9659.
2. **Supply and fibre size (NUMERICAL).** On `23 <= h <= 73` (26 attained offsets) the four offsets with
   exactly two strong classes (`35, 39, 47, 63`; `m(h)` of THM-4884 (2); shares `0.278, 0.350, 0.208, 0.156`) are
   among the five poorest; the third poorest, `67` (`9/34 = 0.265`), has six. Pearson correlation of the covered
   share of `[h, 3h/2]` with the number of strong classes: `0.68` there, `0.66` on `15 <= h <= 73`, `0.41` on all
   attained `3 <= h <= 73` (the offsets `9, 17, 19` also have two classes, with shares `0.60, 0.44, 0.50`);
   `rooted_spectrum/source_supply.out`, section (6). The number of strong classes does not determine the share;
   but by THM-5129 (A2), (A3) a share bounded below at an offset `h` needs `|Delta(h)| >= c h`, and in that sense
   uniform supply at an offset is a statement about multiplicities in the strong spectrum (THM-4810,
   THM-4884 (2), HYP-9453 (3)).
3. **Sources are 2-adic limit points (EXACT).** `D + 2^t U -> D` in `Z_2` (the 2-adic reading of the offset is
   THM-5126 (A); the consequence drawn here is immediate), so every attained number is a limit of other attained
   numbers and the spectrum has no 2-adically isolated point; the grid `(d, r, s)` is the 2-adic polar coordinate
   around 1. "Reusable" means: the limit point stays after any finite use. The head is the one place where the
   2-adic ray does not extrapolate.
4. **A sequence factorisation and an Euler product (EXACT, formal).** Rooted tournaments are sequences of strong
   components, so by THM-5129 (A1), as formal power series in a variable `t` marking the order,

       sum over (B, y) of t^|B| H(B)^(-a) U_B(y)^(-b)  =  (1 - P(a + b))^(-1) * R(a, b) * (1 - Q(a, b))^(-1),

   the sum over isomorphism classes of rooted tournaments, with `P(c) = sum over strong C of t^|C| H(C)^(-c)`,
   `R = sum over strong rooted (C, y) of t^|C| H^(-a) U^(-b)`, `Q = sum over strong C of t^|C| H^(-a) Pi_2^(-b)`
   (the one-vertex tournament included in all three). This is a product of two geometric series over sequences
   of strong components, not an Euler product over primes: strong tournaments are the irreducibles of this free
   monoid (the dictionary strong components <-> primes is kind-pasteur's, SOURCES_AND_DOMINATION.md section 6, in
   which the abundancy index is the local factor; Session 24 has the series at a fixed count). The rational
   primes carry the Euler product of THM-5128, proved for squarefree moduli (the squarefree-level product). The
   first governs how many sources there are, the second how the F1 sources fall into residue classes. No use is
   made of the identity here.
5. **A local-global problem with two exceptions (ANALOGY; HYP-9661 (D)).** Orbit of a digit semigroup, bilinear
   form, and a congruence gap at each fixed modulus: the inputs of the thin-orbit density theorems are present in
   THM-5128 only in this weaker form. Those theorems (Bourgain, Kontorovich: Zaremba, Ann. of Math. 180 (2014)
   137-196; Apollonian curvatures, Invent. Math. 196 (2014) 589-650; CITED) use a gap uniform in the modulus
   (for F1: HYP-9660 (B), OPEN) and a mean multiplicity that is a power of `N` (for F1 it is `0.748 ln N`); the
   minor-arc bound is missing here, and such a bound would give density one at most. The thin-orbit / Zaremba
   analogy for the path spectrum is kind-pasteur S11's (seam_gadget_thickness_spectrum_20261104.md), taken up in
   the S82 note section 7; added here: the comparison of ingredients and the Apollonian counter-example. The
   warning that comes with the analogy is that such problems can fail globally for reasons no congruence sees: the
   local-global conjecture is false for Apollonian packings (Haag, Kertzer, Rickards, Stange, Ann. of Math. 200
   (2024) 749-770: reciprocity obstructions) and for continued fractions over a finite alphabet (Rickards, Stange,
   Duke Math. J. 174 (2025) 3197-3244); CITED, details in HYP-9661 (D).

### 4.5 Salvaged census data (helper runs; FINITE-EXACT, not audited)

`rooted_spectrum/`: all rooted tournaments of orders `<= 9` (and aggregated order 10, 161,870,336 rows); the
rooted census of orders `<= 9` is klein's (THM-4865 (e), (h)), recomputed and matched (`q2.out`). No violation of
`U >= D`, `U = E (mod 2)`, `s + e <= H`, and of nine further inequalities (seven up to converse) listed in
`census10.out`; the composition law over strong components; twelve gadget identities in eight invariants
(`gadgets.out`); least `H` of strong tournaments `3, 5, 9, 15, 25, 45, 75, 125` for orders `3..10` (Busch's
minimum, after Moon: THM-4744; recomputed). `census_k2/`: the one- and two-root-change censuses (HYP-9660 (A)).
`g20_local/`: residue histograms of all `4^n` two-apex words, `n <= 16` (HYP-9660 (D)). `bridges/`: the checks
behind section 3.

## 5. Where this leaves the proof

1. *The local densities of F1 at fixed moduli are determined, and they are not where the difficulty lies.* Local
   densities exist modulo every squarefree odd modulus and every power of 2, multiply over the primes, and are
   positive in every odd class; odd prime-power moduli are numerical only, and nothing is uniform in the modulus
   (HYP-9660 (B), OPEN). They are a thin-interface effect: for a gluing of seam rank `r` with equidistributed
   states the bias is `p^(-r)` (THM-5128 (A)). That no congruence excludes the two missing values was KNOWN
   (Session 25; Session193 modular_sync Theorem 1; Session206 exact-modular-image); THM-5128 adds the size of the
   local densities of F1.
2. *What F1 does at feasible sizes.* Mean multiplicity `0.748 ln N + O(1)` (`+0.67` on the windows below). The
   count is over-dispersed: index of dispersion `Var(r)/mean(r)` about `2.8` on the sheet `s = 0` and `3.1` on the
   sheet `s >= 1`, stable from `2^18` to `2^27`; relative variance `Var/mean^2` from `0.56 / 0.51` on
   `[2^18, 2^21)` down to `0.45 / 0.39` on `[2^24, 2^27)` (`audit/recompute-5128/a7_misses.out`; the constant
   `0.75` of the registered tolerance was a conservative bound, not a measurement). Miss rates `1.38 %`, `0.82 %`,
   `0.62 %`, `0.35 %` on `[2^18, 2^21)`, `[2^21, 2^23)`, `[2^23, 2^24)`, `[2^24, 2^27)` (and `0.087 %` on
   `[2^31, 2^32)`, THM-4875): a decay like a small power of `N`. The misses sit where the class factors are
   smallest: `N = 1 (mod 3)`, the only class without a boost at 3, is missed about three times as often (S78;
   ratio `2.5` to `3.7` on these windows), and integers with `N (N + 1)` free of small primes are the lightest
   (HYP-9657 (B)). The squarefree-level product of THM-5128 (D3) has infimum `8/pi^2 = 0.811` of the mean,
   approached when `N (N + 1)` has no small odd prime factor; at feasible sizes the `1/ln X` term makes the
   contrast larger: on the `N` with `N (N + 1)` coprime to `255255` (8.7 % of the odd `N`) the mean multiplicity
   is `0.72` to `0.73` of the overall mean on the windows from `2^18` to `2^27` (limit `0.821`), and these
   integers carry about a quarter of the misses (`22 %`, `25 %`, `28 %`). The `0.809` of the S82 audit is a
   different statistic, the word-pair-weighted mean on the `N` coprime to `255255` (reproduced: `0.80` to `0.81`),
   whose limit by (D3) is `0.870`: its agreement with `8/pi^2` is a finite-size coincidence
   (`audit/recompute-5128/a13_lightest.out`; `audit/prior-art-wording/a2_coprime_255255_22.out`,
   `a2_coprime_255255_25.out`).
3. *What a density-one statement for F1 still needs is unchanged:* a collision (second-moment) bound, the
   minor-arc side. Such a bound gives at most a density-one statement; completeness needs a pointwise lower
   bound, which no bound of this kind gives (S82 audit M3; THM-4886, NOT CLAIMED) and which F1 alone cannot
   satisfy (it has misses at every scale tested; HYP-9480 (D)). What this session adds to that program is the
   residue-class side for any bilinear join family at a fixed prime (a count of divisible states) and, for one
   state and each fixed odd prime, equidistribution in the regions cut by the form, with a power saving uniform in
   the region but not in the modulus (THM-5128 (B4)). Uniformity in the modulus is HYP-9660 (B), OPEN (for the
   count on full blocks a modulus-uniform bound is THM-4886 (A)); odd prime-power moduli are numerical only.
4. *The extension lane* is section 4: the number of sources on a line is at least `N/64` for every odd `N` in
   `[129, 16384)` at rooted order `<= 10` (THM-5129 (F), a finite-range count), far above the `0.75 ln N`
   representations of F1; linear growth for all `N` is HYP-9661 (A) and, in klein's form, HYP-9454 (1). Explicit
   grammars are a thin part of what exists.

## 6. Principles, named (triggers and counter-indications)

- **Confinement is the only first-order local bias of a bilinear gluing.** Trigger: a residue-class profile of a
  join family looks irregular in `p`. Action: do not fit a formula in `p`; list the states that vanish modulo `p`
  and weigh them (THM-5128 (A)): exactly one class is biased and its excess is that mass. Counter-indications: at
  the slowest moduli (`7`, `31`, `3`) the second-order term is visible at every feasible size; a
  balanced state below `p^2` does not imply a boost (`13`, `19`); modulo a prime power the factor depends on the
  valuation (numerical only).
- **Top digits fix the size, bottom digits fix the ratio, middle digits mix.** For the one-apex state, `y = 2h + 1`
  is read off the top digits, `x/y = 1 - phi(h)/2` off the bottom digits (binary radical inverse), and the residues
  modulo `p` are mixed by the digits in between at the rate of the twisted transfer operator. Trigger: a statement is
  known in a letter-count model and wanted by value. Counter-indication: the by-value limit need not equal the
  letter-count limit; here the shape is the same and the constants are not.
- **On a hyperbola, small states never leave: finite-size effects decay like `1 / ln X`.** A fixed state keeps a
  share of the mass of order `1/ln X`, so the profile at any computable size is dominated by the first few
  divisible states. Trigger: an "excess that decays slowly with N". Action: compute the first-order term from the
  small states before extrapolating.
- **Thin interfaces carry visible local factors.** For a rank-`r` seam with equidistributed states the bias
  is `p^(-r)` (THM-5128 (A)), smaller for larger `r`: a thicker seam still has class factors `1 - p^(-r)` and
  `1 + p^(1-r) - p^(-r)`, only closer to 1. That the grammars `G(k, w)` behave this way is conjectural
  (HYP-9660 (D)). That no congruence explains the two missing values 7 and 21 was KNOWN (section 4.3); the local
  densities of F1 give the sizes.
- **The base prime is exact.** At 2 the digit recursion is unitary up to `sqrt 2` (odd characters): use
  Rudin-Shapiro-type identities there, not mixing estimates.

## 7. Reproduction

    cd 04-computation/experiments/s91_f1_local_densities_20261110
    python3 twisted_spectrum.py 16        > twisted_spectrum.out         # rates, checked against brute force
    python3 twisted_spectrum_2adic.py 18  > twisted_spectrum_2adic.out   # the prime 2: unitary, quadratic form
    python3 mod4_symbol.py                > mod4_symbol.out              # corank of the band form, period 6
    python3 box_identity.py               > box_identity.out             # exact class law (rationals), boxes
    python3 divisible_states.py 24 300    > divisible_states.out         # p-divisible states, limit law of x/y
    python3 constants.py                  > constants.out                # dilogarithm constants, sheet shares
    python3 first_order_exact.py 18 20 22 24 26 > first_order_exact.out  # the two forms of the first-order law
    python3 first_order_exact.py 26 28    > first_order_exact_28.out     # the same at X = 2^28 (slow, memory)
    python3 ap_constants.py               > ap_constants.out             # A_p; reads first_order_exact*.out
    python3 regions_check.py              > regions_check.out            # regions lemma (B4)(ii), random regions
    python3 window_split.py 21 18         > window_split_insample_21_18.out
    python3 window_split.py 23 21         > window_split_insample_23_21.out
    python3 window_split.py 27 24 write-prereg prereg_window_27_24.json   # predictions only
    python3 window_split.py 27 24 prereg_window_27_24.json > window_split_prereg_27_24.out   # about 90 s, 1.5 GB

Every script exits nonzero on a failed check. numba is required for first_order_exact.py, window_split.py,
rates_exact/rates_scan.py and blind_proof/check_blind_proof.py, and for the scripts of `rooted_spectrum/` listed
below except line_density_classes.py and head_values_check.py.

Rates: `rates_exact/rates_by_modulus.txt` (262 moduli: every odd `m <= 521`, and 1023, 1143), from
`rates_exact/rates_scan.py M_LO M_HI OUT.npz` (scan_003_101.npz, scan_103_521.npz with scan_103_521.out,
scan_1023.npz, scan_1143.npz); `m = 1039` is in `bridges/b3_f1_data.out`. The second derivation of THM-5128
(C)-(D): `blind_proof/BLIND_DERIVATION.md`, with `blind_proof/check_blind_proof.py` and its output
`blind_proof/check_blind_proof.out`. The audit's own scripts and outputs: `audit/<auditor>/` next to this file.

Sources (section 4; `rooted_spectrum/`, run from that directory):

    python3 census.py                     > census.out                   # rooted census, orders <= 9 (writes rooted_le9.npz, codes9.npy)
    python3 source_supply.py              > source_supply.out            # THM-5129 (A)-(D); (6) strong classes
    python3 source_supply_levels78.py     > source_supply_levels78.out   # THM-5129 (B'): 129, 257 (needs codes9.npy)
    python3 line_density_classes.py       > line_density_classes.out     # line counts by residue class (exploratory)
    python3 head_parity_check.py          > head_parity_check.out        # THM-5129 (E)
    python3 head_values_check.py          > head_values_check.out        # HYP-9661 (C): heads (needs rooted_le9.npz)
    python3 source_witnesses.py           > source_witnesses.out         # 9-grid sources, levels 2-9, literal recount
    python3 line_density_order10.py 8     > line_density_order9.out      # in sample
    python3 line_density_order10.py 9     > line_density_order10.out     # pre-registered run (needs codes9.npy)

## 8. Next obligations

1. Done: independent audit of THM-5128 and THM-5129, 2026-11-11 (seven auditors; record `audit/AUDIT_REPORT.md`):
   PASS WITH CORRECTIONS, all applied; MISTAKE-761.
2. HYP-9660 (A): the first-order term for two root changes, then a pre-registered window above `2^27`.
3. HYP-9660 (B): the rate scan is complete for the odd moduli up to 521 (and 1023, 1039, 1143): nothing exceeds the
   rate `0.932039` of 7 (`rates_exact/rates_by_modulus.txt`; 1039 in `bridges/b3_f1_data.out`). Extend it, and
   prove a uniform gap: this is the uniformity in the modulus that THM-5128 does not have.
4. HYP-9661 (A): a bivariate local law for `(H(B), U_B(y))`; first test, the saturated shares at rooted order 11
   on `[2^13, 2^15)` by sampling, and the share with `D <= U <= 1.05 D`.
5. The minor-arc side for F1: a Type II bound for the bilinear form using the gap of the twisted operators
   (continues THM-4886 and the collision balance of THM-4890). It would give density one at most, not
   completeness.
6. Exact supply at the offsets `513, 1025, 2049` (levels 9-11): the strong tournaments with at most that many
   paths have order up to 12, 13, 15 (`f(13) = 625`, `f(14) = 1125`, `f(16) = 3125`, by the formula quoted in
   THM-5129 (A2)); the extension method of source_supply_levels78.py applies as long as the lists of strong
   tournaments with few paths stay small.
7. HYP-9661 (C): heads at non-sink roots with `D <= N/3`. Below 382,688 they are the tail values (HYP-9661 (C),
   answer added when the corrections were verified); OPEN beyond.
8. Prime-power moduli (raised by the audit): write the argument of (C), (D2) of THM-5128 for the modulus `p^k`
   (with (B4) for that modulus). The by-value law there is NUMERICAL (9, 25, 27); a proof would make
   `(8/pi^2) sigma(n)/n` the singular series of a sheet.
9. Whether the windows `[h, 3h/2]` fill up for large `h` (OPEN; share `0.50` at `37 <= h <= 73`, `0.75` and `0.72`
   at 129 and 257).
