# AUDIT of the tilted-ensemble packet (THM-5127 / HYP-9659)

**Object audited.** `05-knowledge/results/tilted_ensemble_local_limit_20261110/README.md`
(md5 `eb7928936ec776530ab9a752ca36b02f`, worktree HEAD `84f80f1fe2`), its eleven scripts and `out/`. Also read for
consistency: `01-canon/theorems/THM-5127-*.md`, `05-knowledge/hypotheses/HYP-9659-*.md`, `PREREGISTRATION.md`.
**Date.** 2026-11-10. **Auditor.** Independent adversarial subagent. All programs in this folder were written from
scratch; no packet code was imported or copied (the packet's own scripts were run only to compare their output with
`out/`: `mallows.py`, `mallows2.py`, `random_alternating.py`, `analyze_hist.py`, `verify_cut_mod4.py` reproduce their
files byte for byte, and `ozj.py 30000` was run once, see m12).

## Verdict: PASS WITH CORRECTIONS

- Every stated theorem is mathematically correct as a statement: Theorem A (1)-(3), Theorem B, (4.1)/(4.2), Theorem C
  (1)-(4), Lemma C4, Proposition 5.1, Proposition 6.1 (the last after the hypothesis is worded as in HYP-9659 (C)).
  Every proof was checked line by line and every proof step that can be executed was executed (section 3 below).
- Every exact number of sections 3, 4.1, 5 and 6 was reproduced with independent code, including the complete
  order-10 histogram and the complete two-dimensional census with bases of order `<= 10` (the list of 48
  exceptional `h` with its 185 missing pairs is identical; the (W2) count 57 and its largest exception 363 agree).
  Two numbers differ in the last digit for identifiable reasons (m8, m11).
- Corrections required: **five MAJOR** (two credit errors, one implication that does not follow, one mis-stated
  hypothesis, one statement attributed to a theorem that does not contain it) and **24 minor**.
- The packet does not claim the conjecture, and nothing in it uses a hypothesis as if proved except the items
  flagged in M3 and m16.

---

## 1. MAJOR findings

### M1. Theorem C (1) is not new, and (4.2) is a formula of THM-4118 (credit)

*Where.* README outcome item 3; section 1 ("Closest proved mechanism"); section 4 ((4.2), Theorem C item 1, its
proof, "What `mu_R` is"); THM-5127 status block ("NEW here: ... the parity and strict positivity of the pair term").

*What is wrong.* All three ingredients are prior work in the repository, one of them in a stronger form.

- (4.2) is literally THM-4118 eq. (4) (codex, 2026-08-25): `F_T(S) = H + sum_(i in S) A_i - sum_({i,j} in S) K_ij`
  with `K_ij = Q_ij + Q_ji` (`= mu_R(i,j)`).
- Evenness of `mu`: THM-4097 section 3, eq. (6a) (`F({i,j}) - F({i}) - F({j}) + F(empty) = -2 w_ij`, "the left side
  is even, so `w_ij` is integral", `w = mu/2`) - the same Redei second-difference argument as the README's - and
  THM-4118 eq. (9) ("every generator ... is even").
- Positivity, in the stronger form `mu_R(i,j) >= 2 H(R - {i,j}) >= 2`: Session 136,
  `05-knowledge/results/session136_full_tournament_spectrum_20261102/ear_geometry/positive_cut_geometry.md`, eq. (3)
  (2026-11-02): "`w_ij >= H(T-{i,j}) >= 1` ... pair positivity is strict, regardless of strongness"; restated the same
  day as this packet in THM-5122 ("Capacity": `W_ij >= 2H(C-{i,j}) > 0`). The README cites THM-5122 in section 4.2
  but not for this. (Audit check: `mu >= 2 H(R - {i,j})` holds on all 10,650 (R, pair) cases of orders 2-5, with
  equality in 2,002 of them; `p_proofchecks.out`.)

What is new in Theorem C is (2), (3), (4) and Lemma C4; no earlier occurrence of the `F_2` quadratic-form reading, of
a rank bound for `H mod 4`, or of the comparison with random alternating matrices was found (searches in section 5;
the nearest prior statements on `H mod 4`, THM-466 and THM-478, are not cited either, see m24).
The README's insertion proof of `mu >= 1` is correct (executed on 1,490,330 cases, section 3) and is a second proof
that does not use the odd-cycle formula.

*Replacement text.*

- Outcome item 3: replace "where `mu_R(i,j) >= 2` is the (always even) pair term of the known cut quadratic" by
  "where `mu_R(i,j)` is the pair term of the known cut quadratic (THM-4118 eq. (4)), known to be even (THM-4097
  eq. (6a)) and at least `2 H(R - {i,j}) >= 2` (Session 136 `positive_cut_geometry.md` eq. (3); THM-5122)".
- Section 1, append to "Closest proved mechanism": "THM-4118 states the quadratic in the form (4.2) and the evenness
  of its coefficients; Session 136 (`session136_full_tournament_spectrum_20261102/ear_geometry/positive_cut_geometry.md`)
  proves `mu_R(i,j) >= 2 H(R - {i,j})`. New in Theorem C: (2), (3), (4) and Lemma C4."
- After (4.2): "((4.2) is THM-4118 eq. (4) with `A_i = a_i`, `K_ij = mu_R(i,j)`; in THM-4097's notation `mu = 2w`.)"
- Theorem C item 1: "1. (KNOWN.) `mu_R(i,j)` is even (THM-4097 eq. (6a); THM-4118 eq. (9)) and
  `mu_R(i,j) >= 2 H(R - {i,j}) >= 2` for all `i != j` (Session 136 `positive_cut_geometry.md` eq. (3); THM-5122
  "Capacity")."
- Proof of (1), first words: "(1) is known; the following is a second proof of `mu >= 2` that does not use the
  odd-cycle formula."
- THM-5127 status: delete "the parity and strict positivity of the pair term;" from "NEW here" and add under
  "KNOWN, credited": "the form (4.2) of the cut quadratic and the evenness of its coefficients (THM-4097 eq. (6a),
  THM-4118 eqs. (4), (9)); strict positivity `mu_R(i,j) >= 2 H(R - {i,j})` (Session 136 positive_cut_geometry.md
  eq. (3); THM-5122)".

### M2. Theorem A: the all-`n` algebra is THM-4757 section 5 + THM-4761, not THM-071 (credit and citation)

*Where.* README section 3, first remark after Theorem A; section 5, bullet "small frequencies"; section 7 item 1;
THM-5127 "KNOWN, credited".

*What is wrong.* THM-071 (i) is cited as "the same statement in Walsh language", but that file's status is "PROVED
at n=5,7; CONJECTURED at general n". The statement for every `n` is in canon elsewhere and is not cited:

- THM-4757 section 5 ("Remarks (same session, after audit 2)", PROVED, 2026-10-22):
  `sum_L (n-L)! N_L t^L = (1-t)^(n-1) F(T, (1+t)/(1-t))` and `[r^(n-1-L)] W(r) = 2^(-L) (n-L)! N_L`
  (`F(T,x) = sum_k a_k x^k` is THM-4756's forward-step polynomial);
- THM-4761 (1) (opus S54, 2026-10-22): `kappa_L = N_L/(n)_L` is the degree-`L` Walsh chaos of `kappa = H 2^(n-1)/n!`,
  for every `n`.

At `t = rho` these two give `2^(n-1) (N_rho H)(T) = sum_L (n-L)! N_L rho^L = (1-rho)^(n-1) F(T,(1+rho)/(1-rho))`,
which is Theorem A (1), and the second display is the README's remark about the coefficients of `W_T`. (Verified on
all labelled tournaments of orders 3-5: `a_prior_art_links.out`.) So "what is new is the reading, not the algebra"
is right, but the algebra being read is THM-4757/4761; and the same THM-4761 is the actual Gaussian-regime result
(exact second moment of every chaos, `Var kappa = 2/n - 14/(3 n^3) - ...`, and `sqrt(n)(kappa - 1) => N(0,2)`
conditional on Moran's cited theorem), which section 5 attributes to THM-589 and THM-071.

*Replacement text.*

- First remark after Theorem A: "- In Walsh language (1)-(2) read `N_rho H = sum_j rho^(2j) H^(=2j)`, so the
  coefficient of `r^(n-1-2j)` in `W_T` is `2^(n-1-2j)` times the level-`2j` Walsh part. For every `n` this is already
  in canon: THM-4757 section 5 (`sum_L (n-L)! N_L t^L = (1-t)^(n-1) F(T,(1+t)/(1-t))`,
  `[r^(n-1-L)] W(r) = 2^(-L)(n-L)! N_L`) with THM-4761 (1) (`N_L/(n)_L` is the degree-`L` Walsh chaos of
  `H 2^(n-1)/n!`); put `t = rho`. THM-071 (i) is the case `n = 5, 7` (its general case is marked CONJECTURED in that
  file and follows from THM-4757/4761, or from Theorem A). The gamma-coordinates are the triangular change of
  variable `rho^2 = 1 - 4w`: `gamma_j(T) = (-4)^j sum_{i >= j} C(i,j) H^(=2i)(T)`."
- Section 5 bullet: "- small frequencies at `p = 1/2`: THM-4761 (exact second moment of every Walsh chaos,
  `Var(H)/E[H]^2 = 2/n - 14/(3 n^3) + ...`, and `sqrt(n)(H/E[H] - 1) => N(0, 2)`, CONDITIONAL on Moran's CITED
  central limit theorem for cyclic triangles) and THM-589 (`THM-589-metagraph-H-variance-closed-form-nud-succession`);
  no variance formula or limit law is known for `p != 1/2`;"
- Section 7 item 1, "Source:": "THM-059/062, THM-4757 section 5, THM-4761 (THM-071 for n = 5, 7)."
- THM-5127 "KNOWN, credited", first item: add "for all n: THM-4757 section 5 and THM-4761 (1); THM-071 covers
  n = 5, 7 only".

### M3. HYP-9659 (B) as stated does not imply exponential decay (implication does not follow)

*Where.* README outcome item 5 ("(B) the random-alternating law for `beta_R`, which implies exponential decay of the
mod-4 discrepancy"); HYP-9659 (B) ("within `o(1)` in total variation ... In particular
`E_R[2^(-rank(beta_R)/2)] = O(2^(-m/2))`").

*What is wrong.* For a function with values in `[0, 1]`, closeness `o(1)` in total variation gives only
`E_R[2^(-rank/2)] <= E_random[2^(-rank/2)] + o(1) = o(1)`. The low-rank tail carries weight up to 1 per unit of
probability and is exactly where the census deviates from the model: rank-2 classes are over-represented by factors
7.2, 66 and 145 at orders 8, 9, 10 (2, 8, 7 classes against probabilities `4.0e-5`, `6.3e-7`, `5.0e-9`; with
labelled weights the factors are 5.0, 37, 57). So (B) as stated gives "the mod-4 discrepancy tends to 0", not
`O(2^(-m/2))`.

*Replacement text.*

- README outcome item 5: "(B) the random-alternating law for `beta_R` in the quantitative form
  `E_R[2^(-rank(beta_R)/2)] <= C 2^(-m/2)`, which gives exponential decay of the mod-4 discrepancy by Theorem C (4)
  (agreement of the rank laws in total variation alone gives only `o(1)`)".
- HYP-9659 (B): replace the sentence "In particular ... is `O(2^(-m/2))`." by "**(B') Quantitative form.** There is
  an absolute constant `C` with `E_R[2^(-rank(beta_R)/2)] <= C 2^(-m/2)` (for a uniformly random alternating matrix
  the left side is about `1.62 * 2^(-m/2)` for even `m` and `1.65 * 2^(-m/2)` for odd `m`). By THM-5127 (C) (4),
  (B') gives a mod-4 discrepancy `O(2^(-m/2))` at
  order `m+1`. (B) does not imply (B'): total variation `o(1)` only gives `o(1)`. Evidence for (B') (audit, labelled
  weights `m!/|Aut|`): `E_R[2^(-rank/2)] = 0.10239, 0.07348, 0.05077` at `m = 8, 9, 10` against
  `0.10135, 0.07271, 0.05072`."

(For calibration: the true labelled discrepancy decays faster than either bound, see section 4, "new exact data".)

### M4. Proposition 6.1: the hypothesis as written cannot hold (mis-stated hypothesis)

*Where.* README section 6, Proposition 6.1 ("If (W1) holds for every odd `h >= 237`").

*What is wrong.* (W1) as defined two paragraphs earlier contains "`|T| <= 10`" and "for some `T` with `H(T) = h`".
With `|T| <= 10` it is false for every `h > 15,745` (the largest count at order 10), so the proposition as written
has an unsatisfiable hypothesis; and "for some `T` with `H(T) = h`" presupposes that `h` is a count, which is what
the induction supplies. HYP-9659 (C) has the right wording. The proof itself is correct (arithmetic checked for every
odd `n` in `[1025, 3*10^6]`: `h >= 341`, `n - 3h in {0, 2, 4}`, `n < 5h`).

*Replacement text.* "**Proposition 6.1 (PROVED implication).** Suppose that for every odd `h >= 237` that is a
Hamiltonian path count and every odd `n` with `3h <= n < 5h` there is a tournament `T` (of any order) with
`H(T) = h` and a one-vertex extension with `H(T + z) = n` (HYP-9659 (C); this is (W1) without the restriction
`|T| <= 10`). Then every odd number other than 7 and 21 is a Hamiltonian path count."

### M5. "THM-4876 shows the witnesses must change with `n`" is not in THM-4876

*Where.* README section 6, last paragraph; section 1 ("Corrected near miss ... windows cannot be filled from *one*
realiser"); HYP-9659 (C) "Limits".

*What is wrong.* What is PROVED is that no realiser has an extension value in `(h, (3h+1)/2)` (Session214 forced
shuffle, i.e. HYP-9469 (1); THM-4876 (2) and (7) are the cases "increment 2" and "sole exit"), hence that
`[h, 2h-1]` cannot be filled from one realiser; THM-4876 (5) records as FINITE-EXACT (order `<= 8`) that no realiser
fills `[2h+1, 3h]`. Nothing in THM-4876 concerns `[3h, 5h)`; a single base of order 10 has 1,024 out-sets, so
counting does not exclude it for `h <= 1023`.

The statement is nevertheless true in the census (audit, `d_single_realiser.c`, all 208,767 bases of orders 2-10
with `H < 1024`): the only `h` for which one base fills `[3h, 5h)` is `h = 1`; for `h >= 237` the best single base
covers 107 of the 243 odd numbers (44.0%, `h = 243`, order 10); the best for `h >= 9` is 20 of 23 (`h = 23`).

*Replacement text.*

- Section 6: replace "THM-4876 shows the witnesses must change with `n`;" by "No single base fills a window in the
  census: with bases of order `<= 10` the best realiser of any `h >= 237` covers at most 44% of the odd numbers in
  `[3h, 5h)` (FINITE-EXACT, `audit/d_single_realiser.out`), so the witnesses change with `n`. (PROVED is only that
  `(h, (3h+1)/2)` contains no extension value, Session214 / HYP-9469 (1), with THM-4876 (2), (7) as special cases;
  THM-4876 (5) records FINITE-EXACT, through order 8, that no realiser fills `[2h+1, 3h]`.)"
- Section 1: replace "windows cannot be filled from *one* realiser" by "the windows `[h, 2h-1]` (PROVED) and
  `[2h+1, 3h]` (FINITE-EXACT through order 8) cannot be filled from *one* realiser".
- HYP-9659 (C) "Limits": replace "By THM-4876 the witnesses cannot come from one realiser" by "In the census no single
  realiser covers more than 44% of `[3h, 5h)` for `h >= 237` (FINITE-EXACT, bases of order `<= 10`)".

---

## 2. Minor findings

**m1. Outcome item 2.** "So every real target has a tournament order and a tilt" -> "So every real target `>= 1` has
a tournament order `n >= 3` and a tilt".

**m2. HYP-9659 (A) / Proposition 5.1: `n_0 >= 4` is needed.** `p_n` is undefined at `n = 2` (`E_p[H_2] = 1` for all
`p`) and `p_3 = 1`, where `sd_p(H_3) = 0` and `c/sd_p` is meaningless. Add "`n_0 >= 4`" to the hypothesis ("There are
`n_0 >= 4` and `c > 0` ..."); nothing else changes. With that the proof is correct: `mu_n` is constant for
`n = 1, 2`, strictly increasing from `n = 2`, and for `mu_(n-1) < h <= mu_n` the tilt lies in `[1/2, p_n)`.

**m3. Lemma C4, odd-prime specialisation: missing hypothesis.** "At `theta = a/q` with `q` an odd prime" -> "At
`theta = a/q` with `q` an odd prime not dividing `a`" (for `q | a` the left side is 1 and the claimed bound is
false). The last sentence ("governed by the rank of the pair matrix `mu_R` modulo the denominator") should read
"bounded in terms of the pair matrix `mu_R` reduced modulo the denominator (`mu_R/2` modulo 2 at `theta = 1/4`)".
The lemma and both specialisations hold (480 numerical evaluations of both sides, 0 violations).

**m4. Bridging-tilt remark (NUMERICAL).** (i) With the README's normalisation `A_n(t) = sum_k A(n,k) t^k` (degree
`n-1`) the standard asymptotic is `A_n(t)/n! ~ (1/t)((t-1)/ln t)^(n+1)`; the factor `1/t` is missing (ratio exactly
2.000000 at `t = 1/2` and 1.250000 at `t = 0.8`, for `n = 20` and `60`). The Gaussian exponent `-(4/3) n (p-1/2)^2`
is unaffected. (ii) "`(p_n - 1/2) sqrt(n / ln n)` rises slowly from 0.71 to 0.80" -> "is 0.7355 at `n = 4`, falls
to 0.7105 at `n = 6` and then rises slowly to 0.7954 at `n = 160`". The seven values of `p_n` are right.

**m5. Colliding IDs cited without slug (AGENTS.md rule).** "THM-589" -> `THM-589-metagraph-H-variance-closed-form-nud-succession`
(`E[H^2] = n! W(n)/4^(n-1)`; the other file, `THM-589-cvH-even-run-closed-form-and-2-over-n-rate`, gives the same
numbers and the rate). "THM-1370" -> `THM-1370-h-spectrum-omits-7-21-all-n` (the other THM-1370 is about the
Jacobian conjecture). THM-5127's `related` list already has both slugs; the README has neither.

**m6. Wrong notation attributed to THM-4876 (4).** "`mu = B + C` in THM-4876 (4)'s notation, where the flip
increment is `C - B`": THM-4876 (4) has no `B`, `C` (there `B(.,.)` is the *signed* convolution). Replace by
"`mu_R(i,j) = U(L_i, R_j) + U(L_j, R_i)` in THM-4876 (4)'s notation (`adj(i,j) + adj'(j,i)` in THM-013's), where the
flip increment is `U(L_j, R_i) - U(L_i, R_j)`", and in the next sentence "that `B + C` is even" by "that the two
counts have the same parity". (Both readings of `mu` verified by brute force, orders 2-5.)

**m7. Section 4.2: bit order unstated; "for every `d >= 2`" is provable, not only checked.** The values are right
only if the word is written from the source of the chain to its sink (H = 25/27, 115/117, 487/489, 1999/2001,
8095/8097, 32575/32577, 130687/130689 for `d = 2..8`, `b = 0/1`, actual tournaments on up to 19 vertices); read from
the sink the same words give 11/9, 23/17, 47/33, 95/65, ... . "Apex beats the chain vertices marked 1" alone is
therefore ambiguous (it is the repository convention of Session 133/136 and THM-4876 (6), but the README does not
say so). Replace "(apex beats the chain vertices marked 1)" by "(the word is read along the transitive chain from
its source to its sink, and the apex beats exactly the chain vertices marked 1; THM-4876 (6))". Replace "for every
`d >= 2` (...; checked to `d = 11` with the recurrence of `g`)" by "for every `d >= 2` (PROVED: by THM-4876 (6),
`H = 1 + 2N` with `N = 2^(2d) - 3*2^(d-2) - 1` for `b = 0`, and setting `b = 1` removes `sum_(k=2..d) c(k) = 2^(d-1) - 1`
and adds `c(d+1) = 2^(d-1)`; 25 and 27, 115 and 117, 487 and 489, ...)". The packet contains no script or output
for the `d <= 11` check.

**m8. Section 4.1, comparison table.** (i) Random "lower" at `m = 9` is `0.0013` (exactly `0.0013492`), not
`0.0014` (double rounding of `0.00135`). (ii) The census column is over isomorphism classes, the model is uniform on
labelled matrices; say so. Labelled fractions (audit, weights `m!/|Aut|`): `0.4205 / 0.5505 / 0.0290` at `m = 8` and
`0.8286 / 0.1693 / 0.0021` at `m = 9`; largest deviation 0.80 and 1.08 points (classes: 0.91 and 1.16), so "within
1.2 percentage points" stands. (iii) Section 7 item 3, "over-represented by a factor 1.4-7" -> "the bin
`rank <= m - 4` is over-represented by 1.43 (order 8) and 1.72 (order 9); rank-2 classes by 7.2 and 66".

**m9. Section 5 regime table, second row: false as worded.** "at least `(3/2)^k` words per value at length `k`" ->
"on average at least `(3/2)^k` words per value at length `k`" (S90 section 5 item 3 says "on average"; pointwise
it is false, 29 values have no word).

**m10. "near-Poisson".** The rms normalised second difference 1.36 / 1.45 is a variance of 1.86 / 2.09 times
Poisson, and it grows from order 9 to order 10. In outcome item 4 and in section 5 replace "near-Poisson" by
"of Poisson order (variance about twice Poisson: 1.9 at order 9, 2.1 at order 10)". "No visible singular series" is
fair, but the excess of multiples of 3 (4.9% at order 9, 1.2% at order 10) should be named there as the only
residue effect seen.

**m11. `out/window_le10.out`, band `[3h/2, 2h)`.** The printed `0.7423` counts, for `h = 1 mod 4`, the point
`(3h-1)/2 < 3h/2` (the script starts at `floor(3h/2)`), which is always empty by the forced floor. The density of
odd `n` in `[3h/2, 2h)` is `0.7446`. The README's "74%" is unaffected. The parenthetical "for `h >= 100` the latter
is 50% with bases of order `<= 9` and 74% with order `<= 10`" has no file in `out/`; the audit finds `0.5054` and
`0.7427`.

**m12. `ozj.py`.** README section 9 says "up to 30,000"; the command and `out/ozj_3000.out` are for 3,000. Running
the packet's `ozj.py 30000` gives the same 29 omissions (largest 2,647). Either change the text to 3,000 or ship
the 30,000 output.

**m13. Reproduction gaps.** No command or generator for `out/betarank_labelled_2_7.out` (`betarank.c` reads
gentourng lines only); no merge tool for the eight order-10 parts of `hist.c` or for OR-ing the `pairs.c` bitsets;
`ext_census` (section 8) and the `d <= 11` check (section 4.2) are not in the packet.

**m14. Section 5, "Theorem B supplies the canonical ensemble in which a given target is typical".** Theorem B makes
the target the *mean*; typicality needs concentration under the tilt, which is not proved for `p != 1/2`. Replace
"is typical" by "is the mean".

**m15. Section 5, "are the Gaussian regime at `p = 1/2`".** THM-589 and THM-071 contain no limit law; see M2 for the
replacement (THM-4761, conditional on Moran).

**m16. Section 5, "section 4.1 shows the 2-adic rank is that of a random alternating matrix".** A FINITE-EXACT
comparison worded as a fact. Replace by "section 4.1 shows that at orders 8 and 9 the rank distribution of `beta_R`
agrees with that of a random alternating matrix to within 1.2 points in three bins (HYP-9659 (B))".

**m17. Section 7 item 6, "so a proof of the same shape cannot exist".** Unproved impossibility, and THM-5122 (same
day) gives additive (affine) cubes inside the exchange cube of one core (affine criterion, capacity `k <= n_-(W)`,
a literal two-dimensional witness); THM-5124 is about presence cubes only. Replace by "so the join-additivity used
there is not available on presence cubes (THM-5124) or on the out-set cube of one new vertex (`mu >= 2`); additive
cubes do exist in the exchange cube of a fixed core (THM-5122). This is a typed CONTRAST, not an impossibility
theorem."

**m18. Section 7 item 5.** "There the representation count is critical" is THM-4028
(`THM-4028-sun-two-four-six-eight-average-order-criticality`, cumulative count `~ V X^(25/24)`), not THM-4026/4027;
add the citation. The counterexample is `896,315,812,331,399 = 8.96 * 10^14`.

**m19. Theorem B, credit for the first equality.** `a_k(TT_n) = A(n,k)` is THM-062 ("Transitive tournament"), and
`sum_k A(n,k)(r+1/2)^(n-1-k)(r-1/2)^k` is THM-059's master polynomial `F_(n-1)(r)` ("the W-polynomial of all linear
orderings"). So `E_p[H_n] = F_(n-1)(r)/(2r)^(n-1)`, and for `n = 2k+1`
`sum_j gamma_(n,j) w^j = sum_j b_(k,j) (2j+1)! w^j (1-4w)^(k-j)` with THM-059's central factorial triangle (verified
for `n = 3, 5, ..., 15`). New in Theorem B: the gamma form, its positivity reading and the strict monotonicity. Add
one sentence to the proof: "(the first equality is THM-062's transitive case, i.e. THM-059's `F_(n-1)`)".

**m20. A connection the packet missed (not an error).** THM-4756 (`F(T,x) = H(T) A_n(x) mod 3`) gives
`gamma_j(T) = H(T) gamma_(n,j) (mod 3)` for every `j`, i.e. `N_rho H(T) = H(T) E_p[H_n] (mod 3)` coefficientwise in
`w`: every tournament is `H` copies of the tilted transitive ensemble modulo 3 (verified: all labelled tournaments of
orders 3-5 and samples at 6, 7).

**m21. Section 8, "no realiser of 13 extends to 27".** True for every base of order `<= 10` (audit); type it
FINITE-EXACT and give the range.

**m22. Outcome, last line, and regime table.** "The largest verified range is still `2^35` (THM-4885 (D))": that
file records the `2^35` run as a single run not independently recomputed (`2^30` is). Add "(one run; independently
reproduced to `2^30`)".

**m23. Section 1, "Canonical hostile example".** "density about one half, so no single solvable family can be
complete" does not follow (one family of density one half says nothing about other families). Replace by "density
about one half (THM-4876 (6): 50,006 of the 100,000 odd numbers below 200,000), so this solvable family is far from
complete".

**m24. The closest prior results on `H mod 4` are not cited; they give Theorem C (2) in two lines with an explicit
polar form (credit, not an error).** THM-466 (`THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts`)
and THM-478 (`THM-478-digit-one-degree-law`: digit 1 of `H` is `c_odd mod 2`, by Grinberg-Stanley's
`H = 1 + 2 c_odd (mod 4)`, and has `F_2`-degree at most `2 floor((n-1)/2)` on the arc cube) are the repository's
statements about `H mod 4` as an `F_2` polynomial. An odd cycle through a new vertex `v` uses exactly two arcs at
`v`, so on the star of `v` the degree is 2:
`q_R(X) = sum_{i in X, j notin X} P_ij (mod 2)`, `P_ij` = number of directed paths of `R` from `i` to `j` on an even
number of vertices, hence `beta_R(i,j) = P_ij + P_ji (mod 2)` (verified on all labelled tournaments of orders 2-5
and seeded samples of orders 6-8, every out-set: `b_digit1_link.out`; consistent with Session 136's
`w_ij = R_ij + R_ji` reduced mod 2). Add to section 1: "Closest prior statements on `H mod 4`: THM-466 and THM-478
(digit 1 of `H` is `c_odd mod 2`, of `F_2`-degree at most `2 floor((n-1)/2)` on the arc cube). Theorem C (2) is
their restriction to one vertex star, where the degree is 2; explicitly `beta_R(i,j) = P_ij + P_ji mod 2` with
`P_ij` the number of `i -> j` paths of `R` on an even number of vertices."

---

## 3. Proofs: what was checked and how

| item | result | how |
|---|---|---|
| Theorem A (1) | correct | by hand; brute-force noise average over all `2^C(n,2)` outcomes, exact rationals, all labelled `T` with `n <= 5` (1,098) and 12 random at `n = 6`, `p = 2/3, 1/5, 1/2, 1` |
| Theorem A (2), `abs(r) >= 1/2` | correct | `p = 1/2 + 1/(4r)` lies in `[0,1]` iff `abs(r) >= 1/2`; identity checked at `r = 3/2, -2, 1/2, 7/5`; `W` also computed as a direct product over permutations (`n = 4`) |
| Theorem A (3), `Z`-basis | correct | the forms `(xy)^j (x+y)^(d-2j)` are unitriangular against the monomial symmetric forms; elimination has zero remainder and integer output on every case; `gamma_0 = H`, `sum gamma_j 4^-j = n!/2^(n-1)` exact |
| Walsh remark, `gamma_j = (-4)^j sum C(i,j) H^(=2i)` | correct | exact Walsh transform, all labelled `T`, `n = 3, 4, 5` |
| Theorem B, `gamma_(n,1) = 2^n - 2n`, monotonicity | correct | exhaustive `n = 3..7` as polynomial identities in `p`; `gamma >= 0`, `gamma_1`, endpoints and monotonicity on a grid for `n <= 40`; Foata-Schutzenberger is CITED (not re-proved) |
| (4.1), (4.2) | correct | direct DP on the `(m+1)`-vertex tournament for every out-set: all classes of orders 2-8, all labelled of orders 2-6; one direct control per class at order 9 (191,536) and per base in the window census (208,767), 0 failures |
| Theorem C (1), insertion argument | correct (known, M1) | construction executed literally for every labelled `R` of order 2-6, every pair and every Hamiltonian path `Q` of `R - {i,j}`: 1,490,330 cases, 0 invalid; includes `m = 2` (empty `Q`), `c = r` (6,032 cases) and `c = 0` |
| Theorem C (2) | correct | `q(X) = sum q({i}) + sum beta_ij` checked for every `X`, same universes |
| Theorem C (3), incl. "0 exactly when `q` is non-zero on the radical" | correct | radical computed for every tournament of the same universes: 0 failures |
| Theorem C (4) | correct | `R = T - v` and `X = N^+(v)` are independent and uniform; the identity `P(1) - P(3) = E_R[sign(H(R)) E_X(-1)^q]` reproduced exactly at orders 3-10 by two routes (vertex cubes of order `m` against class weights at order `m+1`) |
| Lemma C4 | correct | Cauchy-Schwarz and expansion redone by hand; 480 numerical evaluations (random `R` on 5-7 vertices, random disjoint `U, V'`, 8 frequencies), 0 violations; `theta = 1/4` gives exactly `2^(-rank of the block)`; odd-prime case holds when `q` does not divide `a` (m3) |
| Proposition 5.1 | correct with `n_0 >= 4` (m2) | by hand |
| Proposition 6.1 | proof correct; hypothesis mis-stated (M4) | arithmetic for all odd `n` in `[1025, 3*10^6]` |

## 4. Recomputation table (packet / audit)

| # | quantity | packet | audit (independent code) |
|---|---|---|---|
| a | gamma-vectors of `E_p[H_n]`, `n = 3..7` | (1,2), (1,8), (1,22,16), (1,52,136), (1,114,720,272) | identical (exhaustive, polynomial identity) |
| a | `E_p[H_n^2]` in `w`, `n = 3..7` | `1+8w`; `1+40w+16w^2`; `1+152w+544w^2+56w^3+48w^4`; `1+552w+6288w^2+6944w^3+3072w^4`; `1+2040w+55760w^2+189848w^3+199464w^4+72640w^5+4128w^6` | identical; value at `w = 1/4` equals `n! W(n)/4^(n-1)`, `W = 8, 32, 158, 928, 6350` |
| a | labelled residues mod 8 (1,3,5,7), orders 5, 6, 7 | .352 .234 .352 .063; .303 .176 .322 .199; .295 .246 .257 .202 | .3516 .2344 .3516 .0625; .3027 .1758 .3223 .1992; .2948 .2462 .2573 .2017 |
| a | `p_n`, `n = 4, 6, 10, 16, 40, 80, 160` | 0.933, 0.888, 0.844, 0.804, 0.731, 0.683, 0.642 | 0.93301, 0.88828, 0.84406, 0.80440, 0.73125, 0.68270, 0.64166 |
| a | `out/lab_3..7.txt` | - | auditor's joint distributions identical |
| b | parity / positivity / Gauss-sum failures, classes of orders 2-9 | 0 / 0 | 0 / 0 / 0 (min `mu` = 2 at every order; max 2, 2, 6, 12, 32, 120, 432, 2016) |
| b | rank table, order 6 | rank 4: 30, rank 6: 26 | 30 (13 zero-bias + 17), 26 |
| b | order 7 | 74, 382 | 74 (32 + 42), 382 (136 + 246) |
| b | order 8 | 2, 205, 3795, 2878 | 2, 205 (160 + 45), 3795 (2532 + 1263), 2878 |
| b | order 9 | 8, 436, 32535, 158557 | 8, 436 (318 + 118), 32535 (26088 + 6447), 158557 (74825 + 83732) |
| b | class mean imbalance, orders 6-9 | 0.1339, 0.0905, 0.0509, 0.0317 | 0.133929, 0.090461, 0.050872, 0.031705 |
| b | labelled mean imbalance, `m = 2..7` | 0.5, 0.5, 0.40625, 0.25, 0.12793, 0.08835 | identical (direct labelled enumeration `m <= 6`; class weights `m <= 7`) |
| b | labelled discrepancy, orders 3-7 | 1/2, 1/2, 13/32, 1/4, 853/8192 | identical (labelled enumeration; vertex cubes; class weights) |
| b | random alternating, `m = 8` / `m = 9` | 0.4205, 0.5585, 0.0210 / 0.8394, 0.1593, 0.0014 | 0.4205, 0.5585, 0.0210 / 0.8394, 0.1593, **0.0013**; closed formula checked against brute force for `m <= 6` |
| b | order 10 ranks (`out/betarank_classes_10_summary.out`, added during the audit; not in the README) | 7, 668, 210562, 5444657, 4077162; no odd `mu`, min 2; P1, P2 confirmed | identical (`b_rankonly.c`); fractions 0.41890 / 0.55940 / 0.02170 lie in the pre-registered intervals. Trivia: the file prints 0.4188985 for 0.41889844, and PREREGISTRATION.md quotes "random: 0.5591" for 0.5590477 |
| c | order 9: classes; solid odd range; first gap; window; least count; rms; residues | 191,536; 2,881; 2,883; [757, 2066]; 108; 1.36; 4.9% | 191,536; 2,881; 2,883; [757.46, 2066.52]; 108 (at `h = 799`); 1.365; 4.9% (`0 mod 3`); histogram identical to `out/hist_9.txt` |
| c | order 10 | 9,733,056; 14,683; [3963, 10192]; 1,504; 1.45; 1.2% | 9,733,056; 14,683 (first gap 14,685); [3963.19, 10192.97]; 1,504 (at `h = 3967`); 1.447; 1.2%; histogram identical to `out/hist_10.txt` |
| d | realised odd `h < 1024`, bases `<= 10` | 510 (all but 7, 21) | 510 |
| d | (W1) exceptions | 48 values, largest 235, 185 pairs listed | identical list |
| d | (W2) exceptions | 57 values, largest 363 | 57, largest 363 |
| d | `n_0(h)/h`, `237 <= h <= 1023` | 1.765 / 2.234 / 2.933 | 1.765 (`h = 345`) / 2.234 / 2.933 (`h = 541`) |
| d | least ratio | 1.5058 | 1.5058 = 387/257, base of order 10; equals the forced value `(3h+3)/(2h)` for `h = 1 mod 4`; no pair below `(3h+1)/2` for any `h < 1024` |
| d | density `[2h,3h)`, `[3h/2,2h)`, `h >= 237` | 0.98927, 0.7423 | 0.98927, **0.7446** (m11) |
| d | out-degree-one density `[3h,5h)`, `[2h,3h)` | 0.7993, 0.8641 | 0.7993, 0.8641 |
| d | `[3h/2,2h)`, `h >= 100`, orders `<= 9` / `<= 10` | 50% / 74% (no file) | 0.5054 / 0.7427 |
| e | section 4.2, `d = 2..8` | `2^(2d+1) - 3*2^(d-1) -+ 1` | confirmed on actual tournaments, source-first reading only (m7) |

**Sub-statement derived by the auditor for bases of order `<= 8`** (7,411 classes, all out-sets): 320 odd `h < 1024`
are realised (every odd `h <= 609` except 7, 21); the window `[3h,5h)` is complete for only 56 of them; for
`h <= 45` the odd `n` in `[3h, 5h)` that are not extension values are exactly the packet's list plus the single pair
`(39, 183)` (supplied by a base of order 9), and the windows of `h = 1, 29, 33, 45` are complete; all 185 pairs the
packet lists as missing are absent at every base order `<= 10`, in particular `(3, 11)`, `(5, 21)`, `(9, 35)`.

**Explicit witnesses** (gentourng line of `T`, out-set of `z` as 0-based vertices; each recounted by an independent
Python DP, `d_witness_verify.out`). Five pairs drawn with `random.Random(20261110)`:

| `h` | `n` | `T` (order 9) | `z` beats |
|---|---|---|---|
| 313 | 1117 | `111111111111011111100111111101111110` | 1 2 4 6 7 8 |
| 381 | 1877 | `111111111111001111110111111101111101` | 2 3 4 6 7 |
| 357 | 1747 | `111111111111011111101111101111111101` | 1 2 4 5 7 |
| 331 | 1319 | `111100011111111111111111111111100101` | 0 4 6 7 8 |
| 283 | 1183 | `111111111111010111101111111111111111` | 1 5 6 |

Further: `(237, 711)`, `(365, 3283)`, `(1023, 5113)`, `(1023, 9205)`, `(235, 837)`, `(39, 183)` at order 9, and the
least-ratio pair `(257, 387)` with `T = 111111101111111111111111111111111111111111110` (order 10), `z` beats `{9}`
(`d_witnesses_extra.txt`).

**New exact data from the audit (not in the packet).**

- Exact labelled mod-4 discrepancy `P(H = 1 mod 4) - P(H = 3 mod 4)`: order 8: `1361/32768 = 0.04153`; order 9:
  `63103/8388608 = 0.007523`; order 10: `650471/268435456 = 0.002423` (each by two routes). The bounds of
  Theorem C (4) are `0.08835`, `0.05028`, `0.03155` (labelled mean imbalance at `m = 7, 8, 9`), so the bound is off
  by factors 2.1, 6.7, 13: most of the decay comes from sign cancellation across `R`, which Theorem C (4) discards.
- Labelled `E_R[2^(-rank(beta_R)/2)]` at `m = 5..10`: `0.28906, 0.19202, 0.14499, 0.10239, 0.07348, 0.05077`.
- Labelled residues mod 8 at order 10: `0.25082, 0.24901, 0.25039, 0.24978`.
- Best single-base coverage of `[3h, 5h)`: see M5.

## 5. Citations and prior-art searches

Every cited theorem and hypothesis file exists (THM-059, 062, 064, 071, 589 (two files), 1370 (two files), 4026,
4027, 4097, 4114, 4875, 4876, 4885, 4887, 4890, 4904, 4960, 5040, 5055, 5120-5125; HYP-9453, 9494, 9657). Statements
attributed to them were compared with the files:

- Correct as cited: THM-059 (definition of `W`), THM-062 (`a_k`, palindromy), THM-064 (`W(i/2)`, i.e. `rho = -i`),
  THM-1370 (`h-spectrum` file: floors `5 < 9`, `15 < 25`), THM-4026/4027, THM-4097 section 2 and THM-4114 section 6
  ((4.1), submodularity), THM-4876 (2), THM-4885 (A), (D), THM-4887, THM-4890, THM-4904 / HYP-9494 (29 omissions to
  `10^9`), THM-4960, THM-5040/5055, THM-5122, THM-5124 (3), THM-4875 (C), HYP-9453 (2), HYP-9657, Session214 (floor
  `(3h+1)/2`), S90 section 5 item 3, META-PATTERNS heading.
- Wrong or incomplete: THM-071 (M2), THM-4876 for the window (M5), THM-4876 (4) notation (m6), THM-589 and THM-1370
  without slug (m5), S90's "on average" dropped (m9), THM-4028 missing (m18).
- Missing credit: THM-4118, THM-4097 (6a), Session 136 `positive_cut_geometry.md`, THM-5122 "Capacity" (M1);
  THM-4757 section 5, THM-4761, THM-4756 (M2, m20); THM-059 `F_(n-1)` (m19); THM-466, THM-478 (m24).

Searches (`grep -ril --include='*.md'` over `01-canon/theorems`, `05-knowledge/hypotheses`, `05-knowledge/results`):

- "noise operator", "noise semigroup", "Ornstein-Uhlenbeck", "heat semigroup": no earlier use for `H` or `W(r)`
  (HYP-2618 uses "noise" for the hard-core measure at activity 2, a different object). The noise *reading* is new;
  the identity is THM-4757/4761 (M2).
- "tilted", "biased (random) tournament", "p-biased", "gamma-vector", "gamma-positiv", "Foata": no earlier tilted
  ensemble or gamma form of the mean. Theorem B's gamma form and monotonicity are new here (ingredients: m19).
- "quadratic form over F_2", "polar form", "alternating form/matrix", "Arf", "H mod 4", "mod-4 law": no earlier
  statement that `(H(R+v_X) - H(R))/2 mod 2` is an `F_2` quadratic form or of a rank bound for `H mod 4`. Nearest:
  THM-466 and THM-478 (digit 1 of `H` as an `F_2` polynomial of degree `2 floor((n-1)/2)` on the whole arc cube; C
  (2) is its restriction to a vertex star, m24); THM-1980 (`H mod 4` as a spectral invariant); THM-4114 section 7
  (rational and odd-characteristic rank of `K = mu`, for another purpose). Theorem C (2)-(4) as stated and Lemma C4
  are new here.
- "random alternating": only Selmer/Sha packets (unrelated).
- "[3h, 5h)", "deletion pair": no earlier all-out-sets window census (Session 191's "strong deletion pairs" is a
  minimum-increment study; HYP-9453 (2) is the near-join form, credited).

## 6. Overclaims and typing (complete list)

M3 (implication), M5 ("shows ... must"), m9 (pointwise for average), m10 ("near-Poisson"), m14 ("typical"), m15
("are the Gaussian regime"), m16 ("shows the 2-adic rank is"), m17 ("cannot exist"), m23 ("so no single solvable
family can be complete"), m7 ("for every `d >= 2`" backed only by a finite check in the packet; provable), m21
(untyped finite fact), m22 (single-run range). Section 8 ("routes tried") is untyped narrative and, apart from
`(13, 27)`, is not supported by anything in the packet. Analogies in
section 7 are all typed (IDENTITY, REFORMULATION, MODEL, ANALOGY, CAUTION, CONTRAST, SHARED PATTERN). No number in
the README disagrees with `out/` except: the text "up to 30,000" against `ozj_3000.out` (m12), "0.0014" against the
exact `0.0013492` (m8), and the two parenthetical percentages that have no file (m11).

## 7. Not checked

- The range `2^35` (THM-4885 (D)) and the `10^9` censuses of THM-4904 / HYP-9494 were not recomputed.
- Foata-Schutzenberger, Dickson, Gilmer-Kopparty, Berkowitz, Moran and Chan-Pak were not consulted (no literature
  lookup was made); they are standard and are used as CITED.
- `out/mallows.out`'s table for the tilted ensemble at `n = 7` (`sd`, `cv`, `P(1 mod 4)` for `p != 1/2`) and the
  values of `p_n` not quoted in the README were not recomputed.
- Section 8's "routes tried" other than `(13, 27)`; the one-apex values for `d = 9, 10, 11` by brute force (closed
  form and proof instead, m7).
- No evidence was examined for HYP-9659 (A) beyond the class histograms.
- The labelled mean imbalance at `m = 10` (needs all 1,024 out-sets for 9.7 million classes) was not computed.

## 8. Addendum: the revised working tree (read as a diff, not re-audited in full)

While this report was being finalised the session applied M1-M5 and m1-m22 to the working tree (README md5
`91afec53f3699673766fc0a43c867eb8`, uncommitted at the time of writing; THM-5127 and HYP-9659 likewise) and added
`gen_labelled.py`, `merge_parts.py`, `run_window.sh`, `flip_family.py`, `out/ozj_30000.out`, the order-10 rank files
and the merged bitsets. The diff was read; the applied texts agree with this report, the corrected
`analyze_window.py` now prints `0.7446`, and `flip_family.out` agrees with `e_oneapex.out` for `d = 2..8`. Still open
in that tree:

- **R1.** m23 and m24 were added to this report after the session read it and are not applied; the header and
  THM-5127 say "22 minor" (now 24).
- **R2. New paragraph "One picture" (section 7), two inaccuracies.** "an ear is a local field with repulsive pair
  couplings (...), which is why presence cubes are never additive (THM-5124) and why no realiser has a full window
  (THM-4876)". (i) The ear cube (out-sets of one new vertex) and the presence cube (vertex selection) are different
  cubes with couplings of opposite sign (THM-4114 section 6 says the two "use different coordinates"); THM-5124's
  rigidity is not a consequence of `mu >= 2`, and it allows additive presence faces of dimension 2 (the
  `c, 3c, 3c, 5c` configuration), forbidding dimension 3. (ii) The second clause repeats M5. Replace by: "an ear is
  a local field with repulsive pair couplings (`+a_i`, `-mu_ij`, THM-4114; by Theorem C (1), known, no coupling
  vanishes), so no two-dimensional face of an ear cube is additive; presence cubes are a different cube, with
  non-negative interactions, and have additive faces only in dimension at most 2 (THM-5124); in the census no
  realiser has a full window `[3h, 5h)` (FINITE-EXACT; PROVED only for `[h, 2h-1]`, THM-4876 (2) and Session214);".
- **R3. `0.5591` -> `0.5590`.** The random rank-8 probability at `m = 10` is `0.5590477`; "0.5591" (double rounding
  of `0.55905`) appears in the new order-10 row of the README table, in HYP-9659 (B) "Evidence", in THM-5127 and in
  `PREREGISTRATION.md` (the last is a committed pre-registration: annotate, do not rewrite). No conclusion changes
  (the pre-registered interval is `[0.550, 0.575]`).
- THM-5127 now mentions an "exact ear-cube Weyl identity" in opus S90's note: it is not in the S90 file of this
  worktree's HEAD but is in `origin/main` as last fetched (S90 note, "The tool at hand"), where the inequality is
  credited to THM-5127 (C4); consistent, no action.

## 9. Rerun

    sh audit/run_all.sh            # everything; 8 min 52 s wall on this machine (4 parallel order-10 processes)
    sh audit/run_all.sh quick      # skips every order-10 pass
    # binaries and tables go to $W (default /tmp/claude-audit-ensemble); set GENTOURNG=/path/to/gentourng if needed

The complete script was run end to end once after all programs were final (exit 0); the `*.out` files in this
folder are from that run (`b_digit1_link.out` from a separate run of the same program).

Programs: `a_thmAB.py`, `a_labjoint.c` + `a_labjoint_analyze.py`, `a_prior_art_links.py` (Theorems A, B);
`p_proofchecks.py` (proof steps, Lemma C4, Propositions, alternating matrices); `b_thmC.c`, `b_classdisc.c`,
`b_rankonly.c`, `b_digit1_link.py` (Theorem C, section 4.1, m24); `d_window.c` + `c_hist_analyze.py` + `d_window_analyze.py` +
`d_spotchecks.py` (sections 5, 6); `d_witness.c` + `d_witness_verify.py`; `d_single_realiser.c` +
`d_single_analyze.py` (M5); `e_oneapex.c` (section 4.2). Outputs: the `*.out` files beside them.
