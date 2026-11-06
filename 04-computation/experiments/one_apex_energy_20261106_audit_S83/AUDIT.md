# Independent audit of THM-4887 (one-apex energy at scale 2^(n/3); X^(2/3) counting bound)

Auditor: independent adversarial audit for opus-2026-11-06-S83, 2026-11-06. Object audited:
`01-canon/theorems/THM-4887-one-apex-energy-at-scale-two-to-the-n-over-three.md` (181 lines, as committed in a2de7874e0;
line numbers below refer to it), plus its section 6 and summary bullets in
`05-knowledge/results/one_apex_exponential_sums_bilinear_weyl_bound_20261106.md` ("note"). All code here was written
from scratch. The session's scripts were read but not imported or run.

## Verdict: PASS WITH CORRECTIONS

The mathematics of (D1)-(D4) is correct. I re-derived every step and every constant, and found no gap:
- the two-block identity;
- same-h separation;
- the strip lemma, including the Dirichlet reduction;
- the small-denominator count and the harmonic sum;
- the absorption 3CBP = 3P^2;
- the (D2) bound 1.39n + 19 and the 0.32 remainder (two printed rows round the wrong way, see m1);
- sharpness;
- the (D3) constants;
- the (D4) Fejer argument and the constant 3.92.

Independent code finds 0 violations everywhere. The lemma was tested with the strongest admissible Dirichlet choice,
including the (D2) splits and the boundary case B = 8.

**Re-derivation of the identity from the recursion** (independent of Session 193). By induction on m:
- m = 1 is the recursion itself.
- Appending a bit b to the low block, l = 2l' + b:
  - b = 1: both sides are unchanged, since g(2l'+1) = g(l'), C_(m+1)(2l'+1) = C_m(l'), and [2l'+1 < 2^m] = [l' < 2^(m-1)].
  - b = 0: g(2V) = g(V) + V + (V mod 2) with V = h 2^m + l', where V mod 2 = l' mod 2. This matches
    g(2l') = g(l') + l' + (l' mod 2) together with C_(m+1)(2l') = 2^m + C_m(l').
- The identity holds for all h >= 0, m >= 1 and 0 <= l < 2^m.

**Edge cases asked about** (none breaks):
- **Small n.** The (D2) splits (12,4), where P = 8C exactly and B = 8, and (13,4), (14,4) pass on all pairs. n = 12 is
  the tightest case: (D1) gives 33.17 against 35.68, and the remainder 0.31995 is just under the 0.32 claimed.
- **h' < h versus h' > h.** The lemma uses only h' >= P and alpha = h/h' in (1/2, 2). Every ordered pair was tested.
- **Parity dependence of gamma_h.** It enters only through [l < 2^(m-1)]. Both the same-h step and the window W only use
  gamma in [0, C] (in fact [0, C - 1]). Both parities were checked exhaustively.
- **c' outside [0, C).** Dropping that constraint only enlarges the count, so the reduction to
  #{c : ||alpha c - delta|| <= t} is valid.
- **Non-unique Dirichlet pair.** The proof works for any admissible (a, b). The audit tests the choice with the largest
  reduced denominator, which is the strongest form.

The corrections are about scope and wording. Three statements in the title or the Statement block claim more than is
proved (M1-M3). There are also rounding, typing, citation and reproducibility fixes, listed under MINOR.

## MAJOR (would change a statement; none is a mathematical error in the proved inequalities)

**M1. "Sharp up to a factor O(n)" is claimed uniformly in Delta, but holds only at Delta of order 2^(n/3).**
Where: title (line 3), (D2) line 52, status line 6, note lines 37 and 294.
- The universal floor is 2^(2n-2)(Delta+1)/(2^n+Delta+1), which is about 2^(n-2)(Delta+1). The ratio of
  (1.39n+19)2^(4n/3) to this floor is about 4(1.39n+19)2^(n/3)/(Delta+1), which is O(n) only when Delta+1 is at least
  c 2^(n/3).
  - At Delta = 2^floor(n/3) the ratio of upper to lower bound is 135, 176 and 230 for n = 12, 13, 14
    (`d_constants.out`), which is fine.
  - At Delta = 0 the bound exceeds the observed E_n = 3.25|X_n| by about 0.6(1.39n + 19)2^(n/3): a factor 390 at n = 12
    and about 1.1e4 at n = 25 (`c_energy.out`, column E_n(0)/2^(4n/3)).
  - It also exceeds, by the same order, every multiset's floor 2^(n-2) and the floor E_n >= 2|X_n| - 2^(n/2) - 1 (the
    THM-4886 audit's U4, from Session 16's involution).
- So it is false that "the pairs ... number at most (1.39n+19)2^(4n/3), which is sharp up to O(n)" for 0 <= Delta <=
  2^floor(n/3).
- **Fix.** Say "sharp up to O(n) at Delta = 2^floor(n/3) (and, by U1, at every dyadic Delta in
  [2^floor(n/3), 2^floor((n-4)/2)]); not sharp up to O(n) for Delta = o(2^(n/3)), in particular not for E_n = E_n(0)".

**M2. (D4) is misread: "on this window, |S_n| has square-root size in mean square, up to O(n)" is false in the standard
sense.**
Where: lines 65-66; note line 299, "square-root cancellation in mean square on that window".
- The inequality I_n := integral over |theta| <= 2^(-t-2) of |S_n|^2 <= 3.92(1.39n+19)2^n is correct.
- But the window has length 2^(-t-1), about 2^(-n/3). The mean square over the window is therefore about 2^(4n/3), so
  the RMS is about 2^(2n/3), not 2^(n/2). This is forced, not an artefact:
  - For |theta| <= 2^(-n-3), every phase theta g(v) lies within pi/4 of 0, so |S_n| >= 2^(n-1)/sqrt 2.
  - That central peak alone gives J_n >= 2^(n-5).
  - Hence the window mean square is at least 2^(n+t-4), and the RMS is at least 2^(2n/3 - 7/3), for every n.
- Exact values (`c_energy.out`, n = 12..22) are I_n = 0.31·2^n and J_n = 0.0618·2^n. So I_n/J_n = 5.0-5.2: the window's
  integral is the trivial peak up to a factor 5.
- **Correct reading.** The window contributes O(n 2^n) to Parseval's integral of |S_n|^2 over [0,1), which equals E_n. It
  has no L^2 mass beyond the central peak, up to O(n). Since I_n >= 2^(n-5), (D4) is itself sharp up to O(n). Say this
  instead.

**M3. "Positive density needs ... / is now reduced to the scales below 2^(n/3)" is unproved, and as stated
unsupported.**
Where:
- title (line 3): "Positive density needs the same equidistribution below the scale 2^(n/3)";
- lines 165-167: "This is exactly the input linear energy E_n = O(2^n) needs";
- note line 40: "Positive density now needs equidistribution of g only below the scale 2^(n/3)";
- note line 331: "Positive density is now reduced to the scales between 1 and 2^(n/3)".

The problems:
- (a) No necessity is proved.
  - Positive density of the spectrum could come from other families: joins, banks, products.
  - For the one-apex family, linear energy is sufficient, not necessary. The note itself says "Linear energy is
    sufficient, not necessary".
- (b) (D1) holds for every admissible m, but THM-4887 draws its consequence only at m = floor(n/3).
  - The claim that the scales above 2^(n/3) are "handled" needs the scale-invariant reading of (D1) (U1 below):
    E_n(Delta) << O(n - 2m)(Delta+1)2^n for all floor(n/3) <= m <= floor((n-4)/2). THM-4887 does not state it. It is a
    two-line corollary.
  - In any case, linear energy E_n(0) = O(2^n) is a statement at scale 0. Larger-scale bounds do not imply it, because
    E_n(0) <= E_n(Delta) runs the wrong way. Logically, nothing has been "reduced".
- **Fix.** Say: "(D1) gives E_n(Delta) << O(n)(Delta+1)2^n for every Delta in [2^floor(n/3), 2^floor((n-4)/2)] (U1).
  Positive density by this route needs that inequality at Delta = 0, which is OPEN. Heuristically, the method uses no
  arithmetic of g below 2^(n/3)."
  - Drop "needs" and "reduced to" from the title.

## MINOR (wording, constants, citations, reproducibility)

**m1. Two (D2) table rows are not upper bounds as printed (lines 134-138).** Exact values at the (D2) split, from
`d_constants.out` (50-digit mpmath):

| r | exact main part | remainder (bound) | printed row | upper bound? |
|---|---|---|---|---|
| 0 | 16.211117 + 1.3862944 n | <= 0.319948 (n = 12) | 16.53 + 1.386 n | no: understates by 0.0011 + 0.00029 n |
| 1 | 17.051628 + 1.1003026 n | <= 0.305519 (n = 13) | 17.36 + 1.101 n | yes |
| 2 | 18.430311 + 0.8733107 n | <= 0.283427 (n = 14) | 18.71 + 0.873 n | no |

- Lines 131-132 say every row includes remainders "at most 0.32 for n >= 12". The rows for r = 1, 2 actually use the
  n = 13 and n = 14 remainders.
- The final 1.39n + 19 holds for all n >= 12, with margin at least 2.51. This was checked analytically, since the slopes
  are below 1.39 and the remainders decrease, and numerically for n <= 3000.
- Fix: print 16.54 + 1.3863 n and 18.72 + 0.8734 n.

**m2. Typing (line 5).** "PROVED (new): the two-block identity (a special case of Session 193's concatenation law)"
contradicts itself.
- The identity is Session 193 `modular_sync/PROOF.md` eq. (1), line 37, with u = the word of h and v = the m-letter word
  of l. Type it INHERITED.
- Same-shape identities already exist:
  - for the tail function Q in Session 206 (`session206/parameter-capacity/audit/REVIEW.md`,
    Q(Mv+r) = Q(v) + Q(r) + v rev_k(r) + topbit_k(r)(v mod 2));
  - the one-block version is Session 16's binary-product identity (`BINARY_PRODUCT_DENSITY.md` eq. (4)).
- 1 <= g(l) <= l + 1 is in Session 158 `binary_bank/NOTE.md` (its accumulator bound 2 <= A(m) <= m + 1) and in the S83
  note, section 5.
- "Its sharpness up to O(n)" is the textbook Cauchy-Schwarz floor for any multiset. Type it as elementary, not "new".

**m3. The comparison in the status (lines 7-11; note lines 320-324) is correct but greatly understated.**

Verified (`e_primes.out`):

| quantity | value |
|---|---|
| pi(2^32) - pi(2^31), segmented sieve and sympy | 98,182,656 |
| C(N, 31) | 6.887e213 > 10^213 |
| (2^32)^31 | 2^992 < 2^1000 |
| X^(2/3) at 2^1000 | 10^200.687 |
| the corollary itself at 2^1000 | 10^196.726 |

- The logic needs only THM-4875 (D), THM-1862 and unique factorisation.
- "Finitely many certified values generate only polylogarithmically many products" is correct: a finite set S in
  [2, infinity) has at most (1 + log2 X)^|S| products up to X. But |S| is about 2·10^8, so "asymptotic" here means
  astronomically large X:
  - Distinct primes from [2^31, 2^32], times an odd certified cofactor below 2^31, already beat
    X^(2/3)/(13 ln X + 121) at every grid point below 2^3456 (`e_primes.out`).
  - Multisets of certified primes taken from all dyadic ranges beat it for every real X < 2^14,751 (`f_smooth.out`).
    This is a rigorous construction, comparing the left endpoint against the right endpoint of each dyadic interval. The
    crossover of this construction is near 2^14,780.
- Fix: replace "below 2^1000" by "for every X < 2^14,751".

**m4. The checks section (lines 173-177) disagrees with the committed script and output.**
- "Max ratio 0.565" is only the W = C maximum, attained at (11, 2). `d_lemma.out` itself reports 0.635 at (13, 3) with
  W = 2C; the note, line 315, correctly says 0.635.
- The "lemma totals" and "actual-shift totals" are W = C quantities that `d_lemma.py` does not print. `h_totals.py`
  reproduces them exactly:

  | (n, m) | lemma total, W = C | actual-shift total, W = C | lemma total, W = 2C | actual-shift total, W = 2C |
  |---|---|---|---|---|
  | (10, 2) | 120,050 | 2,822 | 156,626 | 5,198 |
  | (13, 3) | 2,027,274 | 44,838 | 2,615,946 | 83,322 |

  State "W = C" and print them.
- None of the four tested (n, m) is a (D2) split. In particular the boundary B = 8 used at n = 12, the split (12, 4),
  was untested. This audit closes that gap: see Reproduction.

**m5. Lines 180-181.**
- "E_n/|X_n| rises from 2.93 to 3.26": it peaks at 3.2607 (n = 22), then falls to 3.2526 (n = 26) and 3.2473 (n = 28,
  HYP-9479).
- "Loose by a factor of about 100": the ratio is 53 at n = 12 and 141 at n = 26.

**m6. Line 163 (and note line 327): "a lattice term O(1) per pair".**
- Per pair, the lattice term is 3C/beta + 7, which is up to 3C + 7 when beta = 1.
- On average over pairs it is at most 7 + 3 + (C^2/P)(12 ln B + 31.74) + O(C ln B/P), using the step-4 bound.
- At the (D2) split C^2/P = 2^(1-r), so the average is O(ln B) = O(n), not O(1). This log is exactly the source of the
  O(n) factor.
- The note's "area term C^2/P" should read 4C^2/P. This is immaterial.

**m7. Line 75 and line 79.**
- "l_i the bits of l" must fix the indexing, l = sum l_i 2^i (least significant bit first). Otherwise C_m is not the
  stated reversal.
- Line 79 derives eq. (1) from eq. (2) with w empty. Session 193's eq. (1) states it directly; cite (1).

**m8. Note line 37 says the bound holds "for Delta <= 2^(n/3)".** THM-4887 proves it only for
Delta <= 2^floor(n/3).
- At n = 13, the range (2^4, 2^(13/3)] is not covered: (D1) with m = 5 needs P >= 8C, which fails.
- Align the note with the theorem.

**m9. Citations to add.**
- (a) **Classical lemmas.**
  - The counting lemma of step 3 is the classical spacing lemma behind Weyl's inequality: if |alpha - a'/beta| <=
    beta^(-2), then alpha c mod 1 over beta consecutive c is (1/beta)-spaced up to O(1/beta). See for example the proof
    of Vaughan, *The Hardy-Littlewood Method*, Lemma 2.2.
  - It is also literally THM-4886's step-4 counting lemma, from the same session. Cite both.
  - The count of small denominators is a standard Dirichlet/Farey dissection.
  - The whole argument is an approximate-collision count (lattice points in a thin strip) for the perturbed product
    set {g(h) + hc}. This is in the spirit of the multiplication-table and shifted-product-energy problems; one
    sentence of credit is enough.
- (b) **Upper side.** Session 16 `BINARY_PRODUCT_DENSITY.md` (2) and Session 26 `ASYMPTOTIC_REVIEW.md` show the
  one-apex family has at most (ln 2/2)X + O(sqrt X) values up to X. Add both to `related` next to `CUT_ADDENDUM.md`.
- (c) **(D3).** THM-4886's audited floor E_n >= 2|X_n| - 2^(n/2) - 1 belongs next to "the trivial bound is 2^(2n-2)".
- (d) **Literature.** Chan-Pak, TCS 1015 (2024), section 6.5, Conjecture 6.2, records only the conjecture and Royle's
  verification to 80,557 (MathOverflow 232751). A brief web search found no lower bound on A(X) in the literature.
  Record that.

**m10. Process note.** The ID THM-4887 was first reserved by kind-pasteur-S13 (d36f57c903, 13 s before f5154129b8).
48ba576e1f released it to S83, and the S83 branch already contains that commit. No action needed.

## UPGRADES

**U1. State the scale-invariant consequence of (D1) at general m (recommended; it is what M3's "scales above 2^(n/3)"
needs).**
- Statement. (D1) already holds for every m with P >= 8C. Divided by 2^(n+m), for every m with
  (n-1)/3 <= m <= (n-4)/2 and 0 <= Delta <= 2^m, it reads

      E_n(Delta) <= [17.87 + 6 ln2 (n-2m-1) + 2.5·2^(n-3m) + 2^(-m-1)(4 + 3(n-2m-1)ln 2)] 2^(n+m).

- Consequences:
  - Uniformly over floor(n/3) <= m <= floor((n-4)/2), the bracket is at most 8.3(n - 2m) (n = 12..400). So
    E_n(Delta) <= 17(n - 2m)(Delta + 1)2^n for 2^(m-1) < Delta <= 2^m: linear smoothed energy, up to O(n - 2 log2 Delta),
    at every dyadic scale from 2^(n/3) to 2^(n/2)/4.
  - At the top scale m = floor((n-4)/2) the bracket is an absolute constant:
    - at most 39.9 for every n >= 12, with the maximum at n = 13;
    - at most 34.9 for n >= 20;
    - tending to 30.35 for even n and to 34.51 for odd n.
- Exact values: E_n(2^m)/2^(n+m) = 0.61-0.63 for n = 18..24 and every admissible m, against (D1) values of 30-57
  (`g_scales.out`, `d_constants.out`).

**U2. A slightly sharper perturbation bound.**
- gamma_h(l) <= 2^m - 1, not 2^m. g(l) <= l + 1, with equality exactly at l = 0 and at l = 2^j for j >= 1
  (checked for l < 2^12).
- So W = C - 1 + Delta. This is cosmetic (`a_identity.out`).

**U3. Better constants for the energy itself.**
- With the Delta = 0 constants (8P^2, 33.74), E_n <= (1.39n + 16)2^(4n/3) for every n >= 12, with margin at least 1.01.
- Hence the distinct values number at least 2^(2n/3)/(5.56n + 64), and A(X) >= X^(2/3)/(13 ln X + 102). This is
  marginal.

**U4. Novelty statement.**
- The exponent 2/3 appears to be the first power-law lower bound for A(X) beyond 1/2. This holds both in the repository
  and, after a brief search, in the literature.
- The repository search covered:
  - `git grep` on origin/main for A(X), sqrt(X)/4, X^(2/3), X^(3/4), "distinct ... lower bound";
  - sessions 149-151, 154-158 and 190-206;
  - the Session 16/26 packets, THM-4885 and THM-4886.
  It found only Session 150's sqrt(X)/4 and upper bounds:
  - Session 16/26: at most ln 2/2;
  - Session 158: at most 47/120 for the strong family;
  - Session 197 opposed-binary: O(X^(2/3)), for a different family.
- The status could say "new (no earlier power bound in the repository; none found in the literature)", while keeping
  m3's honest range.

## Reproduction

Folder: `04-computation/experiments/one_apex_energy_20261106_audit_S83/`. Run with `python3 <script>`. Everything
finished with 0 violations.

**`a_identity.py` -> `a_identity.out`: the two-block identity.**
- g(v) = H(one-apex tournament) by brute-force subset DP: all 2,047 words of length <= 10, leading zeros included.
- The two-block identity, re-derived by induction on m (m = 1 base, then appending one bit), and tested:
  - exhaustively for m <= 8 and h < 2^9 (261,120 cases);
  - on 100,000 random big-integer cases with m <= 64 and h < 2^200.
- gamma lies in [0, 2^m - 1].
- 1 <= g(v) <= v + 1 for all v < 2^20.
- Session 193 eq. (1), checked on words by brute force: 1,922 cases.

**`b_strip.py` -> `b_strip_small.out`, `b_strip_big.out`: the strip lemma and the energy decomposition.**
- Worst-case D, for both windows W = C and W = 2C.
- beta = the **largest** admissible reduced Dirichlet denominator, which gives the strongest form of the lemma. The
  min-b rule was also checked.
- The step-4 counts N_beta and the sum of 1/beta.
- Same-h gaps.
- Exact decomposition E_n(Delta) = PC + sum of actual pair counts, with actual <= Z*, and comparison with (D1).
- All ordered pairs for (10,2), (12,4), (13,4), (14,4), (13,3), (15,5), (16,5).
- Random pairs for (18,6) (299,008 pairs) and (21,7) (98,304 pairs).
- Results:
  - 0 violations in every setting;
  - max Z*/bound 0.635, at (13,3) with W = 2C; the sampled maxima are 0.566 and 0.569;
  - max N_beta/bound 0.354;
  - the sum of 1/beta is at most 0.23 of its bound.

**`c_energy.py` -> `c_energy.out`: energies, sharpness and (D4).**
- E_n(0) and E_n(2^floor(n/3)) for n = 12..25:
  - E(Delta)/2^(4n/3) lies in [0.39, 0.67], against (1.39n + 19) = 35.7-53.8;
  - every value is at least the sharpness floor.
- Distinct values are far above 2^(2n/3)/(5.56n + 76).
- (D4) for n = 12..22:
  - the window integral is computed exactly through the autocorrelation of the value histogram: I_n/2^n = 0.31,
    against 140-194;
  - the central peak J_n/2^n = 0.0618;
  - the Fejer step is checked as well.
- Agrees to 4 decimals with the session's `d_energy.out` and `d_lemma.out`.

**`d_constants.py` -> `d_constants.out`: every constant at 50 digits.**
- The (D2) rows and the remainder.
- The harmonic sum (valid, margin 2.11).
- The sharpness floor (minimum 0.15749, at least 1/8).
- (D3): 12.733 ln X + 120.64, at most 13 ln X + 121, minimum ratio 1.023.
- The crossover with sqrt(X)/4, at X = 2^41.8.
- (D4): 3.9168.
- U1 and U3.

**Comparison claims and U1.**
- `e_primes.py` -> `e_primes.out`: the prime count, C(N, 31) and the distinct-prime construction (m3).
- `f_smooth.py` -> `f_smooth.out`: the rigorous certified-prime-product construction, valid for every X < 2^14,751.
- `g_scales.py` -> `g_scales.out`: U1 against exact E_n(2^m).
- `h_totals.py` -> `h_totals.out`: reproduces the unprinted checks-section totals (m4).
