# Independent audit of S83: one-apex exponential sums (THM-4886, HYP-9479)

**Auditor.** Independent adversarial audit, 2026-11-06 (Claude Opus 5.5 subagent). The session's code and prose were not trusted. Every number below comes from this folder's own scripts.

**Object.** Checkpoint 8acbc16ebc:
- `05-knowledge/results/one_apex_exponential_sums_bilinear_weyl_bound_20261106.md` (the note);
- `01-canon/theorems/THM-4886-one-apex-exponential-sums-bilinear-weyl-bound.md`;
- `05-knowledge/hypotheses/HYP-9479-...md` (pre-registered in 185c2b3efa);
- `04-computation/experiments/one_apex_fourier_20261106/`.

Inherited source: Session 193 `modular_sync/PROOF.md`, sections 2, 5 and 6.

**Not audited.** THM-4887 (reserved in f5154129b8 at 13:15; the draft was uncommitted during this audit), and the new `d_*.py` files.

## Verdict: PASS WITH CORRECTIONS

3 MAJOR issues, 11 MINOR issues, 5 upgrades.

- **(A) Theorem A is correct.**
  - Every step was re-derived.
  - Every inequality in the proof chain holds on exact data: 3,892 cases, 0 violations.
  - It holds against exact sums for every admissible q: 1,186 instances, 0 violations.
- **(A1)-(A3): the statements are true, but the written proofs are defective (MAJOR-2).**
  - (A1): the prescribed split does not exist at the top of its range.
  - (A2): the prescribed L = floor(log2 q) gives only 12.70 > 12.
  - "Analytic constant about 8.5" is false.
  - Repaired recipes with analytic constants 5.37 and 11.30 are given in U3.
  - With the optimal split the true constants are 4.223 and 9.543 (n <= 200, every n).
- **The THM-4886 title states a false bound (MAJOR-1).** The (A2) clause lacks the lower limit q >= 2^((n-6)/3). At theta = 1/7 it fails for every n >= 37.
- **(B) and (C) are correct.**
  - (B) is for dyadic blocks; full ranges need a cylinder remark (m4).
  - (C): the irrational case follows from (A1). The rational case is Session 193 Theorem 3 plus a correct transfer. The limiting Fourier coefficient vanishes off (1/2)Z for every modulus, including powers of 2.
- **(D) reproduced exactly with an independent engine.** HYP-9479 P1-P3 hold at n = 26 and 28. The pre-registration timing is valid (section 7).
- **(F) is correct.** The values are at most 2^n, and the lower density is at least 1/(2C). The route can never give more than 1/4, because C >= 2 is forced (U4).
- **Credit is missing (MAJOR-3):**
  - the classical Weyl-Vinogradov Type II method;
  - Session 205's scope correction of S82 (an ancestor of this checkpoint) and the S82 audit's M3, which make section 5's point first;
  - S82 section 7's exact rational-point rates;
  - Session 16's converse involution.
- **Upgrade U1 (proved here) supersedes Theorem A.**
  - The two-block split of Session 193's eq. (1) gives the classical Type II bound
    |S_n| <= 2^(n-1) (6/q + 7·2^(a-n) + 8(1 + ln q) 2^(-a) + 8q(1 + ln q) 2^(-n))^(1/2).
  - On the minor arcs this saves 4 sqrt(n) 2^(-eps n/2), versus 12·2^(-eps n/4).
  - For large q it gives exponent 1/2 instead of 1/4.
  - It is smaller than Theorem A in 1,099 of the 1,186 exact instances.

## 1. Re-derivation of the proofs

### The concatenation law and the split

**H = g(A).** Checked against an explicit tournament: the chain i -> j iff i < j, the apex beating i iff w_i = 1, and Hamiltonian paths counted by a subset DP.
- All 1,023 words of length 0..9 agree (`a1`). The opposite chain orientation fails 930 times.

**Law (1) and law (2).**
- Law (1) is re-derived: each 0 at position j of v adds 2^(j-1) A(u), plus p(u) when j = 1.
- Law (2) follows by applying (1) twice with A(uv) = 2^L A(u) + A(v) and p(uv) = p(v).
- Both are checked exhaustively for |u|, |w| <= 4, |v| <= 4 (28,830 triples), and on 3,000 random triples of length up to 40. There are 0 failures.

**The split facts.** Exhaustive for a, l, L <= 8:
- with p(u) fixed, A(u) = 2^(a-1) + p + 2j for 0 <= j < 2^(a-2);
- with eps(w) fixed, C(w) = eps + 2i for 0 <= i < 2^(l-1);
- A and C are bijections onto [0, 2^L), with C(v) = 2^L - 1 - A(reverse v);
- each difference d is attained by exactly 2^L - |d| ordered pairs.

### Theorem A, step by step

- **B, Hhat and the four pieces.** The terms B(u, w), Hhat(v) and the decomposition into four pieces S_(p,eps) are correct. The direct sum equals the sum of the four pieces in all 3,892 tested (n, theta, split) cases (`a5`, c0).
- **Cauchy-Schwarz:** |S_(p,eps)|^2 <= K_u K_w Q (c1).
- **Expansion of Q.** The geometric sums satisfy |sum_j e(2 theta d j)| <= min(K, 1/(2||2 theta d||)), since |sin(pi x)| >= 2||x||.
- **The off-diagonal terms.** dA and dC are nonzero off the diagonal.
  - The triangle inequality and the min bound hold (c2, c3).
  - Weighted AM-GM, which is the same as Cauchy-Schwarz over the pairs, gives at most 2^L sqrt(T_u T_w) (c4; max ratio 0.615).
- **Counting lemma.** In a block of q consecutive d, #{d : ||beta d|| <= t} <= 2qt + 3.
  - It is correct. It is in fact true with +2 (U2).
  - Exact integer tests over 660,485 blocks (q <= 120, every r, the extreme perturbations |delta| = 1/q^2 and 1/(2q^2)) give max_j (j - 2q t_j) = 1.9915 (`a4`).
- **Layer cake:** sum_block min(K^2, 1/(4||beta d||^2)) <= 3K^2 + 2qK. This is correct.
- **Block count:** at most 2(D/q + 1) blocks, so T_K <= 2(2^L/q + 1)(3K^2 + 2qK) (c5; max ratio 0.704).
- **Assembly:** 4 K_u K_w 2^(L/2) = 2^(a+l-1+L/2) = 2^(n-1) 2^(-L/2). This is correct.

**Theorem A is PROVED as stated** (n >= 4; a >= 2, L, l >= 1; r/q in lowest terms with |2 theta - r/q| <= 1/q^2).

### (A1)-(A3)

**(A1)** (q from 2 to 2^((n-3)/3), constant 6.5): the statement is true.
- The recipe "L = floor(log2 q), a, l >= log2 q + 2" needs n >= floor(log2 q) + 2 ceil(log2 q) + 4.
  - This fails at the top of the range: at q = 2^((n-3)/3) with an integer exponent, it needs n + 1 letters.
  - It is infeasible at 6,273 of 925,317 grid points (`a9`).
- A balanced split proves the bound with 5.37 (U3).
- The optimal-split supremum is 4.2231 (n <= 200, every n, 1/128 grid in log2 q; `a3`). The session's 4.213 came from n = 0 mod 3 only.

**(A2)** (q from 2^((n-6)/3) to 2^(n-6), constant 12): the statement is true. The optimal-split supremum is 9.5428 (`a3`).
- The written recipe "L = floor(log2 q), a, l about (n - L)/2" gives 12.7008 at n = 104, log2 q = 32.667 (`a9`).
  - This is the bottom of the range for n = 2 mod 3. The session's scan sampled only n = 0 mod 3.
  - The analytic supremum of that recipe is 13.61.
- So the written proof does not reach 12. U3 repairs it with L = floor(log2 q) + 1, giving 11.30.
- The note's "(The analytic constant is about 8.5 up to rounding ...)" is false. `const_check.py` itself prints "claimed 11.3".

**(A3)** follows from (A1) and (A2) as written:
- Dirichlet holds with a real Q;
- reducing r/q keeps q > 2^(eps n);
- for eps n < 6 the claimed bound exceeds 4.24 > 1 and is trivial.

The optimal-split supremum of bound · 2^(eps n/4) is 4.793 (n <= 200, eps in [0.05, 1/3]); the note says 4.6. The restriction eps <= 1/3 is unnecessary: the proof works for eps < 1/2.

It is a genuine power saving |X_n|^(-eps/4). It is nontrivial only once eps n > 4 log2 12 = 14.3.

### (B) readouts

Correct as stated for the dyadic blocks X_n.
- In the split, A(x) = 2^(L+l) A(u) + 2^l A(v) + A(w), and e(x) = [x ends in 0] depends only on w.
- The only new v-dependent term, theta y_R 2^l A(v), is absorbed into Hhat.
- The cross terms carry y_h, so beta = 2 y_h theta.
- The identification with the THM-4885 families checks: (g(v), e(v), 1 + v) with v = A(x); the F1 slices have y = (g(u), e(u), 1 + u).

Scope: see m4. The covering obligation (Session 205, section 3) runs v over full ranges [0, V].

### (C) Fourier uniformity

**Irrational theta.** The convergents q_k of 2 theta tend to infinity. For n >= 12, the bound 2^(n/4) <= 2^((n-3)/3) holds, so (A1) applies. The argument is correct.

**Rational theta = r/s.**
- The transfer x = x' 1^t is checked exhaustively for n <= 14: H(x) = H(x'), and x' is strong. There are 2^(m-2) strong words of length m, so S_n(theta) = sum_(m=2..n) 2^(m-2) E[e(theta H(W_m))] + e(theta).
- Session 193 Theorem 3 gives d_TV -> 0 for every fixed M (proof re-read: first-block conditioning on bank members, bijective residues by eq. (5)).
- The limiting coefficient vanishes off (1/2)Z:
  - for odd s it is a full character sum;
  - for even s the law is uniform on odd residues, so E e(rR/s) = e(r/s)(2/s) sum_(k < s/2) e(rk/(s/2)), which is 0 unless (s/2) | r, i.e. s = 2.
  - This covers moduli with powers of 2.
- Numerically, |S_n|/2^(n-1) -> 0 at 1/4, 1/8, 3/8, 1/16, 5/32, 1/6, 1/12, 5/24, 1/3, 1/7, 1/9 and 1/5 (`a8` (iii)). The slowest is 1/7: 0.0756 at n = 26.
- The exact transfer evaluation (`a11`) reproduces S82's rate of 0.932 per step at 1/7.
- "Equidistributed mod 1 after multiplication by any irrational" follows by Weyl's criterion, since k theta is irrational for k != 0.

### (F) linear energy implies positive density

Correct.
- Cauchy-Schwarz gives D_n >= |X_n|^2/E_n.
- The values are at most 2^n (g(v) <= v + 1, checked for v < 2^18 and proved by induction). "At most 2^n + 1" is true but loose.
- For 2^n <= X < 2^(n+1), the count is at least 2^(n-1)/C >= X/(4C). Against X/2 odd numbers this is lower density at least 1/(2C). Along X = 2^n it is 1/C.
- C = 3.25 gives 0.154.
- "The one-apex values alone have density about 0.5": the audit measures 0.4974 at 2^24, using all word lengths (`a6`). The strong-word floor H >= 2^(m-2) + 1 is verified for m <= 26.
- The Session 16 upper bound ln 2 is correctly cited.

## 2. Reproduction (independent code)

**The engine** (`oa_engine.py`):
- g is computed by the (H, A, p) automaton on bits, not by the session's doubling recursion;
- the phases are 64-bit fixed point (error <= h·2^-65), not float64 theta·h;
- the irrationals are 2^-130-accurate rationals (mpmath);
- the float64 method of the session differs from fixed point by at most 2.8e-9 in the ratios.

| n | phi | sqrt2-1 | pi-3 | e-2 | sqrt3-1 | E_n/\|X_n\| | distinct/\|X_n\| | max multiplicity |
|---|---|---|---|---|---|---|---|---|
| 22 | 1.172 | 0.769 | 0.225 | 0.518 | 1.316 | 3.2607 | 0.3729 | 14 |
| 23 | **2.814** | 0.941 | 1.340 | 1.335 | 0.797 | 3.2580 | 0.3730 | 16 |
| 24 | 1.473 | 1.424 | 0.869 | 0.591 | 1.150 | 3.2574 | 0.3730 | 16 |
| 25 | 0.272 | 0.742 | 0.730 | 1.559 | 1.022 | 3.2555 | 0.3731 | 16 |
| 26 | 1.723 | 0.612 | 0.413 | 0.730 | 1.700 | 3.2526 | 0.3734 | 18 |
| 27 | 1.282 | 1.900 | 1.218 | 0.420 | 1.991 | 3.2502 | 0.3735 | 18 |
| 28 | 1.979 | 1.782 | 1.547 | 0.821 | 1.443 | 3.2473 | 0.3737 | 18 |

The ratio columns are |S_n(theta)|/2^((n-1)/2).

- **n = 26 and 28 match outsample.out digit for digit.** So do all in-sample values for n = 12..24 (`a2`, `a2b`).
- **E_n/|X_n| has been slowly decreasing since n = 22.**
- **P3 reproduced:** 48 cases, Q/diagonal in [0.9971, 1.0686], 0 outside [0.85, 1.2] (`a8` (i)). At the major arc 1/7, Q/diagonal is 2.75-7.46.
- **The analogous two-block second moment** sum_u |F(u)|^2/(K 2^m) at n = 22 off the major arcs (`a8` (ii)):
  - 0.76-1.09 at the balanced split a = m = 11 (the low end is 3/37);
  - 0.73-1.55 at a = 8 and a = 14;
  - 14-415 at 1/7.
- **Theorem A against exact sums at every admissible q.** These are all reduced r/q with |2 theta - r/q| <= 1/q^2, q <= 2^min(n,18), enumerated in exact integer arithmetic. There are 16 theta and n = 12, 14, ..., 26.
  - 1,186 instances: 0 violations of (A) (max exact/bound 0.136) and 0 of U1 (max 0.176) (`a7`).

## 3. MAJOR issues

**MAJOR-1. The THM-4886 title states a false bound.** The title reads "12 · 2^(n-1) (q/2^n)^(1/4) for q <= 2^(n-6)". The lower limit is missing.
- Counterexample: theta = 1/7, so 2 theta = 2/7 and q = 7, exactly on the approximation.
  - The exact |S_n|/2^(n-1) is 3.517e-2 at n = 37, against 12 (7/2^37)^(1/4) = 3.206e-2.
  - The gap grows like (0.932/0.841)^n: 6.1e-6 against 1.8e-11 at n = 160 (`a11`, exact transfer evaluation cross-checked with brute force).
  - At theta = 1/3 it fails from n = 113.
- The body (A2) and the note carry the correct range. The note's outcome "C (q^(-1/2) + (q/2^n)^(1/4))" with C = 12 is fine for every q.
- **Correction.** In the THM-4886 title, replace "12 · 2^(n-1) (q/2^n)^(1/4) for q <= 2^(n-6)" with "12 · 2^(n-1) (q/2^n)^(1/4) for 2^((n-6)/3) <= q <= 2^(n-6)".

**MAJOR-2. The written proofs of (A1) and (A2) do not prove the stated constants, and the stated analytic constant is wrong.** The statements remain true.
- (A1)'s prescribed split is infeasible near q = 2^((n-3)/3).
- (A2)'s prescribed L = floor(log2 q) gives 12.70 > 12. The formula's analytic supremum is 13.61; it is not attained, because frac((n-6)/3) lies in {0, 1/3, 2/3}.
- The note's "about 8.5" contradicts both its own script (11.3) and the truth.
- The supporting scan sampled only n = 0 mod 3 and q = 2^(k/4). It missed the n = 2 mod 3 corner where the recipe fails. That gap is invisible with the optimal split.
- **Correction.** In the note section 1 and in THM-4886 (A1)/(A2), replace the two recipes as follows.
  - "(A1) Take L = floor(log2 q) and a = floor((n - L + 1)/2) or ceil((n - L + 1)/2), whichever makes 2^(3-a) q and 2^(2-l) q closer, with l = n - a - L.
    - With f = log2 q - L, their product is at most 2^(2-f), and their ratio is 1 or 2^(+-1).
    - So |S_n| <= 2^(n-1) q^(-1/2) sup_f 2^(f/2) (1 + 2(2^-f + 1)(9 + (9/sqrt2) 2^(1-f/2) + 2^(2-f))^(1/2))^(1/2) <= 5.37 · 2^(n-1) q^(-1/2).
    - The optimal split gives at most 4.23 for n <= 200."
  - "(A2) Take L = floor(log2 q) + 1, with g = L - log2 q in (0, 1], and the same balanced a.
    - Then s = (2^(3-a) q · 2^(2-l) q)^(1/2) >= 2^((g-1)/2) because q^3 >= 2^(n-6).
    - Also |S_n| / (2^(n-1)(q/2^n)^(1/4)) <= (2^((5-g)/2)(1 + 2(2^g + 1)(9 + (9/sqrt2) s + s^2)^(1/2))/s)^(1/2), which is decreasing in s.
    - So it is at most sup_g (2^(3-g)(1 + 2(2^g + 1)(9 + (9/sqrt2) 2^((g-1)/2) + 2^(g-1))^(1/2)))^(1/2) = 11.30.
    - L = floor(log2 q) does not suffice: it gives 12.70 at n = 104, q = 2^(32.67).
    - The optimal split gives at most 9.55 for n <= 200."
  - Delete "(The analytic constant is about 8.5 up to rounding; the numerical maximum is 9.5.)". Replace "(A3) ... at most 4.6" with "4.8 (n <= 200, every n)".

**MAJOR-3. Credit and prior art.**

*(a) The method is classical.*
- Theorem A is the Weyl-Vinogradov bilinear (Type II) estimate: Cauchy-Schwarz over one variable, geometric sums over a progression in the other.
- The counting lemma and its layer-cake/sorted-distance consequences are the standard lemma bounding sum_x min(N, ||alpha x||^(-1)) when |alpha - a/q| <= q^(-2). See Vaughan, *The Hardy-Littlewood Method* (2nd ed.), ch. 2; Davenport, *Multiplicative Number Theory*, ch. 25.
- What is new is the input: Session 193's exact concatenation law makes H bilinear across blocks.
- Add to the note section 1 and the THM-4886 status: "The argument is the classical Weyl-Vinogradov Type II bilinear estimate (Vaughan ch. 2; Davenport ch. 25); the new input is Session 193's concatenation law." Retype (A) "PROVED (new application of a classical method)".

*(b) Session 205's scope correction is not cited.*
- The source is `05-knowledge/results/session205_apex_scope_audit_20261106/S82_SCOPE_CORRECTION.md`, commit 500713e22b at 11:26. It is an ancestor of 8acbc16ebc.
- Its section 3 already shows that an individual bound |sum_v e(theta(a g(v) + cv))| << V^(1-eta), even uniform in the coefficients, together with a growing bank, does not give the covering inequality. It concludes that "cancellation across bank members, stronger integrated estimates, or a different exact covering mechanism" are needed.
- The S82 audit (`one_apex_join_mass_20261106_audit_S82/AUDIT.md`, M3(i)) made the same point first. A sup bound on the minor arcs gives sup |S| ~ |Y| X^(1-eta), against M_Y <= 2|Y| ln 2. With L^2 control it gives at best density one.
- The note section 5 ("The minor-arc bound (B) controls the cross terms only after multiplication by the families' L^2 norms ...") restates this without credit.
- Add both citations there and in THM-4886 `related`. Phrase (B) as delivering the individual dyadic-block bound that the S82 audit and S205 show to be insufficient on its own.

*(c) S82 section 7 is omitted from the prior-work list.* S82 is the same lane, the previous session.
- It already computes exact twisted transfer matrices at rational theta. The decay per doubling is at most 0.9320 over q <= 120, attained at +-1/7, and is sqrt(3)/2 at 1/3.
- It also cites the reversible-prime literature.
- Add to the "Prior work on the Fourier side" list: "S82 (THM-4885 note section 7): exact rational-point rates by twisted transfer matrices (computational; q <= 120; worst 0.9320 at +-1/7); the reversible-prime sums of Dartyge-Martin-Rivat-Shparlinski-Swaenepoel and Bhowmik-Suzuki have linear phases only."

*(d) "The repository had no proved exponential-sum bound for any of these families" is too strong.*
- Session 193 Theorem 3 proves decay at every fixed rational theta, with the explicit, very slow rate 2(1 - T/2^L)^floor((n-2)/L).
- Replace with "no proved minor-arc (uniform in theta) bound; Session 193 Theorem 3 gives proved decay at each fixed rational theta".

*(e) Session 16's converse involution forces E_n >= 2|X_n| - 2^(n/2) - 1 (U4).*
- Cite it wherever E_n = O(2^n) is discussed: "linear energy" means bounded excess over a forced floor of 2.
- The observed excess is 1.25.

## 4. MINOR issues

- **m1.** `bound_check.py` at n = 24 applies Theorem A to inadmissible q: q = 2889 for phi and q = 3363 for sqrt2-1.
  - Fraction.limit_denominator returns semiconvergents. The values of q^2 |2 theta - r/q| are 1.118 and 1.414.
  - 38 of the 68 rows have bound >= 1 and are vacuous. That leaves 28 valid nontrivial checks.
  - Correction: "0 violations at 68 pairs" should become "0 violations at the 66 admissible (n, theta) pairs, 28 of them nontrivial". Better, cite the audit's 1,186 admissible instances (`a7`), and use convergents in the script.
- **m2.** `const_check.py` scanned only n = 9, 12, ..., 90 and q = 2^(k/4).
  - Its labels disagree: it prints "claimed 11.3"; the note says "8.5"; the theorem says 12.
  - The audit scan covers every n <= 200 on a 1/128 grid; the values are in MAJOR-2.
- **m3.** The counting lemma holds with +2 instead of +3 (U2). The perturbations m delta all have one sign.
- **m4.** (B) and (A) are stated for the dyadic blocks X_n. The covering obligation (S205, section 3) sums v over [0, V].
  - The partial top block splits into at most n cylinders (a fixed prefix with free suffixes). The two-block argument of U1 applies to each cylinder, with u = (prefix, free part).
  - Add one sentence saying so.
- **m5.** "Linear energy is the precise missing lemma for positive density" should read "Linear energy is a sufficient condition (the missing lemma on this route)".
  - Likewise, "It needs a growing bank and the pair correlations ..." should read "A route to density one needs ...".
  - The S82 audit flagged the same "needs is stated for every proof" pattern (its M3(ii)).
  - By U4, C >= 2, so this route gives lower density at most 1/4.
- **m6.** HYP-9479 says "in-sample (n = 14..24) they lie in [0.17, 2.46]". This is for even n only: odd n = 23 gives 2.814 at phi.
  - P1's window [0.05, 20] cannot separate |S_n| ~ 2^(n/2) from 2^(0.55 n) at n = 28.
  - P2's window is +-0.2 around an in-sample plateau, 3.24-3.26 since n = 18.
  - CONFIRMED is accurate, but "the minor-arc sums show square-root cancellation" should read "are consistent with square-root cancellation".
- **m7.** "The vectorised block generator of g agrees with the recursion in 4,600 random checks" has no receipt in the folder. The audit's independent engine reproduces n = 26 and 28 exactly, so the numbers stand.
- **m8.** (A3)'s restriction eps <= 1/3 is unnecessary (eps < 1/2 works). Also state that (A3) is nontrivial only once eps n > 14.3.
- **m9.** (C) irrational: state n >= 12, so that 2^(n/4) <= 2^((n-3)/3).
  - "it is equidistributed modulo every odd modulus" should read "every modulus (on the odd residues for even moduli)".
  - "(new)" should read "(new in the repository; OEIS has no entry for g)".
- **m10.** "All at most 2^n + 1" should read "at most 2^n".
- **m11.** THM-4886 `status` should say "AUDIT: PASS WITH CORRECTIONS (this file)" and carry U1 (below) as the stronger form.

## 5. Upgrades (proved here; checked numerically)

### U1. The two-block Type II bound (stronger than Theorem A)

**Statement.** Let n >= 3, 2 <= a <= n - 1, and r/q in lowest terms with |2 theta - r/q| <= 1/q^2. Then

    |S_n(theta)| <= 2^(n-1) (6/q + 7·2^(a-n) + 8(1 + ln q)·2^(-a) + 8q(1 + ln q)·2^(-n))^(1/2).

With U2 the constants 6 and 7 become 4 and 5.

*Proof.*
1. **Split.** Write x = u y with |u| = a (u_1 = 1) and |y| = m = n - a.
2. **Law (1).** Session 193's law (1) reads H(uy) = H(u) + H(y) - 1 + A(u) C(y) + p(u) eps(y).
3. **Fix the boundary bit.** Fix p = p(u). Then A(u) = 2^(a-1) + p + 2j for 0 <= j < K = 2^(a-2), and
   S_p = sum_u e(theta H(u)) sum_y b_y e(theta A(u) C(y)), with |b_y| = 1.
4. **Cauchy-Schwarz over u.** After expanding the square,
   |S_p|^2 <= K sum_(y,y') |sum_j e(2 theta (C(y) - C(y')) j)| <= K sum_(y,y') min(K, 1/(2||2 theta (C(y) - C(y'))||)).
   Since C is a bijection onto [0, 2^m), this is at most K 2^m [K + sum_(0<|d|<2^m) min(K, 1/(2||2 theta d||))].
5. **Sorted distances.** In a block of q consecutive d, sort the distances t_1 <= ... <= t_q.
   - The counting lemma gives t_j >= (j - 3)/(2q).
   - So the block contributes at most 3K + sum_(j=4..q) q/(j - 3) <= 3K + q(1 + ln q).
   - There are at most 2(2^m/q + 1) blocks.
6. **Assemble.**
   - |S_p|^2 / (K 2^m)^2 <= 2^(-m) + 2(1/q + 2^(-m))(3 + q(1 + ln q)/K).
   - |S_n| <= |S_0| + |S_1| <= 2^(n-1)(...)^(1/2), with K = 2^(a-2) and K 2^m = 2^(n-2).

QED.

**Checks.**
- Every inequality of the chain holds on 420 exact (n, theta, a, p) cases (`a5`, d1-d3).
- Against exact sums: 1,186 admissible instances, 0 violations (`a7`).
- U1 is smaller than Theorem A in 1,099 of them. At n = 26 and the five irrationals the best U1 over admissible q is 0.103-0.166, against 0.289-0.324 for (A). The exact values are at most 3.0e-4.

**Corollaries** (analytic constants; the numerical suprema over n <= 300 are in brackets, `a9`):
- **(U1a)** a = ceil(n/2): |S_n|/2^(n-1) <= (6/q + (17.9 + 8 ln q) 2^(-n/2) + 8q(1 + ln q) 2^(-n))^(1/2).
- **(U1b) Minor arcs.** If 2^(eps n) < q <= 2^((1-eps) n), with 0 < eps <= 1/2 and n >= 7, then |S_n| <= 4 sqrt(n) · 2^(n-1) · 2^(-eps n/2) [2.65 sqrt(n)].
  - The three terms of (U1a) are at most 6·2^(-eps n), (17.9 + 5.55 n) 2^(-eps n) and (8 + 5.55 n) 2^(-eps n). So their sum is at most (31.9 + 11.1 n) 2^(-eps n), and 31.9 + 11.1n <= 16n for n >= 7.
  - Under (A3)'s hypothesis, Dirichlet supplies such a q.
  - This doubles (A3)'s exponent, at the cost of sqrt(n). It is better than 12·2^(-eps n/4) as soon as eps n > 2 log2 n - 6.3.
- **(U1c)** For 2 <= q <= 2^((n-3)/3) and n >= 7: |S_n| <= 3.52 · 2^(n-1) q^(-1/2) [3.37]. The constant tends to sqrt 6 = 2.45, against (A1)'s 6.5.
- **(U1d)** For q >= 2^(n/2) and n >= 7: |S_n| <= 4 · 2^(n-1) (q n/2^n)^(1/2) [2.99]. The exponent is 1/2, against (A2)'s 1/4.
- **(U1e) The best saving the bilinear method gives:** at most 4 sqrt(n) 2^(-n/4) at q = 2^(n/2), by (U1d) [1.96 sqrt(n) 2^(-n/4)]. That is |X_n|^(-1/4) up to sqrt(n), against |X_n|^(-1/6) for (A) at q near 2^(n/3).

**Readouts and (C).**
- (B) carries over verbatim: A(x) = 2^m A(u) + A(y), e(x) depends on y, and beta = 2 y_h theta.
- (C) irrational also follows from (U1a) along the convergents q_k <= 2^(n/2).

**Recommendation.** State U1 as the main theorem, and keep Theorem A as a remark. It would also simplify the note's van der Corput outlook.

### U2. The counting lemma with +2

Write beta d_0 + m beta = (beta d_0 + r m/q) + m delta for 0 <= m < q.
- The grid points run over a translate of (1/q)Z/Z, each once.
- All perturbations m delta lie between 0 and (q-1) delta, so they have one sign.
- A point with ||beta d|| <= t therefore comes from a grid point in an interval of length 2t + (q-1)/q^2.
- Such an interval holds at most q(2t + (q-1)/q^2) + 1 < 2qt + 2 grid points.

Consequences:
- The layer cake becomes 2K^2 + 2qK.
- Theorem A holds with (2 + 2^(3-a) q)(2 + 2^(2-l) q).
- U1 holds with 4/q + 5·2^(a-n).
- Exact tests: max_j (j - 2q t_j) = 1.9915 over 660,485 blocks, and 1.926 over 4,000 random real cases (`a4`).

### U3. Valid recipes for (A1) and (A2)

These are the recipes of MAJOR-2, with analytic constants 5.37 and 11.30. Evaluated (`a9`):
- (A1) balanced: 5.341;
- (A2) with L = floor + 1: 10.734;
- (A2) with the note's L = floor: 12.701;
- the analytic suprema of the three formulas are 5.3696, 11.3019 and 13.6104.

### U4. The forced energy floor

The converse of T(w) is T(w*), with w* the reverse complement. So H(w) = H(w*) (Session 16 `BINARY_PRODUCT_DENSITY.md` (2); Session 191 `one_apex/PROOF.md` section 6.2).

On X_n the map iota(x' 1^t) = (x')* 1^t (x' strong), with iota(1^n) = 1^n, is an H-preserving involution. It has 2^(n/2) fixed points for even n (and at most 2^(n/2) + 1 in general). Hence

    E_n >= 2|X_n| - 2^(n/2) - 1.

This is verified at n = 8, 12, 16, 18 (`a10`). So "linear energy" means bounded excess over 2, and the observed excess is 1.25.

### U5. Calibration: what the minor-arc bounds give for the energy

Apply Parseval with U1b at eps = 1/3:
- the major arcs have measure at most 2·2^(-n/3);
- so E_n <= (16n + 2) 2^(5n/3 - 2);
- so D_n >= 2^(n/3)/(16n + 2).

This is weaker than Session 150's elementary A(X) > sqrt(X)/4. Pointwise minor-arc bounds of this type cannot reach linear energy. The session's THM-4887 draft (not audited here) attacks E_n directly; that is the right direction.

## 6. Prior-art search

### Repository

The search ran on origin/main at 2026-11-06 12:4x, with `git grep` and no truncated listings.
- **Fourier, exponential-sum and bilinear terms:** 14,331 files, mostly old LRC and Collatz material.
- **Co-occurrence with one-apex or spectrum terms** (`--all-match`): 282 files. The relevant ones were read.
- **Fourier-analytic terms inside the Session 190-206 packets:**
  - Session 190: the Walsh/Fourier basis of tournament functions. Unrelated.
  - Session 199 `walsh-source-compiler`: Walsh transforms of signed linear observables. Unrelated.
  - Session 201 `hammersley-source-approximation`: slopes g(u)/u approximate the bit-reversal (Hammersley) grid. This is macroscopic, not unit-scale.

Read in full or in the relevant sections:
- Session 193 `modular_sync` (Theorems 1-4: proved rational-point mixing).
- S82 note section 7 (exact rational rates; HEURISTIC programme), and the S82 audit's M3 (minor-arc sup bounds are insufficient for covering).
- Session 205 `S82_SCOPE_CORRECTION.md` (MAJOR-3(b)).
- THM-4881 (D) and HYP-9474 (ribbon exponential sums: NUMERICAL only).
- kind-pasteur S11 `seam_gadget_thickness_spectrum_20261104.md` (Zaremba and Bourgain-Kontorovich analogy only).
- Session 16 `BINARY_PRODUCT_DENSITY.md` (W(X) = X ln 2 + O(sqrt X); upper density ln 2; the converse involution).
- Session 150 HISTORY (A(X) > sqrt(X)/4, the best all-order lower bound found before the session's THM-4887 draft).
- HYP-9355 (H mod M over random labelled tournaments: a different ensemble, OPEN).
- Session 134 `intervals` (backward-arc counting).
- The S78 note, line 302 (a density-one conjecture for one-apex joins).

No repository file contains a proved minor-arc or uniform-in-theta bound for one-apex sums, an energy computation of E_n, or the bilinear Type II argument. The claim "first proved minor-arc bound" stands for the repository, subject to MAJOR-3.

### Literature

- **Dartyge-Martin-Rivat-Shparlinski-Swaenepoel, "Reversible primes"** (arXiv:2309.11380; J. London Math. Soc. 2024; read in full via pdftotext). Their sums are F_n(alpha, theta) = |B_n|^(-1) sum e(alpha a + theta rev(a)), with linear phases. Their Lemma 2.3 factors these as products of cosines. No product or bilinear phase is treated.
- **Bhowmik-Suzuki** (arXiv:2406.13334 and arXiv:2507.08714; abstracts only): reversed primes in arithmetic progressions and arbitrary bases. These are linear digit problems.
- **OEIS**: no sequence contains 1,1,3,1,5,3,5,1,9,5,9,3,11,5,9,1 (searched 2026-11-06).
- **Web**: no exponential-sum result for n·rev(n), and no resolution of the Chan-Pak question.
- **Drmota-Mauduit** (Weyl sums over integers with digital restrictions): a different structure.
- **Credit owed:** the classical method (MAJOR-3(a)).

## 7. Pre-registration timing (E)

- Commit 185c2b3efa (HYP-9479 with P1-P3, and the in-sample scripts) is dated 2026-11-06 12:20:01 -0700.
- `C:/tmp/s83/outsample.py` was last modified at 12:20:55.62.
- `C:/tmp/s83/outsample.out` was born at 12:20:55.74 and last modified at 12:24:28.58. That is after the commit.
- Both files are byte-identical (after CR stripping) to the copies committed in 8acbc16ebc (12:31:37).
- `git diff 185c2b3efa 8acbc16ebc` on HYP-9479 changes only `status` and "Results". The predictions and windows are untouched.

**The pre-registration is valid.** The caveats on window width are in m6.

## 8. Files (this folder)

**Engines**
- `oa_engine.py`: the independent g, fixed-point exponential sums and the numba energy.
- `oa_bounds.py`: Theorem A (exact minimum over splits) and U1.

**Scripts** (each `aN_*.py` writes `aN_*.out`)

| script | what it checks |
|---|---|
| `a1_identities` | Tournament DP equals g(A); laws (1) and (2) exhaustive; split facts; transfer; oddness; g <= v + 1. |
| `a2_sums_energy`, `a2b_n27_28` | Square-root ratios, energy, distinct values and multiplicities for n = 12..28. |
| `a3_constants` | (A1)-(A3) optimal-split constants, every n <= 200. |
| `a4_lemmas` | Counting lemma (exact) and the first-power lemma. |
| `a5_chain` | Every inequality of both proof chains on exact data. |
| `a6_density` | One-apex value density at 2^k; strong-word floor. |
| `a7_bounds_vs_exact` | (A) and U1 at every admissible q. |
| `a8_p3_and_rational` | P3; two-block second moment; rational decay. |
| `a9_recipes_and_U1` | Written versus valid recipes; analytic suprema; U1 constants. |
| `a10_involution_floor` | U4. |
| `a11_transfer_title_check` | MAJOR-1 counterexample by the exact automaton transfer. |

**Runtime.** n <= 26 runs in minutes. n = 27 and 28 take about 2 GB and 70 s (`NUMBA_NUM_THREADS=2..4` was used because a 2^35 coverage job was running).

**Suggested ledger entry** (for the session to file): MISTAKE (S83).
- A canon title dropped the lower range limit of a corollary. The resulting bound is false at fixed q.
- The corollary proof recipes were not checked over all residues of n mod 3.
- The same-day Session 205 scope correction (an ancestor commit) and the classical Type II method were not credited.
- The simpler two-block split was missed. It gives a strictly stronger bound.
