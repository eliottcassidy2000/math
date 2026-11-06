# Representation mass of the one-apex join families: the pi^2/12 law, the constant 48.1976, and what a proof of the {7,21} spectrum would still need

**Session:** opus-2026-11-06-S82. Worktree `C:/tmp/wt-s82`. Canon: THM-4885. Pre-registration: HYP-9478 (commit
09b7a56a98). Independent audit: PASS WITH CORRECTIONS (section 9); corrections applied, MISTAKE-725.

**Owner prompt (verbatim).** "spend another session trying to prove the compelteness of the {7,21} spectrum"

**Target.** Chan-Pak Conjecture 6.2 (TCS 1015 (2024) 114776, arXiv:2308.10214): every odd number except 7 and 21 is the
number H(T) of Hamiltonian paths of a tournament. THM-1370 proves that 7 and 21 never occur. The conjecture is OPEN.
A web search on 2026-11-06 found no published resolution. MathOverflow 232751 records all odd values to 80,555 except
7 and 21.

**Outcome.** No all-order proof. The session proved the first-moment (mass) laws of the known join families, confirmed
them out of sample, and measured the fluctuations a proof through these families would have to control.
- **(A) One-apex readouts (PROVED; extends a known law).**
  - Statement: a g(v) - b e(v) + c(1 + v) <= X on (X/a) ln(1 + a/c) + O(sqrt X) values of v.
  - Credit: the case (a, b, c) = (1, 1, 1) is Session 16's audited W(X) = X ln 2 + O(sqrt X)
    (`feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md` (1), 2026-10-20; Session 26 (4)). The proof is
    the same binary-product coupling.
- **(B) One-apex x one-apex joins (PROVED, new).** There are (pi^2/24) X ln X + kappa X + o(X) of them up to X, with
  2 kappa = -0.4777464 computed by the audit.
- **(C) Mean multiplicities (PROVED, new).**
  - A one-apex x bank family has mean multiplicity sum_y (2/y_h) ln(1 + y_h/y_R) per odd number: 48.1976 for the
    order-<=7 bank and 211.236 for order <= 8.
  - The one-apex x one-apex family has (pi^2/12) ln N.
  - The bank mass is unbounded in the order: M_m >= 2 (ln 2)^2 (m - 1) (audit).
- **(D) Pre-registered tests.**
  - The F2 mass constant holds to four decimals out of sample, and so does the F1 increment of (pi^2/12) ln 2 per
    doubling. These two test the implementation; the laws were derived first.
  - The fixed order-<=7 bank alone covers every odd N <= 2^30 except 7 and 21.
  - Coverage to 2^35 is in section 6.
- **(E) What is missing.**
  - The counts are strongly non-Poisson. F1 misses a density about N^-0.35 of the odd numbers, where an independent
    model with its mean would miss N^-0.82.
  - So a proof through these join families needs a pointwise lower bound for their representation counts. That is a
    Goldbach-type statement (S73, HYP-9457 (c)), not a numeration with free digits.
  - Exponential-sum methods for g would give density-one statements (compare THM-4881 (D) and HYP-9474 for the ribbon
    family). Pointwise completeness needs more.

## 0. Inherited objects

**One-apex tournaments.**
- T(v) is the transitive chain plus an apex that beats the chain vertex i iff bit i of v is 1.
- Its rooted triple at the apex is (g(v), e(v), 1 + v), with g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2) and
  e(v) = [v = 0 or v even].
- Credits: session25; session179 spectrum/REPORT.md (1); session124's f; Session191 one_apex/PROOF.md; THM-4876 (6);
  THM-4875 (D).
- **Closed form** (THM-4875 audit, section 2 (D); Session197 ONE_CORE_BINARY_INVERSE (4)), checked for all v < 2^16:

      g(v) = 1 + #10(v) + sum over the zero bits i of v of floor(v / 2^(i+1)).

- **Weighted inversions.** g - 1 = 2 sum_(i>j) b_i (1 - b_j) w(i - j), with w(1) = 1 and w(l) = 2^(l-2). This is
  Session 16's weighted-inversion formula (`feedback_word_hamilton_spectrum_20261020/proofs.md` section 2) and
  THM-4876 (6).

**One-arc join** (session29 ADDENDUM; THM-4866 (B)). H = x_h y_h - x_e y_e + x_R y_R.

**Families** (THM-4875 (D), audited; recertified to 2^32 by Session 203 with F0 disabled):
- F1 = {(1 + u)(1 + v) + g(u) g(v) - e(u) e(v) : 0 <= u <= v};
- F2 = {y_h g(v) - y_e e(v) + y_R (1 + v) : v >= 0, y in the 2,732 distinct rooted triples of order <= 7}.

**Representation counts.** For odd N, r1(N) counts the pairs u <= v and r2(N) counts the pairs (y, v).

## 1. The attempt: four routes, and where each stalls

Prior work: Sessions 191-204 (navigation S200-S203) and the S73-S81 notes. Relevant obstructions:
- deletion increments: k_min(h) >= h/4 is conjectured (HYP-9469 (1), NUMERICAL); the cases of increment 2 and of
  sole-exit cuts are proved (THM-4876 (2), (7));
- append-only binary classification (Session 202): physical appends give both 2n - 1 and 2n + 1 only for n = 5, 9, 17;
- bounded-defect sparsity (Session 203): sources with Z-defect <= K have O(log^(2K+1) X) distinct values;
- bounded-offset binary lifting (Session 203): from a fixed seed bank it has only O(log M) outputs;
- the one-port digital-context obstruction (Session 203).

The four routes:
1. **H-affine digit operations** (H -> 2H +- 1, 3H + {0, 2, 4}, applied to every member of a class). These need the
   auxiliary data (end counts, root responses) to be affine in H on the class, i.e. bounded defects. That is sparse by
   Session 203. No new route here.
2. **Erdos-type two-register numerations.**
   - Every positive integer is a sum of numbers 2^a 3^b, none dividing another (Erdos; see Erdos-Lewin, Math. Comp. 65
     (1996), and Blecksmith-McCallum-Selfridge, Amer. Math. Monthly 105 (1998)). A family with
     H = 1 + 2 sum 2^(a_i) 3^(b_i) over free staircases would therefore be complete except where its constraints bite.
   - The stationary right-sewing algebra (Session197 RETAINED_ROOT_SEWING (2)) is upper triangular in (h, x, y, d), with
     diagonal characters (H(B), C(B), 0, H(B)), C = labelled two-bin covers.
   - Its readouts are sums over intervals of blocks of prod_outside H x prod_inside C. A block's ratio C/H is not always
     a power of 2 (32/5 for the strong 4-tournament; 36/5, 100/11, 28/3, 64/5 at order 5), so mixed radices occur.
   - But every interval of the word contributes; there is no free choice of terms. No free-staircase family was found.
3. **Intercyclic linearisation (KNOWN and closed).**
   - By the odd-cycle formula H = I(Omega, 2) (Grinberg-Stanley), H(T) >= 1 + 2 oc(T), with equality iff T is
     intercyclic.
   - The intercyclic tournaments are classified. By Bang-Jensen-Bessy (feedback number <= 2) their spectrum is
     S_1 ∪ {3·2^t + 1 : t >= 2} ∪ {9·2^t - 3 : t >= 1} (Session 16 `intercyclic_rays.md` (1), audited ACCEPT).
   - So the linear law reaches at most (ln 2 / 2) X + O(sqrt X) values up to X (BINARY_PRODUCT_DENSITY (3)). Numerically
     it reaches about half the odd numbers (HYP-9469 (3), THM-4878 (5)).
   - The session's own examples (three transitively inserted apexes; two apexes through one hinge vertex) agree with
     that theorem.
4. **Counting.** Every odd N has many representations in F1 u F2 (section 5). A proof by counting needs the mean
   (sections 2-4) and pointwise control of the fluctuations (section 7). This is the route the session pursued.

## 2. (A) The mass of one one-apex readout (PROVED; extends Session 16/26)

**Theorem A.** Let a, c >= 1 and b >= 0 be integers, and Phi(v) = a g(v) - b e(v) + c (1 + v). Then

    C(X) = #{v >= 0 : Phi(v) <= X} = (X / a) ln(1 + a / c) + O(sqrt(X / c) + (a / c) log X + b / c + 1).

The case (a, b, c) = (1, 1, 1) is Session 16's law W(X) = X ln 2 + O(sqrt X). Indeed Phi(w) = g(2w), and the reduced
one-apex words are the even v >= 2. The proof below is Session 16's binary-product coupling
(BINARY_PRODUCT_DENSITY (4)-(8); Session 193 modular_sync Theorem 4) applied to a general readout.

*Proof.*
1. **Linearisation.** For v >= 1 with top bit n, put rho(v) = sum over the zero bits i < n of 2^-(i+1).
   - By the closed form, delta(v) = g(v) - v rho(v) = 1 + #10(v) - (sum of the fractional parts of v / 2^(i+1) over
     those zero bits).
   - Each maximal run of zero bits is topped by exactly one "10", and its fractional parts sum to less than 1. Hence
     1 <= delta(v) <= 1 + #10(v) <= 1 + ceil(n/2) (audit U6).
   - Therefore |Phi(v) - v (c + a rho(v))| <= Delta = a (log2 X + 2) + b + c.
2. **Two independent coordinates.** For v in [2^n, 2^(n+1)), write v = 2^n (1 + t).
   - rho(v) depends only on the n low bits. Bit i has weight 2^-(i+1), so the lowest floor(n/2) bits fix rho to within
     2^-floor(n/2).
   - The remaining ceil(n/2) bits fix t to within 2^-ceil(n/2).
   - The two bit sets are disjoint. So as v runs over the block, the pair (top part, bottom part) runs over the full
     2^ceil(n/2) x 2^floor(n/2) grid, one v per cell, and (t, rho(v)) lies in its cell's rectangle. The point set
     {(t, rho)} is a Hammersley set (van der Corput radical inverse of the complemented low bits).
3. **Counting a block.** The region {(t, s) in [0,1)^2 : 2^n (1 + t)(c + a s) <= Y} is bounded by a decreasing curve,
   which meets O(2^(n/2)) cells. So the block contributes 2^n Area_n(Y) + O(2^(n/2)) values v with
   v (c + a rho(v)) <= Y.
4. **Summing the blocks** up to n = log2(Y / c) (beyond that v c > Y): the error is O(sqrt(Y / c)), and
   sum_n 2^n Area_n(Y) = integral over s in [0, 1] of (Y / (c + a s) - 1)_+ ds = (Y / a) ln(1 + a / c) - 1 for
   Y >= a + c.
   - This uses that w = 2^n (1 + t) runs exactly once over [1, infinity) as n and t vary, with measure 2^n dt = dw.
5. **Sandwich.** {Phi <= X} lies between the sets for Y = X - Delta and Y = X + Delta. Their counts differ by at most
   2 Delta / c, plus the grid error. QED.

In every application b <= a, so the b/c term is harmless: e(u) <= 1 <= g(u), and y_e <= y_h for every bank triple.

## 3. (B) One-apex x one-apex joins: the pi^2/12 law (PROVED, new)

**Theorem B.** #{0 <= u <= v : F1(u, v) <= X} = (pi^2 / 24) X ln X + kappa X + o(X), with
kappa = kappa_1 - c' = -0.2388732 (audit U1).

*Proof.*
1. For fixed u, F1(u, v) is Phi(v) with (a, b, c) = (g(u), e(u), 1 + u). Only u with (1 + u)^2 <= X contribute.
2. Since g(u) <= 1 + u, Theorem A applies with uniform error O(sqrt(X/(1 + u)) + log X). Summed over u <= sqrt X, the
   error is O(X^(3/4)). So

       sum_u C_u(X) = X sum_(u <= sqrt X) phi(t_u) / (1 + u) + O(X^(3/4)),
       phi(t) = ln(1 + t) / t,   t_u = g(u) / (1 + u) = rho(u) + O(log u / u).

3. **Each dyadic block of u.** Over u in [2^m, 2^(m+1)), rho(u) runs exactly once through j / 2^m. The weight 1/(1 + u)
   depends on the top bits and rho(u) on the bottom bits, which are independent as in Theorem A. Hence

       sum_(u in block m) phi(t_u) / (1 + u) = (pi^2/12) ln 2 + O(m 2^(-m/2)),

   since integral_0^1 ln(1 + t)/t dt = -Li_2(-1) = pi^2/12. Numerically the block deviation is about 1.4·2^-m.
4. **The top block.** The last block [2^m, U], with U = floor(sqrt X) - 1 not a power of 2, is handled by the same cell
   argument with the extra boundary t <= U/2^m - 1. Hence sum_(u <= U) phi(t_u)/(1 + u) = (pi^2/12) ln U + kappa_1 + o(1),
   with kappa_1 = 0.24748886 (audit, partial sums to 2^31). This gives (pi^2/24) X ln X + kappa_1 X.
5. **The pairs with v < u.** C_u counts all v >= 0, but the theorem counts only v >= u. Subtract the pairs
   v < u <= sqrt X with F1 <= X.
   - With u = lambda sqrt X and v = mu sqrt X, apply the two-dimensional cell argument to the top and bottom halves of
     u and of v, restricted to u > eps sqrt X.
   - The count is X times c' = integral over 0 < mu < lambda < 1 of P(lambda mu (1 + r s) <= 1), plus o(X), with r, s
     independent uniforms. The audit evaluates c' = 0.4863621.

QED.

The predicted block means (pi^2/12)(k + 2) ln 2 + 2 kappa are 16.624975 at k = 28 and 17.195065 at k = 29. The data give
16.625002 and 17.195083.

## 4. (C) Mean multiplicities

For a realisable triple y, every value y_h g(v) - y_e e(v) + y_R (1 + v) is an H-count, hence odd; so is every F1 value.
So Theorem A summed over a fixed bank Y gives, over the odd N <= X (or in a dyadic block):

    mean r_Y = M_Y + O_Y(X^(-1/2)),   M_Y = sum over y in Y of (2 / y_h) ln(1 + y_h / y_R).

Theorem B gives a mean r1 of (pi^2/12) ln X + O(1) over the odd N <= X. The block mean is
(pi^2/12)(k + 2) ln 2 + 2 kappa + o(1), so it gains (pi^2/12) ln 2 = 0.570091 per doubling.

**Bank masses** (distinct rooted triples of tournaments of order <= m; the audit regenerated all of them independently):

| m | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| distinct triples | 1 | 2 | 4 | 12 | 42 | 254 | 2,732 | 51,268 |
| M (one-apex x bank) | 1.386 | 2.197 | 3.106 | 5.134 | 8.380 | 16.954 | 48.198 | 211.236 |

- **m = 1 is the one-apex family itself.** Its mass 2 ln 2 = 1.386 is KNOWN (Session 16/26). It gives about 2.8
  representations per attained value, since it reaches about half of the odd numbers (HYP-9469 (3), THM-4878 (5); the
  audit measured 0.4966-0.4989 per block up to 2^30).
- **The mass is unbounded in the order** (PROVED, audit U2).
  - The order-<=m bank contains the one-apex triples (g(v), e(v), 1 + v) for v < 2^(m-1).
  - Each contributes 2 phi(t_v)/(1 + v) >= 2 ln 2 / (1 + v).
  - Hence M_m >= 2 ln 2 H_(2^(m-1)) >= 2 (ln 2)^2 (m - 1).
- **HEURISTIC.** Using every rooted tournament of order <= m(N) with m(N) growing, the expected number of
  representations of N grows much faster than ln N. Theorem A's error is not uniform once y_R is comparable with N, so
  this is a heuristic, not a theorem. At the level of means, completeness beyond small N has a large margin.

## 5. (D) Data and pre-registered tests

**Pre-registered** (HYP-9478, commit 09b7a56a98; `chanpak_join_multiplicity_20261106.out`, LOGX = 30, reproduced line by
line by the audit's independent sieve):

| k | mean r1 | increment | mean r2 | min r2 | r2 <= 10 | r2 <= 12 | r2 <= 20 | mean r2 by N mod 3 |
|---|---|---|---|---|---|---|---|---|
| 27 | 16.0549 | 0.5701 | 48.1976 | 6 | 4.87e-6 | 2.85e-5 | 3.42e-3 | 60.08, 41.56, 42.95 |
| 28 | 16.6250 | 0.5701 | 48.1976 | 5 | 4.40e-6 | 2.64e-5 | 3.35e-3 | 60.05, 41.62, 42.92 |
| 29 | 17.1951 | 0.5701 | 48.1976 | 4 | 4.62e-6 | 2.57e-5 | 3.33e-3 | 60.11, 41.62, 42.87 |

- **P1-P4 CONFIRMED.**
  - The mean of r2 sits in [48.1966, 48.1986], and the increments in [0.555, 0.585].
  - F2 alone misses only 7 and 21 up to 2^30.
  - The block minimum of r2 is >= 4.
  - The r2 <= 12 and r2 <= 20 tails are below their k = 27 values.
- **Caveats** (audit m7).
  - The program that produced the output was committed after the pre-registration commit, at the checkpoint. The
    audit's independent code gives identical output.
  - P1 and P2 follow from Theorems A and B, which were derived before the commit. So they test the implementation, not
    the laws. P3 and P4 are the genuine out-of-sample tests.

**Structure of the fluctuations** (exploratory; `chanpak_mass_exploratory_20261106/`; the audit's
`aud82_sieve_26_mods2.out` and `aud82_localfactor_28.out`):
- **2-adically flat.** The mean r2 by N mod 8 or 16 is 48.20 in every class (k = 25, 27). This extends S78's F1
  residue profile ("2-adically balanced").
- **A bounded local factor at every small odd prime**, the N ≡ 0 class being over-represented:
  - r2 by N mod 3 is (60.5, 41.3, 42.8) at k = 25 and (60.1, 41.6, 42.9) at k = 27-29;
  - r2 by N mod 5 is 56.28 against 45.9-46.7, and by N mod 7 it is 52.64 against 47.2-48.0 (k = 25);
  - r1 by N mod 3 is (18.1, 12.5, 14.1) at k = 25, and going to mod 9, 27 or 81 does not lower the worst class (12.47,
    12.45, 12.45);
  - over N coprime to 3·5·7·11·13·17 the mean r2 is 0.813 of the overall mean, and over prime N it is 0.796 (k = 26,
    audit U3). The deficit on nonzero classes falls roughly like 1/p^2, consistent with a convergent singular series.
  - These extend S73's 3-adic singular series (HYP-9457 P6; THM-4866 (D)) and S78's F1 residue profile.
  - **The 3-adic factor of F2, predicted** (after the audit; `local_factor_mod3_prediction.py`).
    - Assume g(v) mod 3 is asymptotically uniform and independent of v mod 3, of the parity of v and of rho (Session 193
      modular_sync Theorems 3-4).
    - Then a triple spreads its mass evenly over the three classes unless 3 | y_h and 3 | y_R. In that case
      Phi = -y_e e(v) (mod 3): the odd v (rho < 1/2) land in class 0 and the even v (rho >= 1/2) in class -y_e.
    - The predicted class means are (60.443, 41.411, 42.738). The observed means oscillate around these with decaying
      amplitude: class 0 is 70.8, 64.3, 59.4, 58.2, 59.3, 60.8, 61.3, 60.9, 60.3, 60.0 at k = 10, 12, ..., 28.
    - So the 3-adic local factor of F2 comes entirely from the triples with 3 | y_h and 3 | y_R, split by the parity
      of v.
- **Far from Poisson.**
  - F1 misses 0.0074 of the class N ≡ 1 mod 3 at k = 25, where a Poisson count with mean 12.5 would miss 4e-6.
  - Overall, F1's miss density falls like N^-0.35 (ratio 1.57 per doubling, THM-4875 (D)), against N^-(pi^2/12) =
    N^-0.82 for a Poisson count with mean (pi^2/12) ln N.
  - F2's lower tail P(r2 <= 10) = 4.6e-6 is about 10^5 times the Poisson value 2.7e-11 for mean 48.2.
- **The tail thins slowly and is sampled ever deeper.** For F2 with the fixed bank:
  - the fraction with r2 <= 20 peaks at 4.5e-3 (k = 22) and is 3.3e-3 at k = 29;
  - the block minimum falls 6, 5, 6, 5, 4 over k = 25..29;
  - a crude extrapolation of the tail puts P(r2 = 0) between 1e-14 and 1e-12 per number (HYP-9478's pre-registered
    range). So a fixed bank may well fail somewhere beyond 2^35; P6 tests this to 2^35.
- **Larger banks concentrate.** With the order-<=8 bank (51,268 triples; to 2^24), the mean is 211.2 and the block
  minimum is at least 102 for k = 17..23, against 4 out of a mean of 48 for order <= 7. The audit confirmed this
  independently; F2 with that bank also misses only 7 and 21 to 2^24.

## 6. Coverage to 2^35 (pre-registered HYP-9478 P5-P7)

The pre-registered program `chanpak_f1f2_cover34_20261106.c` makes one random write into a 2 GB bitmap per value. At
2^35 it is memory-bound: after 67 minutes it had not finished F1, and it was stopped before producing any output.

`chanpak_f1f2_cover35_bucketed_20261106.c` computes the same families and output with bucketed writes; it was committed
in 82376badc0 before its run. At 2^28 and 2^30 it reproduces the pre-registered program's and the audit's outputs
exactly.

**Status at close-out (2026-11-06 11:50): PENDING.** The bucketed run started at 11:25 and is expected to finish around 12:45-13:25. Following the pre-registration, P5-P7 are UNDECIDED at close. The result will be added here and in HYP-9478 as a dated addendum. Until then the certified range stays 2^32 (THM-4875 (D), recertified by Session 203), plus F2 alone to 2^30 (P3).

## 7. What a proof through these families would need (HEURISTIC)

- **The shape of the problem.** A proof built on one-arc joins of one-apex tournaments is a lower-bound problem for
  representation counts, not a numeration with free digits.
  - The representations of N are solutions of y_h g(v) - y_e e(v) + y_R (1 + v) = N over a bank, or of the bilinear F1
    equation.
  - Their expected number is M_Y, or (pi^2/12) ln N for F1. It is unbounded when the bank grows (section 4).
  - This Goldbach-type reading is S73's (HYP-9457 (c)). The covering is collective (mac-mini, HYP-9474), and the
    exponential-sum view with major arcs is THM-4881 (D), for the ribbon family.
- **The circle-method form.**
  - Write r_Y(N) = sum_y #{v : Phi_y(v) = N} and expand in additive characters.
  - The major arcs would give M_Y S(N) for a local factor S(N). The data suggest that S is 2-adically trivial and has a
    bounded factor at every small odd prime (section 5).
  - A power-saving bound sum over v in [V, 2V) of e(theta (a g(v) + c v)) << V^(1 - eta) off the major arcs, with
    suitable L^2 control of the minor arcs, would at best give the expected count for almost every N (density one), as
    in Bourgain-Kontorovich.
  - A pointwise lower bound needs control of the minor-arc integral itself, which a family with one free parameter per
    triple does not provide (kind-pasteur S11, `seam_gadget_thickness_spectrum_20261104.md` line 262).
- **Where the major arcs are** (exploratory: S(theta) = sum over v in [V, 2V) of e(theta g(v)), V = 2^12..2^21).
  - theta = 1/2 is trivial: g is always odd (as THM-4881 (D) notes for the ribbon family).
  - The largest remaining sums are at theta = +-1/7 (mod 1/2). There |S|/V goes from 0.199 (V = 2^12) to 0.103
    (V = 2^21), roughly V^-0.11.
  - At exactly theta = 1/3 it goes from 0.096 to 0.025. At 1/21 and 1/49 it is about 0.01 and 0.002 at 2^21, and at 1/5
    it is below 1e-4 from 2^20.
  - g mod 7 carries a bias that rotates with the scale: the low residues are {0, 1} at k = 10, {1, 2} at k = 14,
    {2, 3} at k = 18 and {3, 4} at k = 22, with frequencies 0.79-1.15 of uniform at k = 22. The value set S_1 has a
    7-adic bias too (THM-4878 (5): 0.455-0.581 by residue below 2^21).
  - Session 193 `modular_sync/PROOF.md` Theorem 3 proves that H mod M tends to uniform for uniform strong one-apex
    words, for every fixed M, with an explicit but very slow rate. So the 7-adic factor vanishes eventually. The
    computation below gives the exact rate.
- **Exact rates at rational points** (computed; `chanpak_mass_exploratory_20261106/twisted*.py`).
  - Reading the bits after the top one from the most significant end, g mod q is computed by an automaton on the states
    (P mod q, last bit). A zero bit adds the prefix value P plus [last bit = 1], and g = 1 + the total.
  - So for theta = a/q the sum over v in [2^k, 2^(k+1)) is exactly e(a/q) 1^T T_(q,a)^k delta_(1,1), where T_(q,a) is an
    explicit 2q x 2q twisted transfer matrix. The FFT values agree; for example |S(1/7)|/V = 0.1989 at k = 12.
  - The decay per doubling is the spectral radius divided by 2. Over all a/q with q <= 120 and theta not in (1/2)Z, it is
    at most 0.9320, attained only at theta = +-1/7 (mod 1/2) (`twisted_scan_q120.out`).
  - The next worst are 11/93 (0.9137), 1/31 (0.9002), 1/15 (0.8922) and 1/63 (0.8895); the rate is sqrt(3)/2 = 0.866 at
    1/3. All are Mersenne-related denominators, on which the powers of 2 cycle quickly.
  - This is a power saving V^-0.10 at every rational point of height <= 120. Irrational theta is not covered.
- **Why this is hard.** v rho(v) is essentially a product a·rev(a) of an integer and its bit reversal.
  - The nearest analytic results concern sums with the bit reversal entering linearly:
    Dartyge-Martin-Rivat-Shparlinski-Swaenepoel, "Reversible primes", arXiv:2309.11380, estimate sums of
    e(alpha a + theta rev(a)); Bhowmik-Suzuki, arXiv:2507.08714, give the arbitrary-base version.
  - To our knowledge the product form is not covered.
  - The analogous question for continued fractions with bounded partial quotients is Zaremba's conjecture: density one
    is known (Bourgain-Kontorovich), and the full statement is open. The analogy is kind-pasteur S11's
    (`seam_gadget_thickness_spectrum_20261104.md` lines 45-50, 153-156).
  - A density-one theorem for the H-spectrum looks like the realistic next milestone. (Opinion:) pointwise completeness
    through these families is at least as hard.
- **What is not available.**
  - A fixed bank: P6 tests whether it ever fails.
  - Unbounded F1 multiplicity alone: F1 has misses of positive, decaying density.
  - Bounded-register digit compilers: Sessions 202-203.

## 8. Typing (after the audit)

| claim | type |
|---|---|
| (A) one-readout mass (X/a) ln(1 + a/c) + O(sqrt X) | PROVED; extends Session 16/26's (1,1,1) law by the same coupling |
| (B) F1 count (pi^2/24) X ln X + kappa X + o(X); kappa = -0.2388732 | PROVED (new); kappa computed (audit U1) |
| (C) mean multiplicities M_Y; F1 gains (pi^2/12) ln 2 per doubling; M_m >= 2 (ln 2)^2 (m - 1) | PROVED (new; the last by the audit) |
| m = 1 mass 2 ln 2 | KNOWN (Session 16/26) |
| bank masses for orders 1-8 | FINITE-EXACT (computed; independently regenerated) |
| (D) mass laws to four decimals at 2^28-2^30; F2 alone covers every odd N <= 2^30 except 7, 21 | FINITE-EXACT (pre-registered HYP-9478 P1-P4; independently reproduced) |
| non-Poisson tails, local factors, order-8 concentration, rotating 7-adic bias | exploratory data |
| twisted transfer identity at rational theta; spectral radii for q <= 120 | exact identity; radii computed |
| coverage of F1 u F2 to 2^35 | section 6 |
| H >= 1 + 2 oc(T), equality iff intercyclic; intercyclic spectrum | KNOWN (Grinberg-Stanley; Session 16) |
| section 7 programme | HEURISTIC |
| Chan-Pak completeness | OPEN |

## 9. Audit

`04-computation/experiments/one_apex_join_mass_20261106_audit_S82/AUDIT.md`: **PASS WITH CORRECTIONS.**

The audit reproduced with its own code:
- the banks of orders 1-8, two ways for order <= 7;
- every line k = 10..29 of the 2^30 output;
- F2-alone coverage to 2^30;
- the order-8 bank results to 2^24;
- the one-arc join law and the closed form.

It confirmed that the 2^30 output postdates the pre-registration.

**Majors (all applied).**
- **M1.** Theorem A, its coupling and the one-apex mass were Session 16's (2026-10-20) and Session 26's. (A) is
  retyped as an extension.
- **M2.** Route 3 was closed by Session 16's audited intercyclic theorem. The note's "not classified" was wrong.
- **M3.** Section 7 overclaimed: the stated bound gives at best density one. The local factor is not 3-adic only.
  The framing is S73's, mac-mini's and kind-pasteur S11's.

**Minors (all applied).**
- the b/c term in the uniform error;
- the two gaps in Theorem B;
- oddness wording;
- a mislabelled block;
- the tail-range inconsistency;
- the pre-registration caveats;
- the order-8 recipe;
- HYP-9469 (1) is conjectural;
- the heuristic typing;
- literature hedges;
- credits;
- the van der Corput/Hammersley naming;
- the HYP in-sample "5 to 8" (k = 19 has 10);
- the exact 1/3 values;
- housekeeping.

**Upgrades (adopted).**
- U1: kappa explicit.
- U2: M_m unbounded.
- U3: local factors over primes.
- U4: one-apex density per block to 2^30.
- U5: order-8 F2 to 2^24.
- U6: 1 <= delta <= 1 + #10(v).

## 10. Reproduction

All programs are in `04-computation/experiments/`; the bank is `rooted_triple_bank_order7_20261104.txt` (THM-4875).

    gcc -O3 -march=native -fopenmp -o mult chanpak_join_multiplicity_20261106.c -lm && ./mult 30 rooted_triple_bank_order7_20261104.txt
    gcc -O3 -march=native -fopenmp -o cover chanpak_f1f2_cover35_bucketed_20261106.c && ./cover 35 rooted_triple_bank_order7_20261104.txt

The order-8 bank comes from the THM-4875 generator, keeping the distinct triples:

    gcc -O3 -o bank rooted_triple_bank_20261104.c && ./bank 8 | awk '{print $2, $3, $4}' | sort -u > bank8.txt

It has 51,268 lines. The audit's `aud82_bank_le8.txt` is an independent copy.
