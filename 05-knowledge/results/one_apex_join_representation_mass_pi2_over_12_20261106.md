# Representation mass of the one-apex join families: the pi^2/12 law, the constant 48.1976, and what a proof of the {7,21} spectrum still needs

**Session:** opus-2026-11-06-S82. Worktree `C:/tmp/wt-s82`. Canon: THM-4885. Pre-registration: HYP-9478 (commit
09b7a56a98).

**Owner prompt (verbatim).** "spend another session trying to prove the compelteness of the {7,21} spectrum"

**Target.** Chan-Pak Conjecture 6.2 (TCS 1015 (2024) 114776, arXiv:2308.10214): every odd number except 7 and 21 is the
number H(T) of Hamiltonian paths of a tournament. THM-1370 proves that 7 and 21 never occur. The conjecture is OPEN.
A web search on 2026-11-06 found no published resolution. MathOverflow 232751 records all odd values to 80,555 except
7 and 21.

**Outcome.** No all-order proof. The session measured what any proof built from the known join families would have to
control, and proved the first-moment laws behind them.
- **(A, B) Mass laws (PROVED).**
  - A one-apex readout a g(v) - b e(v) + c(1 + v) takes at most X on (X/a) ln(1 + a/c) + O(sqrt X) values of v.
  - The one-apex x one-apex joins number (pi^2/24) X ln X + kappa X + o(X) up to X.
- **(C) Mean multiplicities (PROVED).**
  - A one-apex x bank family has mean multiplicity sum_y (2/y_h) ln(1 + y_h/y_R) per odd number: 48.1976 for the
    order-<=7 bank and 211.236 for order <= 8.
  - The one-apex x one-apex family has (pi^2/12) ln N.
- **(D) Pre-registered tests (CONFIRMED).**
  - Both laws hold to four decimals out of sample.
  - The fixed order-<=7 bank alone already covers every odd N <= 2^30 except 7 and 21.
  - Coverage of the two families together to 2^35 is reported in section 6.
- **(E) What is missing.**
  - The counts are strongly non-Poisson. The one-apex x one-apex family misses a density about N^-0.35 of the odd
    numbers, where an independent model with its mean would miss N^-0.82.
  - So completeness needs a pointwise lower bound for the representation count, i.e. a circle-method or thin-orbit
    estimate of Zaremba type, not a numeration with free digits. The missing analytic input is an exponential-sum bound
    for the one-apex function g (section 7).

## 0. Inherited objects

**One-apex tournaments.**
- T(v) is the transitive chain plus an apex that beats the chain vertex i iff bit i of v is 1.
- Its rooted triple at the apex is (g(v), e(v), 1 + v), with g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2) and
  e(v) = [v = 0 or v even].
- Credits: session25; session179 spectrum/REPORT.md (1); session124's f; Session191 one_apex/PROOF.md; THM-4876 (6);
  THM-4875 (D).
- **Closed form** (THM-4875 audit, section 2 (D); Session197 ONE_CORE_BINARY_INVERSE (4)):

      g(v) = 1 + #10(v) + sum over the zero bits i of v of floor(v / 2^(i+1)),

  checked for all v < 2^16 here.

**One-arc join** (session29 ADDENDUM; THM-4866 (B)). H = x_h y_h - x_e y_e + x_R y_R.

**Families** (THM-4875 (D), audited; recertified to 2^32 by Session 203 with F0 disabled):
- F1 = {(1 + u)(1 + v) + g(u) g(v) - e(u) e(v) : 0 <= u <= v};
- F2 = {y_h g(v) - y_e e(v) + y_R (1 + v) : v >= 0, y in the 2,732 distinct rooted triples of order <= 7}.

**Representation counts.** For odd N, r1(N) counts the pairs u <= v and r2(N) counts the pairs (y, v).

## 1. The attempt: four routes, and where each stalls

Prior work: Sessions 191-204 (navigation S200-S203) and the S73-S81 notes. Recent obstructions:
- the k_min(h) >= h/4 deletion increment (HYP-9469, THM-4876);
- append-only binary classification (Session 202): physical appends give both 2n - 1 and 2n + 1 only for n = 5, 9, 17;
- bounded-defect sparsity (Session 203): sources with Z-defect <= K have O(log^(2K+1) X) distinct values;
- bounded-offset binary lifting (Session 203): from a fixed seed bank it has only O(log M) outputs;
- the one-port digital-context obstruction (Session 203).

The four routes:
1. **H-affine digit operations** (H -> 2H +- 1, 3H + {0, 2, 4}, applied to every member of a class). These need the
   auxiliary data (end counts, root responses) to be affine in H on the class, i.e. bounded defects. That is sparse by
   Session 203, and impossible as small increments by THM-4876. No new route here.
2. **Erdos-type two-register numerations.** Every integer is a sum of numbers 2^a 3^b, none dividing another (Erdos).
   A family with H = 1 + 2 sum 2^(a_i) 3^(b_i) over free staircases would therefore be complete except where its
   constraints bite.
   - The stationary right-sewing algebra (Session197 RETAINED_ROOT_SEWING (2)) is upper triangular in (h, x, y, d), with
     diagonal characters (H(B), C(B), 0, H(B)), C = labelled two-bin covers.
   - Its readouts are sums over intervals of blocks of prod_outside H x prod_inside C. A block's ratio C/H is not
     always a power of 2 (32/5 for the strong 4-tournament; 36/5, 100/11, 28/3, 64/5 at order 5), so mixed radices
     occur.
   - But every interval of the word contributes; there is no free choice of terms. No free-staircase family was found.
3. **Intercyclic linearisation.** By the odd-cycle formula H = I(Omega, 2), H(T) >= 1 + 2 oc(T) (oc = number of odd
   directed cycles).
   - Equality holds iff no two odd cycles are disjoint, iff T has no two vertex-disjoint cycles (a cycle's vertex set
     induces a strong subtournament, which contains a 3-cycle).
   - If the cycles share a common vertex, then T minus that vertex is acyclic, so T is one-apex.
   - Intercyclic tournaments without a common vertex were not classified. The examples examined (three transitively
     inserted apexes; two apexes forced through one hinge vertex) have boundedly many cycle types, so their counts form
     sparse sets.
   - So the exact linear law H = 1 + 2 oc has so far yielded only the one-apex half of the odd numbers.
4. **Counting.** Every odd N has many representations in F1 u F2 (section 5). A proof by counting needs the mean
   (sections 2-4) and pointwise control of the fluctuations (section 7). This is the route the session pursued.

## 2. (A) The mass of one one-apex readout (PROVED)

**Theorem A.** Let a, c >= 1 and b >= 0 be integers, and Phi(v) = a g(v) - b e(v) + c (1 + v). Then

    C(X) = #{v >= 0 : Phi(v) <= X} = (X / a) ln(1 + a / c) + O(sqrt(X / c) + (a / c) log X + 1).

*Proof.*
1. **Linearisation.** For v >= 1 with top bit n, put rho(v) = sum over zero bits i < n of 2^-(i+1). By the closed form,
   g(v) = v rho(v) + delta(v), with |delta(v)| <= 2n + 1: each floor loses less than 1, and #10(v) <= n. Hence
   |Phi(v) - v (c + a rho(v))| <= Delta = a (2 log2 X + 1) + b + c.
2. **Two independent coordinates.** For v in [2^n, 2^(n+1)), write v = 2^n (1 + t).
   - rho(v) depends only on the n low bits. Bit i has weight 2^-(i+1), so the lowest floor(n/2) bits fix rho to within
     2^-floor(n/2).
   - The remaining ceil(n/2) bits fix t to within 2^-ceil(n/2).
   - The two bit sets are disjoint, so as v runs over the block, the pair (top part, bottom part) runs over the full
     product grid of 2^ceil(n/2) x 2^floor(n/2) cells, one v per cell. The point (t, rho(v)) lies in its cell's
     rectangle.
3. **Counting a block.** The region {(t, s) in [0,1)^2 : 2^n (1 + t)(c + a s) <= Y} is bounded by a monotone curve,
   which meets O(2^(n/2)) cells. So the block contributes 2^n Area_n(Y) + O(2^(n/2)) values v with
   v (c + a rho(v)) <= Y.
4. **Summing the blocks** up to n = log2(Y / c) (beyond that v c > Y): the error is O(sqrt(Y / c)), and
   sum_n 2^n Area_n(Y) = integral over s in [0, 1] of (Y / (c + a s) - 1)_+ ds = (Y / a) ln(1 + a / c) - 1 for
   Y >= a + c.
   - This uses that w = 2^n (1 + t) runs exactly once over [1, infinity) as n and t vary, with measure 2^n dt = dw.
5. **Sandwich.** {Phi <= X} lies between the sets for Y = X - Delta and Y = X + Delta. Their counts differ by at most
   (2 Delta / a) ln(1 + a / c) <= 2 Delta / c, plus the grid error. QED.

The proof uses only that rho is a bit-reversal of the low bits. Rho(v) runs exactly once through the grid j / 2^n as
the n low bits vary.

## 3. (B) One-apex x one-apex joins: the pi^2/12 law (PROVED)

**Theorem B.** #{0 <= u <= v : F1(u, v) <= X} = (pi^2 / 24) X ln X + kappa X + o(X) for a constant kappa.

*Proof.*
1. For fixed u, F1(u, v) is Phi(v) with (a, b, c) = (g(u), e(u), 1 + u). Only u with (1 + u)^2 <= X contribute.
2. Theorem A applies with uniform error O(sqrt(X/(1 + u)) + log X), since a / c = g(u)/(1 + u) <= 1. Summed over
   u <= sqrt X, the error is O(X^(3/4)). So

       sum_u C_u(X) = X sum_(u <= sqrt X) phi(t_u) / (1 + u) + O(X^(3/4)),
       phi(t) = ln(1 + t) / t,   t_u = g(u) / (1 + u) = rho(u) + O(log u / u).

3. **Each dyadic block of u.** Over u in [2^m, 2^(m+1)), rho(u) runs exactly once through j / 2^m. The weight 1/(1 + u)
   depends on the top bits and rho(u) on the bottom bits, which are independent as in Theorem A. Hence

       sum_(u in block m) phi(t_u) / (1 + u) = (integral_0^1 phi) ln 2 + O(m 2^(-m/2)),   integral_0^1 ln(1 + t)/t dt = pi^2/12.

4. The errors are summable, so sum_(u <= U) phi(t_u)/(1 + u) = (pi^2/12) ln U + kappa_1 + o(1). This gives
   (pi^2/24) X ln X + kappa_1 X.
5. **The pairs with v < u.** C_u counts all v >= 0; the theorem counts only v >= u. Subtract the pairs v < u <= sqrt X
   with F1 <= X.
   - With u = lambda sqrt X and v = mu sqrt X, their count is X times the measure of
     {mu < lambda : lambda mu (1 + rho rho') <= 1}, averaged over independent uniform rho, rho', plus o(X).
   - This is a constant times X.

QED.

## 4. (C) Mean multiplicities

Every value is an H-count, hence odd. So Theorem A summed over a bank Y gives, over the odd N <= X (or in a dyadic
block):

    mean r_Y = M_Y + O(X^(-1/2)),   M_Y = sum over y in Y of (2 / y_h) ln(1 + y_h / y_R).

Theorem B gives a mean r1 of (pi^2/12) ln X + O(1) over the odd N <= X. The block mean is
(pi^2/12)(k + 2) ln 2 + 2 kappa + o(1), so it gains (pi^2/12) ln 2 = 0.570091 per doubling. The data give
2 kappa = -0.4772 (k = 20) and -0.4777 (k = 29).

**The bank mass grows quickly with the order** (distinct rooted triples of tournaments of order <= m):

| m | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| distinct triples | 1 | 2 | 4 | 12 | 42 | 254 | 2,732 | 51,268 |
| M (one-apex x bank) | 1.386 | 2.197 | 3.106 | 5.134 | 8.380 | 16.954 | 48.198 | 211.236 |

- m = 1 is the one-apex family itself: 2 ln 2 = 1.386, about 2.8 representations per attained value, since it reaches
  about half of the odd numbers.
- Using every rooted tournament of order <= m(N) with m(N) growing, the expected number of representations of N tends
  to infinity much faster than ln N.
- At the level of means, completeness beyond small N has an enormous margin. The difficulty is entirely in the
  fluctuations.

## 5. (D) Data and pre-registered tests

**Pre-registered** (HYP-9478, commit 09b7a56a98; `chanpak_join_multiplicity_20261106.out`, LOGX = 30).

| k | mean r1 | increment | mean r2 | min r2 | r2 <= 10 | r2 <= 12 | r2 <= 20 | mean r2 by N mod 3 |
|---|---|---|---|---|---|---|---|---|
| 27 | 16.0549 | 0.5701 | 48.1976 | 6 | 4.87e-6 | 2.85e-5 | 3.42e-3 | 60.08, 41.56, 42.95 |
| 28 | 16.6250 | 0.5701 | 48.1976 | 5 | 4.40e-6 | 2.64e-5 | 3.35e-3 | 60.05, 41.62, 42.92 |
| 29 | 17.1951 | 0.5701 | 48.1976 | 4 | 4.62e-6 | 2.57e-5 | 3.33e-3 | 60.11, 41.62, 42.87 |

- **P1-P4 CONFIRMED.** The mean of r2 sits in [48.1966, 48.1986] and the increments in [0.555, 0.585]. F2 alone misses
  only 7 and 21 up to 2^30. The block minimum of r2 is >= 4, and the r2 <= 12 and r2 <= 20 tails are below their k = 27
  values.
- **Convergence of the means.** The F2 mean sits at 48.1976 from k = 24. The F1 increment is within 0.0001 of
  (pi^2/12) ln 2 from k = 23.

**Structure of the fluctuations** (exploratory, in sample unless marked):
- **2-adically flat.** The mean r2 by N mod 8 or 16 is 48.20 in every class (k = 25, 27).
- **3-adic local factors, bounded.** By N mod 3 the means are r2 about (60.1, 41.6, 42.9) and r1 about (18.1, 12.5,
  14.1) at k = 25. Going to mod 9, 27 or 81 does not lower the worst class (r1 stays near 12.45). So the 3-adic
  deficit is a constant-factor local density, like a singular series, not an unbounded depletion.
- **Far from Poisson.**
  - F1 misses about 0.0074 of the class N = 1 mod 3 at k = 25, where a Poisson count with mean 12.5 would miss 4e-6.
    Overall, F1's miss density falls like N^-0.35 (ratio 1.57 per doubling, THM-4875 (D)), against N^-(pi^2/12) =
    N^-0.82 for a Poisson count with mean (pi^2/12) ln N.
  - F2's lower tail P(r2 <= 10) = 4.6e-6 is about 10^5 times the Poisson value 3e-11 for mean 48.2.
- **The tail thins slowly and is sampled ever deeper.** For F2 with the fixed bank, the fraction with r2 <= 20 peaks at
  4.5e-3 (k = 22) and is 3.3e-3 at k = 29. The block minimum falls 6, 5, 6, 5, 4 over k = 25..29. A crude
  extrapolation of the tail puts P(r2 = 0) somewhere between 1e-14 and 1e-11 per number. So a fixed bank may well fail
  somewhere beyond 2^35. P6 tests this to 2^35.
- **Larger banks concentrate.** With the order-<=8 bank (51,268 triples; exploratory, to 2^24), the mean is 211.2 and
  the block minimum is at least 102 for k = 17..23, against 4 out of a mean of 48 for order <= 7.

## 6. Coverage to 2^35 (pre-registered HYP-9478 P5-P7)

(Filled when `chanpak_f1f2_cover34_20261106` finishes at LOGX = 35.)

## 7. What a proof of completeness needs

- **The shape of the problem.** The data and (A)-(C) say that a proof built on one-arc joins of one-apex tournaments
  is a lower-bound problem for representation counts, not a numeration with free digits.
  - The representations of N are solutions of y_h g(v) - y_e e(v) + y_R (1 + v) = N over a bank, or of the bilinear
    F1 equation.
  - Their expected number is M (or (pi^2/12) ln N). It grows without bound when the bank is allowed to grow with N.
- **The circle-method form.** Write r_Y(N) = sum_y #{v : Phi_y(v) = N} and expand in additive characters. The major arcs
  give M_Y S(N) for a local factor S(N). The data show S is 2-adically trivial and 3-adically bounded (section 5). A
  proof needs, for theta away from rationals with small denominator,

      sum over v in [V, 2V) of e(theta (a g(v) + c v)) << V^(1 - eta),

  uniformly in a <= c, together with a bank growing with N so that M_Y -> infinity. Then r_Y(N) > 0 for all large N.
- **Why this is hard.** g is 2-regular but quadratic in the bits: g - 1 = 2 sum_(i>j) b_i (1 - b_j) w(i - j), with
  w(1) = 1 and w(l) = 2^(l-2). So v rho(v) is a bit-reversal product. Mauduit-Rivat-type methods handle sums of digits
  and Rudin-Shapiro-type block counts, but not yet this long-range form. The analogous question for continued fractions
  with bounded partial quotients is Zaremba's conjecture: density one is known (Bourgain-Kontorovich), and the full
  statement is open. A density-one theorem for the H-spectrum looks like the realistic next milestone. Pointwise
  completeness is at least as hard as the Zaremba-type step.
- **What is not available.**
  - A fixed bank: P6 tests whether it ever fails.
  - Unbounded F1 multiplicity alone: F1 has misses of positive, decaying density.
  - Bounded-register digit compilers: Sessions 202-203.

## 8. Typing

| claim | type |
|---|---|
| (A) one-readout mass (X/a) ln(1 + a/c) + O(sqrt X) | PROVED (new) |
| (B) F1 count (pi^2/24) X ln X + kappa X + o(X) | PROVED (new) |
| (C) mean multiplicities M_Y; F1 gains (pi^2/12) ln 2 per doubling | PROVED (corollary) |
| bank masses for orders 1-8 | FINITE-EXACT (computed) |
| (D) mass laws to four decimals at 2^28-2^30; F2 alone covers every odd N <= 2^30 except 7, 21 | FINITE-EXACT (pre-registered HYP-9478 P1-P4) |
| non-Poisson tails, 3-adic local factors, order-8 concentration | exploratory data |
| coverage of F1 u F2 to 2^35 (section 6) | see section 6 |
| H >= 1 + 2 oc(T), equality iff no two disjoint cycles | elementary (odd-cycle formula) |
| Chan-Pak completeness | OPEN |

## 9. Audit

(Filled after the independent audit.)

## 10. Reproduction

All programs are in `04-computation/experiments/`; the bank is `rooted_triple_bank_order7_20261104.txt` (THM-4875).

    gcc -O3 -march=native -fopenmp -o mult chanpak_join_multiplicity_20261106.c -lm && ./mult 30 rooted_triple_bank_order7_20261104.txt
    gcc -O3 -march=native -fopenmp -o cover chanpak_f1f2_cover34_20261106.c && ./cover 35 rooted_triple_bank_order7_20261104.txt
