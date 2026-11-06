# Independent audit of opus-2026-11-06-S82 (THM-4885, HYP-9478)

**Auditor:** independent adversarial audit subagent, 2026-11-06.

**Audited version.**
- Pre-registration commit 09b7a56a98 (HYP-9478 predictions, `chanpak_f1f2_cover34_20261106.c`) and the pre-audit
  checkpoint 5bde18fade (note, THM-4885, HYP-9478 results, `chanpak_join_multiplicity_20261106.c/.out`).
- Note line numbers refer to the **checkpoint** file `05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md`
  (256 lines). During the audit the session added, uncommitted, a 12-line bullet "Where the major arcs are" after note
  line 221 and a folder `04-computation/experiments/chanpak_mass_exploratory_20261106/`; these are commented on separately
  (section 4, m15).
- THM line numbers refer to `01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md`, HYP line numbers
  to `05-knowledge/hypotheses/HYP-9478-one-apex-join-mass-laws-and-coverage-to-2-35.md`.

**Scope rule.** Nothing outside this folder was modified. Nothing was committed. Executables were deleted after use.
All programs here share no code with the session's programs (nor with the S78 auditor's).

## Verdict: PASS WITH CORRECTIONS

**The mathematics is correct and every number I checked reproduces.**
- (A) and (B) are correct proofs (with two small gaps that are easily filled, m1-m2). I computed Theorem B's constant
  explicitly: 2 kappa = -0.477746, which matches the 2^30 data to about 2e-5 (upgrade U1).
- (C): all eight bank sizes and all eight bank masses were regenerated from scratch (identical sets of triples).
- (D): an independent sieve to 2^30 (own g, own split law, own bank, own sieve) reproduces **every** line k = 10..29 of
  `chanpak_join_multiplicity_20261106.out` exactly. F2 alone misses exactly {7, 21} up to 2^30; the block minima of r2
  are 5 (k = 28) and 4 (k = 29).
- (E): the pre-registration is sound: the 2^30 output was written 4.5 minutes after the HYP-9478 commit.
- **Not audited:** HYP-9478 P5-P7 and note section 6 (coverage to 2^35). The session's run (`cover34.exe`,
  `C:/tmp/s82/cover35.out`) was still running and its output was empty throughout this audit. P5-P7 remain PENDING.

**Three MAJOR corrections are needed** (all three concern credit or the reach of the conclusions, not a false theorem):
- **M1.** Theorem A, the binary-product mechanism behind it, its special case (a, b, c) = (1, 1, 1) and the one-apex
  mass 2 ln 2 are prior, audited repository results (Session 16, 2026-10-20; Session 26, 2026-10-22), not cited.
  (A) is a generalisation by the same method, not "PROVED (new)".
- **M2.** Route 3 (intercyclic linearisation, note lines 77-85) is closed by an audited theorem the note does not cite:
  the intercyclic spectrum is exactly S_1 ∪ {3·2^t + 1} ∪ {9·2^t − 3}, and it misses at least 1 − ln 2 of the odd
  numbers (Session 16). "Intercyclic tournaments without a common vertex were not classified" is false as a statement
  about the repository.
- **M3.** Section 7 and outcome (E) overclaim. The stated exponential-sum bound plus a growing bank does **not** imply
  r_Y(N) > 0 for all large N (it is the "almost all, not all" situation the repository already recorded); "the data show
  S is 2-adically trivial and 3-adically bounded" omits 5- and 7-adic factors of comparable size that the session's own
  exploratory output contains; and the whole Goldbach / circle-method / Zaremba framing is uncredited prior work (S73,
  mac-mini ribbon invariant, kind-pasteur S11).

## 1. What I reproduced, and how

| claim | auditor's method (files in this folder) | result |
|---|---|---|
| one-apex triple (g(v), e(v), 1 + v) | `aud82_onepex.c` (1): own construction of T(v) (chain c_i -> c_j iff i > j; apex -> c_i iff bit i = 1), subset DP for (h, e_x, R_x) | 8,191 words (all L <= 12, v < 2^L, leading zeros included): 0 failures (R checked for L <= 8) |
| closed form of g; g <= v + 1 | `aud82_onepex.c` (0): own closed form g(v) = 1 + sum over zero bits i < top of (floor(v/2^(i+1)) + bit_(i+1)(v)) vs the recursion | 0 failures for all v <= 2^24 |
| linearisation error delta(v) = g(v) − v rho(v) | `aud82_onepex.c` (2), exact per top bit n | delta in [1, 9.44] for n <= 23 (max about 0.33 n + 1.8); the note's bound abs(delta) <= 2n + 1 holds with room |
| one-arc join law, oddness | `aud82_onepex.c` (3): build the join, count H by DP | 2,000 random one-apex x (B, b) of order <= 7 and 2,000 random F1 joins: 0 failures; all values odd |
| bank sizes and masses, orders 1-8 | `aud82_bank.c`: (i) all 2,131,019 labelled tournaments of orders 1-7 rooted at vertex 0; (ii) independently, isomorphism classes by own brute-force canonical form (Heap's algorithm over all n!), class counts 1, 1, 2, 4, 12, 56, 456 (A000568), then every rooted tournament of order n as (rep of order n − 1) + new vertex in all 2^(n−1) ways; order 8 = 456 x 128 = 58,368 rooted extensions | distinct triples 1, 2, 4, 12, 42, 254, 2,732, 51,268; M = 1.386294361, 2.197224577, 3.105609800, 5.134391326, 8.379469122, 16.954223268, 48.197612382, **211.236115928**; (i) = (ii) at order <= 7; the order-7 set equals `rooted_triple_bank_order7_20261104.txt` and the order-8 set equals the session's `C:/tmp/s82/bank8.txt` exactly (as sets, CR stripped) (`aud82_bank.out`) |
| (D) to 2^30 | `aud82_sieve.c`: own g (closed form on v < 2^15, own split law beyond, checked on all v <= 2^24, 10^7 random v and the endpoints), (u, v-chunk) work items, 16-bit atomic counters with a saturation counter, oddness counter; 4 threads, 24 min | 11,246,843,046 F1 pairs, 33,502,684,574 F2 evaluations, 0 even values, 0 saturations; every k = 10..29 line of the session's `.out` reproduced exactly (means, increments, minima, #r2 = 0, tails, mod-3 means); **F2 alone misses exactly {7, 21} on [1, 2^30]**; max r1 = 111, max r2 = 210 (k = 29); F1-only misses k = 28, 29: 241,881 and 378,695, equal to S78's `aud_cover_32.out` (`aud82_sieve_30.out`) |
| order-8 bank, exploratory claim | same sieve, order-<=8 bank, F2 only, to 2^24 | block means 211.19-211.235 for k = 17..23; block minima 120, 118, 110, 112, 109, 103, 102: "at least 102" confirmed; F2 (order <= 8) alone misses only 7, 21 to 2^24 (`aud82_sieve_bank8_24.out`) |
| Theorem B constant | `aud82_kappa.c` + `aud82_kappa.py` (see U1) | 2 kappa = −0.4777464; predicted block means 16.624975 (k = 28) and 17.195065 (k = 29) vs data 16.625002 and 17.195083 |
| section 5 residue claims | `aud82_sieve_30.out` (k = 25, 27), `aud82_sieve_26_mods2.out` (mod 5, 7), `aud82_localfactor_28.out` (coprime classes, primes) | see section 2.5 and U3 |
| one-apex value density | `aud82_onepex_density.c` to 2^30 | reproduces THM-4878 (5) (521,915 of 1,048,576 below 2^21, all 18 block miss counts) and HYP-9469 (3) (50,006 of 100,000 below 200,000); attained fraction per block 0.4966-0.4989 for k = 17..29; mean multiplicity per odd N 1.38629 = 2 ln 2 for k >= 20 |
| C/H ratios (route 2) | `aud82_cratio.py`, all labelled tournaments of order <= 6 | 32/5 (strong 4-tournament, H = 5, C = 32); order 5: 36/5, 100/11, 28/3 (strong) and 64/5 (not strong): all as stated |
| retained-root sewing algebra | read Session197 `RETAINED_ROOT_SEWING.md` eq. (2) | the matrix in (h, x, y, d) is upper triangular with diagonal (H, C, 0, H): correct |
| uncommitted section-7 bullet | `aud82_expsum.c` (exact rationals, V = 2^10..2^24) | 1/7: 0.1989 (2^12), 0.1023 (2^21); 1/21: 0.0100, 1/49: 0.0020 at 2^21; 1/5 < 1e-4 from 2^20; g mod 7 low classes {0,1} k = 10, {1,2} k = 14, {2,3} k = 18, {3,4} k = 22, frequencies 0.79-1.15 at k = 22: confirmed. At exactly 1/3: 0.0964 at 2^12 (the note says 0.12), 0.0245 at 2^21 |

**Pre-registration timing** (file times in `C:/tmp/s82`, all −0700; commit times from `git log`):
- HYP-9478 was committed at 10:09:03 (author) / 10:09:07 (committer), 09b7a56a98, which is on origin/main.
- `multiplicity.exe` was built at 10:09:34, and `multiplicity30.out` was last written at 10:13:39. It is byte-identical to the
  committed `chanpak_join_multiplicity_20261106.out`.
- The in-sample run `repcount28.out` (09:56:33) stops at k = 27, as the HYP's in-sample section says.
- The predictions P1-P7 are textually unchanged between 09b7a56a98 and 5bde18fade (`git diff -M`).
- So the out-of-sample output postdates the pre-registration. The two caveats are m7.

## 2. Verification of the proofs

### 2.1 Theorem A (note lines 91-117): correct

- **Step 1.** From the closed form, delta(v) = g(v) − v rho(v) = 1 + #10(v) − (sum of the fractional parts of v/2^(i+1)
  over the zero bits i < n). The stated |delta| <= 2n + 1 holds. It can be sharpened (U6):
  - For a maximal run of zero bits at positions p..q, the bits below i are zero down to p, so
    frac(v/2^(i+1)) < 2^(p−i−1). The run therefore contributes less than 1 to the sum.
  - Each such run is topped by exactly one "10".
  - Hence **1 <= delta(v) <= 1 + #10(v) <= 1 + ceil(n/2)**, with delta = 1 iff v has no zero bit below its top bit.
  - Exact values (`aud82_onepex.out`): delta is in [1, 9.44] for n <= 23.
- **Step 2.** The bit split is right:
  - the low floor(n/2) bits fix rho within 2^(−floor(n/2)), because the top-half zero bits contribute at most
    2^(−floor(n/2)) − 2^(−n) < 2^(−floor(n/2));
  - the high ceil(n/2) bits fix t within 2^(−ceil(n/2));
  - the map v -> cell is a bijection onto the 2^ceil(n/2) x 2^floor(n/2) grid, and (t, rho(v)) lies in the half-open
    cell.
- **Step 3.** The region boundary s = (Y/(2^n (1 + t)) − c)/a is decreasing in t, so it meets at most p + q − 1 cells.
  Interior cells are counted exactly. Hence count = 2^n Area_n(Y) + O(2^(n/2)). Only about log2(2(1 + a/c)) blocks have
  a boundary inside the unit square, and the geometric sum gives O(sqrt(Y/c)).
- **Steps 4-5.** Correct, including the substitution w = 2^n (1 + t) and (Y/a) ln(1 + a/c) − 1 for Y >= a + c.
- **v = 0 and small v** contribute O(1).
- **Gap (minor, m1).** The sandwich error is 2 Delta/c with Delta = a(2 log2 X + 1) + b + c. So the uniform error term is
  O(sqrt(X/c) + (a/c) log X + b/c + 1). The note (line 93) drops b/c, and THM line 43 drops b/c and the "+1".
  - For fixed (a, b, c) this is immaterial.
  - For arbitrary b (THM: "b >= 0") the uniform form is false as written: a = c = 1, b = X gives about X extra values.
  - In every application b <= a (e(u) <= 1 <= g(u); y_e <= y_h for all 51,268 bank triples, checked), so nothing
    downstream changes.

### 2.2 Theorem B (note lines 121-144): correct; one gap is filled numerically

- **Steps 1-2.** Correct. g(u) <= 1 + u (checked to 2^24; it also follows from R >= h) gives a/c <= 1, b <= 1, hence a
  per-u error O(sqrt(X/(1 + u)) + log X) and a total of O(X^(3/4)).
- **Step 3.** Correct. The point set (u/2^m, rho(u)) is a Hammersley set (bit-reversed van der Corput; see m13).
  - The claimed block error O(m 2^(−m/2)) is conservative. The actual block deviations
    B(m) = sum over the block of phi(t_u)/(1 + u) − (pi^2/12) ln 2 decay like 2^(−m): B(29) = −2.53e-9, B(30) = −1.30e-9
    (`aud82_kappa.out`).
- **Gap 1 (m2).** Step 4 needs the partial block [2^m, U] when U = floor(sqrt X) − 1 is not a power of 2. It follows from
  the same cell argument with the extra boundary t <= U/2^m − 1.
- **Gap 2 (m2).** Step 5's "plus o(X)" for the pairs v < u needs the two-dimensional version of the cell argument: four
  coordinates (top and bottom halves of u and of v), restricted to u > eps sqrt X.
- **Both gaps are routine.** U1 confirms the whole structure numerically, including the v < u term c' = 0.486362.
- **The pi^2/12 step.** integral_0^1 ln(1 + t)/t dt = −Li_2(−1) = pi^2/12 is classical.

### 2.3 (C) (note lines 146-169): correct, with three wording issues

- **mean r_Y = M_Y + O(X^(−1/2)).** This holds for a fixed bank; the constant depends on Y.
- **Oddness.** "All the values in (A) and (B)" (THM line 53) should read "all the values for realisable triples": (A) is
  stated for arbitrary integers (m4). I checked oddness of every produced value: 0 even values in 4.5e10 evaluations.
- **The block-mean formula** (pi^2/12)(k + 2) ln 2 + 2 kappa is right: [Count(2X) − Count(X)]/(X/2).
- **Table, lines 159-162.** All 16 entries reproduced (section 1).
- **Line 164 (one-apex mass).** The mass 2 ln 2 = 1.386 and "about 2.8 representations per attained value" are right:
  per block the attained fraction is 0.4966-0.4989 for k = 17..29 and the mean multiplicity is 1.38629. But both the mass
  and the half-density are prior results (M1, m12).
- **Lines 166-169 ("the expected number of representations of N tends to infinity much faster than ln N").** This is a
  heuristic: Theorem A's error is not uniform once y_R is comparable with N. Its weaker part, M_m -> infinity, is now
  proved (U2).

### 2.4 Section 1 (the four routes)

- **Route 2** (lines 68-76) is right.
  - Session197 eq. (2) is upper triangular with diagonal (H, C, 0, H).
  - C/H = 32/5 for the strong 4-tournament, and 36/5, 100/11, 28/3, 64/5 occur at order 5 (64/5 at the non-strong
    padding of the 4-tournament).
  - The interval readout structure follows from unrolling the triangular product.
  - Erdős's 3-smooth fact is classical (Erdős-Lewin, Math. Comp. 65 (1996); Blecksmith-McCallum-Selfridge, Amer. Math.
    Monthly 105 (1998)); add a citation.
- **Route 3** (lines 77-85): the mathematics stated is right, but the route is closed by prior art (M2).
  - H >= 1 + 2 oc(T) by the Grinberg-Stanley odd-cycle formula, with equality iff there are no two vertex-disjoint odd
    cycles, iff there are no two vertex-disjoint cycles (every cycle's vertex set induces a strong subtournament, which
    contains a 3-cycle).
  - A common vertex implies one-apex.
- **Route 1 and the obstruction list** (lines 57-67). Session 202 (`append-binary/PROOF.md`: parents 5, 9, 17) and
  Session 203 (`bounded-defect-sparsity` O_K(log^(2K+1) X); `binary-lift-obstruction` Theorem 3, O(log M);
  `digital-context-obstruction`) are cited correctly.
  - "The k_min(h) >= h/4 deletion increment (HYP-9469, THM-4876)" is presented as an obstruction. It is HYP-9469 (1), a
    NUMERICAL conjecture; only special cases are proved (THM-4876 (2), (7)) (m9).

### 2.5 Section 5 (exploratory) checked against the auditor's data

| statement (line) | auditor | verdict |
|---|---|---|
| 2-adically flat, r2 = 48.20 in every class mod 8, 16 (188) | k = 25: [48.191, 48.206] mod 8, [48.184, 48.212] mod 16; k = 27: [48.196, 48.200], [48.195, 48.202] | confirmed |
| r2 by N mod 3 "about (60.1, 41.6, 42.9) ... at k = 25" (189) | k = 25: 60.531, 41.262, 42.800; (60.1, 41.6, 42.9) are the k = 27-29 values | mislabelled (m5) |
| r1 by N mod 3 about (18.1, 12.5, 14.1) at k = 25 (189-190) | 18.106, 12.535, 14.103 | confirmed |
| worst r1 class mod 9, 27, 81 stays near 12.45 (190) | 12.467, 12.451, 12.446 (k = 25) | confirmed |
| F1 misses 0.0074 of N = 1 mod 3 at k = 25; Poisson(12.5) misses 4e-6 (193) | 0.00741; e^(−12.535) = 3.6e-6 | confirmed |
| F1 miss density ~ N^(−0.35), ratio 1.57; Poisson N^(−0.82) (194-195) | ratios 1.551-1.572 for k = 24..29; log2(1.57/2) = −0.349 | confirmed |
| P(r2 <= 10) = 4.6e-6 vs Poisson 3e-11 (196) | 4.616e-6 (k = 29); P(Poisson(48.1976) <= 10) = 2.73e-11 | confirmed |
| r2 <= 20 fraction peaks 4.5e-3 at k = 22, 3.3e-3 at k = 29; minima 6, 5, 6, 5, 4 for k = 25..29 (197-198) | 4.527e-3; 3.334e-3; 6, 5, 6, 5, 4 | confirmed |
| P(r2 = 0) "between 1e-14 and 1e-11" (199) | HYP-9478 line 62 says "between 1e-14 and 1e-12" | inconsistent (m6) |
| order-8 bank: mean 211.2, block minimum >= 102 for k = 17..23 (201-202) | 211.19-211.235; minima 120 ... 102 | confirmed |
| "the 3-adic deficit is a constant-factor local density ... not an unbounded depletion" (190-191) | true for 3 alone, but 5, 7, ... add comparable factors; over all primes the depletion is bounded (about 0.80) | incomplete (M3 (iii), U3) |

## 3. MAJOR issues

### M1. Theorem A, its mechanism and the one-apex mass are prior, audited results; (A) is typed "PROVED (new)"

**Prior art (all on origin/main before S82, not cited by S82):**
- **Session 16, `05-knowledge/results/feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md`** (commit
  cb3e06b22e, 2026-10-20; independently audited, `density-audit/AUDIT.md` line 3, "Status: ACCEPT"):
  - line 19, W(X) = #{reduced one-apex words with H <= X} = X log 2 + O(sqrt X), at all cutoffs;
  - lines 41-77, the exact binary-product identity H(b) − N(b) R(b)/2^m = 1 + A(b) − B(b), with R(b) the
    reversed-complement number, i.e. g(v) = v rho(v) + O(log v);
  - lines 79-126, the disjoint-prefix coupling: the first s free digits of U_m and of V_m use disjoint index sets, and the
    coupled U, V are independent uniforms. This is exactly the note's "two independent coordinates" (step 2);
  - lines 148-194, the dyadic-phase cancellation, i.e. the integral in step 4;
  - lines 21-35, |S_1 ∩ [1, X]| <= (ln 2/2) X + O(sqrt X).
- **Session 26,
  `05-knowledge/results/session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md`**
  lines 90-124 (2026-10-22, reviewed in `session26-cut-review/REVIEW.md` line 22): the same product law and
  W(X) = (log 2) X + O(sqrt X), eq. (4).
- **Later restatements:**
  - `session133_fractran_spectrum_refinement_20261101/SYNTHESIS.md` line 66;
  - `session193_triangulated_proof_20261105/modular_sync/PROOF.md` line 185 and Theorem 4 (line 154), the joint
    size/residue law (H/2^n, H mod M) => (UV, R_M);
  - `metagraph_spectrum_torus_session16_20261020.md` line 37.
- **The weighted-inversion formula** of note line 222 (g − 1 = 2 sum b_i (1 − b_j) w(i − j), w(1) = 1, w(l) = 2^(l−2))
  is Session 16 `proofs.md` section 2 (lines 17-40) and THM-4876 (6) (line 107: N(A) = sum c(b − a), c(1) = 1,
  c(d) = 2^(d−2)). Note section 0 (lines 41-45) credits only the THM-4875 audit and Session197.

**Consequence.**
- Theorem A for (a, b, c) = (1, 1, 1) is Session 16's (1): Phi(w) = g(w) − e(w) + 1 + w = g(2w), and the reduced words
  are the even v >= 2.
- The general case is a direct extension by the same coupling and the same integral.
- The note's own m = 1 column (line 164, "2 ln 2 = 1.386") is therefore KNOWN.

**Still new.** Nothing in the repository has (as far as I could find; searches below):
- the (a, b, c) generalisation;
- the pi^2/24 bilinear law (B);
- the mass formula M_Y for a bank and the eight masses.

I searched for "ln(1+", "log(1+", "pi^2/12", "pi^2/24", "zeta(2)/2", "0.57009", "48.197", "211.236" (all tracked files),
"representation mass", "mean multiplicity" and "expected number of representations".

**Retype.**
- (A): "PROVED; extends the audited one-apex word law W(X) = X ln 2 + O(sqrt X) (Session 16 BINARY_PRODUCT_DENSITY (1);
  Session 26 (4)) to every readout a g − b e + c(1 + v), by the same binary-product coupling".
- (B), (C): "PROVED (new)".
- The m = 1 mass: KNOWN.
- Typing table lines 237 and 239, THM status line 5, note lines 15-16 and 116-117.

### M2. Route 3 is closed by an audited theorem; "Intercyclic tournaments without a common vertex were not classified" is false

**Prior art.** `05-knowledge/results/feedback_word_hamilton_spectrum_20261020/intercyclic_rays.md` (commit 422c69c14b,
2026-10-20; audit `audit/INTERCYCLIC_AUDIT.md` line 3: "Status: ACCEPT"):
- line 9: **S_ic = S_1 ∪ {3·2^t + 1 : t >= 2} ∪ {9·2^t − 3 : t >= 1}**, the complete Hamiltonian-path spectrum of the
  intercyclic tournaments;
- line 17: proved via Bang-Jensen-Bessy, *Small cycle transversals in tournaments*, Thm 1.8 (feedback number <= 2), and
  a complete reduced-word grammar;
- lines 72-76: H = 1 + 2(c_x + c_y + c_xy), which is the note's "exact linear law H = 1 + 2 oc", with the Grinberg-Stanley
  / Irving-Omar citation;
- line 181: "the proposed H = 1 + 2·(number of odd cycles) intercyclic strategy cannot establish full odd-spectrum
  coverage";
- `BINARY_PRODUCT_DENSITY.md` (3) and line 35: at least 1 − ln 2 = 30.7% of the odd numbers lie outside S_ic at every
  cutoff;
- `metagraph_spectrum_torus_session16_20261020.md` lines 11-13 and 41-43 restate it.

**What the note says instead** (lines 82-85):
- "Intercyclic tournaments without a common vertex were not classified."
- "The examples examined ... have boundedly many cycle types, so their counts form sparse sets."
- "So the exact linear law H = 1 + 2 oc has so far yielded only the one-apex half."

The intercyclic tournaments that are not one-apex contribute exactly the two dyadic rays, O(log X) values: this is the
sparsity the session observed in its examples.
THM-4876 (6) (line 107) also already has "every odd cycle passes through u, so H(X) = 1 + 2 N(A)", and notes that a second
monotone centre adds nothing.

**Retype.**
- Route 3: "KNOWN and closed (Session 16 intercyclic_rays (1)-(2); BINARY_PRODUCT_DENSITY (3))".
- Typing row line 244: "H >= 1 + 2 oc(T), equality iff intercyclic: KNOWN (Grinberg-Stanley; Session 16)".
- Replacement text: section 6.

### M3. Section 7 and outcome (E) overclaim, misdescribe the local factor, and leave the framing uncredited

**(i) The sufficiency claim is not valid.**
- Note line 221 says: "A proof needs ... sum over v in [V, 2V) of e(theta(a g(v) + c v)) << V^(1 − eta) ... together with a
  bank growing with N so that M_Y -> infinity. **Then r_Y(N) > 0 for all large N.**"
- r_Y(N) = integral of S(theta) e(−N theta), where S has about X M_Y / 2 terms and the main term is about M_Y.
- A sup bound on the minor arcs bounds their contribution only by sup |S| (or by integral |S|). For a family with one
  free variable per triple, that is far larger than M_Y however the bank grows: sup |S| ~ |Y| X^(1 − eta) against
  M_Y <= 2|Y| ln 2.
- Power-saving minor-arc bounds for such "binary" problems give, at best (with L^2 / Parseval control of the minor arcs,
  which the stated sup bound alone does not even supply for small eta), asymptotics for almost all N (density one), not
  for every N. This is precisely the Bourgain-Kontorovich situation the note invokes at line 225.
- The repository already says so: `seam_gadget_thickness_spectrum_20261104.md` (kind-pasteur S11, audited) line 262,
  "Such a theorem would at best give 'almost all' odd N, not all", and lines 154-155.
- The note contradicts itself: line 226 says only density one is realistic.

**(ii) "Needs" is stated for every proof.**
- Lines 29-31 say: "So completeness needs a pointwise lower bound for the representation count ... The missing analytic
  input is an exponential-sum bound for the one-apex function g".
- The data concern the two families F1, F2. They say nothing about proofs through other constructions; the section 7
  header (line 210) correctly restricts to "a proof built on one-arc joins of one-apex tournaments".
- Even for those families, the bound named would not suffice, by (i).

**(iii) The local factor is misdescribed.**
- Line 216 says "The data show S is 2-adically trivial and 3-adically bounded". Data can only suggest. More importantly,
  the factor is not purely 3-adic.
- At k = 25 (`aud82_sieve_26_mods2.out`):

  | | N ≡ 0 | other classes |
  |---|---|---|
  | r2 mod 5 | 56.28 | 45.94, 46.16, 45.95, 46.65 |
  | r2 mod 7 | 52.64 | 47.22-47.99 |
  | r1 mod 5 | 17.71 | 13.98-14.90 (mean 14.91) |
  | r1 mod 7 | 17.07 | 14.05-15.24 |

- The 5- and 7-adic excesses (r2: x1.17, x1.09) are of the same order as the 3-adic one (x1.26).
- The session's own exploratory output contains the mod-5 and mod-7 means, identical to mine
  (`chanpak_mass_exploratory_20261106/r2resid_26.out` lines 7-8), but the note reports only 2- and 3-adic ones.
- Cumulatively (U3), N coprime to 3·5·7·11·13·17 has mean r2 = 0.813 of the overall mean, and prime N has 0.796 (k = 26).
  So the depletion stays bounded, as the note suggests, but it is not 3-adic alone.

**(iv) Unsupported comparisons.**
- "Pointwise completeness is at least as hard as the Zaremba-type step" (line 227) is an opinion stated as fact.
- "Mauduit-Rivat-type methods ... not yet this long-range form" (lines 223-224) is an unreferenced literature claim (m11).

**(v) The framing is prior work and is not credited.** The note's only credit is the blanket "the S73-S81 notes" (line 57).
- **The Goldbach-type reading, the comet and the 3-adic singular series** are S73's (opus-2026-11-04):
  - HYP-9457 (c), lines 3 and 106-109: "As a Goldbach-type statement it would follow from a lower bound on representation
    counts of the kind the comet shows";
  - P6, lines 73-98: the 3-adic singular series, ratio 1.37;
  - THM-4866 (D), lines 106-112;
  - `single_cycle_floor_crevasses_lorentzian_joins_20261104.md` section 5, lines 227-276, including the generic local
    factor 1 + p/(p^3 − 1).
- **The F1 residue profile is S78's.** `old_mark_bound_lorentzian_gadgets_one_apex_join_coverage_20261104.md` lines
  305-312 found "2-adically balanced ... Their deficit is 3-adic", and line 302 states the density-one conjecture for F1.
- **The collective-covering, local-limit and exponential-sum framing is mac-mini's.**
  - HYP-9474 lines 22-26: "The covering is therefore collective".
  - HYP-9474 (2), line 37: minor arcs for S_n(theta) = sum e(theta H(u)).
  - THM-4881 (D), lines 89-106: exponential sums of H with major arcs at denominators 2·3^i·5^j, and "1/2 is a major arc
    since H is odd" (the audit correction in `ribbon_reachable_set_audit_20261104.md` line 99). The uncommitted bullet's
    "theta = 1/2 is trivial" repeats this.
- **The Zaremba / Bourgain-Kontorovich analogy is kind-pasteur S11's:** `seam_gadget_thickness_spectrum_20261104.md`
  lines 45-50 and 153-156.

**Retype.**
- Section 7 is HEURISTIC, with credits.
- Replace the sufficiency sentence and (E) as in section 6.
- Replace line 216 with the measured local factors (U3).

## 4. MINOR issues

1. **Uniform error of (A)** (note line 93; THM line 43). Add "+ b/c + 1": O(sqrt(X/c) + (a/c) log X + b/c + 1). As stated,
   THM line 43 fails for large b (section 2.1). It is harmless for b <= a, which covers every application.
2. **Theorem B proof gaps** (note lines 136-142).
   - Add the partial-block argument for U not a power of 2.
   - Add one sentence for the two-dimensional cell argument behind "plus o(X)" for v < u.
   - Optionally give kappa explicitly (U1).
3. **Block error** (line 134). O(m 2^(−m/2)) is valid but loose; numerically the deviation is about 1.4·2^(−m)
   (`aud82_kappa.out`). No correction is needed.
4. **Oddness wording** (note line 148; THM line 53). Write "every value y_h g(v) − y_e e(v) + y_R(1 + v) with y realisable,
   and every F1 value, is an H-count, hence odd". (A) itself is for arbitrary integers.
5. **Mislabelled block** (line 189). "r2 about (60.1, 41.6, 42.9) ... at k = 25" should read "(60.5, 41.3, 42.8) at
   k = 25; (60.1, 41.6, 42.9) at k = 27-29".
6. **Inconsistent tail extrapolation.** Note line 199 says P(r2 = 0) is in [1e-14, 1e-11]; HYP-9478 line 62 says
   [1e-14, 1e-12]. Use one range.
7. **Pre-registration caveats.**
   - (a) 09b7a56a98 froze the predictions and `chanpak_f1f2_cover34_20261106.c`, but not the program that produced the
     P1-P4 output: `multiplicity.exe` was built 27 s after the commit, from a source first committed in 5bde18fade.
     This is harmless, since the auditor's independent code gives identical output, but it should be stated.
   - (b) P1 and P2 follow from Theorems A and B, which HYP-9478 line 22 says were derived before the commit. So they test
     the implementation, not the mass law. P3 and P4 are the genuine out-of-sample tests.
   - (c) The exploratory programs before the commit (`r2tail.exe` 10:06:52, `r2resid.exe` 10:03:13) printed to the
     console. Their range cannot be checked from files, but everything reported in-sample stops at k = 27.
8. **The order-8 bank is not reproducible from the note.** The note and THM give no recipe for the 51,268-triple bank.
   - The session's `C:/tmp/s82/bankgen.exe` has no source in the scratch folder. Its 58,653-line raw output matches
     `rooted_triple_bank_20261104.c` run with NMAX = 8 (one line per (class, vertex), orders 1-8).
   - Add that command and the dedup step to section 10, and commit the bank or the generator invocation.
   - The auditor's `aud82_bank_le8.txt` is an independent copy.
9. **HYP-9469 (1) is cited as an obstruction** (line 58). "k_min(h) >= h/4" is a NUMERICAL conjecture. Write
   "(conjectured, HYP-9469 (1); proved for increment 2 and for sole-exit cuts, THM-4876 (2), (7))".
10. **Heuristic typed as fact** (lines 166-169). "tends to infinity much faster than ln N" and "The difficulty is entirely
    in the fluctuations" should be typed HEURISTIC. M_m -> infinity itself is now proved (U2).
11. **Literature claim** (lines 223-224). The nearest analytic results concern exponential sums with the binary digit
    reversal, which is linear in rev(a): Dartyge-Martin-Rivat-Shparlinski-Swaenepoel, "Reversible primes",
    arXiv:2309.11380, estimates sums of e(alpha a + theta rev(a)); Bhowmik-Suzuki, arXiv:2507.08714, gives the
    arbitrary-base version. g(v) ≈ v rho(v) is a product a·rev(a). Cite these and hedge "not yet" as "to our knowledge".
12. **Further credits.**
    - The "half of the odd numbers" (lines 85 and 164-165) is HYP-9469 (3) (NUMERICAL, line 40), THM-4876 (6) (line 107)
      and THM-4878 (5) (FINITE-EXACT to 2^21, line 101). The proved upper bound is <= ln 2 (Session 16 (2)). The auditor
      extends the per-block attained fraction to 2^30: 0.4966-0.4989 for k = 17..29 (`aud82_onepex_density.out`).
    - Section 5's 2-adic and 3-adic observations extend S78's F1-miss residue profile (M3 (v)) and S73's singular series.
13. **Classical objects.**
    - rho(v) is the base-2 van der Corput radical inverse of the complemented low bits, and {(v/2^n, rho)} is a
      Hammersley point set. Its discrepancy O(n/2^n) on boxes explains the 2^(−m) block deviations of m3. Name it in the
      proof of (A).
    - integral_0^1 ln(1 + t)/t dt = −Li_2(−1) = pi^2/12.
14. **HYP-9478 in-sample line 39.** "The block minimum of r2 is 5 to 8 for k = 18..27" is wrong at k = 19, where the minimum
    is 10 (both the session's `.out` and the auditor's). Write "5 to 10".
15. **Uncommitted section-7 bullet** (working tree, after line 221). It is correct as data (section 1), with three fixes:
    - (a) "theta = 1/3 ... 0.12 (V = 2^12)" is the value at the nearest FFT grid point j/(8V); at exactly 1/3 it is
      0.0964 at 2^12 and 0.0245 at 2^21. The session's later `twisted.out` agrees.
    - (b) Credit THM-4878 (5) for the mod-7 bias of the one-apex value set: 0.455-0.581 by residue below 2^21.
    - (c) Credit Session 193 `modular_sync/PROOF.md` Theorem 3 (line 138). For uniform strong one-apex words,
      H mod M tends to uniform on the admissible residues for every fixed M, with the explicit but extremely slow rate
      (1 − T/2^L)^floor((n − 2)/L); for M = 7, L = 42. So the "slowly decaying" 7-adic factor is proved to vanish
      eventually, and the observed V^(−0.11) is consistent with that theorem.
16. **Housekeeping.** THM-4885 status still says "AUDIT: pending" (line 12); note sections 6 and 9 are placeholders.
    THM (D) and the HYP results should cite the auditor's independent 2^30 reproduction.

## 5. Upgrades proved or computed by the auditor

### U1. The constant of Theorem B: 2 kappa = −0.4777464

From the proof's own decomposition, kappa = kappa_1 − c', where:
- kappa_1 = lim_U [sum over u < U of phi(t_u)/(1 + u) − (pi^2/12) ln U], with phi(t) = ln(1 + t)/t and t_u = g(u)/(1 + u);
- c' = integral_0^1 dlambda integral_0^lambda dmu P(lambda mu (1 + r s) <= 1), with r, s independent uniforms.

**Computation.**
- The partial sums D(m) at U = 2^m are exact up to floating point for m <= 31 (`aud82_kappa.out`): D(31) = 0.247488858178.
  The block deviations decay like 2^(−m), so kappa_1 = 0.24748886 (+-2e-9).
- With G(p) = P(rs <= 1/p − 1) = 1 for p <= 1/2 and w(1 − ln w) with w = 1/p − 1 otherwise, and Gi(q) = integral_0^q G,
  c' = integral_0^1 Gi(q)/(2q) dq = 0.4863620732. A direct two-dimensional quadrature agrees to 3e-11
  (`aud82_kappa_py.out`).
- Hence kappa = −0.2388732 and 2 kappa = −0.4777464.

**Test against the data.** The predicted block means (pi^2/12)(k + 2) ln 2 + 2 kappa are:

| k | predicted | data |
|---|---|---|
| 25 | 14.91470 | 14.9148 |
| 28 | 16.624975 | 16.625002 |
| 29 | 17.195065 | 17.195083 |

The data are exact, from the cumulative pair counts in `aud82_sieve_30.out`. This independently confirms Theorem B's
whole structure, including the v < u correction, which contributes −c' X. The note's "2 kappa = −0.4777 (k = 29)" becomes
a theorem-level constant.

### U2. The bank mass is unbounded in the order (PROVED)

- The bank of order <= m contains the 2^(m−1) distinct apex triples (g(v), e(v), 1 + v), v < 2^(m−1), of the one-apex
  tournaments T(v), since T(v) has order <= m.
- Since phi(t) = ln(1 + t)/t >= ln 2 on (0, 1] and t_v <= 1, each such triple contributes
  (2/g(v)) ln(1 + g(v)/(1 + v)) = 2 phi(t_v)/(1 + v) >= 2 ln 2/(1 + v).
- Hence M_m >= 2 ln 2 · H_(2^(m−1)) >= 2 (ln 2)^2 (m − 1).
- More precisely, by U1's sums the one-apex part alone is (pi^2/6)(m − 1) ln 2 + 2 D(m − 1); for m = 7 that is 7.37 of the
  48.20.
- This proves the note's "It grows without bound when the bank is allowed to grow" (line 214); it is not a heuristic.

### U3. Local factors at every small prime; the depletion over primes is bounded (FINITE-EXACT data, k <= 27)

- **Mod 5 and 7** (`aud82_sieve_26_mods2.out`): see the table in M3 (iii). Multiples of 5 and of 7 are over-represented
  in r1 and r2.
- **Over N avoiding small primes** (`aud82_localfactor_28.out`), mean r2 at k = 26 (r1 in brackets) as a fraction of
  the overall mean:

  | N restricted to | r2 | r1 |
  |---|---|---|
  | coprime to 3 | 0.875 | 0.892 |
  | coprime to 15 | 0.842 | 0.845 |
  | coprime to 105 | 0.831 | 0.825 |
  | coprime to 1155 | 0.826 | 0.819 |
  | coprime to 255255 | 0.813 | 0.809 |
  | primes | 0.796 | 0.804 |

- The per-prime deficits on nonzero classes fall roughly like 1/p^2 (3: 0.128, 5: 0.042, 7: 0.015). This is consistent
  with a convergent singular-series product. So the "bounded local factor" picture survives, but with factors at every
  odd prime.
- Primes have the lowest local density. The block minimum of r2 over primes is 6-7 for k = 25..27 (overall 5-6).

### U4. One-apex value set to 2^30

- The attained fraction of the odd numbers per dyadic block is 0.4976, 0.4991, 0.4966, 0.4974, 0.4973, 0.4975, 0.4973,
  0.4982, 0.4980, 0.4984, 0.4985, 0.4988, 0.4989 for k = 17..29.
- Mean multiplicity is 1.38629 = 2 ln 2 for k >= 20 (`aud82_onepex_density.out`). This extends THM-4878 (5) and
  HYP-9469 (3) from 2^21 to 2^30.
- The half-density remains conjectural; the proved upper bound is ln 2.

### U5. The F2 family with the order-<=8 bank alone covers every odd N <= 2^24 except 7 and 21

FINITE-EXACT, auditor's run (`aud82_sieve_bank8_24.out`): the block minima of r2 are >= 102 for k = 17..23.

### U6. A sharp form of the linearisation (PROVED)

For v >= 1 with top bit n,

    1 <= delta(v) = g(v) − v rho(v) <= 1 + #10(v) <= 1 + ceil(n/2),

with delta = 1 iff v = 2^(n+1) − 1.

*Proof.* delta = 1 + #10(v) − sum over the zero bits i < n of frac(v/2^(i+1)).
- Group the zero bits into maximal runs p..q. Each run is topped by a 1 at q + 1 <= n, so it contains exactly one "10".
- For i in the run, v mod 2^i = v mod 2^p < 2^p, so frac(v/2^(i+1)) < 2^(p−i−1). The run's total is therefore < 1.
- Hence 0 <= sum of fracs < #10(v) whenever a zero bit exists. ∎

This is Session 16's bound (5), |H − NR/2^m| <= 1 + m/2, made one-sided. It replaces the note's |delta| <= 2n + 1 in
step 1. It matches the computed range [1, 9.44] for n <= 23.

## 6. Exact suggested corrections

1. **Note lines 15-16, 89, 116-117; THM status line 5 and title.** After (A) add: "(A) extends the audited one-apex word
   law W(X) = X ln 2 + O(sqrt X) (Session 16 `feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md` (1),
   2026-10-20; Session 26 APEX_CLASSIFICATION_AND_SPECTRAL_LAW (4)), which is its case (a, b, c) = (1, 1, 1); the proof
   is the same binary-product coupling (BINARY_PRODUCT_DENSITY (4)-(8)) applied to a g − b e + c(1 + v)."
   - Typing row 237: "PROVED (extends Session 16/26)".
   - Line 164: "2 ln 2 = 1.386 (KNOWN, Session 16/26)".
2. **Note section 0, line 45.** Add: "This is Session 16's weighted-inversion formula (`proofs.md` section 2) and
   THM-4876 (6)."
3. **Note lines 77-85 (route 3).** Replace with:
   "3. **Intercyclic linearisation (KNOWN and closed).** By the odd-cycle formula H = I(Omega, 2) (Grinberg-Stanley),
   H(T) >= 1 + 2 oc(T), with equality iff T is intercyclic. The intercyclic tournaments are classified: by
   Bang-Jensen-Bessy (feedback number <= 2) their spectrum is S_1 ∪ {3·2^t + 1 : t >= 2} ∪ {9·2^t − 3 : t >= 1}
   (Session 16 `intercyclic_rays.md` (1), audited ACCEPT). So the linear law reaches at most (ln 2/2) X + O(sqrt X)
   values up to X (BINARY_PRODUCT_DENSITY (3)), and numerically about half the odd numbers (HYP-9469 (3), THM-4878 (5)).
   The session's examples (three transitive apexes; two apexes through a hinge) agree with that theorem."
   - Typing row 244: "KNOWN (Grinberg-Stanley; Session 16)".
4. **Note line 221.** Replace "Then r_Y(N) > 0 for all large N." with:
   "Such a bound (with suitable L^2 control of the minor arcs) would at best give the expected count M_Y S(N) for almost
   every N (density one), as in Bourgain-Kontorovich. A pointwise lower bound needs control of the minor-arc integral
   itself, which a one-parameter-per-triple family does not provide (compare kind-pasteur S11,
   `seam_gadget_thickness_spectrum_20261104.md` line 262)."
5. **Note lines 29-31 (E).** Replace with:
   "So a proof through these join families needs a pointwise lower bound for their representation counts, a Goldbach-type
   statement (S73, HYP-9457 (c)), not a numeration with free digits. Exponential-sum methods for g would give density-one
   statements (cf. THM-4881 (D), HYP-9474 for the ribbon family); pointwise completeness needs more."
6. **Note line 216.** Replace "The data show S is 2-adically trivial and 3-adically bounded (section 5)" with:
   "The data suggest that S is 2-adically trivial and has a bounded local factor at every small odd prime. At k = 25, r2 by
   N mod 3, 5 and 7 has its N ≡ 0 class at x1.26, x1.17 and x1.09 of the mean, and over prime N the mean is 0.80 of the
   overall mean (audit U3)."
7. **Note line 227.** "Pointwise completeness is at least as hard as the Zaremba-type step" becomes "(opinion)". Credit
   kind-pasteur S11 for the Zaremba analogy.
8. **Note section 5.** Fix line 189 (m5) and line 199 or HYP line 62 (m6). Add credits to S78 lines 305-312 and to S73's
   singular series. Add the mod-5 and mod-7 means.
9. **THM line 43.** Use "O(sqrt(X/c) + (a/c) log X + b/c + 1)". **THM line 53.** "For realisable triples ... hence odd."
10. **THM (B), or note line 155.** Add "kappa = kappa_1 − c' = −0.2388732 (audit U1)".
11. **Note section 10.** Add the order-8 bank recipe (m8). **HYP-9478 line 39.** Fix "5 to 8" to "5 to 10" (m14).
12. **MISTAKE ledger.** Recommend a MISTAKE entry: Session 16/26's audited W(X) law and intercyclic theorem were not found
    before (A) was typed new and route 3 was called unclassified. The searches that would have caught them are
    `git grep -n "X log 2"` and `git grep -n -i intercyclic`.

## 7. Files (this folder)

| file | purpose |
|---|---|
| `aud82_bank.c`, `aud82_bank.out`, `aud82_bank_le2.txt` ... `aud82_bank_le8.txt` | independent banks of orders <= 2..8 (labelled mode to 7, class mode to 8); counts and masses |
| `aud82_onepex.c`, `aud82_onepex.out` | one-apex triple, closed form, delta range, one-arc join law with direct counts |
| `aud82_sieve.c`; `aud82_sieve_22.out`, `aud82_sieve_30.out` (+ `.time`), `aud82_sieve_bank8_24.out`, `aud82_sieve_26_mods2.out` | independent r1/r2 sieve; the 2^30 reproduction; order-8 bank; mod 5/7 class means (build with `-DMODS2`; the `#ifdef MODS2` blocks were added after the 2^30 run and do not change the default build) |
| `aud82_kappa.c`, `aud82_kappa.out`, `aud82_kappa.py`, `aud82_kappa_py.out` | U1: kappa_1 partial sums to 2^31, c' by quadrature, predicted block means |
| `aud82_localfactor.c`, `aud82_localfactor_28.out` | U3: means over N coprime to 3, 15, ..., 255255 and over primes |
| `aud82_onepex_density.c`, `aud82_onepex_density.out` | U4: one-apex value density and mass to 2^30 |
| `aud82_cratio.py`, `aud82_cratio.out` | route 2: C/H ratios for orders <= 6 |
| `aud82_expsum.c`, `aud82_expsum.out` | uncommitted section-7 bullet: exponential sums at a/q and g mod 7 |

**Reproduction** (from this folder, gcc 15 with OpenMP, Python 3.13 with mpmath):

    gcc -O3 -march=native -fopenmp -o aud82_bank aud82_bank.c -lm && ./aud82_bank 7 8
    gcc -O3 -march=native -o aud82_onepex aud82_onepex.c -lm && ./aud82_onepex 12 24 2000
    gcc -O3 -march=native -fopenmp -o aud82_sieve aud82_sieve.c -lm && ./aud82_sieve 30 aud82_bank_le7.txt 1 25
    ./aud82_sieve 24 aud82_bank_le8.txt 0 21
    gcc -O3 -march=native -fopenmp -DMODS2 -o aud82_sieve_mods2 aud82_sieve.c -lm && ./aud82_sieve_mods2 26 aud82_bank_le7.txt 1 23
    gcc -O3 -march=native -fopenmp -o aud82_kappa aud82_kappa.c -lm && ./aud82_kappa 31 && python aud82_kappa.py
    gcc -O3 -march=native -fopenmp -o aud82_localfactor aud82_localfactor.c -lm && ./aud82_localfactor 28 aud82_bank_le7.txt
    gcc -O3 -march=native -o aud82_onepex_density aud82_onepex_density.c -lm && ./aud82_onepex_density 30
    python aud82_cratio.py 6
    gcc -O3 -march=native -o aud82_expsum aud82_expsum.c -lm && ./aud82_expsum 24

**Web checks (2026-11-06).**
- MathOverflow 232751 (bof, Royle): all odd values to 80,555 except 7 and 21 (cited by arXiv:2308.04512 and Grinberg's
  Rédei-Berge notes). The note's line 10 is correct.
- No published resolution of Chan-Pak Conjecture 6.2 was found.
- Schweser-Stiebitz-Toft (arXiv:2510.10659) concerns parity, not the spectrum.
- No literature result on the density or multiplicity of Hamiltonian path counts was found.
