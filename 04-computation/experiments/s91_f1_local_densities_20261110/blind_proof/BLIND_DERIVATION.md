VERDICT. The statement is true, and I prove a sharper form: the 1/ln X term is exhibited and the remainder after it is power-saving. [PROVED by hand below; not machine-checked.] Nothing was refuted. One inaccuracy in the justification of (L1) is corrected in Lemma 2.

NOTATION. B_n = {2^(n-1) <= h < 2^n}, B_0 = {0}; rho = x/y; M = x1 x2 + 2^s y1 y2 = N + [s = 0]; c' = M mod p, so the distinguished class is c' = 0 on both sheets; Y = X/2^s; chi(h) = [p | g(h)] - 1/p^2; w_s(h) = (1/x) ln((x + 2^s y)/(x/2 + 2^s y)); w^0 = w_0; w^1 = sum_{s>=1} w_s = (1/x) ln(1 + x/(2y)) (telescoping). T_s(X) = #{(h1, h2), h_i >= 0 : M <= X}; T^0(X) = T_0(X + 1); T^1 = sum_{s>=1} T_s.

THEOREM.
(A) T^sigma(X) = kappa_sigma X ln X + O(X), with
- kappa_0 = pi^2/12 + 2 Li2(-1/2) - Li2(-1/4) = 0.1615389173,
- kappa_1 = Li2(-1/4) - Li2(-1/2) = 0.2125139092,
- kappa_0 + kappa_1 = (pi^2/6 + 2 Li2(-1/2))/2 = 0.3740528265, shares 0.43186124 and 0.56813876.

For one s: T_s(X) = I_s X ln X + O((s+1) X 2^-s), I_s = l_s - 2 l_(s+1) + l_(s+2), l_j = -Li2(-2^-j).

(B) For each odd prime p there are theta_p in (0, 1/2) and reals A_p^0, A_p^1 such that for every class c

    #{t in T^sigma(X) : N = c mod p} = mu_c T^sigma(X) + ([c = c_sigma] - 1/p) 2 A_p^sigma X + O_p(X^(1 - theta_p/2)),

mu_c = 1/p + (1/p^2)([c = c_sigma] - 1/p), A_p^sigma = sum_{h>=0} chi(h) w^sigma(h) (convergent in the natural order). Hence the density tends to mu_c, and density - mu_c = ([c = c_sigma] - 1/p)(2 A_p^sigma/kappa_sigma)/ln X + O(1/ln^2 X).

LEMMA 1 (L2; uniform law of x/y). [PROVED; FINITE-EXACT for h < 2^24]
(a) q(h) = sum over one-bits i of ceil(H_i/2), H_i = floor(h/2^(i+1)); 0 <= q <= h/2. So y/2 < x <= y, with x = y iff h = 0 or 2^j.
(b) For h in B_n: |rho - (1 - phi/2)| < n/(2h+1) (L2 holds with n instead of n+1).
(c) On B_n, phi takes each value (2i+1)/2^n exactly once, so rho is within n 2^-n of an exactly equidistributed grid of (1/2, 1].

Proof. (a) Appending a digit d shifts the old bits and adds d ceil(h/2); induct. q(2h) = q(h) <= h/2 and q(2h+1) <= h/2 + (h+1)/2. (b) 2(2h+1)(rho - 1 + phi/2) = sum over one-bits of [(2h+1)/2^(i+1) - 4 ceil(H_i/2)]. The first term lies in (2H_i + 1, 2H_i + 2), the second is 2H_i or 2H_i + 2, so each bracket is in (-1, 2). (c) The top bit is fixed, the others free. QED

LEMMA 2 (L1). [PROVED; constants not explicit] For odd m and k != 0, rho(T_k) < 2. Fix rho_m with max(1/2, max_k rho(T_k)/2) < rho_m < 1 (strictly above, because of possible Jordan blocks). Then ||T_k^r|| <= C_m (2 rho_m)^r in the max-row-sum norm (Gelfand), and every cell sum with r free digits is at most C_m (2 rho_m)^r, for any prescribed top and bottom words and any l.

Proof. Spectral radius: if T_k v = lambda v with |lambda| = 2, take sigma with |v(sigma)| maximal. Then 2|v(sigma)| = |v(sigma0) + omega v(sigma1)| forces both successors to be maximal and v(sigma0) = omega_sigma v(sigma1). Every state is reachable (for 2^r >= 2m the endpoints (2^r eta + D mod m, D mod 2) exhaust the states by CRT), so this holds at every state. The states (eta, 0) and (eta, 1) have the same successors, with phases differing by zeta^(k/2). Hence zeta^(k/2) = 1 and k = 0.

Uniformity: the cell sum is e_0^T M_a T_k^r M_b v_l, with M_a, M_b products of one-digit matrices. A one-digit matrix has exactly one non-zero entry, of modulus 1, in each ROW. Its columns have 0 or 2, so "per row and column" in the task is inaccurate; only rows are needed. So e_0^T M_a is a unimodular multiple of a coordinate vector, M_b v_l has unimodular entries, and |u^T A w| <= ||u||_1 ||A||_(inf->inf) ||w||_inf. QED

LEMMA 3 (cells mod p). [PROVED] For a cell C of B_n with r free digits and every v in F_p^2: |#{h in C : (x,y) = v} - 2^r/p^2| <= C'_p (2 rho_p)^r, C'_p = max(C_p, p/2).

Proof. alpha x + beta y = (alpha + beta)(2h+1) - 2 alpha q. If alpha != 0, apply Lemma 2 with k = -2 alpha. If alpha = 0 != beta, the sum over the free digits is geometric with ratio e_p(2 beta 2^j) != 1, so at most 1/sin(pi/p) <= p/2. Invert. QED

LEMMA 4 (regions). [PROVED] Let b >= a > 0, X > 0, Omega = {h : a x + b y <= X}, I an interval, n >= 1. Constants do not depend on a, b, X, I.
(i) #(B_n ∩ I ∩ Omega) = 2^n area{(u,v) in [1/2,1) x [0,1) : 2^n u in I, 2^(n+1) u (a(1 - v/2) + b) <= X} + O(2^(n/2)).
(ii) #{h in B_n ∩ I ∩ Omega : (x,y) = v} = p^-2 #(B_n ∩ I ∩ Omega) + O_p(2^(n(1-theta))), theta = theta_p = ln(1/rho_p)/ln(2/rho_p^2).
(iii) If |W| <= W_0 on B_n and W oscillates by O(W_0 2^-j) on each cell below, then sum over B_n ∩ I of chi W = O_p(W_0 2^(n(1-theta))).

Proof. Use cells with j top and j bottom digits, j <= n/2.
- The top word A gives 2^(n-j+1) A < y < 2^(n-j+1)(A+1).
- The bottom word gives phi in (Phi, Phi + 2^-j), so by Lemma 1 rho lies in (R - 2^(-j-1) - eps, R + eps), with R = 1 - Phi/2 >= 1/2 + 2^(-j-1) and eps = n 2^-n.
- A cell is inner if 2^(n-j+1)(A+1)(a(R + eps) + b) <= X (it lies in Omega), outer if 2^(n-j+1) A (a(R - 2^(-j-1) - eps) + b) >= X (it misses Omega).
- For a fixed bottom word the remaining A lie in an open interval (P - 1, P(1 + 2^(-j-1) + 2 eps)) with P < 2^j; this uses a <= b. There are fewer than 2.5 + 2^(j+1) eps < 5 of them.
- So at most 4·2^j cells are cut by Omega and 2·2^j by I: at most 6·2^(n-j) points.

(i) Take j = floor(n/2); an inner (outer) cell has its rectangle inside (outside) the region and holds exactly 2^n times its area. (ii) Take j = floor(theta n); Lemma 3 on at most 2^(2j-1) inner cells costs C'_p 2^n rho_p^(n-2j) = O(2^(n(1-theta))), because rho_p^(1-2 theta) = 2^-theta. (iii) The same, freezing W on each cell. QED

LEMMA 5 (one free state). [PROVED; this is THM-4885 Theorem A] G(X; h1) := #{h2 > h1 : M <= X} = Gamma_X(h1) + O(sqrt(Y/y1) + 1), where Gamma_X(h1) = int_0^1 (Z(v) - max(h1,1))^+ dv and Z(v) = X/(2 y1 (rho1 (1 - v/2) + 2^s)).
- G = 0 unless y1 < sqrt Y, since M > y1 y2 (2^s + 1/4).
- int_0^1 Z = X w_s(h1), so Gamma_X = X w_s - Xi_X with 0 <= Xi_X(h1) <= max(h1, 1).

Proof. Lemma 4(i) with a = x1 <= b = 2^s y1, summed over n; t = 2^n u runs once over [1, infinity) with 2^n du = dt. QED

LEMMA 6 (class law for a fixed smaller state; (L3) with D1 a point mass). [PROVED] If p does not divide g(h1), then (x2, y2) -> x1 x2 + 2^s y1 y2 is a non-zero linear form on F_p^2 and each value has p preimages. Lemma 4(ii) gives #{h2 > h1 : M <= X, M = c'} = G(X; h1)/p + O_p((Y/y1)^(1-theta)). If p | g(h1), then M = 0 mod p for every h2.

PROOF OF (A).
- M is symmetric, so T_s = 2 sum_{h1} G(X; h1) + diagonal, and the diagonal is at most sqrt Y + 1.
- By Lemma 5, T_s = 2 sum_{y1 < sqrt Y} Gamma_X(h1) + O(Y^(3/4)) = 2X sum_{y1 < sqrt Y} w_s(h1) + O(Y).
- Write w_s = F_s(rho)/y with F_s(rho) = (1/rho) ln((rho + 2^s)/(rho/2 + 2^s)) <= 2^(-s-1), Lipschitz. With the cells j = floor(n/2), sum over B_n of w_s = ln 2 · int_(1/2)^1 F_s + O(2^-s 2^(-n/2)), also for a partial block.
- int_(1/2)^1 F_s = [Li2(-rho/2^(s+1)) - Li2(-rho/2^s)]_(1/2)^1 = I_s <= 2^(-s-2). So T_s = I_s X ln Y + O(Y).
- I_0 = kappa_0 and sum_{s>=1} I_s = l_1 - l_2 = kappa_1. Terms with 2^s > X vanish; the shift X+1 on sheet 0 costs O(ln X). QED

PROOF OF (B).
Step 1. Let D_s(X) be the number of (h1, h2) in T_s(X), h1 != h2, whose smaller state is p-divisible. Summing Lemma 6 over y1 < sqrt Y:

    T_s(X; c') = T_s/p + ([c' = 0] - 1/p) D_s + O_p(Y^(1 - theta/2)),

since sum_{y1 < sqrt Y} (Y/y1)^(1-theta) = O(Y^(1-theta/2)) and the diagonal is O(sqrt Y).

Step 2. D_s - T_s/p^2 = 2 sum_{h1} chi(h1) G(X; h1) + O(sqrt Y) = 2 sum chi Gamma_X + O(Y^(3/4)), by Lemma 5, uniformly in p.

Step 3. On B_n both w_s (W_0 = 2^(-s-n)) and Xi_X = int min(Z, max(h1,1)) (W_0 = 2^n) satisfy Lemma 4(iii).
- Block sums of chi w_s are O(2^-s 2^(-theta n)). So A_(p,s) = sum chi w_s converges, |A_(p,s)| = O_p(2^-s), and its tail beyond sqrt Y is O(2^-s Y^(-theta/2)).
- sum_{y1 < sqrt Y} chi Xi_X = O(Y^(1-theta/2)).
- Hence D_s = T_s/p^2 + 2 A_(p,s) X + O_p(Y^(1-theta/2)).

Step 4. Sheet 0 is s = 0 with X+1. For sheet 1 sum over 1 <= s <= log2 X: constants are uniform in s, sum_s (X/2^s)^(1-theta/2) = O(X^(1-theta/2)), and sum_{s>=1} A_(p,s) = A_p^1 (absolutely convergent by blocks). QED

WHERE THE REQUESTED ERROR TERMS GO.
- Two orders and symmetry: each pair is counted through its smaller state; diagonal O(sqrt Y).
- Small h1: not an error. With the larger state equidistributed (remainder (Y/y1)^(1-theta) per h1), the small divisible states produce exactly the term 2 A_p X.
- Boundary cells: at most 6·2^(n-j) points per block.
- Rounding: none in cells (exactly 2^(n-2j) points each); only O(2^(n/2)) in the area count.
- Large s: T_s = 0 for 2^s > X, and everything decays geometrically in s.
- Order of limits: none is needed, because the cell depth j = floor(theta n) grows with the block. With J fixed (as planned, and as in Session206), Lemma 4 leaves a relative error O(2^-J) + O_J(rho_p^n), and states with fewer than 2J + L digits carry a share O((J+L)/ln X) of the triples. Then X goes to infinity first, J and L afterwards, and no rate comes out.
- (L3): I checked the general statement; Cauchy-Schwarz gives |Rhat| <= (p-1)^2 (p+1) p^-2 E1 E2 <= (p-1) E1 E2. Only its point-mass case is used here.

RATE.
- Confirmed, conditionally on A_p^sigma != 0: the convergence is then exactly of order 1/ln X, with coefficient ([c = c_sigma] - 1/p) 2 A_p^sigma/kappa_sigma. A_p != 0 is NUMERICAL only, because C_p is not explicit.
- The session's "small divisible states" are the positive terms of A_p (h = 1 contributes (8/9) w(1) at p = 3). A_p is the regularised sum and is negative for p = 13 and 31, where the distinguished class approaches its limit from below.
- The remainder is power-saving after subtracting 2 A_p X/T(X), or equivalently with the measured share F = D/T, or with the area sum 2 sum chi Gamma_X/T (error O(X^(-1/4)/ln X), uniform in p). After subtracting only the leading 1/ln X term it is O(1/ln^2 X).
- Proved exponents theta_p/2 (any smaller number): 0.073 (p = 3), 0.116 (p = 5), 0.042 (p = 7), from the rates 0.8660, 0.7418, 0.9320 [VERIFIED].
- In the session's variables the excess is E = p F/(1 - F). F and Delta = 2f - 1/p^2 differ by the discrepancy of the larger state, which is power-saving.

NUMERICAL CHECKS. Script: C:/tmp/wt-s91/04-computation/experiments/s91_f1_local_densities_20261110/blind_proof/check_blind_proof.py. `python3 check_blind_proof.py 24 26` runs in about 75 s and ends with ALL CHECKS PASSED. No output file is kept.
- Constants [VERIFIED, mpmath 30 digits]: closed forms against quadrature; pair normalisation sum (s+1) I_s = pi^2/12, shares 0.196408 and 0.803592.
- Lemma 1 [FINITE-EXACT, 1 <= h < 2^24]: max of |rho - (1 - phi/2)|(2h+1)/n is 0.75.
- Lemma 2 norm inequality [VERIFIED]: all cells with n = 16, J = 3, m in {3, 5, 7, 9, 11, 13}.
- Lemma 4 boundary count [VERIFIED]: 1,320 random regions, at most 2 cut top words per bottom word (proved <= 4).
- Census [FINITE-EXACT, all triples with N <= 2^24, both sheets, also s = 1 and s = 2 alone]:
  - (T - kappa X ln X)/X is stable at +0.0478 (sheet 0) and -0.0891 (sheet 1) for X = 2^16..2^24 [NUMERICAL constants].
  - |T - 2 sum Gamma - diagonal| <= 0.05 X^(3/4), and |D - 2 sum_{h1 divisible} Gamma_X| <= 0.004 X^(3/4), at X = 2^24.

Sheet 0 at X = 2^24 (rem = max_c |T(c') - T/p - ([c'=0] - 1/p) D|/T):

| p | density of c_sigma | mu | F = D/T | area prediction of F | rem |
|---|---|---|---|---|---|
| 3 | 0.432935 | 0.407407 | 0.146125 | 0.146135 | 4.7e-3 |
| 5 | 0.249309 | 0.232000 | 0.062043 | 0.062047 | 3.3e-4 |
| 7 | 0.175210 | 0.160350 | 0.012615 | 0.012620 | 2.2e-2 |
| 13 | 0.079717 | 0.082385 | 0.003557 | 0.003557 | 4.9e-4 |

For p = 7 the remainder is flat from 2^16 to 2^24, consistent with the tiny proved exponent; with the two-sided Delta it is 1.0e-2.

A_p by blocks, h < 2^26 [NUMERICAL]:

| p | A_p^0 | A_p^1 | note |
|---|---|---|---|
| 3 | 0.0557 | 0.0789 | oscillates by 0.002 |
| 5 | 0.02942 | 0.04250 | |
| 7 | about 0.006 | about 0.007 | not converged |
| 11 | 0.0110 | 0.0145 | |
| 13 | -0.00377 | -0.00488 | |
| 17 | 0.01355 | 0.01931 | |
| 31 | -0.00194 | -0.00260 | |
| 43 | 0.00235 | 0.00323 | |

ALREADY IN THE REPOSITORY (not new).
- Lemma 1, Lemma 4(i), Lemma 5 and the summation over the smaller state up to sqrt X are THM-4885 Theorems A and B, in the coordinates u = 2h: `05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md`, sections 2-3 (Hammersley cells, (X/a) ln(1 + a/c)). New here: the split by s with closed forms for kappa_0, kappa_1; the session already has the values (CONTEXT item 5).
- Lemma 6 without a rate is Session206: `session206/weighted-parameter-capacity/WEIGHTED_CAPACITY.md` sections 2-3 and `session206/root-residue-capacity/ROOT_RESIDUE_CAPACITY.md` (endpoint cylinders; R U - d V uniform on the multiples of gcd(R, d, p); limits "X, then epsilon, then L"; "No uniform rate in p, s, or a growing bank is asserted"). New here: the uniform power-saving rate, which is what allows summing over smaller states up to sqrt X, and the constant A_p.
- Twisted matrices and rates for the count alone: the same S82 note, section 7 (0.9320 at 1/7). Mixing and size/residue coupling for one word: `05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md`, Theorems 3-4.
- The census limit is recorded as open in THM-4916 / HYP-9506 (A) according to CONTEXT; I did not open those files.

DELICATE POINTS.
1. Lemma 2: C_m is not explicit and rho_m must be strictly above the spectral rate. So (B) has no explicit constant, and the sign and non-vanishing of A_p are numerical.
2. Lemma 4: the bound on cut cells needs a <= b and j <= n/2. Inner and outer cells must be defined with the true ranges of y and rho (margin n 2^-n), not with the model 1 - phi/2.
3. Step 3 near the diagonal y1 ~ sqrt Y: the oscillation of Xi_X = int min(Z, h1) on a cell is O(2^(n-j)) because Z matters only where Z <= 2^(n+1).
4. Step 1 is dominated by balanced pairs h1 ~ h2 ~ sqrt Y, where the larger state has half the digits; this is why the exponent is theta/2. For slow primes the saving is invisible at census scale; the theorem is asymptotic.
5. A_p converges only conditionally (by blocks or in the natural order; sum |chi| w diverges). On sheet 1 the s- and h-sums are exchanged using the block bounds O(2^-s 2^(-theta n)).
6. Conventions: h = 0 is a state (x = y = 1); excluding it changes A_p by w(0)/p^2 but not the limit. On sheet 0, c = c' - 1 and X becomes X + 1.
7. Not proved here: the O(X) terms of T^sigma as constants (same method, values numerical), and a two-sided version with both states in cells, which would give the session's Delta = 2f - 1/p^2 with a bilinear remainder and a better exponent.