# One-apex exponential sums: a Type II bound, Fourier uniformity, the energy at scale 2^(n/3), and the target for positive density

**Session:** opus-2026-11-06-S83. Worktree `C:/tmp/wt-s83`. Canon: THM-4886, THM-4887. Pre-registration: HYP-9479
(commit 185c2b3efa). Audits: section 8.

**Owner prompt (verbatim).** "start another session working toward a density-one theorem for the spectrum"

**Target.** Chan-Pak (TCS 1015 (2024) 114776) asks that the Hamiltonian path counts H(T) of tournaments cover every odd
number except 7 and 21. That is OPEN. A density-one theorem (almost every odd number is a count) is a weaker milestone.

S82 (THM-4885, section 7) named two analytic inputs a circle-method proof would need:
- **minor-arc bounds** for the exponential sums of the one-apex counts, known only numerically;
- **pair correlations** (energy) of the join families.

The S82 audit (M3) and Session 205 (`session205_apex_scope_audit_20261106/S82_SCOPE_CORRECTION.md`, section 3) then
showed that individual minor-arc bounds, even with a growing bank, do not give the covering inequality.

**Outcome.**
- **(A) The first minor-arc bound for these families, uniform in theta.**
  - Session 193's exact concatenation law makes H bilinear across two blocks. The classical Weyl-Vinogradov Type II
    estimate then gives |sum_x e(theta H(x))| <= 2^(n-1) (6/q + 7 · 2^(a-n) + 8(1 + ln q) 2^(-a) + 8q(1 + ln q) 2^(-n))^(1/2)
    whenever |2 theta - r/q| <= 1/q^2.
  - On the minor arcs this saves 4 sqrt(n) 2^(-eps n/2).
  - This two-block form is the independent audit's upgrade (U1). It supersedes the session's three-block version (A').
- **(B)** The same holds for every one-apex readout y_h H - y_e e + y_R (1 + A) of the join families.
- **(C) Fourier uniformity.** sum_x e(theta H(x)) = o(2^n) for every real theta outside (1/2)Z.
- **(D) Pre-registered and CONFIRMED out of sample** (HYP-9479):
  - consistent with square-root cancellation at irrational theta;
  - the energy E_n (pairs of words with equal count) stays at 3.25 |X_n|;
  - the bilinear second moment equals its diagonal within 7%.
- **(E) What remains.**
  - *Positive density.* Linear energy E_n = O(2^n) is sufficient. The converse involution forces E_n >= 2 |X_n| - 2^(n/2) - 1.
    By Parseval, the minor-arc bounds give only E_n << n 2^(5n/3).
  - *Density one.* A route needs pair correlations of the bank families.
- **(F) The energy down to the scale 2^(n/3) (THM-4887).**
  - A two-block identity g(h 2^m + l) = g(h) + h C_m(l) + O(2^m) and a lattice-point count give
    #{(v, v') : |g(v) - g(v')| <= Delta} <= (1.39 n + 19) 2^(4n/3) for Delta <= 2^(n/3). This is sharp up to O(n).
  - Hence E_n = O(n 2^(4n/3)), and the number of Hamiltonian path counts up to X is at least X^(2/3)/(13 ln X + 121).
    The previous all-X exponent was 1/2.
  - Positive density now needs equidistribution of g only below the scale 2^(n/3).

## 0. Objects and inherited facts

**Words and counts** (Session 193 `modular_sync/PROOF.md` section 2).
- For a binary word w, T(w) is the increasing chain plus an apex beating the i-th chain vertex iff w_i = 1. H(w) is its
  Hamiltonian path count; H(empty) = 1.
- A(w) is the binary integer of w (most significant letter first), and C(w) = sum_j (1 - w_j) 2^(j-1) (the first letter
  has weight 1).
- p(w) is the last bit, and eps(w) = [w nonempty and begins in 0].
- H(w) = g(A(w)), with g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2).

**Concatenation laws** (Session 193 eqs. (1), (2); checked here on 20,000 random triples, and by the audit exhaustively
for short words). For words u, y:

    H(uy) = H(u) + H(y) - 1 + A(u) C(y) + p(u) eps(y).                                                   (1)

For |v| = L:

    H(uvw) = H(u) + H(v) + H(w) - 2 + A(u) C(v) + (2^L A(u) + A(v)) C(w) + p(u) eps(v) + p(v) eps(w).    (2)

**The sum.** X_n = words of length n beginning with 1 (so A runs over [2^(n-1), 2^n)), and
S_n(theta) = sum_(x in X_n) e(theta H(x)).

**Prior work on the Fourier side.**
- Session 193's Theorems 3-4 give the size/residue coupling and the exponential mixing of H mod M for every fixed M.
  Theorem 3 is a proved decay at each fixed rational theta, with the explicit, slow rate 2(1 - T/2^L)^floor((n-2)/L).
- S82 (THM-4885 note section 7) computed the exact rational-point rates with twisted transfer matrices (q <= 120). The
  worst is 0.9320 per doubling, at +-1/7. It also cites the reversible-prime sums of Dartyge-Martin-Rivat-Shparlinski-
  Swaenepoel and Bhowmik-Suzuki, whose phases are linear in the digits.
- The binary-product law H ~ A R / 2^m and the converse involution H(w) = H(w*) (w* the reverse complement) are Session
  16/26's (`BINARY_PRODUCT_DENSITY.md`).
- For the ribbon family, mac-mini measured minor-arc decay numerically and conjectured an integer-scale local limit law
  (HYP-9474; THM-4881 (D); open).
- The "at best almost all" limitation of circle methods is kind-pasteur S11's. Its sharpening for banks of readouts is
  the S82 audit's M3 and Session 205's scope correction.
- The repository had no uniform-in-theta (minor-arc) bound for any of these families. The audit confirmed this with a
  full-tree search: no minor-arc or uniform bound, no energy computation, no bilinear argument.
- **Method.** The argument is the classical Weyl-Vinogradov Type II (bilinear) estimate: Cauchy-Schwarz over one block,
  geometric sums over a progression in the other, and the standard bound for sum min(N, 1/||alpha x||) when
  |alpha - a/q| <= 1/q^2 (Vaughan, *The Hardy-Littlewood Method*, 2nd ed., ch. 2; Davenport, *Multiplicative Number
  Theory*, ch. 25). The new input is the concatenation law, which makes H bilinear.

## 1. (A) The Type II bound (PROVED)

**Theorem A** (the audit's U1). Let n >= 3, 2 <= a <= n - 1, and let r/q be in lowest terms with
|2 theta - r/q| <= 1/q^2. Then

    |S_n(theta)| <= 2^(n-1) (6/q + 7 · 2^(a-n) + 8 (1 + ln q) 2^(-a) + 8 q (1 + ln q) 2^(-n))^(1/2).

*Proof.*
1. **Split.** Write x = u y with |u| = a (u begins with 1) and |y| = m = n - a.
   - Fix p = p(u). Then A(u) runs over 2^(a-1) + p + 2j for 0 <= j < K = 2^(a-2).
   - By (1), S_p = sum_u e(theta H(u)) sum_y b_y e(theta A(u) C(y)), with b_y = e(theta (H(y) - 1 + p eps(y))).
2. **Cauchy-Schwarz over u.** |S_p|^2 <= K sum_(y,y') |sum_j e(2 theta j (C(y) - C(y')))|.
   This is at most K sum_(y,y') min(K, 1/(2 ||2 theta (C(y) - C(y'))||)).
   - C is a bijection onto [0, 2^m). So the diagonal contributes K^2 2^m.
   - Each nonzero difference d occurs at most 2^m times.
3. **Counting lemma.** Suppose |beta - r/q| <= 1/q^2. Then in any block of q consecutive integers d,
   #{d : ||beta d|| <= t} <= 2qt + 3.
   - Write beta d = beta d_0 + r m/q + m delta, with |m delta| < 1/q.
   - The terms r m/q run over all classes mod q.
   - So the sorted distances in a block satisfy t_j >= (j - 3)/(2q), and the block contributes at most
     3K + sum_(j=4..q) q/(j - 3) <= 3K + q(1 + ln q).
   - At most 2(2^m/q + 1) blocks cover 0 < |d| < 2^m.
4. **Assemble.**
   - |S_p|^2 / (K 2^m)^2 <= 2^(-m) + 2(1/q + 2^(-m))(3 + q(1 + ln q)/K).
   - Then |S_n| <= |S_0| + |S_1|, with K 2^m = 2^(n-2).

QED.

- The perturbations m delta all have one sign, so the counting lemma holds with +2 (audit U2). This turns 6 and 7 into
  4 and 5.
- Every inequality of the chain was checked on 420 exact cases, and the bound against exact sums at every admissible
  q: 1,186 instances, 0 violations (audit `a5`, `a7`).

**Corollaries** (n >= 7, a = ceil(n/2); numerical suprema over n <= 300 in brackets, audit `a9`).
- **(A1)** |S_n| <= 3.52 · 2^(n-1) q^(-1/2) for 2 <= q <= 2^((n-3)/3) [3.37].
- **(A2)** |S_n| <= 4 · 2^(n-1) (q n/2^n)^(1/2) for q >= 2^(n/2) [2.99].
- **(A3) Minor arcs.** If 2^(eps n) < q <= 2^((1-eps) n) with 0 < eps <= 1/2, then
  |S_n| <= 4 sqrt(n) 2^(n-1) 2^(-eps n/2) [2.65 sqrt(n)].
  - The three terms of the bound total at most (31.9 + 11.1 n) 2^(-eps n), which is at most 16 n 2^(-eps n).
  - Suppose no r/q with q <= 2^(eps n) has |2 theta - r/q| <= 1/(q 2^((1-eps) n)). Then Dirichlet supplies such a q.
- The best saving the method gives is 4 sqrt(n) 2^(-n/4), at q = 2^(n/2): |X_n|^(-1/4) up to sqrt(n).

### 1'. (A') The three-block version (the session's theorem; superseded)

The session first proved a three-block bound from (2). It is correct, but (A) is smaller in 1,099 of the 1,186 exact
instances and doubles the minor-arc exponent.

For n = a + L + l (a >= 2, L, l >= 1) and r/q as above:

    |S_n(theta)| <= 2^(n-1) 2^(-L/2) ( 1 + 2 (2^L/q + 1) sqrt((3 + 2^(3-a) q)(3 + 2^(2-l) q)) )^(1/2).

*Proof.*
1. Fix p = p(u) and eps = eps(w).
   - A(u) runs over 2^(a-1) + p + 2j.
   - C(w) runs over eps + 2i.
   - H(uvw) = B(u, w) + Hhat(v) + A(u) C(v) + A(v) C(w).
2. Cauchy-Schwarz over (u, w) gives |S_(p,eps)|^2 <= K_u K_w Q.
3. Expanding Q, the u- and w-sums are geometric. Weighted AM-GM bounds the off-diagonal part by 2^L sqrt(T_u T_w), where
   T_K = sum over 0 < |d| < 2^L of min(K, 1/(2 ||2 theta d||))^2.
4. The counting lemma and a layer-cake estimate give T_K <= 2(2^L/q + 1)(3K^2 + 2qK).

**Corollaries.** The recipes below were repaired by the audit (MAJOR-2). The pre-audit recipes did not reach the stated
constants.
- **(A'1)** 6.5 · 2^(n-1) q^(-1/2) for 2 <= q <= 2^((n-3)/3).
  - Take L = floor(log2 q) and a balanced split a = floor or ceil of (n - L + 1)/2.
  - The analytic constant is 5.37, and the optimal split gives 4.23 for n <= 200.
- **(A'2)** 12 · 2^(n-1) (q/2^n)^(1/4) for 2^((n-6)/3) <= q <= 2^(n-6).
  - Take L = floor(log2 q) + 1 and the same balanced a.
  - The analytic constant is 11.30, and the optimal split gives 9.55.
  - L = floor(log2 q) reaches only 12.70, at n = 104.
  - The lower limit on q is necessary. At theta = 1/7 the bound fails for every n >= 37: 3.517e-2 against 3.206e-2.
- **(A'3)** The minor-arc form 12 · 2^(n-1) 2^(-eps n/4) holds for eps < 1/2. It is nontrivial only once eps n > 14.3,
  and the optimal-split supremum of bound · 2^(eps n/4) is 4.8 (n <= 200, every n).

**Checks of (A').**
- Against exact sums, `bound_check.py`: 0 violations at the 66 admissible (n, theta) pairs, 28 of them nontrivial.
  - At n = 24, limit_denominator returned two inadmissible semiconvergents.
  - The audit's 1,186 admissible instances also give 0 violations.
- The counting lemma: 0 violations in 3,000 random cases (max lhs/rhs 0.815).
- Q <= 2^L K_u K_w (1 + R): 0 violations in 100 exact cases.
- The actual Q / diagonal is 0.92-1.15 off the major arcs and up to 4.4 at theta = 1/7 (`core_check.py`).
- `const_check.py` scanned only n = 0 mod 3 and missed the n = 2 mod 3 corner where the old recipe fails.

## 2. (B) One-apex readouts (PROVED)

Let y = (y_h, y_e, y_R) be a rooted triple and Phi_y(x) = y_h H(x) - y_e e(x) + y_R (1 + A(x)), with e(x) = [x ends in 0].
Then the sum over x in X_n of e(theta Phi_y(x)) satisfies (A) and (A'), with 2 theta replaced by 2 y_h theta.

*Proof.*
- In the split, A(x) = 2^m A(u) + A(y) (two blocks) or 2^(L+l) A(u) + 2^l A(v) + A(w) (three blocks).
- e(x) depends only on the last block.
- The new terms that depend on the inner variable have modulus 1.
- The cross terms acquire the factor y_h. QED.

These are the one-apex x bank families F2 of THM-4885; for y = (g(u), e(u), 1 + u) they are the slices of F1.

**Full ranges.** For v over [0, V], split the partial top block into at most n cylinders (a fixed prefix with free
suffixes), and apply the two-block argument to each cylinder.

**These are individual bounds.** The S82 audit (M3) and Session 205 (section 3) show that such bounds do not give the
covering inequality, even uniformly in the coefficients and with a growing bank. A triangle bound sums errors
proportional to the lengths raised to 1 - eta. Each bank member contributes only a constant-order mean
representation count. So enlarging the bank does not repair the scaling.

## 3. (C) Fourier uniformity of the one-apex counts (PROVED)

For every real theta outside (1/2)Z, S_n(theta) / 2^(n-1) -> 0. Since every H is odd, S_n(theta + 1/2) = -S_n(theta), so
|S_n| = 2^(n-1) on (1/2)Z.

*Proof.*
- **Irrational theta.**
  - 2 theta has convergents r_k/q_k with q_k -> infinity and |2 theta - r_k/q_k| <= 1/q_k^2.
  - For n >= 12, take the largest q_k <= 2^(n/4) <= 2^((n-3)/3). It tends to infinity with n.
  - Then (A1) gives |S_n(theta)| / 2^(n-1) <= 3.52 q_k^(-1/2) -> 0.
- **Rational theta = r/s outside (1/2)Z.**
  - Session 193 Theorem 3: for uniform strong words W_m (first letter 1, last letter 0), H(W_m) mod s tends to the
    uniform law on the parity-admissible residues.
  - The Fourier coefficient of that law at r/s vanishes off (1/2)Z, for every modulus s:
    - for odd s it is a full character sum;
    - for even s it is a sum over the odd residues, which vanishes unless s = 2.
  - Every x in X_n is x' 1^t with x' strong of length m = n - t >= 2, or x = 1^n, and H(x) = H(x').
  - So S_n(theta) = sum_(m=2..n) 2^(m-2) E[e(theta H(W_m))] + e(theta), and the weights 2^(m-n-1) concentrate on m near
    n.

QED.

So the multiset {H(x) : x in X_n} has no large Fourier coefficients except at 0 and 1/2:
- it is equidistributed modulo every modulus, on the odd residues for even moduli (Session 193);
- it is equidistributed modulo 1 after multiplication by any irrational number (new in the repository; OEIS has no entry
  for g);
- the minor arcs have the explicit power saving (A3).

Numerically, |S_n|/2^(n-1) -> 0 at twelve rationals; the slowest is 1/7, at 0.0756 for n = 26 (audit `a8`).

## 4. (D) Pre-registered tests (HYP-9479, commit 185c2b3efa; `one_apex_fourier_20261106/outsample.out`)

The table gives |S_n|/2^((n-1)/2) at five irrational theta, then E_n/|X_n| and distinct values/|X_n|.

| n | phi | sqrt2-1 | pi-3 | e-2 | sqrt3-1 | E_n / \|X_n\| | distinct / \|X_n\| |
|---|---|---|---|---|---|---|---|
| 26 | 1.723 | 0.612 | 0.413 | 0.730 | 1.700 | 3.2526 | 0.3734 |
| 28 | 1.979 | 1.782 | 1.547 | 0.821 | 1.443 | 3.2473 | 0.3737 |

- **(P1) CONFIRMED.**
  - All ten square-root ratios lie in [0.05, 20].
  - In sample (n = 14..24) they lie in [0.17, 2.46] for even n. The odd n = 23 gives 2.814 (audit m6).
  - The minor-arc sums are consistent with square-root cancellation, far beyond the proved saving.
  - The window [0.05, 20] cannot separate 2^(n/2) from 2^(0.55 n) at n = 28.
- **(P2) CONFIRMED.**
  - The energy ratio is 3.2526 and 3.2473, against 2.93, 3.06, 3.19, 3.23, 3.24, 3.26, 3.26 for n = 12..24.
  - It has been decreasing slowly since n = 22 (audit: 3.2607 at n = 22, 3.2502 at n = 27).
- **(P3) CONFIRMED.**
  - Setup: the splits (9,5,9), (10,5,10) and (9,7,9), each at the four theta and the four boundary-bit pairs.
  - Q / diagonal lies in [0.9971, 1.0686] in all 48 cases.
- **Reproduction.** The audit's independent engine reproduces the n = 26 and n = 28 rows digit for digit, with
  fixed-point phases and an automaton for g. The "4,600 random checks" of the block generator have no receipt in the
  folder; the audit's reproduction supersedes them.
- **Pre-registration timing.** Verified by the audit (section 7): outsample.out was born after the commit, and the
  predictions are unedited.

## 5. (E) What this gives toward density one, and what it does not

**Positive density from linear energy** (elementary, CONDITIONAL).
- Suppose E_n <= C · 2^(n-1) for all n >= n_0. By Cauchy-Schwarz, the words of length n then give at least
  |X_n|^2 / E_n >= 2^(n-1)/C distinct counts, all at most 2^n (since g(v) <= v + 1).
  - For 2^n <= X < 2^(n+1) this is at least X/(4C) counts up to X.
  - That is lower density at least 1/(2C) among the odd numbers, and 1/C along X = 2^n.
- **The forced floor.** The converse involution x' 1^t -> (x')* 1^t preserves H (Session 16's (2); Session 191
  `one_apex/PROOF.md` section 6.2). Its fixed points number at most 2^(n/2) + 1. Hence E_n >= 2 |X_n| - 2^(n/2) - 1
  (audit U4; checked at n = 8, 12, 16, 18).
  - So C >= 2, and this route gives lower density at most 1/4.
  - "Linear energy" means bounded excess over the floor 2. The observed excess is 1.25 (C = 3.25 gives 0.154).
- **Comparison.** The one-apex values alone have density about 0.5 (audit: 0.4974 at 2^24) and at most ln 2
  (Session 16).
- **Linear energy is sufficient, not necessary.** It is the missing lemma on this route.
  - By Parseval, E_n is the integral of |S_n|^2, so it is square-root cancellation in mean square.
  - Pointwise minor-arc bounds do not reach it: with (A3) at eps = 1/3, Parseval gives only E_n <= (16 n + 2) 2^(5n/3 - 2),
    i.e. 2^(n/3)/(16 n + 2) distinct values (audit U5). That is weaker than Session 150's sqrt(X)/4.
  - THM-4887 (section 6) bounds the energy directly in physical space, and gets 2^(4n/3) up to O(n).

**Density one** (S82's programme).
- A route to density one needs a growing bank and the pair correlations of the bank families. The second moment of r_Y
  around its local mean must be o(mean^2).
- The S82 audit (M3) and Session 205 (section 3) showed this first: individual bounds such as (B), uniform or not, do
  not suffice. One needs cancellation across bank members, integrated estimates, or another exact covering mechanism.
- (A) is the minor-arc half of the circle method. The other half (linear energy, pair correlations) is open.

**Why the bound is weak and how it could improve.**
- The proof bounds the off-diagonal terms in absolute value. The data (P3) show they nearly cancel: Q equals its diagonal
  to within 7% off the major arcs.
- The best saving of the bilinear method is about |X_n|^(-1/4).
- To go further one must use the phases e(theta (H(y) - H(y'))) of the inner block. Each H(y) again obeys (1), so a
  van der Corput-type iteration on the inner block is the natural next step.

## 6. (F) The energy at the scale 2^(n/3) and an X^(2/3) counting bound (PROVED, THM-4887)

(A) controls S_n pointwise, but the quantity positive density needs is a mean square. A physical-space argument bounds
the energy and its smoothed versions down to the scale 2^(n/3), and is sharp there.

**Two-block identity** ((1) with the low block written as an integer). For h >= 1, m >= 1 and 0 <= l < 2^m:

    g(h 2^m + l) = g(h) + h C_m(l) + gamma_h(l),

- C_m(l) = sum_(i<m) (1 - l_i) 2^(m-1-i) is the complement of l reversed, a bijection of [0, 2^m).
- gamma_h(l) = g(l) - 1 + (h mod 2)[l < 2^(m-1)] lies in [0, 2^m].

So a length-n count is a product h·c of the top block h and the reversed complement c of the low block, shifted by
g(h) and by a perturbation of size at most C = 2^m.

**Theorem F (THM-4887).** Let E_n(Delta) = #{(v, v') in X_n^2 : |g(v) - g(v')| <= Delta}.
- **Smoothed energy.** For n >= 12 and 0 <= Delta <= 2^floor(n/3), E_n(Delta) <= (1.39 n + 19) 2^(4n/3).
  - This is sharp up to O(n): any 2^(n-1) values in [1, 2^n] have E_n(Delta) >= 2^(2n-2)(Delta+1)/(2^n+Delta+1).
- **Energy.** E_n <= (1.39 n + 19) 2^(4n/3). The counts on X_n take at least 2^(2n/3)/(5.56 n + 76) values.
- **Counting.** The number A(X) of Hamiltonian path counts up to X satisfies A(X) >= X^(2/3)/(13 ln X + 121) for
  X >= 4096.
- **Fourier form.** The integral of |S_n|^2 over |theta| <= 2^(-floor(n/3)-2) is at most 3.92 (1.39 n + 19) 2^n. That is
  square-root cancellation in mean square on that window, up to O(n). (A) cannot see this.

*Proof sketch* (full proof in THM-4887).
- **Same top block.** Within one top block h in [P, 2P) (P = 2^(n-m-1) >= 8C), the C values are separated by at least
  P - C > C.
- **Different top blocks.** For h != h', a near-collision forces |h c - h' c' - D| <= 2C with D = g(h') - g(h).
  - For each c there is at most one c'.
  - Let beta be the reduced Dirichlet denominator of h/h' with parameter P/C.
  - The counting lemma over blocks of beta consecutive c bounds the strip count by 4C^2/P + 3C/beta + 7, for every
    shift D.
- **Small denominators are rare.** At most P(beta + 1)(4C/beta + 1) pairs have reduced denominator beta.
- **Sum.** E_n(Delta) <= PC + 10 P^2 + P C^2 (35.74 + 12 ln(P/C)) + 3 P C (1 + ln(P/C)). Take m = floor(n/3).

**Checks** (`d_identity.py`, `d_lemma.py`, `d_energy.py`).
- **Identity:** 200,000 random cases.
- **Strip lemma:** checked with the worst-case shift on all ordered pairs, for four (n, m) and both windows (C and 2C):
  0 violations, max ratio 0.635.
- **Small-denominator counts and harmonic sum:** within their bounds.
- **Energies:** E_n(2^floor(n/3)) / 2^(4n/3) lies in [0.39, 0.67] for n = 12..26, against the proved 36 to 55.

**What it means.**
- **The exponent.** The previous all-X exponent was 1/2 (Session 150's injective family, A(X) > sqrt(X)/4).
  - For any X one can write down, the certified range [1, 2^32] (THM-4875 (D)) and order-join multiplicativity
    (THM-1862) give far more: more than 10^213 counts below 2^1000, against X^(2/3) = 10^200.7.
  - Asymptotically, though, finitely many certified values generate only polylogarithmically many products, so the
    exponent 2/3 is the content.
- **Why 2/3.**
  - The proof knows only that gamma lies in [0, C], so a perturbation of size C forces the window up to C.
  - The strip count then carries an area term C^2/P and a lattice term O(1) per pair. Balancing P^2 against P C^2
    places C at 2^(n/3).
  - Below that scale one needs arithmetic information on g: where D = g(h') - g(h) falls relative to the lattice
    {h c - h' c'}. That is the same input linear energy needs.
- **Positive density** is now reduced to the scales between 1 and 2^(n/3). Linear energy is E_n(Delta) << (Delta + 1) 2^n
  at Delta = 0; Theorem F proves that inequality at Delta = 2^(n/3), up to O(n).

## 7. Typing

| claim | type |
|---|---|
| (A) the Type II bound, (A1)-(A3) | PROVED (new application of the classical Type II method; the audit's U1) |
| (A') the three-block bound and its corollaries | PROVED (session; recipes repaired by the audit) |
| (B) the same for the one-apex readouts (bank families) | PROVED (new) |
| (C) Fourier uniformity off (1/2)Z | PROVED (irrational theta new in the repository; rational theta by Session 193 Theorem 3) |
| (D) square-root cancellation, energy 3.25 \|X_n\|, diagonal dominance | FINITE-EXACT / numerical (pre-registered HYP-9479, confirmed) |
| E_n >= 2 \|X_n\| - 2^(n/2) - 1 | PROVED (the audit's U4, from Session 16's involution) |
| linear energy implies positive density (at most 1/4) | elementary (conditional) |
| (F) energy at scale 2^(n/3), A(X) >= X^(2/3)/(13 ln X + 121) | PROVED (new; THM-4887) |
| linear energy; pair correlations; density one; Chan-Pak | OPEN |

## 8. Audits

**THM-4886 and HYP-9479:** `04-computation/experiments/one_apex_fourier_20261106_audit_S83/AUDIT.md`. Verdict PASS WITH
CORRECTIONS: 3 major, 11 minor, 5 upgrades, all applied (MISTAKE-726).

- **Confirmed.**
  - Theorem A' was re-derived and every inequality of its proof checked on exact data: 3,892 cases.
  - There are 0 violations against exact sums at every admissible q (1,186 instances).
  - (B) and (C) are correct; (C) covers moduli with powers of 2.
  - The conditional positive-density claim is correct.
  - n = 26 and 28 were reproduced digit for digit by an independent engine.
  - The pre-registration timing is valid.
- **MAJOR-1.** The pre-audit THM-4886 title dropped the lower limit q >= 2^((n-6)/3) of the (q/2^n)^(1/4) corollary. That
  bound fails at theta = 1/7 for every n >= 37. Fixed.
- **MAJOR-2.** The written recipes for (A'1) and (A'2) did not reach 6.5 and 12. The A'1 split did not exist at the
  top of its range, and the A'2 recipe gave 12.70 at n = 104; the "analytic constant about 8.5" was false.
  The statements are true. The audit's recipes (5.37 and 11.30) replace them.
- **MAJOR-3.** Missing credit. Now added:
  - the classical Type II method;
  - the S82 audit's M3 and Session 205's scope correction, for section 5's point;
  - S82 section 7's rational rates;
  - Session 193 Theorem 3 as a proved fixed-theta decay;
  - Session 16's involution.
- **Minors applied.**
  - m1: the admissible-pair count of `bound_check.py`.
  - m2: the n mod 3 scan.
  - m3: the +2 counting lemma.
  - m4: cylinders for full ranges.
  - m5: "sufficient condition"; "a route to density one".
  - m6: odd n in sample; "consistent with".
  - m7: the missing receipt.
  - m8: (A'3) holds for eps < 1/2 and is nontrivial once eps n > 14.3.
  - m9: n >= 12; every modulus; "new in the repository".
  - m10: at most 2^n.
  - m11: status line.
- **Upgrades adopted.**
  - U1, the two-block Type II bound, is now Theorem A.
  - U2: the +2 counting lemma.
  - U3: the repaired recipes.
  - U4: the energy floor.
  - U5: Parseval calibration.

**THM-4887:** `04-computation/experiments/one_apex_energy_20261106_audit_S83/AUDIT.md` (pending at this revision).

## 9. Reproduction

From `04-computation/experiments/one_apex_fourier_20261106/`:

    python3 concat_check.py      # Session 193's concatenation laws on random triples
    python3 bound_check.py       # Theorem A' against exact sums, n = 12..24
    python3 core_check.py        # the counting lemma and Q <= diagonal (1 + R)
    python3 const_check.py       # constants of (A'1)-(A'3) (n = 0 mod 3 only; see the audit's a3, a9)
    python3 insample.py          # square-root ratios and energy, n <= 24
    python3 outsample.py         # HYP-9479 P1-P3 (n = 26, 28; three splits)
    python3 d_identity.py        # THM-4887: the two-block identity
    python3 d_lemma.py           # THM-4887: strip lemma, small denominators, energies (n = 10..13)
    python3 d_energy.py          # THM-4887: E_n(0), E_n(2^floor(n/3)) against the bound, n = 12..26

The audits' independent scripts are in their folders (section 8). The Type II bound (A) is checked by the audit's
`a5_chain.py`, `a7_bounds_vs_exact.py` and `a9_recipes_and_U1.py`.
