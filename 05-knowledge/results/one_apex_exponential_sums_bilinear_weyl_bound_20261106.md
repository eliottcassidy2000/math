# One-apex exponential sums: a bilinear Weyl bound, Fourier uniformity, and the energy target for positive density

**Session:** opus-2026-11-06-S83. Worktree `C:/tmp/wt-s83`. Canon: THM-4886. Pre-registration: HYP-9479 (commit
185c2b3efa).

**Owner prompt (verbatim).** "start another session working toward a density-one theorem for the spectrum"

**Target.** Chan-Pak (TCS 1015 (2024) 114776) asks that the Hamiltonian path counts H(T) of tournaments cover every odd
number except 7 and 21. That is OPEN. A density-one theorem (almost every odd number is a count) is a weaker milestone.
S82 (THM-4885, section 7) identified the analytic inputs a circle-method proof would need:
- **minor-arc bounds** for the exponential sums of the one-apex counts, known only numerically;
- **pair correlations** (energy) of the join families.

**Outcome.**
- **(A) The first proved minor-arc bound for these families.** Session 193's exact concatenation law turns the middle
  block of a one-apex word into a bilinear form. Cauchy-Schwarz and geometric sums then give
  |sum_x e(theta H(x))| <= 2^(n-1) C (q^(-1/2) + (q/2^n)^(1/4)) whenever |2 theta - r/q| <= 1/q^2. So there is a
  power saving 2^(-eps n / 4) on the minor arcs.
- **(B)** The same bound holds for every one-apex readout y_h H - y_e e + y_R (1 + A) of the join families.
- **(C) Fourier uniformity.** sum_x e(theta H(x)) = o(2^n) for every real theta outside (1/2)Z.
- **(D) Pre-registered and CONFIRMED out of sample** (HYP-9479):
  - square-root cancellation at irrational theta;
  - the energy E_n (pairs of words with equal count) stays at 3.25 |X_n|;
  - the bilinear second moment equals its diagonal to within 7%.
- **(E) What remains.**
  - For POSITIVE density: E_n = O(2^n), i.e. linear energy (observed: 3.25). It gives positive lower density of the
    spectrum by Cauchy-Schwarz.
  - For density ONE: pair correlations of the bank families. Neither follows from (A); both are pair-correlation
    statements.

## 0. Objects and inherited facts

**Words and counts** (Session 193 `modular_sync/PROOF.md` section 2).
- For a binary word w, T(w) is the increasing chain plus an apex beating the i-th chain vertex iff w_i = 1. H(w) is its
  Hamiltonian path count; H(empty) = 1.
- A(w) is the binary integer of w (most significant letter first), and C(w) = sum_j (1 - w_j) 2^(j-1) (the first letter
  has weight 1).
- p(w) is the last bit, and eps(w) = [w nonempty and begins in 0].
- H(w) = g(A(w)), with g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2).

**Concatenation law** (Session 193 eq. (2), checked here on 20,000 random triples). For |v| = L:

    H(uvw) = H(u) + H(v) + H(w) - 2 + A(u) C(v) + (2^L A(u) + A(v)) C(w) + p(u) eps(v) + p(v) eps(w).

**The sum.** X_n = words of length n beginning with 1 (so A runs over [2^(n-1), 2^n)), and
S_n(theta) = sum_(x in X_n) e(theta H(x)).

**Prior work on the Fourier side.**
- The size/residue coupling and the exponential mixing of H mod M for every fixed M are Session 193's (Theorems 3-4).
- The binary-product law H ~ A R / 2^m is Session 16/26's.
- For the ribbon family, mac-mini measured minor-arc decay numerically and conjectured an integer-scale local limit law
  (HYP-9474; THM-4881 (D); open).
- The "at best almost all" limitation of circle methods is kind-pasteur S11's.
- The repository had no proved exponential-sum bound for any of these families (`git grep` for exponential sum, Fourier,
  Weyl, minor arc over sessions 190-205 and the S16 packets).

## 1. (A) The bilinear Weyl bound (PROVED)

**Theorem A.** Let n >= 4, and n = a + L + l with a >= 2, L >= 1, l >= 1. Let theta be real and r/q a fraction in
lowest terms with |2 theta - r/q| <= 1/q^2. Then

    |S_n(theta)| <= 2^(n-1) 2^(-L/2) ( 1 + 2 (2^L/q + 1) sqrt((3 + 2^(3-a) q)(3 + 2^(2-l) q)) )^(1/2).

*Proof.*
1. **Split.** Write x = u v w with |u| = a (u begins with 1), |v| = L, |w| = l. This is a bijection with X_n.
   - Fix p = p(u) and eps = eps(w). Then A(u) runs over the progression 2^(a-1) + p + 2j (0 <= j < K_u = 2^(a-2)).
   - C(w) runs over eps + 2i (0 <= i < K_w = 2^(l-1)).
   - By the concatenation law, H(uvw) = B(u, w) + Hhat(v) + A(u) C(v) + A(v) C(w), with
     B = H(u) + H(w) - 2 + 2^L A(u) C(w) and Hhat(v) = H(v) + p eps(v) + p(v) eps.
2. **Cauchy-Schwarz over (u, w).** S_n is the sum of the four pieces S_(p,eps), and
   |S_(p,eps)|^2 <= K_u K_w Q, where

       Q = sum_(u,w) |F(u, w)|^2,   F(u, w) = sum_v e(theta (Hhat(v) + A(u) C(v) + A(v) C(w))).

3. **Expand.** Q = sum over pairs v, v' of e(theta (Hhat(v) - Hhat(v'))) [sum_u e(theta A(u) dC)] [sum_w e(theta C(w) dA)],
   with dC = C(v) - C(v') and dA = A(v) - A(v').
   - The progressions give |sum_u| <= min(K_u, 1/(2 ||2 theta dC||)) and |sum_w| <= min(K_w, 1/(2 ||2 theta dA||)).
   - The diagonal v = v' contributes 2^L K_u K_w.
   - Off the diagonal, dA != 0 and dC != 0, because A and C are bijections onto [0, 2^L). Each nonzero difference d
     occurs 2^L - |d| times.
   - Using 2xy <= lambda x^2 + y^2/lambda and optimising lambda, the off-diagonal part is at most 2^L sqrt(T_u T_w),
     where T_K = sum over 0 < |d| < 2^L of min(K, 1/(2 ||2 theta d||))^2.
4. **Counting lemma.** If |beta - r/q| <= 1/q^2, then in any block of q consecutive integers d,
   #{d : ||beta d|| <= t} <= 2qt + 3. (beta d = beta d_0 + r m/q + m delta with |m delta| < 1/q, and r m/q runs over all
   classes mod q.)
   - Layer-cake integration gives sum_block min(K^2, 1/(4 ||beta d||^2)) <= 3K^2 + 2qK.
   - At most 2(D/q + 1) blocks cover 0 < |d| < D. Hence T_K <= 2(2^L/q + 1)(3K^2 + 2qK).
5. **Assemble.**
   - Q <= 2^L K_u K_w (1 + R), with R = 2 (2^L/q + 1) sqrt((3 + 2q/K_u)(3 + 2q/K_w)).
   - So |S_n| <= 4 K_u K_w 2^(L/2) sqrt(1 + R), and 4 K_u K_w 2^(L/2) = 2^(n-1) 2^(-L/2).

QED.

**Checks.**
- Theorem A against exact sums: 0 violations at 68 (n, theta) pairs, n = 12..24, rational and irrational theta.
- The counting lemma: 0 violations in 3,000 random cases (max lhs/rhs 0.815).
- The inequality Q <= 2^L K_u K_w (1 + R): 0 violations in 100 exact cases. The actual Q / diagonal is 0.92-1.15 off the
  major arcs and up to 4.4 at theta = 1/7 (`core_check.py`). The proof's off-diagonal bound is therefore very loose.

**Corollaries** (absolute constants; the analytic constants are confirmed by evaluating the exact bound,
`const_check.py`).
- **(A1)** If 2 <= q <= 2^((n-3)/3), take L = floor(log2 q) and a, l >= log2 q + 2. Then |S_n(theta)| <= 6.5 · 2^(n-1) q^(-1/2).
  (The numerical maximum of the exact bound times sqrt q is 4.2.)
- **(A2)** If 2^((n-6)/3) <= q <= 2^(n-6), take L = floor(log2 q) and a, l about (n - L)/2. Then
  |S_n(theta)| <= 12 · 2^(n-1) (q/2^n)^(1/4). (The analytic constant is about 8.5 up to rounding; the numerical
  maximum is 9.5.)
- **(A3) Minor arcs.** Fix 0 < eps <= 1/3. Suppose no fraction r/q with q <= 2^(eps n) satisfies
  |2 theta - r/q| <= 1/(q 2^((1-eps) n)). Then |S_n(theta)| <= 12 · 2^(n-1) 2^(-eps n/4).
  - Proof: Dirichlet gives q <= 2^((1-eps) n) with |2 theta - r/q| <= 1/q^2. By hypothesis q > 2^(eps n). Apply (A1) or
    (A2).
  - Evaluating the exact bound over n = 30..120 gives at most 4.6 · 2^(-eps n/4).

## 2. (B) One-apex readouts (PROVED)

For a rooted triple y = (y_h, y_e, y_R) and Phi_y(x) = y_h H(x) - y_e e(x) + y_R (1 + A(x)), with e(x) = [x ends in 0]:
the sum over x in X_n of e(theta Phi_y(x)) satisfies Theorem A and (A1)-(A3) with 2 theta replaced by 2 y_h theta.

*Proof.*
- In the split, A(x) = 2^(L+l) A(u) + 2^l A(v) + A(w), and e(x) depends only on w.
- The only v-dependent addition is the phase theta y_R 2^l A(v), which has modulus 1.
- The cross terms acquire the factor y_h. QED.

These are the one-apex x bank families F2 (and, for y = (g(u), e(u), 1 + u), the slices of F1) of THM-4885.

## 3. (C) Fourier uniformity of the one-apex counts (PROVED)

For every real theta outside (1/2)Z, S_n(theta) / 2^(n-1) -> 0. Since every H is odd, S_n(theta + 1/2) = -S_n(theta), so
|S_n| = 2^(n-1) on (1/2)Z.

*Proof.*
- **Irrational theta.** 2 theta has convergents r_k/q_k with q_k -> infinity and |2 theta - r_k/q_k| <= 1/q_k^2. For each
  n take the largest q_k <= 2^(n/4); it tends to infinity with n. Then (A1) gives |S_n(theta)| / 2^(n-1) <= 6.5 q_k^(-1/2)
  -> 0.
- **Rational theta = r/s outside (1/2)Z.**
  - Session 193 Theorem 3: for uniform strong words W_m (first letter 1, last letter 0), H(W_m) mod s tends to the
    uniform law on the parity-admissible residues. The Fourier coefficient of that law at r/s is 0 unless 2r/s is an
    integer.
  - Every x in X_n is x' 1^t with x' strong of length m = n - t >= 2, or x = 1^n, and H(x) = H(x').
  - So S_n(theta) = sum_(m=2..n) 2^(m-2) E[e(theta H(W_m))] + e(theta), and the weights 2^(m-n-1) concentrate on m near
    n.

QED.

So the multiset {H(x) : x in X_n} has no large Fourier coefficients except at 0 and 1/2:
- it is equidistributed modulo every odd modulus (Session 193);
- it is equidistributed modulo 1 after multiplication by any irrational number (new);
- the minor arcs have the explicit power saving (A3).

## 4. (D) Pre-registered tests (HYP-9479, commit 185c2b3efa; `one_apex_fourier_20261106/outsample.out`)

| n | |S_n|/2^((n-1)/2) at phi, sqrt2-1, pi-3, e-2, sqrt3-1 | E_n / |X_n| | distinct values / |X_n| |
|---|---|---|---|
| 26 | 1.723, 0.612, 0.413, 0.730, 1.700 | 3.2526 | 0.3734 |
| 28 | 1.979, 1.782, 1.547, 0.821, 1.443 | 3.2473 | 0.3737 |

- **(P1) CONFIRMED.** All ten square-root ratios lie in [0.05, 20], and in-sample (n = 14..24) they lie in [0.17, 2.46]:
  the minor-arc sums show square-root cancellation, far beyond the proved saving.
- **(P2) CONFIRMED.** The energy ratio is 3.2526 and 3.2473, against 2.93, 3.06, 3.19, 3.23, 3.24, 3.26, 3.26 for
  n = 12..24.
- **(P3) CONFIRMED.** At the splits (9,5,9), (10,5,10), (9,7,9), the four theta and the four boundary-bit pairs,
  Q / diagonal lies in [0.9971, 1.0686] in all 48 cases.
- The block generator of g used for n = 26, 28 agrees with the recursion in 4,600 random checks.

## 5. (E) What this gives toward density one, and what it does not

**Positive density from linear energy** (elementary, CONDITIONAL).
- Suppose E_n <= C · 2^(n-1) for all n >= n_0. By Cauchy-Schwarz, the words of length n then give at least
  |X_n|^2 / E_n >= 2^(n-1)/C distinct counts, all at most 2^n + 1 (since g(v) <= v + 1).
- So the spectrum has at least X/(4C) elements up to X for large X: positive lower density among the odd numbers,
  at least 1/(2C).
- The observed C = 3.25 would give density at least 0.15. The one-apex values alone have density about 0.5, and at most
  ln 2 (Session 16).
- Linear energy is the precise missing lemma for positive density. By Parseval, E_n is the integral of |S_n(theta)|^2 over
  theta, so it is equivalent to square-root cancellation in mean square. The proved saving (A3) gives only
  E_n <= 2^(2n(1 - eps/4)) + (major arcs).

**Density one** (S82's programme).
- It needs a growing bank and the pair correlations of the bank families: the second moment of r_Y around its local mean
  must be o(mean^2).
- The minor-arc bound (B) controls the cross terms only after multiplication by the families' L^2 norms (their
  energies). So density one still rests on energy-type statements for every family in the bank.
- (A) is the minor-arc half of the circle method; the other half (linear energy, pair correlations) is open.

**Why the bound is weak and how it could improve.**
- The proof bounds the off-diagonal terms in absolute value. The data (P3) show they nearly cancel: Q equals its diagonal
  to within 7% off the major arcs.
- That cancellation comes from the phases e(theta (Hhat(v) - Hhat(v'))). Hhat(v) itself has a bilinear decomposition, so
  iterating the method on the middle block (van der Corput style) is the natural next step toward the observed square-root
  cancellation. If it reached L^2 square-root cancellation, it would prove linear energy, hence positive density.

## 6. Typing

| claim | type |
|---|---|
| (A) the bilinear Weyl bound; (A1)-(A3) | PROVED (new; from Session 193's concatenation law) |
| (B) the same for the one-apex readouts (bank families) | PROVED (new) |
| (C) Fourier uniformity off (1/2)Z | PROVED (irrational theta new; rational theta by Session 193 Theorem 3) |
| (D) square-root cancellation, energy 3.25 |X_n|, diagonal dominance | FINITE-EXACT / numerical (pre-registered HYP-9479, confirmed) |
| linear energy implies positive density | elementary (conditional) |
| linear energy; pair correlations; density one; Chan-Pak | OPEN |

## 7. Audit

(Filled after the independent audit.)

## 8. Reproduction

From `04-computation/experiments/one_apex_fourier_20261106/`:

    python3 concat_check.py      # Session 193's concatenation laws on random triples
    python3 bound_check.py       # Theorem A against exact sums, n = 12..24
    python3 core_check.py        # the counting lemma and Q <= diagonal (1 + R)
    python3 const_check.py       # constants of (A1)-(A3)
    python3 insample.py          # square-root ratios and energy, n <= 24
    python3 outsample.py         # HYP-9479 P1-P3 (n = 26, 28; three splits)
