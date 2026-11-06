---
id: THM-4886
title: "A Type II (bilinear) bound for one-apex exponential sums. (A) If |2 theta - r/q| <= 1/q^2 with r/q in lowest terms, then for every split n = a + m (2 <= a <= n - 1), |sum over words x of length n beginning with 1 of e(theta H(x))| <= 2^(n-1) (6/q + 7 · 2^(a-n) + 8 (1 + ln q) 2^(-a) + 8 q (1 + ln q) 2^(-n))^(1/2): the classical Weyl-Vinogradov Type II estimate applied to Session 193's concatenation law H(uy) = H(u) + H(y) - 1 + A(u) C(y) + p(u) eps(y) (found by the independent audit). Hence, for n >= 7, 3.52 · 2^(n-1) q^(-1/2) for q <= 2^((n-3)/3), 4 · 2^(n-1) (q n / 2^n)^(1/2) for q >= 2^(n/2), and 4 sqrt(n) 2^(n-1) 2^(-eps n/2) when 2^(eps n) < q <= 2^((1-eps) n). (A') The session's three-block version: 2^(n-1) 2^(-L/2) (1 + 2 (2^L/q + 1) sqrt((3 + 2^(3-a) q)(3 + 2^(2-l) q)))^(1/2), with 6.5 · 2^(n-1) q^(-1/2) for q <= 2^((n-3)/3) and 12 · 2^(n-1) (q/2^n)^(1/4) for 2^((n-6)/3) <= q <= 2^(n-6). (B) The same for every one-apex readout y_h H - y_e e + y_R (1 + A). (C) Fourier uniformity: the sum is o(2^n) for every theta outside (1/2)Z. Linear energy E_n = O(2^n) (forced floor 2 |X_n|; observed 3.25 |X_n|) would give positive density; density one and Chan-Pak OPEN"
status: >
  PROVED (new application of a classical method): (A), (A'), (B), (C). The argument is the classical Weyl-Vinogradov
  Type II bilinear estimate (Vaughan, The Hardy-Littlewood Method, 2nd ed., ch. 2; Davenport, Multiplicative Number
  Theory, ch. 25); the new input is Session 193's exact concatenation law (modular_sync/PROOF.md eqs. (1), (2)).
  (A) is the audit's upgrade U1 (two-block split, eq. (1)). It is smaller than (A') in 1,099 of 1,186 exact admissible
  instances and doubles the minor-arc exponent. (A') is the session's three-block version (eq. (2)).
  (C) is new in the repository for irrational theta; for rational theta it is Session 193's Theorem 3 transferred to all
  words.
  Prior work:
  - Session 193 Theorem 3 proves decay at each fixed rational theta (rate 2(1 - T/2^L)^floor((n-2)/L)).
  - S82 (THM-4885 note section 7) computed the exact rational rates (worst 0.9320 per doubling at +-1/7, q <= 120).
  - No uniform-in-theta (minor-arc) bound existed before.
  FINITE-EXACT / numerical, pre-registered as HYP-9479 (commit 185c2b3efa) and CONFIRMED:
  - consistent with square-root cancellation at five irrational theta to n = 28;
  - E_n / |X_n| = 3.25 at n = 26, 28;
  - Q / diagonal in [0.997, 1.069] at unseen splits.
  NOT CLAIMED: linear energy, positive density, density one, Chan-Pak (all OPEN).
  - Individual minor-arc bounds, even with a growing bank, do not give the covering inequality (S82 audit M3;
    Session 205 S82_SCOPE_CORRECTION.md section 3).
  - By Parseval, (A) gives only E_n <= (16 n + 2) 2^(5n/3 - 2) (audit U5); THM-4887 bounds E_n directly.
  AUDIT: PASS WITH CORRECTIONS (3 major, 11 minor, 5 upgrades; MISTAKE-726),
  04-computation/experiments/one_apex_fourier_20261106_audit_S83/AUDIT.md.
session: opus-2026-11-06-S83
source: 05-knowledge/results/one_apex_exponential_sums_bilinear_weyl_bound_20261106.md
related:
  - THM-4887
  - THM-4885
  - THM-4881
  - THM-4875
  - HYP-9479
  - HYP-9474
  - HYP-9478
  - 05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md
  - 05-knowledge/results/feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md
  - 05-knowledge/results/seam_gadget_thickness_spectrum_20261104.md
  - 05-knowledge/results/session205_apex_scope_audit_20261106/S82_SCOPE_CORRECTION.md
  - 04-computation/experiments/one_apex_join_mass_20261106_audit_S82/AUDIT.md
verification:
  - 04-computation/experiments/one_apex_fourier_20261106/concat_check.py
  - 04-computation/experiments/one_apex_fourier_20261106/bound_check.py
  - 04-computation/experiments/one_apex_fourier_20261106/core_check.py
  - 04-computation/experiments/one_apex_fourier_20261106/const_check.py
  - 04-computation/experiments/one_apex_fourier_20261106/insample.py
  - 04-computation/experiments/one_apex_fourier_20261106/outsample.py
  - 04-computation/experiments/one_apex_fourier_20261106/outsample.out
  - 04-computation/experiments/one_apex_fourier_20261106_audit_S83/a5_chain.py
  - 04-computation/experiments/one_apex_fourier_20261106_audit_S83/a7_bounds_vs_exact.py
  - 04-computation/experiments/one_apex_fourier_20261106_audit_S83/a9_recipes_and_U1.py
  - 04-computation/experiments/one_apex_fourier_20261106_audit_S83/a10_involution_floor.py
---

# THM-4886 - a Type II bound and Fourier uniformity for one-apex counts

**Notation** (Session 193).
- X_n is the set of binary words of length n beginning with 1, and H(x) = g(A(x)) is the one-apex Hamiltonian path count.
- A(w) is the binary value of w, and C(w) = sum_j (1 - w_j) 2^(j-1).
- p(w) is the last bit, and eps(w) = [w begins in 0].
- S_n(theta) = sum_(x in X_n) e(theta H(x)).

**Concatenation laws** (Session 193 eqs. (1), (2)).
- For words u, y: H(uy) = H(u) + H(y) - 1 + A(u) C(y) + p(u) eps(y).
- For |v| = L: H(uvw) = H(u) + H(v) + H(w) - 2 + A(u) C(v) + (2^L A(u) + A(v)) C(w) + p(u) eps(v) + p(v) eps(w).

## (A) The Type II bound (found by the independent audit, U1)

For n >= 3, 2 <= a <= n - 1, and r/q in lowest terms with |2 theta - r/q| <= 1/q^2:

    |S_n(theta)| <= 2^(n-1) (6/q + 7 · 2^(a-n) + 8 (1 + ln q) 2^(-a) + 8 q (1 + ln q) 2^(-n))^(1/2).

*Proof.*
1. **Split.** Write x = u y with |u| = a and |y| = m = n - a. Fix p = p(u). Then A(u) runs over 2^(a-1) + p + 2j for
   0 <= j < K = 2^(a-2). By eq. (1), S_p = sum_u e(theta H(u)) sum_y b_y e(theta A(u) C(y)), with |b_y| = 1.
2. **Cauchy-Schwarz over u.** |S_p|^2 <= K sum_(y,y') min(K, 1/(2 ||2 theta (C(y) - C(y'))||)).
   - C is a bijection onto [0, 2^m).
   - So this is at most K 2^m (K + sum over 0 < |d| < 2^m of min(K, 1/(2 ||2 theta d||))).
3. **Counting lemma.** In a block of q consecutive d, at most 2qt + 3 have ||2 theta d|| <= t.
   - So the sorted distances satisfy t_j >= (j - 3)/(2q).
   - A block therefore contributes at most 3K + q(1 + ln q).
   - At most 2(2^m/q + 1) blocks cover 0 < |d| < 2^m.
4. **Assemble.** |S_p|^2 / (K 2^m)^2 <= 2^(-m) + 2(1/q + 2^(-m))(3 + q(1 + ln q)/K). Then |S_n| <= |S_0| + |S_1|,
   with K 2^m = 2^(n-2).

QED.

The counting lemma holds with +2, because the perturbations have one sign (audit U2). This turns 6 and 7 into 4 and 5.

**Corollaries** (n >= 7, a = ceil(n/2); numerical suprema over n <= 300 in brackets, audit `a9`):
- **(A1)** |S_n| <= 3.52 · 2^(n-1) q^(-1/2) for 2 <= q <= 2^((n-3)/3) [3.37].
- **(A2)** |S_n| <= 4 · 2^(n-1) (q n / 2^n)^(1/2) for q >= 2^(n/2) [2.99].
- **(A3) Minor arcs.** If 2^(eps n) < q <= 2^((1-eps) n) with 0 < eps <= 1/2, then
  |S_n| <= 4 sqrt(n) 2^(n-1) 2^(-eps n/2).
  - Suppose no r/q with q <= 2^(eps n) has |2 theta - r/q| <= 1/(q 2^((1-eps) n)). Then Dirichlet supplies such a q.
- The best saving this method gives is about |X_n|^(-1/4), at q = 2^(n/2).

## (A') The three-block version (the session's theorem; superseded by (A))

For n = a + L + l (a >= 2, L, l >= 1) and r/q as above:

    |S_n(theta)| <= 2^(n-1) 2^(-L/2) (1 + 2 (2^L/q + 1) sqrt((3 + 2^(3-a) q)(3 + 2^(2-l) q)))^(1/2).

*Proof.*
- Cauchy-Schwarz over (u, w) in eq. (2).
- Geometric sums over the progressions A(u) and C(w).
- The same counting lemma for the off-diagonal pairs of middle words.

**Corollaries** (recipes repaired by the audit, MAJOR-2):
- 6.5 · 2^(n-1) q^(-1/2) for 2 <= q <= 2^((n-3)/3).
  - Recipe: a balanced split with L = floor(log2 q).
  - The analytic constant is 5.37; the optimal split gives 4.23 for n <= 200.
- 12 · 2^(n-1) (q/2^n)^(1/4) for 2^((n-6)/3) <= q <= 2^(n-6).
  - Recipe: L = floor(log2 q) + 1.
  - The analytic constant is 11.30; the optimal split gives 9.55 for n <= 200.
  - The lower limit on q is necessary: at theta = 1/7 the bound fails for every n >= 37.

## (B) Readouts

The sum of e(theta (y_h H(x) - y_e [x ends in 0] + y_R (1 + A(x)))) over X_n satisfies (A) and (A'), with 2 theta
replaced by 2 y_h theta. These are the one-apex x bank families of THM-4885.
- **Full ranges.** For v over [0, V], split the partial top block into at most n cylinders (a fixed prefix with free
  suffixes), and apply the two-block argument to each.
- **These are individual bounds.** The S82 audit (M3) and Session 205 (S82_SCOPE_CORRECTION.md section 3) show that
  such bounds do not give the covering inequality. This holds even when they are uniform in the coefficients and the
  bank grows.

## (C) Fourier uniformity

S_n(theta) = o(2^n) for every real theta outside (1/2)Z, and S_n(theta + 1/2) = -S_n(theta).
- **Irrational theta:** (A1) along the convergents of 2 theta (n >= 12).
- **Rational theta:** Session 193 Theorem 3 (mixing of H mod M on strong words), transferred through x = x' 1^t.
  The limiting Fourier coefficient vanishes off (1/2)Z for every modulus, including powers of 2.

## Toward density one (not claimed)

- **The energy has a forced floor.** The converse involution x' 1^t -> (x')* 1^t, where (x')* is the reverse
  complement (Session 16), preserves H. So E_n >= 2 |X_n| - 2^(n/2) - 1 (audit U4).
  - Linear energy E_n = O(2^n) therefore means bounded excess over 2. The observed value is 3.25 |X_n| at n <= 28.
  - Linear energy would give positive lower density, at most 1/4 by this route.
- **Linear energy is a sufficient condition**, the missing lemma on this route. THM-4887 proves E_n = O(n 2^(4n/3)).
- **A route to density one** needs, in addition, pair correlations of the bank families.

**Proofs.** Source note, sections 1-3, and section 5 of the audit report.
