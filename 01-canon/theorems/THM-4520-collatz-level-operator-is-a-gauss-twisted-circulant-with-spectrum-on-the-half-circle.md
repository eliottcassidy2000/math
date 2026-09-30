---
id: THM-4520
title: "The level-n Syracuse frequency recursion is a weighted circulant digraph on the cyclic unit group Z/L_n (L_n = 2 3^(n-1), jumps -a with weights 2^-a) composed with the Gauss-sum diagonal omega_n(k) = e(2^k/3^n); because {2^k} is the full set of units, the characteristic polynomial of T_n = (S/2)(I - S/2)^-1 diag(omega_n) is (2^L - 1) lambda^L + lambda^(L/2) - 1, so every eigenvalue has modulus 1/2 -+ 2^(-L/2)/(2L) (two circles), while the singular values are |1/(2 e(xi/L) - 1)| filling [1/3, 1] on the circle |w - 1/3| = 2/3; omega_n has two-level autocorrelation (L, -L/2 at +-L/3, 0 elsewhere) and the mean square of the one-step Gauss sums over the cycle is 1/3: the Parseval rate 3^(-1/2) of the Syracuse law is a non-normal transient of operators of spectral radius 1/2"
status: >
  PROVED (elementary: the circulant/diagonal factorisation, the eigenvalue
  equation through the cyclotomic polynomial Phi_(3^n)(x) = x^L + x^(L/2) + 1, the
  explicit eigenvectors, the singular values, the autocorrelation and the
  mean-square Gauss-sum identity via Ramanujan sums) + VERIFIED (numerically for
  n = 2..7: eigenvalue moduli [0.490904, 0.511945] = |mu_-+|^(1/3) at n = 2,
  [0.499946, 0.500054] at n = 3, 0.500000 for n >= 4; characteristic-polynomial
  residuals 1e-16..1e-155 (scaled); eigenvector identity 1e-12..1e-16; autocorrelation
  exact to 1e-13; mean |G_n(2^k)|^2 = 0.333333333 for n >= 5 with A = 40) +
  FINITE-EXACT (transients). NOT independently audited. Setting: L = L_n =
  2 3^(n-1), S the cyclic shift (Sv)(k) = v(k-1) on C^L, D = diag(omega_n),
  omega_n(k) = e(2^k/3^n) (all primitive 3^n-th roots of unity, once each);
  the untruncated recursion f_n = T_n lift f_(n-1) with T_n = sum_(a>=1) 2^-a S^a D
  = (S/2)(I - S/2)^-1 D (S20/S22/S23 closed family; THM-4519).
  (1) Eigenvalues: T_n v = lambda v iff (2 lambda)^L = prod_k (lambda + omega_k) =
  Phi_(3^n)(-lambda) = lambda^L - lambda^(L/2) + 1, i.e. (2^L - 1) lambda^L +
  lambda^(L/2) - 1 = 0: lambda^(L/2) = mu_-+ with (2^L - 1) mu^2 + mu - 1 = 0, so
  the L eigenvalues are lambda = mu_-+^(2/L) e(2j/L), j < L/2, of moduli
  |mu_-+|^(2/L) = 1/2 (1 -+ 2^(-L/2)/L + O(2^-L)); the spectral radius is 1/2 +
  2^(-L/2)/(2L) + O(2^-L/L) and the geometric mean of the moduli is (2^L - 1)^(-1/L).
  Eigenvectors: with w(k-1) = 2 lambda/(omega_k + lambda) w(k), v = S w/(2 lambda).
  (2) Singular values: sigma(T_n) = sigma(C) = |Ghat(xi)| = 1/|2 e(xi/L) - 1|,
  xi in Z/L, ranging over [1/3, 1]; the eigenvalues of the circulant C alone are
  Ghat(xi) = 1/(2 e(xi/L) - 1), on the circle |w - 1/3| = 2/3 (through 1 and -1/3).
  (3) omega_n has Fourier transform tau(bar psi_xi) (Gauss sums, |.| = 3^(n/2) on
  primitive xi, 0 on 3 | xi, n >= 2), hence autocorrelation R(0) = L, R(+-L/3) = -L/2,
  R(d) = 0 otherwise: a two-level near-perfect sequence on the cyclic group.
  (4) (1/L) sum_k |G_n(2^k)|^2 = sum_a 4^-a (1 + O(2^(-L/3))) = 1/3 + O(2^(-L/3)) for
  the untruncated kernel (Ramanujan sums), and exactly sum_(a<=A) 4^-a when
  A < 2 3^(n-2): the rms one-step Gauss sum over the units is the Parseval rate 3^(-1/2);
  equivalently ||T_n 1||^2/||1||^2 = 1/3 + O(2^(-L/3)).
  (5) Transient (FINITE-EXACT): ||T_n^m 1||/||T_n^(m-1) 1|| = 0.577, 0.576, 0.573,
  0.568, 0.571, 0.566, 0.563, 0.577, 0.581, 0.566, 0.562 at n = 6 (m = 2..12) and
  ||T_n^8||^(1/8) = 0.61, 0.66, 0.70, 0.72 at n = 3, 4, 5, 6: the contraction at the
  Parseval rate persists for many iterations of one level operator although its
  spectral radius is 1/2 (non-normal, sigma_max = 1). The true cocycle f_n = T_n lift
  f_(n-1) reproduces the S20 profile exactly (max |f_n| = M(n): 0.37792, 0.25224,
  0.17700, 0.12927, 0.09611, 0.07587 at n = 2..7) with rms |f_n| 3^(n/2) = 0.836 constant.
  Reading: circulants enter Collatz exactly where a cyclic symmetry is exact -- the
  fair split of a cycle necklace (THM-4515, Z/j) and the unit group (Z/L_n) here; the
  Syracuse step is a weighted jump on the circulant digraph Cay(Z/L_n, {-1,-2,...}) with
  a Gauss-sum gauge field, and H1 (HYP-9166) concerns the product of the DIFFERENT
  operators T_n (a cocycle), to which the spectral radius 1/2 does not apply.
  NOT claimed: any bound on the cocycle; H1; a spectral gap.
source: collatz-necklace-20260929 session (mac-mini), part 3 (2026-09-30), owner directive "circulant graphs and their connection to Collatz"
depends_on:
  - 01-canon/theorems/THM-4519-frequency-one-coefficient-carries-the-renewal-series.md (the Gauss-sum multiplier structure; the recursion)
related:
  - 05-knowledge/results/collatz_three_mirrors_20260929.md (S23 Theorem 1: the character transform of the closed family is tau(bar psi) E[psi(Y_n)]; Proposition 2: the full-period one-pole recursion)
  - 05-knowledge/results/collatz_five_mirrors_20260929.md (section 2d, Lemma G, the one-step Gauss sums; section 3, the Dalfo-Fiol-Reyes circulant mirror)
  - 01-canon/theorems/THM-4515-fair-consecutive-splits-of-cycle-necklaces-are-circulant.md (the other circulant in Collatz)
  - 05-knowledge/results/collatz_circulant_20260930_circulants_lucas_cubic_monotile.md (the session note)
note: 05-knowledge/results/collatz_circulant_20260930_circulants_lucas_cubic_monotile.md
scripts:
  - 04-computation/experiments/collatz_circulant_20260930_spectrum_exact.py
  - 04-computation/experiments/collatz_circulant_20260930_circulant_spectrum.py
  - 04-computation/experiments/collatz_circulant_20260930_transient.py
outputs: 05-knowledge/results/collatz_circulant_20260930_{spectrum_exact, circulant_spectrum_n7, transient}.out
script_sha256: a10099fd4d7e6647fd515f92673e86fd8de9553035af6c0d5f8464ac9fff3e30 (spectrum_exact), f1735577bcd25b47ee0575b2617b1d7a3fdb46ad31ab049535da6822b4f0037e (circulant_spectrum), de0d8964c6eee28d1d32eb6119cd69bc09bada77ed6bcc21451b8aea636ca77c (transient)
hash_basis: raw LF bytes
audit: NOT independently audited; the proof is a five-line computation once the eigenvalue equation is written, and every identity is checked numerically at n <= 7.
---

# THM-4520 -- the Collatz level operator is a Gauss-twisted circulant with spectrum on the half-circle

**PROVED + VERIFIED; not independently audited.** Session note: [collatz_circulant_20260930_circulants_lucas_cubic_monotile](../../05-knowledge/results/collatz_circulant_20260930_circulants_lucas_cubic_monotile.md), section 1.

## 1. The operator

Let `L = L_n = 2 3^{n-1}`, let `S` be the cyclic shift `(Sv)(k) = v(k-1)` on `C^L`, and `D = diag(omega_n)` with `omega_n(k) = e(2^k/3^n)`. Since `2` generates `(Z/3^n)^x`, `{omega_n(k)}` is the set of all primitive `3^n`-th roots of unity, each once. The frequency recursion of the closed family (S20; THM-4519) is `f_n = T_n lift f_{n-1}` with

`T_n = sum_{a >= 1} 2^{-a} S^a D = (S/2)(I - S/2)^{-1} D`

(the Neumann series converges, `||S/2|| = 1/2`). `T_n` is the adjacency operator of the weighted circulant digraph `Cay(Z/L, {-1, -2, ...})` with weights `2^{-a}`, composed with the unitary "gauge" `D`.

## 2. Spectrum

**Theorem (eigenvalues).** `lambda` is an eigenvalue of `T_n` iff `(2 lambda)^L = prod_k (lambda + omega_n(k))`. Since `prod_{u unit} (x - e(u/3^n)) = Phi_{3^n}(x) = x^L + x^{L/2} + 1`, the right side is `Phi_{3^n}(-lambda) = lambda^L - lambda^{L/2} + 1` (`L/2 = 3^{n-1}` odd, `L` even), so the characteristic equation is

`(2^L - 1) lambda^L + lambda^{L/2} - 1 = 0`, i.e. `lambda^{L/2} = mu_±`, `(2^L - 1) mu^2 + mu - 1 = 0`.

*Proof.* Put `w = (I - S/2)^{-1} D v`; then `T_n v = (S/2) w = lambda v` gives `v = S w/(2 lambda)` and `D S w/(2 lambda) = (I - S/2) w`, i.e. `(D/lambda + I)(S w)/2 = w`, i.e. `w(k-1) = 2 lambda/(omega_n(k) + lambda) w(k)`. Around the cycle the product of the ratios is `1`: `prod_k 2 lambda/(omega_n(k) + lambda) = 1`. Conversely each root `lambda` (none equals `-omega_n(k)`, since the roots have modulus near `1/2`, and none is `0`) gives the eigenvector `w`, `v`. The polynomial has degree `L`, so these are all the eigenvalues. ∎

Hence the eigenvalues are `mu_±^{2/L} e(2j/L)`, `0 <= j < L/2`, on two circles of radii `|mu_±|^{2/L} = (1/2)(1 ∓ 2^{-L/2}/L + O(2^{-L}))`: **every eigenvalue has modulus `1/2` up to `2^{-L/2}/(2L)`**, the spectral radius is `1/2 + 2^{-L/2}/(2L) + O(2^{-L}/L)`, and the geometric mean of the moduli is `(2^L - 1)^{-1/L}` (Jensen's formula for `|det C| = prod |Ghat(xi)|` gives the same). Verified: moduli `[0.490904, 0.511945]` at `n = 2` (`|mu_∓|^{1/3}` with `mu = (-1 ± sqrt 253)/126`), `[0.499946, 0.500054]` at `n = 3`, `0.500000` at `n = 4..7`; characteristic-polynomial residuals `10^{-16} .. 10^{-155}`; eigenvector identity to `10^{-12}`.

**Singular values.** `D` is unitary, so the singular values of `T_n` are those of the circulant `C = (S/2)(I - S/2)^{-1}`: `|Ghat(xi)| = 1/|2 e(xi/L) - 1|`, `xi in Z/L`, filling `[1/3, 1]`; the eigenvalues of `C` are `Ghat(xi) = 1/(2 e(xi/L) - 1)`, which lie on the circle `|w - 1/3| = 2/3` (through `1` at `xi = 0` and `-1/3` at `xi = L/2`). So `T_n` is strongly non-normal: `sigma_max = 1`, `sigma_min = 1/3`, spectral radius `1/2`.

## 3. The gauge field is a near-perfect sequence

The Fourier transform of `omega_n` on `Z/L` is the Gauss sum `tau(bar psi_xi)` of the multiplicative character `psi_xi(2^k) = e(k xi/L)` (S23, Theorem 1 of the three-mirrors note; THM-4519): `|hat omega_n(xi)| = 3^{n/2}` for `3 ∤ xi`, `0` for `3 | xi` (`n >= 2`). Hence the periodic autocorrelation `R(d) = sum_k omega_n(k) bar omega_n(k+d)` is `L` at `d = 0`, `-L/2` at `d = ± L/3` (the shift by `L/3` multiplies `2^k` by an element of order three, `2^{L/3} = 1 + 3^{n-1} c`, so `omega_n(k + L/3) = omega_n(k) e(2^k c/3)`), and `0` at every other `d`: a two-level near-perfect sequence on the cyclic group, the exponential analogue of the Legendre/Paley sequences behind Paley's Hadamard matrices. Verified exactly (`10^{-13}`) for `n <= 7`.

## 4. The Parseval rate as a transient

`(1/L) sum_k |G_n(2^k)|^2 = 1/3 + O(2^{-L/3})` for the one-step Gauss sums `G_n(2^k) = (T_n 1)(k)` (Ramanujan sums: `sum_u |G_n(u)|^2 = sum_{a,a'} 2^{-a-a'} c_{3^n}(2^{-a} - 2^{-a'})`, and `c_{3^n}(t)` is `L`, `-3^{n-1}`, `0` according as `3^n | t`, `3^{n-1} || t`, otherwise; exactly `sum_{a <= A} 4^{-a}` for the kernel truncated at `A < 2 3^{n-2}`). So the first step from the constant vector contracts the energy by exactly `1/3` — the Parseval rate `3^{-1/2}` — and the rms of the one-step Gauss sum over the units is `3^{-1/2}` (S22 measured the mean modulus `0.543` and the sup `-> 1`). Numerically the per-step ratios `||T_n^m 1||/||T_n^{m-1} 1||` stay at `0.56–0.58` for `m <= 12` at `n = 4..6` and `||T_n^8||^{1/8} = 0.61, 0.66, 0.70, 0.72` at `n = 3..6`: the contraction at the Parseval rate is a non-normal transient that lengthens with `L`, while the spectral radius is `1/2`. The actual Syracuse cocycle `f_n = T_n lift f_{n-1}` reproduces the S20 profile (`max |f_n| = M(n)` to five digits at `n <= 7`) with `rms |f_n| 3^{n/2} = 0.836` constant — every level is a fresh transient of a new operator, and H1 (HYP-9166) is a statement about this cocycle, not about any spectral radius.

## 5. Boundary

Nothing is claimed about the cocycle beyond the exact first-step identity; no spectral gap for the products exists in this formulation. The two circulants of Collatz (THM-4515 on `Z/j` for fair necklace splits; this one on `Z/L_n`) are the two places where a cyclic symmetry is exact.
