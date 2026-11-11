---
id: THM-5128
title: "By-value local densities of the one-arc join family F1. (A) Exact class law of a bilinear gluing: for independent states u ~ D1, v ~ D2 in (Z/p)^r and a nondegenerate bilinear form B, P(B(u,v) = c) = 1/p + (pi1 + pi2 - p^-r)([c = 0] - 1/p) + Rhat(c), pi_i = D_i(0), |Rhat| <= ((p-1)/p) p^(r/2) |D1 - U|_2 |D2 - U|_2: one class is biased and to first order the bias is the mass of the states that vanish modulo p. (C) For F1 (N = x1 x2 + 2^s y1 y2 - [s = 0] over one-apex states (x, y) = (2 psi(h) + 1, 2h + 1), counted by value N <= X) on each sheet: #{N = c mod p} = #T/p + ([c = c_s] - 1/p) D(X) + O(X^(1 - theta_p/2)), D(X) the number of triples whose smaller state is p-divisible (p | gcd(x, y)), c_s = -1 for s = 0 and 0 for s >= 1. (D) #T = kappa_s X ln X + O(X) and D(X) = #T/p^2 + 2 A_p X + O(X^(1 - theta_p/2)): the by-value local factors tend to 1 - 1/p^2 + [c = c_s]/p (an Euler product over the odd primes) with an explicit correction of order 1/ln X; over both sheets the boost 1/p is split (0.43186, 0.56814) between the classes -1 and 0 for triples and (0.19641, 0.80359) for word pairs, not the (1/4, 3/4) of the letter-count model. (E) At the prime 2 the twisted transfer operator of the digit recursion is sqrt(2) times a unitary matrix; F1 has no 2-adic bias; the mod-4 law of the one-apex family is a band form with symbol (1 + z)^2 (1 + z + z^2), of rank periodic with period 6. (F) x = y exactly at h = 2^j (x = y = 2^(j+1) + 1), so the prime factors of 2^m + 1 have early p-divisible states."
status: >
  PRE-AUDIT DRAFT (2026-11-10). (A), (B), (C), (D), (E), (F): PROVED below (elementary; proofs complete in this
  file; (C)-(D) were derived twice, by the session and by an independent blind derivation, and merged).
  The constants A_p and the rates rho_p are VERIFIED / NUMERICAL (floating point). (G): NUMERICAL, with eight
  PRE-REGISTERED tests (279 intervals) on the unseen window [2^24, 2^27), all CONFIRMED (pre-registration commit
  1f97e27805 precedes the output).
  KNOWN before this session, and used: the digit chain on (h mod m, h mod 2), its uniform stationary law and its
  mixing (Session206 root-residue-capacity and weighted-parameter-capacity; Session193 modular_sync Theorems 3-4;
  THM-4886 (C)); the twisted transfer matrices and their rates for the count alone, including 0.9320 at 1/7 and
  sqrt(3)/2 at 1/3 (S82 note section 7); the uniform law of x/y, the Hammersley cells (top digits / bottom digits),
  the count of second states (X/a) ln(1 + a/c) and the summation over the smaller state, with the mean
  (pi^2/12) ln N for unordered word pairs (THM-4885 Theorems A, B); the class law for a fixed smaller state without
  a rate (Session206); in the letter-count model of O/Z/J words with one root change, the two favoured classes 0
  and -1, the confinement mechanism at d = U = 0 (mod p) and the limits L(0) = 1 + (3p-4)/(4p^2),
  L(-1) = 1 - 1/p^2 + 1/(4p) (THM-4916, kind-pasteur-2026-11-07-S15; 3-adic case Session223); the "Fermat boost"
  at census scale as an observation (HYP-9506 (B)); the mod-4 quadratic form of a vertex cube (THM-5127 (C)), of
  which the form in (E) is the case of a chain.
  NEW here: the exact class law (A) with its remainder; the power-saving rate in the regions lemma (B4); the
  by-value statements (C) and (D) with the constants kappa_s, A_p and the sheet split, which THM-4916 lists as NOT
  CLAIMED and HYP-9506 (A), (C) record as OPEN for the root-change census (settled here for one root change: same
  shape as the letter-count law, different split, and an explicit 1/ln X term); the unitary law (E) for every
  power of 2 and the symbol of the mod-4 form; the characterisation (F); the a-priori window predictions (G).
  Chan-Pak remains OPEN; F1 has at most density-one reach and nothing here bears on completeness.
session: opus-2026-11-10-S91
date: 2026-11-10
source: 05-knowledge/results/s91_f1_local_densities_20261110/README.md
scripts:
  - 04-computation/experiments/s91_f1_local_densities_20261110/twisted_spectrum.py, twisted_spectrum_2adic.py, mod4_symbol.py, box_identity.py, divisible_states.py, constants.py, regions_check.py, ap_constants.py, first_order_exact.py (+ .out each; first_order_exact_28.out)
  - 04-computation/experiments/s91_f1_local_densities_20261110/window_split.py, PREREGISTRATION.md, prereg_window_27_24.json, window_split_insample_21_18.out, window_split_insample_23_21.out, window_split_prereg_27_24.out, window_split_extra_24_23.out
  - 04-computation/experiments/s91_f1_local_densities_20261110/blind_proof/check_blind_proof.py (the independent derivation's checks)
related:
  - 01-canon/theorems/THM-4916-*, 05-knowledge/hypotheses/HYP-9506-*, 05-knowledge/results/root_change_local_factors_20261107/ (letter-count local factors, census OPEN)
  - 05-knowledge/results/session223_root_changes_and_local_factors_20261107/, session206/root-residue-capacity/, session206/weighted-parameter-capacity/, 05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md
  - 01-canon/theorems/THM-4875-*, THM-4885-*, THM-4886-*, THM-4878-*, 05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md
  - 05-knowledge/hypotheses/HYP-9657-* (S90; corrected by this session), HYP-9660-* (this session)
  - 01-canon/theorems/THM-5127-* (C) (mod-4 law), THM-5125-* (4b)
---

# THM-5128 - by-value local densities of the one-arc join family

## Setting

`psi(0) = 0`, `psi(2h) = psi(h) + h`, `psi(2h+1) = psi(h) + floor(h/2) + 1`; `q(h) = h - psi(h)`, so
`q(2h) = q(h)`, `q(2h+1) = q(h) + ceil(h/2)`. The one-apex tournament with word `h` has rooted state

    x(h) = 2 psi(h) + 1 = y(h) - 2 q(h)   (its number of Hamiltonian paths),      y(h) = 2h + 1

(canon coordinates: word value `u = 2h`, `g(u) = x`, `1 + u = y`; THM-4875 (D), THM-4885). `rho(h) = x(h)/y(h)`.
The F1 family (one-arc joins of two one-apex tournaments at their apexes) is

    N(h1, h2, s) = M - [s = 0],     M = x(h1) x(h2) + 2^s y(h1) y(h2),        h1, h2, s >= 0,

counted over ORDERED triples. Sheet 0 is `s = 0`, sheet 1 is `s >= 1`; `c_0 = -1`, `c_1 = 0` (so the distinguished
class is `M = 0`). `T_sigma(X)` is the set of triples of the sheet `sigma` with `N <= X`. A state is *p-divisible*
when `p | g(h)`, `g(h) = gcd(x(h), y(h)) = gcd(2h + 1, q(h))`; `chi(h) = [p | g(h)] - 1/p^2`. The canon's ordered
pairs of words `(u, v)` correspond to triples weighted by `s + 1` (census_k2/dict_k1.py). `e(t) = exp(2 pi i t)`;
`B_n` is the set of `h` with exactly `n` binary digits; `l_j = -Li_2(-2^(-j))`.

## Statements

**(A) Class law of a bilinear gluing (exact).** Let `p` be a prime, `r >= 1`, `B` a nondegenerate bilinear form on
`V = (Z/p)^r`, `D1, D2` probability distributions on `V`, `pi_i = D_i(0)`, `U` the uniform distribution. For
independent `u ~ D1`, `v ~ D2` and every `c` in `Z/p`:

    P(B(u, v) = c) = 1/p + (pi_1 + pi_2 - p^(-r)) ([c = 0] - 1/p) + Rhat(c),                                  (1)

where `sum_c Rhat(c) = 0`, `Rhat = 0` identically if `D1` or `D2` is uniform, and

    |Rhat(c)| <= ((p - 1)/p) p^(r/2) |D1 - U|_2 |D2 - U|_2 <= ((p - 1)/p) p^(r/2) E_1 E_2,                       (2)

`|.|_2` the Euclidean norm on functions on `V`, `E_i` the largest modulus of a non-trivial Fourier coefficient of
`D_i`. In words: only the class `0` is biased; to first order its excess is the mass of the states that vanish
modulo `p`; the correction is of second order in the distances from uniformity. If both states are uniform the law
is `1/p + p^(-r) ([c = 0] - 1/p)`. If `D1` is a point mass at `u_0 != 0`, then `P(B = c) = sum of D2 over an affine
hyperplane`, which is `1/p` up to `E_2`. For F1, `r = 2` and `B((a1, b1), (a2, b2)) = a1 a2 + 2^s b1 b2`.

**(B) Digit chain, cells, regions (ingredients; (B1)-(B3) are KNOWN in substance, see the status).** Let `m >= 3`
be odd, `zeta = e(1/m)`. Reading the binary digits of `h` from the top, `(h mod m, h mod 2)` is a Markov state and
`q` an additive functional: with the `2m x 2m` matrices

    T_k[(eta, beta) -> (2 eta, 0)] = 1,       T_k[(eta, beta) -> (2 eta + 1, 1)] = zeta^(k (eta + beta)/2)

(division by 2 modulo `m`) one has `sum_{h < 2^n} zeta^(k q(h) + l h) = e_(0,0)^T T_k^n v_l`, `v_l(eta, beta) = zeta^(l eta)`.
1. `rho(T_k) < 2` for every `k != 0 (mod m)`. Fix `rho_m` with `max(1/2, max_k rho(T_k)/2) < rho_m < 1`.
2. (Cells.) A cell of `B_n` with `r` free digits is the set of `h` in `B_n` with prescribed top digits and
   prescribed bottom digits, `n - r` digits in all. There is `C_m` such that for every cell and every pair of
   residues `(alpha, beta)` modulo `m`:
   `| #{h in cell : (x(h), y(h)) = (alpha, beta) mod m} - 2^r / m^2 | <= C_m (2 rho_m)^r`.
3. `0 <= q(h) <= h/2`, so `y/2 < x <= y`; and `|rho(h) - (1 - phi(h)/2)| < n/(2h + 1)` for `h` in `B_n`, with
   `phi(h) = sum_i h_i 2^(-i-1)` the binary radical inverse. Hence `rho` is asymptotically uniform on `[1/2, 1]`.
4. (Regions.) Let `p` be an odd prime, `theta_p = ln(1/rho_p) / ln(2/rho_p^2)` (so `0 < theta_p < 1/2`). Let
   `b >= a > 0`, `X > 0`, `Omega = {h : a x(h) + b y(h) <= X}`, `I` an interval. Uniformly in `a, b, X, I`:
   (i) `#(B_n ∩ I ∩ Omega) = 2^n area{(u, v) in [1/2, 1) x [0, 1) : 2^n u in I, 2^(n+1) u (a (1 - v/2) + b) <= X} + O(2^(n/2))`;
   (ii) for every residue vector `w`: `#{h in B_n ∩ I ∩ Omega : (x, y) = w mod p} = p^(-2) #(B_n ∩ I ∩ Omega) + O_p(2^(n (1 - theta_p)))`;
   (iii) if `|W| <= W_0` on `B_n` and `W` varies by `O(W_0 2^(-j))` on every cell with `j` top and `j` bottom digits
   prescribed (`j <= n/2`), then `sum over B_n ∩ I of chi W = O_p(W_0 2^(n (1 - theta_p)))`.

Numerical rates `max_k rho(T_k)/2` (VERIFIED, floating point): `0.8660` (m = 3), `0.7418` (5), `0.9320` (7),
`0.8587` (11), `0.8272` (13), `0.8289` (17), `0.8053` (19), `0.7931` (23), `0.9002` (31), `0.8342` (43),
`0.8516` (73), `0.8754` (127); between `0.776` and `0.806` for the other primes up to 89. So
`theta_p / 2 = 0.073, 0.116, 0.042` for `p = 3, 5, 7`.

**(C) First-order by-value law.** Let `p` be an odd prime, `sigma` a sheet, `c` a class, `X >= 2`.
1. (Fixed smaller state.) For `h1, s` let `G = #{h2 > h1 : M <= X}` and `Y = X / 2^s`. If `p` does not divide
   `g(h1)`: `#{h2 > h1 : M <= X, M = c' (mod p)} = G/p + O_p((Y/y1)^(1 - theta_p))`. If `p | g(h1)`: `M = 0 (mod p)`
   for every `h2`.
2. `#{t in T_sigma(X) : N = c (mod p)} = (1/p) #T_sigma(X) + ([c = c_sigma] - 1/p) D_sigma(X) + O_p(X^(1 - theta_p/2))`,
   where `D_sigma(X)` is the number of triples of `T_sigma(X)` with `h1 != h2` whose state of smaller word is
   p-divisible.
3. Equivalently, in shares, with `F_small = D_sigma / #T_sigma` and `f` the share of triples whose first state is
   p-divisible: `share(c) = 1/p + ([c = c_sigma] - 1/p) F_small + O(X^(-theta_p/2))`, and the same with
   `2f - 1/p^2` in place of `F_small` (they differ by the discrepancy of the larger state, `O(X^(-theta_p/2))`).
So, by value, the class distribution of F1 is flat except at the distinguished class, whose excess is the mass of
p-divisible states up to a power-saving error. For a squarefree odd modulus `m` the same holds with the product law:
the class `c` has share `sum over d | m of F^(d)(X) prod_{p | d} [c = c_sigma (mod p)] prod_{p | m/d} (1/p)` up to a
power-saving error, `F^(d)` being the share of triples whose smaller state has `gcd(g, m) = d`.
(VERIFIED, first_order_exact.py, cumulative `X = 2^18, ..., 2^26`: the remainder of the symmetric form, times `p`,
is `0.0003 / 0.0002` at `p = 5`, `0.002` at `17`, `0.004` at `13`, `0.009 / 0.007` at `3`, `0.011 / 0.008` at `31` and
`0.059 / 0.044` at `p = 7` for `X = 2^26` (sheets 0 / 1), and `0.00013`, `0.0013`, `0.002`, `0.014 / 0.008`,
`0.009 / 0.006`, `0.053 / 0.049` at `X = 2^28`: decreasing with `X` on the whole, with oscillations at `3` and `7`; the
one-sided form is less accurate at the slow primes, `0.117` at `p = 7` and `X = 2^26`, because the larger state is
not yet mixed there: its divisible share is still `1.47 / p^2`. The regions lemma (B4)(ii) itself is checked on
random regions in regions_check.py.)

**(D) Mass, limit, and the `1 / ln X` term.**
1. `#T_sigma(X) = kappa_sigma X ln X + O(X)`, `kappa_0 = l_0 - 2 l_1 + l_2 = pi^2/12 + 2 Li_2(-1/2) - Li_2(-1/4) = 0.1615389173`,
   `kappa_1 = l_1 - l_2 = Li_2(-1/4) - Li_2(-1/2) = 0.2125139092`. For one padding `s`:
   `#{(h1, h2) : M <= X} = (l_s - 2 l_(s+1) + l_(s+2)) X ln X + O((s + 1) X 2^(-s))`. So the mean of `r` over odd `N` is
   `g_bar ln N + O(1)`, `g_bar = 2 (kappa_0 + kappa_1) = pi^2/6 + 2 Li_2(-1/2) = 0.7481056530`, and the sheets have
   limiting shares `w_0 = 0.43186124`, `w_1 = 0.56813876`. With the weight `s + 1` (word pairs) the total constant is
   `pi^2/12`, in agreement with THM-4885, and the shares are `w_0 = kappa_0 / (pi^2/12) = 0.19640777`, `w_1 = 0.80359223`.
2. `D_sigma(X) = #T_sigma(X) / p^2 + 2 A_p^sigma X + O_p(X^(1 - theta_p/2))`, where
   `A_p^sigma = sum_{h >= 0} chi(h) w^sigma(h)` (convergent when summed by blocks `B_n`, geometrically),
   `w^0(h) = (1/x) ln(2 (x + y)/(x + 2y))`, `w^1(h) = (1/x) ln(1 + x/(2y))`.
3. Hence, with `mu_c = 1/p + (1/p^2)([c = c_sigma] - 1/p)`:

       #{N = c} / #T_sigma(X) = mu_c + ([c = c_sigma] - 1/p) (2 A_p^sigma / kappa_sigma) / ln X + O(1 / ln^2 X).

   The by-value local factors are `p mu_c = 1 - 1/p^2 + [c = c_sigma]/p`, and for squarefree odd `m` the limiting
   density of a class is the product of these: an Euler product over the odd primes, with singular series
   `S_sigma(N) = prod_{p odd} (1 - 1/p^2 + [p | N + [sigma = 0]]/p) = (8/pi^2) prod_{p | N + [sigma = 0], p odd} (1 + p/(p^2 - 1))`.
   Over both sheets the boost `1/p` is divided between the classes `-1` and `0` in the proportions `(w_0, w_1)`.
4. (NUMERICAL, states below `2^26`; ap_constants.py and, independently, blind_proof/check_blind_proof.py.)
   `A_p^0, A_p^1`: `0.0557, 0.0789` (p = 3, partial sums oscillating by `0.002`); `0.02943, 0.04250` (5); about
   `0.006, 0.007` (7, not converged); `0.01096, 0.01449` (11); `-0.00377, -0.00488` (13); `0.01355, 0.01931` (17);
   `-0.00194, -0.00261` (31); `0.00235, 0.00323` (43). The sign is that of the excess over the generic value at
   feasible sizes: positive when a p-divisible state occurs early, negative (approach from below) at `13` and `31`.
   That `A_p != 0` is numerical. Against the census: `(F_small - 1/p^2) #T / (2X)` at `X = 2^28` is
   `0.0300, 0.0434` (p = 5), `0.0138, 0.0197` (17), `0.0118, 0.0155` (11), `-0.0035, -0.0044` (13),
   `-0.0017, -0.0021` (31), approaching `A_p`; at `p = 7` it is still `-0.005, -0.010` and rising (rate `0.932`).

**(E) The prime 2.** Modulo `2^j` the increment `ceil(h/2)` is a function of `w = h mod 2^(j+1)`, and for a character
`k` of `Z/2^j` the twisted matrix is `T_k[w -> 2w] = 1`, `T_k[w -> 2w + 1] = e(k c(w) / 2^j)`,
`c(w) = (w + (w mod 2))/2`. For every odd `k`, `T_k T_k^* = 2 I`: `T_k / sqrt 2` is unitary, and
`|sum_{h < 2^n} e(k q(h)/2^j)| <= 2^((n + j + 1)/2)`. In particular `q(h) = sum_i h_i (h_(i+1) + h_(i+2)) (mod 2)` in
the binary digits of `h`, a quadratic form of Rudin-Shapiro type, and `(x(h), y(h))` is equidistributed on the
pairs of odd residues modulo `2^j`. Consequently `N` is asymptotically equidistributed on the odd classes modulo
every power of 2, on each sheet: F1 has no local factor at 2.
*The mod-4 imbalance, exactly.* The polar form of `q mod 2` on words with `n` digits is the `n x n` band matrix `B_n`
over `F_2` with ones at distances 1 and 2 from the diagonal, a Toeplitz matrix with symbol
`z + z^2 + z^(-1) + z^(-2) = z^(-2) (1 + z)^2 (1 + z + z^2)`. Its corank is `1, 0, 1, 2, 1, 0` for
`n = 1, 2, 3, 4, 5, 0 (mod 6)`, and `S_n = sum_{h < 2^n} (-1)^q(h)` and
`D_n = #{h < 2^n : x(h) = 1 mod 4} - #{h < 2^n : x(h) = 3 mod 4}` are `0` or `+- 2^(n - rank/2)` (for example
`D_18 = -512`, `D_20 = 1024`, `D_27 = 16384`, `D_n = 0` for `n = 1, 4, 5, 7, 9, 10 (mod 12)`; mod4_symbol.py, exact
for `n <= 48`). This is the chain case of THM-5127 (C) made explicit; the factor `1 + z + z^2` is the polynomial of
THM-5125 (4b).

**(F) Divisible states.** `x(h) = y(h)` if and only if `h = 0` or `h = 2^j`, and then `x = y = 2^(j+1) + 1`. Hence
every prime factor of `2^m + 1` has the p-divisible state `h = 2^(m-1)`: `3` at `h = 1`, `5` at `2`, `17` at `8`,
`11` at `16`, `13` at `32`, `43` at `64`, `257` at `128`, `65537` at `2^15`. The density of p-divisible states is
`1/p^2`, and the density of the `h` with `g(h)` coprime to a given squarefree odd `m` is `prod_{p | m} (1 - 1/p^2)`.
(FINITE-EXACT: first p-divisible states `7: 38`, `19: 218`, `23: 103`, `29: 1435`, `31: 8664`, `37: 1350`,
`41: 266`, `127: 109283`; proportion of `h < 2^24` with `g(h) = 1`: `0.81176`, against `8/pi^2 = 0.81057`; that the
density of `g = 1` equals `8/pi^2` is NOT proved, the tail of the product being uncontrolled.)

**(G) Windows (NUMERICAL; pre-registered).** By (B3) the number of triples with first state `h1` and padding `s`
whose value lies in `[A, B)` is `W(h1, s) = (1/x1) int_{f0}^{f1} w(f) df`, `f0 = x1/2 + 2^s y1`, `f1 = x1 + 2^s y1`,
`w(f) = (B - A)/f` for `f < A`, `B/f` for `A <= f < B`. With `f` the `W`-weighted share of p-divisible first states
and `D = 2f - 1/p^2`, the predicted excess of the distinguished class over the mean of the other classes is
`E = p D / (1 - D)`, and after removal of the p-divisible triples `E_nd = p (f^2 - 1/p^2)/(1 - D)`.
On `[2^24, 2^27)` (unseen; 8.3e8 triples): total mass `14.1405` per odd `N` against `14.1417` predicted; weights of
`h1 = 0, 1, 2, 3, 8, 16, 64, 128, 1024` within `0.006 %`; `E` measured / predicted on sheet 0 and sheet 1:
`p = 5`: `+0.3130 / +0.3120`, `+0.3280 / +0.3284`; `11`: `+0.1690 / +0.1705`, `+0.1707 / +0.1749`;
`17`: `+0.2107 / +0.2081`, `+0.2296 / +0.2288`; `19`: `+0.0510 / +0.0511`, `+0.0495 / +0.0497`;
`29`: `+0.0051 / +0.0032`, `+0.0020 / +0.0007`; `43`: `+0.0885 / +0.0877`, `+0.0930 / +0.0941`;
`257`: `+0.1859 / +0.1815`, `+0.2023 / +0.2044`; `2731`: `+0.063 / +0.061`, `+0.075 / +0.069`;
`65537`: `+0.144 / +0.183`, `+0.168 / +0.206` (standard error `0.022`); controls `65521`, `65539`, `8191`, `2729`,
`2741`, `251`, `263`, `677`, `691` within noise of zero. At the slow primes the second-order term is visible:
`p = 3`: `+0.539 / +0.515`, `+0.544 / +0.532`; `p = 7`: `+0.233 / +0.174`, `+0.230 / +0.175`. All eight
pre-registered tests (279 intervals) CONFIRMED. The same rules also hold on `[2^23, 2^24)` (not pre-registered).

## Proofs

**(A).** Write `Dhat(xi) = sum_v D(v) e(xi . v / p)` and `B(u, v) = (B* u) . v` with `B*` invertible. Then
`P(B = c) = (1/p) sum_t e(-tc/p) Phi(t)`, `Phi(t) = sum_u D1(u) D2hat(t B* u)`, `Phi(0) = 1`. For `t != 0` separate
`u = 0` (it contributes `pi_1`) and put `D1 = p^(-r) + delta_1` on `u != 0`. Since `u -> t B* u` permutes the
non-zero vectors and `sum_xi D2hat(xi) = p^r D2(0)`,

    Phi(t) = pi_1 + p^(-r) (p^r pi_2 - 1) + R(t) = pi_1 + pi_2 - p^(-r) + R(t),       R(t) = sum_{u != 0} delta_1(u) D2hat(t B* u).

Summing over `t` with `sum_{t != 0} e(-tc/p) = p [c = 0] - 1` gives (1) with
`Rhat(c) = (1/p) sum_{t != 0} e(-tc/p) R(t)`. `sum_c Rhat(c) = 0` because `sum_c e(-tc/p) = 0` for `t != 0`. If `D2`
is uniform, `D2hat` vanishes off `0`; if `D1` is uniform, `delta_1 = 0`. By Cauchy-Schwarz and Parseval
(`sum_{xi != 0} |Dhat(xi)|^2 = p^r |D - U|_2^2`), `|R(t)| <= |delta_1|_2 (sum_{xi != 0} |D2hat(xi)|^2)^(1/2) = p^(r/2) |D1 - U|_2 |D2 - U|_2`,
which gives (2); the second form follows from `|D - U|_2^2 <= E^2`. For a point mass at `u_0 != 0` the event
`B(u_0, v) = c` is `v` in an affine hyperplane, and `sum of D2 over it = 1/p + (1/p) sum_{t != 0} e(-tc/p) D2hat(t B* u_0)`.

**(B1).** The entrywise modulus of `T_k` is `T_0`, whose rows sum to `2`, so `rho(T_k) <= 2`. `P = T_0 / 2` is the
transition matrix of an irreducible aperiodic chain: from any state, after `t` steps with `2^t >= 2m`, the states
`(2^t eta + j, j mod 2)`, `0 <= j < 2^t`, are all states. Suppose `T_k g = lambda g` with `|lambda| = 2`, `g != 0`.
Then `|g| <= P |g|` entrywise; pairing with the stationary vector of `P` (positive) gives equality, so `|g|` is
`P`-harmonic, hence constant, and at every state the two terms of `(T_k g)(eta, beta)` have the same argument:

    g(2 eta, 0) = zeta^(k (eta + beta)/2) g(2 eta + 1, 1) = (lambda / 2) g(eta, beta).

The states `(eta, 0)` and `(eta, 1)` have the same two successors. Comparing the two equations gives
`zeta^(k eta / 2) = zeta^(k (eta + 1)/2)`, i.e. `zeta^(k/2) = 1`, i.e. `m | k`. So `rho(T_k) < 2` for `k != 0 (mod m)`.

**(B2).** Let `T^(0)`, `T^(1)` be the parts of `T_k` for the digits `0`, `1`. Each has exactly one non-zero entry, of
modulus `1`, in every row. Hence a product of them maps `e_(0,0)^T` to a unimodular multiple of a coordinate row
vector, and maps a column vector of sup norm at most `1` to one of sup norm at most `1`. The character sum of a
cell is `(e_(0,0)^T T[top]) T_k^r (T[bottom] v_l)`, so its modulus is at most the norm of `T_k^r` as an operator on
`l^infinity`, which is at most `C (2 rho_m)^r` (Gelfand; finitely many `k`; `rho_m` is strictly above the spectral
rates because of possible Jordan blocks). For `k = 0`, `l != 0` the sum over the free digits is geometric with
ratio `zeta^(l 2^j) != 1`, of modulus at most `m/2`. Fourier inversion on `(Z/m)^2` gives the claim for `(q, h)`,
and `(x, y) = (2h + 1 - 2q, 2h + 1)` is an affine bijection of `(Z/m)^2`.

**(B3).** `q(h) = sum over one-bits i of ceil(H_i / 2)`, `H_i = floor(h / 2^(i+1))`, by induction on the recursion;
`q(2h) = q(h) <= h/2` and `q(2h + 1) <= h/2 + (h + 1)/2` give `q <= h/2`. Next
`2 (2h + 1)(rho - 1 + phi/2) = sum over one-bits i of [(2h + 1)/2^(i+1) - 4 ceil(H_i/2)]`; the first term lies in
`(2 H_i + 1, 2 H_i + 2)` and the second is `2 H_i` or `2 H_i + 2`, so every bracket lies in `(-1, 2)` and there are at
most `n` of them. On `B_n` the radical inverse takes each value `(2i + 1)/2^n` once.

**(B4).** Use the cells of `B_n` with `j` top and `j` bottom digits prescribed, `j <= n/2`. The top word `A` (an
integer in `[2^(j-1), 2^j)`) gives `2^(n-j+1) A < y < 2^(n-j+1) (A + 1)`. The bottom word gives `phi` in
`[Phi, Phi + 2^(-j))`, hence by (B3) `rho` in `(R - 2^(-j-1) - eps, R + eps)`, `R = 1 - Phi/2`, `eps = n 2^(-n)`. Since
`a x + b y = y (a rho + b)`, a cell lies inside `Omega` if `2^(n-j+1) (A + 1)(a (R + eps) + b) <= X` and outside it if
`2^(n-j+1) A (a (R - 2^(-j-1) - eps) + b) >= X`. For a fixed bottom word the remaining `A` lie in an interval
`(P - 1, P (1 + 2^(-j-1) + 2 eps))` with `P < 2^j` (here `a <= b` is used), which contains fewer than five integers.
So at most `4 * 2^j` cells are cut by `Omega` and the cells of at most two top words by `I`: at most
`6 * 2^(n-j)` points. (i): take `j = floor(n/2)`; a cell inside (outside) `Omega` has its rectangle inside (outside)
the region and contains exactly `2^n` times its area. (ii): take `j = floor(theta_p n)`; by (B2) the uncut cells
contribute an error at most `C 2^n rho_p^(n - 2j)`, and `rho_p^(1 - 2 theta_p) = 2^(-theta_p)`. (iii): the same,
freezing `W` on each cell.

**(C1).** Here `M = a x2 + b y2` with `a = x1 <= b = 2^s y1`. If `p` does not divide `g(h1)` the linear form
`w -> x1 w_1 + 2^s y1 w_2` on `(Z/p)^2` is non-zero and takes each value on `p` residue vectors; sum (B4)(ii) over
these vectors and over the blocks `B_n` with `I = (h1, infinity)`; `M > y1 y2 2^s` bounds the blocks by
`2^n <= Y/y1`. If `p | g(h1)` both coefficients vanish modulo `p`.

**(C2).** `M` is symmetric in `(h1, h2)`, so each triple with `h1 != h2` is counted through its smaller state, and
those with `h1 = h2` number `O(Y^(1/2))` for each `s`. `G = 0` unless `y1 < Y^(1/2)`. Summing (C1) over the smaller
states, the errors add up to `sum_{y1 < Y^(1/2)} (Y/y1)^(1 - theta) = O(Y^(1 - theta/2))`; summing over
`1 <= s <= log_2 X` on sheet 1 keeps `O(X^(1 - theta/2))`. On sheet 0, `N <= X` is `M <= X + 1`. (C3): by (B4)(ii)
with `w = 0` the larger state is p-divisible in a share `1/p^2` of the triples up to the same error, and
`2f = F_small + F_large`. The squarefree case is the same proof with (B2) for `m` and (C1) applied prime by prime.

**(D1).** By (B4)(i) summed over the blocks (the variable `t = 2^n u` runs once over `[1, infinity)`),
`G = Gamma_X(h1) + O((Y/y1)^(1/2) + 1)`, `Gamma_X(h1) = int_0^1 (Z(v) - max(h1, 1))^+ dv`,
`Z(v) = X / (2 y1 (rho_1 (1 - v/2) + 2^s))`, and `int_0^1 Z = X w_s(h1)` with
`w_s(h) = (1/x) ln((x + 2^s y)/(x/2 + 2^s y))`; so `Gamma_X = X w_s - Xi_X`, `0 <= Xi_X(h1) <= max(h1, 1)`. (This is
THM-4885 Theorem A in the present coordinates.) Hence `#{(h1, h2) : M <= X} = 2X sum_{y1 < Y^(1/2)} w_s(h1) + O(Y)`.
Write `w_s = F_s(rho)/y`, `F_s(rho) = (1/rho) ln((rho + 2^s)/(rho/2 + 2^s)) <= 2^(-s-1)`; by the cells with
`j = floor(n/2)`, `sum over B_n of w_s = ln 2 * int_{1/2}^{1} F_s + O(2^(-s) 2^(-n/2))`, and
`int_{1/2}^{1} F_s = [Li_2(-rho/2^(s+1)) - Li_2(-rho/2^s)]_{1/2}^{1} = l_s - 2 l_(s+1) + l_(s+2)`. The blocks up to
`Y^(1/2)` number `(1/2) log_2 Y`, which gives `(l_s - 2 l_(s+1) + l_(s+2)) X ln Y + O(Y)`. Summing over `s >= 1`
telescopes to `l_1 - l_2`; with the weight `s + 1`, `sum_{s >= 0} (s + 1)(l_s - 2 l_(s+1) + l_(s+2)) = l_0 = pi^2/12`.

**(D2).** `D_s - T_s/p^2 = 2 sum_{h1} chi(h1) G + O(Y^(1/2)) = 2 sum chi Gamma_X + O(Y^(3/4))`. On `B_n` the weights
`w_s` (`W_0 = 2^(-s-n)`) and `Xi_X` (`W_0 = 2^n`) satisfy the hypothesis of (B4)(iii) (they are smooth functions of
`y` and `rho`, and `Xi_X = int_0^1 min(Z, max(h1, 1))`). So the block sums of `chi w_s` are `O(2^(-s) 2^(-theta n))`:
the series `A_(p,s) = sum chi w_s` converges, its tail beyond `Y^(1/2)` is `O(2^(-s) Y^(-theta/2))`, and
`sum_{y1 < Y^(1/2)} chi Xi_X = O(Y^(1 - theta/2))`. Summing over `s` gives the statement with
`w^1 = sum_{s >= 1} w_s = (1/x) ln(1 + x/(2y))`. **(D3)** follows from (C2), (D1), (D2); the product law from the
squarefree case of (C).

**(E).** The rows `w` and `w + 2^j` of `T_k` have the same two non-zero columns, with entries `(1, phi)` and
`(1, (-1)^k phi)` because `c(w + 2^j) = c(w) + 2^(j-1)`; for odd `k` they are orthogonal, rows with different
column pairs have disjoint supports, and every row has norm `sqrt 2`. So `T_k T_k^* = 2I`, and
`|e_0^T T_k^n 1| <= 2^(n/2) 2^((j+1)/2)`. An even `k` is a character of `Z/2^(j-1)`. `ceil(m/2) = m_0 + m_1 (mod 2)`
in the binary digits of `m` gives the quadratic form. Joint equidistribution with the last `i` digits of `h`: on the
state space `h mod 2^max(j+1, i)` let `v` be a function of the last `i` digits; `T_k^i v` depends on the state only
modulo `2^(j+1)` (the phases do, and the last `i` digits are the appended ones), so the remaining `n - i` steps are
governed by the unitary matrix above and the sum is `O(2^(n/2))`. Then `x = y - 2q` and `y = 2h + 1` are
independent uniform odd residues in the limit, `x1 x2` is a uniform odd residue, and adding `2^s y1 y2` (sheet 1) or
`y1 y2 - 1` (sheet 0; the sum of two independent uniform odd residues is a uniform even residue) leaves `N` uniform
on the odd classes.
*Corank of `B_n`.* A vector `v` is in the kernel iff, extended by `v_(-1) = v_0 = v_(n+1) = v_(n+2) = 0`, it
satisfies `v_(i-2) + v_(i-1) + v_(i+1) + v_(i+2) = 0` for `1 <= i <= n`. The solutions of this recurrence with
`v_(-1) = v_0 = 0` form a space of dimension 2; its characteristic polynomial `(1 + z)^2 (1 + z + z^2)` divides
`z^6 - 1` over `F_2`, so every solution has period 6 and the two end conditions depend only on `n` modulo 6; they
hold for the whole space exactly when `n + 1 = -1` and `n + 2 = 0 (mod 6)`. The values `1, 0, 1, 2, 1, 0` are read
off the two basic solutions. A quadratic form of rank `r` on `F_2^n`, plus a linear form, has character sum `0` or
`+- 2^(n - r/2)`.

**(F).** `x = y` iff `q(h) = 0` iff every one-bit of `h` has no one-bit above it. The densities follow from (B2)
with `m` the product of the primes considered.

## Remarks

1. *Relation to the letter-count law.* THM-4916 (one root change) gives `L(0) = 1 - 1/p^2 + 3/(4p)`,
   `L(-1) = 1 - 1/p^2 + 1/(4p)`: the boost `1/p` split `(3/4, 1/4)` by the parity of the tail. By value the split is
   `(w_1, w_0) = (0.80359, 0.19641)` for word pairs. The census limit therefore differs from the letter-count limit
   already for one root change. This bears on HYP-9506 (A), which conjectures equality of the two limits for two
   root changes; see HYP-9660.
2. *Why the feasible windows look irregular.* The excess at `p` is governed by `F_small - 1/p^2`, which is
   `2 A_p / (kappa ln X)` to first order: far above zero when a p-divisible state is small (`3, 5, 17, 11, 43, 257,
   65537`, the prime factors of `2^m + 1`: the "Fermat boost" of HYP-9506 (B)), and negative when the first
   p-divisible states are late or sparse (`13`, `29`, `31`, `37`, `61`, `73`, `127`: the distinguished class is then
   under-represented, as the S90 audit observed at 73 and 127). The other classes are flat up to the power-saving
   term, which is visible only at slowly mixing primes (`7`: rate `0.932`).
3. *What the law does not say.* (A) is a statement about one prime; the representation count is over-dispersed
   (relative variance about `0.75` on a sheet), so a positive singular series does not give `r(N) > 0`.
