---
id: THM-5128
title: "By-value local densities of the one-arc join family F1. (A) Exact class law of a bilinear gluing: for independent states u ~ D1, v ~ D2 in (Z/p)^r and a nondegenerate bilinear form B, P(B(u,v) = c) = 1/p + (pi1 + pi2 - p^-r)([c = 0] - 1/p) + Rhat(c), pi_i = D_i(0), |Rhat| <= ((p-1)/p) p^(r/2) |D1 - U|_2 |D2 - U|_2: only the class 0 is biased, and to first order its excess over the mean of the other classes is pi1 + pi2 - p^-r, the mass of the states that vanish modulo p less p^-r. (C) For F1 (N = x1 x2 + 2^s y1 y2 - [s = 0] over one-apex states (x, y) = (2 psi(h) + 1, 2h + 1), counted by value N <= X) on each sheet: #{N = c mod p} = #T/p + ([c = c_s] - 1/p) D(X) + O(X^(1 - theta_p/2)), D(X) the number of triples whose smaller state is p-divisible (p | gcd(x, y)), c_s = -1 for s = 0 and 0 for s >= 1. (D) #T = kappa_s X ln X + O(X) and D(X) = #T/p^2 + 2 A_p X + O(X^(1 - theta_p/2)): the by-value class factors modulo an odd prime tend to 1 - 1/p^2 + [c = c_s]/p, with an explicit correction of order 1/ln X, and the limiting density of a class modulo a squarefree odd modulus is their product (the squarefree-level product is a class average, not the singular series; modulo odd prime powers the by-value law is numerical only); over both sheets the boost 1/p is split (0.43186, 0.56814) between the classes -1 and 0 for triples and (0.19641, 0.80359) for word pairs, not the (1/4, 3/4) of the letter-count model. (E) At the prime 2, for odd characters the twisted transfer operator of the digit recursion is sqrt(2) times a unitary matrix; by value every odd class modulo 2^j receives the share 2^(1-j) of the triples of a sheet up to O(X^(7/8)) triples, so F1 has no 2-adic bias; the mod-4 law of the one-apex family is a band form with symbol (1 + z)^2 (1 + z + z^2), of corank periodic with period 6, and its character sums have period 12 up to the factor -64. (F) x = y exactly at h = 0 and h = 2^j (x = y = 2^(j+1) + 1 at h = 2^j), so every prime factor p of 2^m + 1 has the p-divisible state h = 2^(m-1), which is below the generic scale p^2 only for small m (3, 5, 11, 13, 17, 19, 41, 43, 113, 241, 257, 65537, ...)."
status: >
  Independent audit 2026-11-11 (seven auditors; record
  05-knowledge/results/s91_f1_local_densities_20261110/audit/AUDIT_REPORT.md): PASS WITH CORRECTIONS, all applied;
  MISTAKE-761.
  PROVED below (elementary; proofs complete in this file): (A); (B); (C) and (D1)-(D3) at every odd prime and
  every squarefree odd modulus; (E); (F). (C)-(D) were derived twice, by the session and by an independent blind
  derivation (text blind_proof/BLIND_DERIVATION.md, checks blind_proof/check_blind_proof.out), and merged. In (E)
  the by-value statement and its proof, and the period-12 law of the mod-4 sums, were supplied by the independent
  audit (auditor proof-5128-CF); the session had proved the statement on boxes only.
  VERIFIED / NUMERICAL (floating point): the rates rho_m; the constants A_p (A_7 is not determined, sign
  included). NUMERICAL only in this file: the by-value law modulo odd prime powers in (D3) (the factor for independent
  equidistributed states is an exact elementary count; the by-value limit is census data). ADDENDUM 2026-11-11
  (opus-2026-11-11-S92): the by-value law at every odd modulus, prime powers included, is now PROVED in THM-5134
  (the argument announced in (D3), written out); the series (8/pi^2) sigma_1(n)/n is THM-5134 (E).
  (G): NUMERICAL, with eight PRE-REGISTERED tests (279 intervals) on the unseen window [2^24, 2^27), all CONFIRMED
  (pre-registration commit 1f97e27805 precedes the output); the constant 0.75 of the registered tolerance was a
  conservative bound, so the registered "4 se" bands were about 5.2 to 5.5 measured standard errors wide, and the
  registered intervals of 43691 and 241 contain 0.
  SCOPE: the residue-class side of F1 is determined at each fixed odd prime and squarefree odd modulus (and every
  power of 2), with a power saving uniform in the region but not in the modulus; uniformity in the modulus is
  HYP-9660 (B), OPEN; odd prime-power moduli are numerical only.
  KNOWN before this session, and used: the digit chain on (h mod m, h mod 2), its uniform stationary law and its
  mixing (Session206 root-residue-capacity and weighted-parameter-capacity; Session193 modular_sync Theorems 3-4;
  THM-4886 (C)); the twisted transfer matrices and their rates for the count alone: the scan of all a/q with
  q <= 120, the extremality of 1/7 (0.9320), the next-worst fractions 11/93, 1/31, 1/15, 1/63, and sqrt(3)/2 at
  1/3 (S82 note section 7); the uniform law of x/y, the Hammersley cells (top digits / bottom digits),
  the count of second states (X/a) ln(1 + a/c) and the summation over the smaller state, with the mean
  (pi^2/12) ln N for unordered word pairs (THM-4885 Theorems A, B); the class law for a fixed smaller state without
  a rate, that is the two-class by-value law of a fixed p-divisible seed with the weights c_0, c_1 (Session206
  weighted-parameter-capacity; THM-4888 (C); HYP-9480 (A)); for the bank family F2 at p = 3, the first-order
  prediction (the local factor comes entirely from the triples with 3 | y_h and 3 | y_R, split by the parity of
  v), with the weights 2 w^1, 2 w^0 of (D2), and its comparison with the census (S82 note section 5,
  chanpak_mass_exploratory_20261106/local_factor_mod3_prediction.py); in the letter-count model of O/Z/J words
  with one root change, the two favoured classes 0 and -1, the confinement mechanism at d = U = 0 (mod p) and the
  limits L(0) = 1 + (3p-4)/(4p^2), L(-1) = 1 - 1/p^2 + 1/(4p) (THM-4916, kind-pasteur-2026-11-07-S15; 3-adic case
  Session223); the "Fermat boost" at census scale as an observation, and the boost of the primes dividing 2^i + 1
  through the multiplier tails (HYP-9506 (B); THM-4916 (1c)); x = y exactly at h = 0 and h = 2^j (THM-5126 (A),
  kind-pasteur-2026-11-10-S21: the only unit-offset one-apex roots are the fans and the transitive tournament
  rooted at its sink; its delta' is 2 q(h)); the count 2^m + 1 of one back arc (classical, Moon 1972 as cited in
  THM-4748; THM-4861); the 2-adic balance of F1 as an observation (S78 note, residue profile of the F1 misses;
  S82 note section 5); equidistribution of the one-apex count modulo 2^j (Session192 dyadic; Session193
  modular_sync Theorem 3; THM-4886 (C)); the earlier H mod 4 laws (THM-466, THM-478); the mod-4 quadratic form of
  a vertex cube (THM-5127 (C)), of which the form in (E) is the case of a chain; that no congruence excludes any
  odd number, 7 and 21 included (the one-apex family alone attains every residue modulo every odd modulus and every odd
  residue modulo every even modulus: Session 25; Session193 modular_sync Theorem 1; Session206 exact-modular-image).
  NEW here: the exact class law (A) with its remainder, as a statement of the join calculus (the identity is a
  standard character-sum computation and the bound the standard Fourier bound for a nondegenerate bilinear
  kernel); the power-saving rate in the regions lemma (B4); the extension of the rate scan to every odd m <= 521
  and 1023, 1039, 1143; the summation over the smaller state with a power-saving error, for F1 and every odd
  prime, that is the by-value statements (C) and (D) with the constants kappa_s, A_p and the sheet split, which
  THM-4916 lists as NOT CLAIMED and HYP-9506 (A), (C) record as OPEN for the root-change census (settled here for
  one root change: same shape as the letter-count law, different split, and an explicit 1/ln X term); in (E): the
  unitarity T_k T_k^* = 2I for odd characters with the bound 2^((n+j+1)/2), the by-value equidistribution of N
  modulo 2^j with a power saving, the symbol, the corank and the period 12 of the mod-4 form; in (F): the reading
  of the fans as p-divisible states of the first-order law, the density 1/p^2 of p-divisible states and the table
  of first p-divisible states; the a-priori window predictions (G) for F1 (method of the S82 prediction). What
  this adds to the known absence of congruence obstructions is the size of the local densities of F1 at fixed odd
  primes, squarefree odd moduli and powers of 2.
  Chan-Pak remains OPEN; F1 has at most density-one reach and nothing here bears on completeness.
session: opus-2026-11-10-S91
date: 2026-11-10
source: 05-knowledge/results/s91_f1_local_densities_20261110/README.md
scripts:
  - 04-computation/experiments/s91_f1_local_densities_20261110/twisted_spectrum.py, twisted_spectrum_2adic.py, mod4_symbol.py, box_identity.py, divisible_states.py, constants.py, regions_check.py, ap_constants.py, first_order_exact.py (+ .out each; first_order_exact_28.out)
  - 04-computation/experiments/s91_f1_local_densities_20261110/window_split.py, PREREGISTRATION.md, prereg_window_27_24.json, window_split_insample_21_18.out, window_split_insample_23_21.out, window_split_prereg_27_24.out, window_split_extra_24_23.out
  - 04-computation/experiments/s91_f1_local_densities_20261110/rates_exact/rates_by_modulus.txt (receipt for every rate; rates_scan.py), bridges/b3_f1_data.out (the modulus 1039)
  - 04-computation/experiments/s91_f1_local_densities_20261110/blind_proof/BLIND_DERIVATION.md, check_blind_proof.py, check_blind_proof.out (the independent derivation and its checks)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/ (AUDIT_REPORT.md and the auditors' receipts in proof-5128-AB/, proof-5128-CF/, recompute-5128/, numbers-vs-outputs/, citations/, prior-art-wording/)
related:
  - 01-canon/theorems/THM-4916-*, 05-knowledge/hypotheses/HYP-9506-*, 05-knowledge/results/root_change_local_factors_20261107/ (letter-count local factors, census OPEN)
  - 05-knowledge/results/session223_root_changes_and_local_factors_20261107/, session206/root-residue-capacity/, session206/weighted-parameter-capacity/, 05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md
  - 01-canon/theorems/THM-4875-*, THM-4885-*, THM-4886-*, THM-4878-*, 05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md
  - 05-knowledge/hypotheses/HYP-9657-* (S90; corrected by this session), HYP-9660-* (this session)
  - 01-canon/theorems/THM-5127-* (C) (mod-4 law), THM-5125-* (4b)
  - 01-canon/theorems/THM-5126-* (A) (fans), THM-4748-*, THM-4861-* (the count 2^m + 1), THM-4888-* (C), 05-knowledge/hypotheses/HYP-9480-* (A), 01-canon/theorems/THM-466-*, THM-478-* (H mod 4)
  - 05-knowledge/results/old_mark_bound_lorentzian_gadgets_one_apex_join_coverage_20261104.md (S78 note), 05-knowledge/results/session192_surprising_approaches_20261105/dyadic/PROOF.md, session206/exact-modular-image/, 05-knowledge/results/session25_spectrum_wiring_20261021/
  - 04-computation/experiments/chanpak_mass_exploratory_20261106/local_factor_mod3_prediction.py, twisted_scan_q120.out (S82 note sections 5, 7)
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
`B_n` is the set of `h` with exactly `n` binary digits (`2^(n-1) <= h < 2^n` for `n >= 1`, `B_0 = {0}`);
`l_j = -Li_2(-2^(-j))`; `sigma_1(n)` is the sum of the divisors of `n` (`sigma` alone always denotes a sheet).

## Statements

**(A) Class law of a bilinear gluing (exact).** Let `p` be a prime, `r >= 1`, `B` a nondegenerate bilinear form on
`V = (Z/p)^r`, `D1, D2` probability distributions on `V`, `pi_i = D_i(0)`, `U` the uniform distribution. For
independent `u ~ D1`, `v ~ D2` and every `c` in `Z/p`:

    P(B(u, v) = c) = 1/p + (pi_1 + pi_2 - p^(-r)) ([c = 0] - 1/p) + Rhat(c),                                  (1)

where `sum_c Rhat(c) = 0`, `Rhat = 0` identically if `D1` or `D2` is uniform, and

    |Rhat(c)| <= ((p - 1)/p) p^(r/2) |D1 - U|_2 |D2 - U|_2 <= ((p - 1)/p) p^(r/2) E_1 E_2,                       (2)

`|.|_2` the Euclidean norm on functions on `V`, `E_i` the largest modulus of a non-trivial Fourier coefficient of
`D_i`. For `c != 0` the bound also holds with the constant `p^((r-1)/2)`, which is smaller than that of (2) for
`p >= 3` (at `p = 2` it is larger); for `c = 0` and `r = 2` the constant `p - 1` of (2) is attained (two lines
through the origin against their `B`-orthogonals; see the proof). In words: the classes
`c != 0` have equal shares to first order, and the excess of the class `0` over their mean is
`pi_1 + pi_2 - p^(-r)`, the mass of the states that vanish modulo `p` less the generic `p^(-r)`, up to
`(p/(p-1)) Rhat(0)`; the correction is of second order in the distances from uniformity. If both states are
uniform the law is `1/p + p^(-r) ([c = 0] - 1/p)`. If `D1` is a point mass at `u_0 != 0`, then
`P(B = c) = sum of D2 over an affine hyperplane`, which is `1/p` up to `E_2`. For F1 and an odd prime `p`, `r = 2`
and `B((a1, b1), (a2, b2)) = a1 a2 + 2^s b1 b2` (nondegenerate because `2^s` is a unit; at `p = 2`, `s >= 1` the
form has rank 1 and (A) does not apply: see (E)). The identity (1) is a standard character-sum computation and
(2) is the standard Fourier bound for the kernel `e(t B(u, v)/p)`, which is `p^(r/2)` times a unitary matrix
(literature name "Lindsey-type lemma", and for `c != 0` the point-hyperplane incidence bound: reported by the
audit; source not located); what is used below is the reading of (1) for a gluing of two rooted states.

**(B) Digit chain, cells, regions (ingredients; (B1)-(B3) are KNOWN in substance, see the status).** Let `m >= 3`
be odd, `zeta = e(1/m)`. Reading the binary digits of `h` from the top, `(h mod m, h mod 2)` is a Markov state and
`q` an additive functional: with the `2m x 2m` matrices

    T_k[(eta, beta) -> (2 eta, 0)] = 1,       T_k[(eta, beta) -> (2 eta + 1, 1)] = zeta^(k (eta + beta)/2)

(division by 2 modulo `m`) one has `sum_{h < 2^n} zeta^(k q(h) + l h) = e_(0,0)^T T_k^n v_l`, `v_l(eta, beta) = zeta^(l eta)`.
1. `rho(T_k) < 2` for every `k != 0 (mod m)`. Fix `rho_m` with `max(1/2, max_k rho(T_k)/2) < rho_m < 1`.
2. (Cells.) A cell of `B_n` with `r` free digits is the set of `h` in `B_n` with prescribed top digits (the leading
   digit 1 included, so at least one) and prescribed bottom digits, `n - r` digits in all; it has `2^r` elements.
   There is `C_m` such that for every cell and every pair of residues `(alpha, beta)` modulo `m`:
   `| #{h in cell : (x(h), y(h)) = (alpha, beta) mod m} - 2^r / m^2 | <= C_m (2 rho_m)^r`.
3. `0 <= q(h) <= h/2`, so `y/2 < x <= y`; and `|rho(h) - (1 - phi(h)/2)| < n/(2h + 1)` for `h` in `B_n`, `n >= 1`,
   with `phi(h) = sum_i h_i 2^(-i-1)` the binary radical inverse. Hence `rho` is asymptotically uniform on
   `[1/2, 1]`. (FINITE-EXACT, audit, `1 <= h < 2^24`: the deviation is one-signed and smaller,
   `0 < rho(h) - (1 - phi(h)/2) <= 3n/(4(2h + 1))`, with equality only at `h = 1`; proof-5128-AB/check_B3.out,
   check_B3_sign.out.)
4. (Regions.) Let `p` be an odd prime, `theta_p = ln(1/rho_p) / ln(2/rho_p^2)` (so `0 < theta_p < 1/3`, because
   `rho_p > 1/2`). Let `b >= a > 0`, `X > 0`, `Omega = {h : a x(h) + b y(h) <= X}`, `I` an interval. Uniformly in
   `a, b, X, I`:
   (i) `#(B_n ∩ I ∩ Omega) = 2^n area{(u, v) in [1/2, 1) x [0, 1) : 2^n u in I, 2^(n+1) u (a (1 - v/2) + b) <= X} + O(2^(n/2))`;
   (ii) for every residue vector `w`: `#{h in B_n ∩ I ∩ Omega : (x, y) = w mod p} = p^(-2) #(B_n ∩ I ∩ Omega) + O_p(2^(n (1 - theta_p)))`;
   (iii) if `|W| <= W_0` on `B_n` and `W` varies by `O(W_0 2^(-j))` on every cell with `j` top and `j` bottom digits
   prescribed (`j <= n/2`), then `sum over B_n ∩ I of chi W = O_p(W_0 2^(n (1 - theta_p)))`.
   (iv) (Odd moduli.) (ii) and (iii) hold, with the same proof, for an odd modulus `m >= 3` in place of `p`: residue
   vectors modulo `m`, exponent `theta_m = ln(1/rho_m) / ln(2/rho_m^2)`, and in (iii), in place of `chi`, any
   `chi_S(h) = [(x(h), y(h)) mod m in S] - #S/m^2`, `S` a set of residue vectors modulo `m`.
   The hypothesis `b >= a` only fixes the constants: for all `a > 0`, `b >= 0` the same proof gives a bounded number
   of cut top words, because `rho >= 1/2`. "Uniformly" refers to the region: the constants depend on the modulus
   (`C_m` comes from Gelfand's formula and is not explicit, and `theta_m` depends on the rate at `m`). Nothing here
   is uniform in the modulus; a uniform gap is HYP-9660 (B), OPEN.

Numerical rates `max_k rho(T_k)/2` (VERIFIED, floating point; the receipt for every rate is
rates_exact/rates_by_modulus.txt, 262 moduli: every odd `m <= 521`, and `1023`, `1143`; `1039` is in
bridges/b3_f1_data.out): `0.8660` (m = 3), `0.7418` (5), `0.9320` (7), `0.8587` (11), `0.8272` (13), `0.8289` (17),
`0.8053` (19), `0.7931` (23), `0.9002` (31), `0.8342` (43), `0.8516` (73), `0.8754` (127); between `0.769` and
`0.806` for the other primes up to 89 (least `0.7691` at 79, largest `0.8060` at 37; twisted_spectrum.out has no
row for 71, 79, 83). Over the whole table no rate exceeds `0.932039`, which is attained exactly at the multiples
of 7. The scan for `q <= 120`, the extremality of `1/7` and the next-worst list `11/93`, `1/31`, `1/15`, `1/63` are
KNOWN (S82 note section 7); new here is the extension to every odd `m <= 521` and `1023`, `1039`, `1143`
(HYP-9660 (B)). `theta_p / 2` can be taken as close as wished to, but below, `0.0733, 0.1157, 0.0422` for
`p = 3, 5, 7`: these are the values at `rho_p` equal to the rate, which (B1) excludes (audit,
proof-5128-AB/check_B1.out; regions_check.py takes the rate plus `0.005`, that is `theta_p = 0.1425, 0.2286,
0.0790`).

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
power-saving error (exponent `theta_m / 2`, by (B4)(iv)), `F^(d)` being the share of triples whose smaller state
has `gcd(g, m) = d`. Odd prime-power moduli are not covered by (C): see (D3).
(VERIFIED, first_order_exact.out and first_order_exact_28.out, cumulative `X = 2^18, ..., 2^28`; every value below
was reproduced by the audit, recompute-5128/a3_analysis_28.out and a12_compare.out, up to the last digit of one
value, noted below. Two forms are compared: the statement (C2) itself, with `F_small` (the smaller state only;
column `R_small` of the outputs), and the two-state form of (C3), with `2f - 1/p^2` (the measured divisible share
of both states; column `R_sym`). For (C2) the remainder times `p`, sheets 0 / 1, is `0.0019 / 0.0021` at `p = 5`,
`0.003 / 0.002` at `17`, `0.004 / 0.006` at `13`, `0.018 / 0.011` at `3`, `0.008 / 0.007` at `31` and
`0.117 / 0.116` at `7` for `X = 2^26`, and `0.00139 / 0.00185`, `0.0016 / 0.0022`, `0.0026 / 0.0026`,
`0.017 / 0.015`, `0.006 / 0.005`, `0.079 / 0.093` at `X = 2^28`; at `p = 5`, sheet 0, it runs
`0.0050, 0.0038, 0.0005, 0.0016, 0.0019, 0.0014` over `X = 2^18, ..., 2^28` (the session's script counts the
diagonal `h1 = h2` in `F_small`; with `D_sigma` as in (C2) the first value is `0.0051`, as in the audit's output,
and no other value quoted here changes): not monotone, and 12 and 18 times the sampling noise `sqrt(p/#T)` at
`2^26` and `2^28` (`0.00016`, `0.00008`). For the two-state form it is
`0.0003 / 0.0002`, `0.002 / 0.002`, `0.004 / 0.004`, `0.009 / 0.007`, `0.011 / 0.008`, `0.059 / 0.044` at
`X = 2^26` and `0.00013 / 0.00014`, `0.0012 / 0.0014`, `0.0026 / 0.0019`, `0.014 / 0.008`, `0.009 / 0.006`,
`0.053 / 0.049` at `X = 2^28`: decreasing with `X` on the whole, with oscillations at `3` and `7`. The two forms
differ by the discrepancy of the larger state, whose divisible share at `X = 2^26` is still `0.987 / p^2` at `5`,
`1.47 / p^2` at `7` (sheet 0; `1.59 / p^2` on sheet 1) and `0.74 / p^2` at `31` (sheet 0; `0.66 / p^2` on sheet 1).
So (C2) is the less accurate form at `5`, `3`, `11` and above all `7`, while at `31` it is the two-state form that is
less accurate on sheet 0 (`0.011` against `0.008`): an unmixed larger state does not decide which form is better.
For squarefree `m` the session committed no check; the audit's are recompute-5128/a11_squarefree.out
(`m = 15, 55, 85, 165`: `m max_c |share - first-order law|` between `0.006` and `0.038` at `X = 2^26`) and
proof-5128-CF/a4_moduli_26.out (`m = 15, 21`). The regions lemma (B4)(ii) itself is checked on random regions in
regions_check.py.)

**(D) Mass, limit, and the `1 / ln X` term.**
1. `#T_sigma(X) = kappa_sigma X ln X + O(X)`, `kappa_0 = l_0 - 2 l_1 + l_2 = pi^2/12 + 2 Li_2(-1/2) - Li_2(-1/4) = 0.1615389173`,
   `kappa_1 = l_1 - l_2 = Li_2(-1/4) - Li_2(-1/2) = 0.2125139092`. For one padding `s`:
   `#{(h1, h2) : M <= X} = (l_s - 2 l_(s+1) + l_(s+2)) X ln X + O((s + 1) X 2^(-s))`. So the mean of `r` over odd `N` is
   `g_bar ln N + O(1)`, `g_bar = 2 (kappa_0 + kappa_1) = pi^2/6 + 2 Li_2(-1/2) = 0.7481056530`, and the sheets have
   limiting shares `w_0 = 0.43186124`, `w_1 = 0.56813876`. With the weight `s + 1` (word pairs) the total constant is
   `pi^2/12`, in agreement with THM-4885, and the shares are `w_0 = kappa_0 / (pi^2/12) = 0.19640777`, `w_1 = 0.80359223`.
2. `D_sigma(X) = #T_sigma(X) / p^2 + 2 A_p^sigma X + O_p(X^(1 - theta_p/2))`, where
   `A_p^sigma = lim_{H -> infinity} sum_{h < H} chi(h) w^sigma(h)`,
   `w^0(h) = (1/x) ln(2 (x + y)/(x + 2y))`, `w^1(h) = (1/x) ln(1 + x/(2y))`. The series is not absolutely convergent
   (`sum_{h < H} |chi(h)| w^sigma(h)` grows like `(2 (p^2 - 1)/p^4) kappa_sigma ln H`); the limit exists because the
   sums over `B_n ∩ [0, H)` are `O_p(2^(-theta_p n))` uniformly in `H`, by (B4)(iii), so that the block sums decay
   geometrically.
3. Hence, with `mu_c = 1/p + (1/p^2)([c = c_sigma] - 1/p)`:

       #{N = c} / #T_sigma(X) = mu_c + ([c = c_sigma] - 1/p) (2 A_p^sigma / kappa_sigma) / ln X + O(1 / ln^2 X).

   The by-value class factors modulo an odd prime are `p mu_c = 1 - 1/p^2 + [c = c_sigma]/p`, and for squarefree odd
   `m` the limiting density of a class modulo `m` is the product of these (PROVED). The product

       S^sf_sigma(N) = prod_{p odd} (1 - 1/p^2 + [p | N + [sigma = 0]]/p)
                     = (8/pi^2) prod_{p | N + [sigma = 0], p odd} (1 + p/(p^2 - 1))

   is the squarefree-level product (a class average): the limit, along the products `m` of the first odd primes, of
   `m` times the limiting density of the class of `N` modulo `m`. It is NOT the product of the p-adic densities at
   `N`: inside the distinguished class modulo `p` the classes modulo `p^2` are not equally likely (data, audit,
   `N <= 2^22`, sheet 1, `3 | N`: shares `0.377, 0.315, 0.309` of `N = 0, 3, 6 (mod 9)`;
   citations/check_singular_series_prime_powers.out). For independent equidistributed states modulo `p^k`
   the value of the form `a1 a2 + u b1 b2` (`u` a unit) has, at a class `c` of valuation `v = v_p(c) < k`, the factor

       p^k P(c) = (1 - p^-2) sigma_1(p^v)/p^v = (1 + 1/p)(1 - p^(-v-1))

   (an exact elementary count, see the proof; enumerated at 9, 25, 27 and 49 by the audit,
   proof-5128-CF/a15_padic_factor_exact.out and citations/check_singular_series_prime_powers.out). It is `1 - 1/p^2`
   at `v = 0`, below the class factor by `p^-3` at `v = 1` and above it for `v >= 2`; its mean over `v >= 1` is the
   class factor. For F1, with `v = v_p(c - c_sigma)`, the by-value limit modulo `p^k`, `k >= 2`, is NUMERICAL: audit
   census at `X = 2^26`, modulo 9, sheet 1: `9 * share = 0.848, 1.229, 1.456` for `v = 0`, `1`, `2 or more`, against
   the limits `0.889, 1.185, 1.296`; the first-order law at 9, with the measured shares `0.1215`, `0.0268` of smaller
   states with `3 || g`, `9 | g`, predicts `0.852, 1.216, 1.457` (proof-5128-CF/a4_moduli_26.out, also for 25 and
   27). It would follow from the argument of (C), (D2) carried out with (B4)(iv) for the modulus `p^k` (a smaller
   state with `p^i || g` spreads its mass uniformly over the classes with `p^i | c - c_sigma`), which is not written
   out here. (Addendum 2026-11-11: it is written out in THM-5134 (A)-(E), for every odd modulus; the statement
   below is THM-5134 (E).) The singular series proper would then be

       S_sigma(N) = prod_{p odd} (1 + 1/p)(1 - p^(-v_p(n)-1)) = (8/pi^2) sigma_1(n)/n,

   `n` the odd part of `N + [sigma = 0]` (mean 1, minimum `8/pi^2`): `8/pi^2` times the abundancy index of the odd
   part. Either product is at least `8/pi^2` of its mean for every odd `N`. Over both sheets the density is
   `w_0 S_0 + w_1 S_1` (with the triple shares `w_0 = 0.43186`, `w_1 = 0.56814` of (D1)), which is not an Euler
   product; at one prime, the boost `1/p` is divided between the classes
   `-1` and `0` in the proportions `(w_0, w_1)`.
4. (NUMERICAL. Partial sums over the states below `2^26`: ap_constants.out, blind_proof/check_blind_proof.out and,
   independently, the audit, proof-5128-CF/a10_ap_26.out. Continued by the audit to `2^28` over all states
   (recompute-5128/a5_Ap.out) and over the p-divisible states to `2^36`, `2^38` or `2^40`
   (recompute-5128/a5b_Ap_extended_*.out).) `A_p^0, A_p^1`:

        p     states below 2^26       audit continuation
        3     0.0556, 0.0789          0.0553, 0.0785 at 2^36; swings of +-0.003 since 2^16 (sheet 0: 0.0597,
                                      0.0539, 0.0556, 0.0563 at 2^16, 2^22, 2^26, 2^28)
        5     0.02943, 0.04250        0.02943, 0.04250 at 2^36
        7     0.0058, 0.0066          0.0006, -0.0002 at 2^40: NOT determined, sign included
        11    0.01096, 0.01449        0.0110, 0.0146 at 2^38
        13    -0.00377, -0.00488      -0.00377, -0.00488 at 2^38
        17    0.01355, 0.01931        0.01357, 0.01933 at 2^38
        31    -0.00194, -0.00261      -0.0019, -0.0026 at 2^40
        43    0.00235, 0.00323        0.00235, 0.00323 at 2^28 (all states; not continued further)

   `A_7` is not determined by the states below `2^26`: the partial sums are `0.0058, 0.0066` there and still falling
   (on sheet 0 they rise from `0.0014` at `2^14` to `0.0080` at `2^21` and fall again); continued over the
   7-divisible states to `2^40` they are `0.0006, -0.0002`, where the block terms change sign
   (recompute-5128/a5b_Ap_extended_p7_40.out). The oscillation decays at the rate `0.932` with a half-period of
   about 18 binary digits (argument of the eigenvalue `0.0282 * 2 pi`, recompute-5128/a1b_rates_detail.out). The
   sign of `A_7` is undetermined.
   At `3` and at the primes where the partial sums have settled (`5, 11, 13, 17, 31, 43`) the sign of `A_p` is that
   of the measured excess `F_small - 1/p^2`: positive at `3, 5, 11, 17, 43`, negative at `13` and `31` (at `31`
   trivially, see below), where by (D3) the distinguished class approaches its limit from below. At `13` the sign
   is negative although the first divisible state `h = 32` is early (the next is `h = 721`): an early divisible
   state does not decide the sign. That `A_p != 0` is numerical.
   Against the census: `(F_small - 1/p^2) #T / (2X)` at `X = 2^28` is `0.0300, 0.0434` (p = 5), `0.0138, 0.0197` (17),
   `0.0118, 0.0155` (11), `-0.0035, -0.0044` (13): about 2 % from `A_p` at 5 and 17 and 6 to 9 % at 11 and 13. By the
   proof of (D2) the census value at `X` is the partial sum of the series over the smaller states with
   `y1 < (X/2^s)^(1/2)`, less a boundary term, so the approach is slow and need not be monotone: the gap decreases
   from `2^26` to `2^28` at 5, 13 and 17, while at 11 the measured value has overshot (`0.0155` against `0.0146` on
   sheet 1) and is not yet returning (exact census values: recompute-5128/a3_analysis_28.out,
   numbers-vs-outputs/verify_fsmall_28.out; the column of ap_constants.out is computed from three-decimal inputs).
   At `p = 31` no smaller state can be divisible below `X = 17329^2 = 3.0e8` (first divisible state `h = 8664`,
   `y = 17329`), so `D = 0` in the whole census range and the value `-0.0017, -0.0021` at `2^28` is `-#T/(2 p^2 X)`:
   it is not a test of `A_31`. At `p = 7` the census value is `-0.005, -0.010` at `2^28` and rising (rate `0.932`); it
   does not point to a target.

**(E) The prime 2.** Modulo `2^j` the increment `ceil(h/2)` is a function of `w = h mod 2^(j+1)`, and for a character
`k` of `Z/2^j` the twisted matrix is `T_k[w -> 2w] = 1`, `T_k[w -> 2w + 1] = e(k c(w) / 2^j)`,
`c(w) = (w + (w mod 2))/2`. For every odd `k`, `T_k T_k^* = 2 I`: `T_k / sqrt 2` is unitary, and
`|sum_{h < 2^n} e(k q(h)/2^j)| <= 2^((n + j + 1)/2)` (for even `k` the identity fails; such a `k` is a character of
`Z/2^(j-1)`). In particular `q(h) = sum_i h_i (h_(i+1) + h_(i+2)) (mod 2)` in
the binary digits of `h`, a quadratic form of Rudin-Shapiro type, and over `h < 2^n` (index order) `(x(h), y(h))` is
equidistributed on the pairs of odd residues modulo `2^j`, up to `O_j(2^(n/2))`. Consequently, on a box
`h1, h2 < 2^n` with `s` fixed, every odd class modulo `2^j` receives `2^(2n+1-j) + O_j(2^(3n/2))` triples.
By value the same is true (statement and proof supplied by the independent audit): for every `j >= 1`, every sheet
and every odd `c`,

    #{t in T_sigma(X) : N = c (mod 2^j)} = 2^(1-j) #T_sigma(X) + O_j(X^(7/8)).

F1 has no local factor at 2. By value the states themselves are NOT equidistributed modulo `2^j`: the bottom digits
of `h` fix both `y mod 2^j` and the size of `rho`. For instance the share of the triples whose larger word is odd
tends to `(Li_2(-1/2) - Li_2(-1/4) + Li_2(-3/8) - Li_2(-3/4)) / kappa_0 = 0.5296` on sheet 0 and to
`(Li_2(-1/4) - Li_2(-3/8)) / kappa_1 = 0.5120` on sheet 1, not to `1/2`. What is uniform is `x` given `y`.
(VERIFIED by the audit, proof-5128-CF/a6_two_adic_by_value_26.out: at `X = 2^26` the largest relative deviation of an
odd class modulo `4, 8, 16, 32` is `0.0002, 0.0003, 0.0005, 0.0013` on sheet 0 and `0.00004, 0.00008, 0.0003,
0.0007` on sheet 1, and the larger word is odd in the shares `0.53345`, `0.51418`; limits a13_parity_limit.out; for
a fixed smaller state, a16_conditional_uniformity.out. The session's own check, twisted_spectrum_2adic.out, is on
the boxes `[0, 2^11)^2` only.)
*The mod-4 imbalance, exactly.* The polar form of `q mod 2` on words with `n` digits is the `n x n` band matrix `Q_n`
over `F_2` with ones at distances 1 and 2 from the diagonal, a Toeplitz matrix with symbol
`z + z^2 + z^(-1) + z^(-2) = z^(-2) (1 + z)^2 (1 + z + z^2)`. Its corank is `1, 0, 1, 2, 1, 0` for
`n = 1, 2, 3, 4, 5, 0 (mod 6)` (the corank is periodic, not the rank), and `S_n = sum_{h < 2^n} (-1)^q(h)` and
`D_n = #{h < 2^n : x(h) = 1 mod 4} - #{h < 2^n : x(h) = 3 mod 4}` are `0` or `+- 2^(n - rank/2)`. The signed transfer
matrix `M` of `q mod 2` on the pairs `(h_0, h_1)` of last digits (`4 x 4`; in (E) and its proof `M` is this matrix,
not the value `N + [s = 0]` of the Setting) satisfies `M M^T = 2I` and `M^12 = -64 I`, so `S_(n+12) = -64 S_n` and
`D_(n+12) = -64 D_n` for every `n >= 0` (the audit's observation, proved below). With
`(D_0, ..., D_11) = (1, 0, 2, 4, 0, 0, 8, 0, -16, 0, 0, -64)` and
`(S_0, ..., S_11) = (1, 2, 2, 0, 0, 0, -8, -16, -16, -32, -64, -64)`: `D_n = 0` exactly for
`n = 1, 4, 5, 7, 9, 10 (mod 12)` and `S_n = 0` exactly for `n = 3, 4, 5 (mod 12)` (for example `D_18 = -512`,
`D_20 = 1024`, `D_27 = 16384`; mod4_symbol.py, exact for `n <= 48`; proof-5128-CF/a8_period12.out, a5_two_adic.out).
This is the chain case of THM-5127 (C) made explicit; the factor `1 + z + z^2` is the polynomial of
THM-5125 (4b).

**(F) Divisible states.** `x(h) = y(h)` if and only if `h = 0` (`x = y = 1`) or `h = 2^j` (`x = y = 2^(j+1) + 1`).
This is KNOWN: the only unit-offset one-apex roots are the fans and the transitive tournament rooted at its sink
(THM-5126 (A); unit offset means `q(h) = 0`), and `2^m + 1` is the count of one back arc (THM-4748, THM-4861).
Hence every prime factor `p` of `2^m + 1` has the p-divisible state `h = 2^(m-1)`:
`3` at `h = 1`, `5` at `2`, `17` at `8`, `11` at `16`, `13` at `32`, `43` at `64`, `257` at `128`, `65537` at `2^15`.
Every prime with `ord_p(2)` even divides some `2^m + 1`, but the state `2^(m-1)` (least `m`) is below the generic
scale `p^2` only for small `m`: among the primes below 300, only for `3, 5, 11, 13, 17, 19, 41, 43, 113, 241, 257`;
for `29` and `61` even the first divisible state (`1435`, `5886`) is beyond `p^2` (audit,
proof-5128-CF/a7_divisible.out). So an early state `h = 2^(m-1)` is a property of these few primes, not of the prime
factors of `2^m + 1` in general. The density of p-divisible states is `1/p^2`, and the density of the `h` with
`g(h)` coprime to a given squarefree odd `m` is `prod_{p | m} (1 - 1/p^2)`.
(FINITE-EXACT: first p-divisible states `7: 38`, `19: 218`, `23: 103`, `29: 1435`, `31: 8664`, `37: 1350`,
`41: 266`, `127: 109283`; proportion of `h < 2^24` with `g(h) = 1`: `0.81176`, against `8/pi^2 = 0.81057`; that the
density of `g = 1` equals `8/pi^2` is NOT proved, the tail of the product being uncontrolled.)

**(G) Windows (NUMERICAL; pre-registered).** By (B3) (the uniform law of `rho`) the number of triples with first
state `h1` and padding `s` whose value lies in `[A, B)` is approximately
`W(h1, s) = (1/x1) int_{f0}^{f1} w(f) df`, `f0 = x1/2 + 2^s y1`, `f1 = x1 + 2^s y1`,
`w(f) = (B - A)/f` for `f < A`, `B/f` for `A <= f < B`. With `f` the `W`-weighted share of p-divisible first states
and `D = 2f - 1/p^2`, the predicted excess of the distinguished class over the mean of the other classes is
`E = p D / (1 - D)`, and after removal of the p-divisible triples `E_nd = p (f^2 - 1/p^2)/(1 - D)`. (The method is
that of the S82 prediction for the bank family F2 at `p = 3`: S82 note section 5.)
On `[2^24, 2^27)` (unseen; 8.3e8 triples): total mass `14.1405` per odd `N` against `14.1417` predicted (per sheet
`6.2379` against `6.2819` and `7.9026` against `7.8598`: the rule is off by `+0.7 %` and `-0.5 %` on the sheets, a
constant offset of `+0.044` and `-0.043` triples per odd `N` seen on every window, which cancels in the total);
weights of `h1 = 0, 1, 2, 3, 8, 16, 64, 128, 1024` within `0.006 %`; `E` measured / predicted on sheet 0 and sheet 1:
`p = 5`: `+0.3130 / +0.3120`, `+0.3280 / +0.3284`; `11`: `+0.1690 / +0.1705`, `+0.1707 / +0.1749`;
`17`: `+0.2107 / +0.2081`, `+0.2296 / +0.2288`; `19`: `+0.0510 / +0.0511`, `+0.0495 / +0.0497`;
`29`: `+0.0051 / +0.0032`, `+0.0020 / +0.0007`; `43`: `+0.0885 / +0.0877`, `+0.0930 / +0.0941`;
`257`: `+0.1859 / +0.1815`, `+0.2023 / +0.2044`; `2731`: `+0.063 / +0.061`, `+0.075 / +0.069`;
`43691`: `+0.018 / +0.061`, `+0.046 / +0.069` (standard error `0.018 / 0.017`: 2.4 and 1.4 standard errors below the
prediction; inside the registered intervals `[-0.046, +0.168]` and `[-0.039, +0.177]`, which contain 0, so this
prime is not resolved; likewise `241`: `+0.012 / +0.014`, `+0.016 / +0.015`);
`65537`: `+0.144 / +0.183`, `+0.168 / +0.206` (standard error `0.022 / 0.021` from the measured dispersion, 1.75 and
1.82 standard errors below the prediction; the registered tolerance used `0.029`); controls `65521`, `65539`,
`8191`, `2729`, `2741`, `251`, `263`, `677`, `691` within noise of zero. At the slow primes the second-order term is
visible: `p = 3`: `+0.539 / +0.515`, `+0.544 / +0.532`; `p = 7`: `+0.233 / +0.174`, `+0.230 / +0.175`. All eight
pre-registered tests (279 intervals) CONFIRMED. The same rules also hold on `[2^23, 2^24)` (not pre-registered).
The registered tolerances used `se = sqrt(0.75 p / #N)`. The constant `0.75` was a conservative bound, not a
measurement: no output shows it, and the measured relative variance on this window is `0.45 / 0.39` (Remark 3), so
the registered "4 se" bands were about `5.2` to `5.5` measured standard errors wide. All 279 items are also inside
the bands recomputed with the measured standard errors, at 4, 3 and 2 standard errors (audit,
recompute-5128/a10_prereg_full.out). The registered files PREREGISTRATION.md and window_split.py are unchanged.

## Proofs

**(A).** Write `Dhat(xi) = sum_v D(v) e(xi . v / p)` and `B(u, v) = (B* u) . v` with `B*` invertible. Then
`P(B = c) = (1/p) sum_t e(-tc/p) Phi(t)`, `Phi(t) = sum_u D1(u) D2hat(t B* u)`, `Phi(0) = 1`. For `t != 0` separate
`u = 0` (it contributes `pi_1`) and put `D1 = p^(-r) + delta_1` on `u != 0`. Since `u -> t B* u` permutes the
non-zero vectors and `sum_xi D2hat(xi) = p^r D2(0)`,

    Phi(t) = pi_1 + p^(-r) (p^r pi_2 - 1) + R(t) = pi_1 + pi_2 - p^(-r) + R(t),       R(t) = sum_{u != 0} delta_1(u) D2hat(t B* u).

Summing over `t` with `sum_{t != 0} e(-tc/p) = p [c = 0] - 1` gives (1) with
`Rhat(c) = (1/p) sum_{t != 0} e(-tc/p) R(t)`. `sum_c Rhat(c) = 0` because `sum_c e(-tc/p) = 0` for `t != 0`. If `D2`
is uniform, `D2hat` vanishes off `0`; if `D1` is uniform, `delta_1 = 0`. By Cauchy-Schwarz and Parseval
(`sum_{xi != 0} |Dhat(xi)|^2 = p^r |D - U|_2^2`),
`|R(t)| <= |delta_1|_2 (sum_{xi != 0} |D2hat(xi)|^2)^(1/2) <= p^(r/2) |D1 - U|_2 |D2 - U|_2`
(the norm of `delta_1` is over `u != 0`, hence at most `|D1 - U|_2`),
which gives (2); the second form follows from `|D - U|_2^2 <= E^2`. For a point mass at `u_0 != 0` the event
`B(u_0, v) = c` is `v` in an affine hyperplane, and `sum of D2 over it = 1/p + (1/p) sum_{t != 0} e(-tc/p) D2hat(t B* u_0)`.
*The constant for `c != 0`.* `Rhat(c) = sum_{w != 0} g(w) K(w)` with `g(w) = delta_1(B*^(-1) w)` and
`K(w) = (1/p) sum_{t != 0} e(-tc/p) D2hat(t w)`. On a punctured line `{s w_0 : s != 0}` put `a_t = D2hat(t w_0)`; then
`K(s w_0) = (1/p) sum_{t != 0} e(-t (c/s)/p) a_t` (`c/s` in `Z/p`), and `c/s` runs over the non-zero residues with
`s`, so by Parseval on `Z/p`, `sum_{s != 0} |K(s w_0)|^2 <= (1/p) sum_{t != 0} |a_t|^2`. Summing over the lines,
`|K|_2^2 <= (1/p) sum_{xi != 0} |D2hat(xi)|^2 = p^(r-1) |D2 - U|_2^2`, and Cauchy-Schwarz gives
`|Rhat(c)| <= p^((r-1)/2) |D1 - U|_2 |D2 - U|_2`.
*Equality in (2) at `c = 0`, `r = 2`.* Let `L != L'` be two lines through the origin, `Lperp`, `L'perp` their
`B`-orthogonals (`Lperp = {v : B(u, v) = 0 for u in L}`), `0 < t <= p^(-2)`, `D1 = U + t (1_L - 1_L')` and
`D2 = U + t (1_Lperp - 1_L'perp)` on the non-zero points. Then `pi_i = p^(-2)`, the cross terms vanish, and
`Rhat(0) = sum_{u, v} delta_1(u) delta_2(v) [B(u, v) = 0] = 2 t^2 (p - 1)^2`, which is the bound
`(p - 1) |D1 - U|_2 |D2 - U|_2 = (p - 1) * 2 t^2 (p - 1)` (exact check at `p = 3, 5, 7, 11`:
proof-5128-AB/check_extra.out; the sharp constants for `r <= 3`: check_A.out).

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
ratio `zeta^(l 2^j) != 1` (`j` the number of prescribed bottom digits), of modulus at most `m/2`. Fourier inversion
on `(Z/m)^2` gives the claim for `(q, h)`, and `(x, y) = (2h + 1 - 2q, 2h + 1)` is an affine bijection of `(Z/m)^2`.

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
(Without `a <= b` the factor is at most `1 + (2^(-j-1) + 2 eps)/(1/2 - eps)`, because `R - 2^(-j-1) >= 1/2`: the
number of cut top words stays bounded for all `a > 0`, `b >= 0`.)
So at most `4 * 2^j` cells are cut by `Omega` and the cells of at most two top words by `I`: at most
`6 * 2^(n-j)` points. (i): take `j = floor(n/2)` (for `n = 1` there is nothing to prove); a cell inside (outside)
`Omega` has its rectangle inside (outside) the region and contains exactly `2^n` times its area. (ii): take
`j = floor(theta_p n)`; if `j = 0` the claim is trivial, the error term being larger than `#B_n`; otherwise by (B2)
the uncut cells contribute an error at most `C 2^n rho_p^(n - 2j)`, and `rho_p^(1 - 2 theta_p) = 2^(-theta_p)`.
(iii): the same, freezing `W` on each cell (again trivial if `j = 0`). (iv): (B2) is stated for every odd `m`, and
nothing else in the proof of (ii) depends on the modulus; in (iii) the sum of `chi_S` over an uncut cell is
`O_m((2 rho_m)^r)`, by (B2) summed over `S`.

**(C1).** Here `M = a x2 + b y2` with `a = x1 <= b = 2^s y1`. If `p` does not divide `g(h1)` the linear form
`w -> x1 w_1 + 2^s y1 w_2` on `(Z/p)^2` is non-zero and takes each value on `p` residue vectors; sum (B4)(ii) over
these vectors and over the blocks `B_n` with `I = (h1, infinity)`; `M > y1 y2 2^s` bounds the blocks by
`2^n <= Y/y1`. If `p | g(h1)` both coefficients vanish modulo `p`.

**(C2).** `M` is symmetric in `(h1, h2)`, so each triple with `h1 != h2` is counted through its smaller state, and
those with `h1 = h2` number `O(Y^(1/2))` for each `s`. `G = 0` unless `y1 < Y^(1/2)`. Summing (C1) over the smaller
states, the errors add up to `sum_{y1 < Y^(1/2)} (Y/y1)^(1 - theta) = O(Y^(1 - theta/2))`; summing over
`1 <= s <= log_2 X` on sheet 1 keeps `O(X^(1 - theta/2))`. On sheet 0, `N <= X` is `M <= X + 1`. (C3): by (B4)(ii)
with `w = 0` the larger state is p-divisible in a share `1/p^2` of the triples up to the same error, and
`2f = F_small + F_large`. The squarefree case is the same proof with (B4)(iv) for `m`: for a smaller state with
`gcd(g(h1), m) = d` the form `w -> x1 w_1 + 2^s y1 w_2` vanishes modulo every `p | d` and is non-zero modulo every
`p | m/d`, so by the Chinese remainder theorem it takes each value that is `0` modulo `d` on `m d` residue vectors
modulo `m`, and no other value. Hence `M` is uniform on the classes modulo `m` that are `0` modulo `d`, up to
`O_m((Y/y1)^(1 - theta_m))`, and the sum over the smaller states is as before.

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
`y` and `rho`, and `Xi_X = int_0^1 min(Z, max(h1, 1))`). So the sums of `chi w_s` over `B_n ∩ I`, for every interval
`I`, are `O(2^(-s) 2^(-theta n))`: the series `A_(p,s) = sum chi w_s` converges when summed over initial segments
(not absolutely: `|chi|` has mean `2 (p^2 - 1)/p^4` and the block sums of `w_s` tend to
`ln 2 * (l_s - 2 l_(s+1) + l_(s+2))`), its tail beyond `Y^(1/2)` is `O(2^(-s) Y^(-theta/2))`, and
`sum_{y1 < Y^(1/2)} chi Xi_X = O(Y^(1 - theta/2))`. Summing over `s` gives the statement with
`w^1 = sum_{s >= 1} w_s = (1/x) ln(1 + x/(2y))`. **(D3)** follows from (C2), (D1), (D2). For squarefree `m`, (D2)
holds with `chi` replaced by `chi_d(h) = [gcd(g(h), m) = d] - c_d`,
`c_d = prod_{p | d} p^(-2) prod_{p | m/d} (1 - p^(-2))` (same proof, with (B4)(iii) in the form (iv): `gcd(g, m) = d`
is a condition on `(x, y)` modulo `m`, met by `c_d m^2` residue vectors). So `F^(d)(X)` tends to `c_d` at the rate
`1/ln X`, and the squarefree case of (C) gives the limit
`sum_{d | m} c_d prod_{p | d} [c = c_sigma (mod p)] prod_{p | m/d} (1/p) = prod_{p | m} mu_c^(p)`, `mu_c^(p)` the
`mu_c` of the prime `p`: the product law.
*Uniform states modulo `p^k`.* Let `(a1, b1)`, `(a2, b2)` be independent and uniform modulo `p^k` and `u` a unit. With
probability `p^(-2i) (1 - p^(-2))` the first state has `p^i || gcd(a1, b1)`, `i < k`, and then
`(a2, b2) -> a1 a2 + u b1 b2` is uniform on the multiples of `p^i` modulo `p^k` (it is `p^i` times a surjective
homomorphism onto `Z/p^k`). So for `v = v_p(c) < k`:
`P(c) = sum_{i <= v} p^(-2i) (1 - p^(-2)) p^(i-k) = p^(-k) (1 - p^(-2)) sigma_1(p^v)/p^v`. The mean of
`(1 + 1/p)(1 - p^(-v-1))` over `v >= 1` (weights `(1 - 1/p) p^(1-v)`) is `1 + 1/p - 1/p^2`, and over `v >= 0` it is 1.

**(E).** The rows `w` and `w + 2^j` of `T_k` have the same two non-zero columns, with entries `(1, omega)` and
`(1, (-1)^k omega)`, `omega = e(k c(w)/2^j)`, because `c(w + 2^j) = c(w) + 2^(j-1)`; for odd `k` they are orthogonal,
rows with different column pairs have disjoint supports, and every row has norm `sqrt 2`. So `T_k T_k^* = 2I`, and
`|e_0^T T_k^n 1| <= 2^(n/2) 2^((j+1)/2)`. For even `k` the two rows coincide, so the identity fails; an even `k` is a
character of `Z/2^(j-1)`. `ceil(H/2) = H_0 + H_1 (mod 2)` in the binary digits of an integer `H`, applied to the
`H_i = floor(h / 2^(i+1))` of the proof of (B3), gives the quadratic form.
*Index order (boxes).* Joint equidistribution with the last `i` digits of `h`: on the state space
`h mod 2^max(j+1, i)` let `v` be a function of the last `i` digits; `T_k^i v` depends on the state only modulo
`2^(j+1)` (the phases do, and the last `i` digits are the appended ones), so the remaining `n - i` steps are
governed by the unitary matrix above and the sum is `O(2^(n/2))`. Hence, over `h < 2^n`, `x = y - 2q` and
`y = 2h + 1` are independent uniform odd residues modulo `2^j` up to `O_j(2^(n/2))`. For fixed `h1`, `x1` being
odd, `x1 x2` is then a uniform odd residue given `y2`, and adding the even number `2^s y1 y2 - [s = 0]` leaves `N`
uniform on the odd classes; summing over `h1 < 2^n` gives the statement on boxes.
*By value (proof supplied by the independent audit, auditor proof-5128-CF).* The case `j = 1` is trivial (`N` is
odd); let `j >= 2`. Fix the smaller state `(h1, s)` and, for the larger word `h2`, a cell of `B_n` with `i` top and
`i` bottom digits prescribed, `j <= i <= n/2`, and `r = n - 2i` free digits. The bottom digits fix `y2 mod 2^j`. A
non-trivial character of `q mod 2^(j-1)` is an odd character `k` of some `Z/2^(j')`, `1 <= j' <= j - 1`, and the
cell sum of `e(k q(h2)/2^(j'))` is `(e_0^T T[top]) T_k^r (T[bottom] 1)`, with the matrices of the modulus `2^(j')`.
The digit matrices have one unimodular entry in every row, so `e_0^T T[top]` is a unimodular multiple of a
coordinate row vector and `T[bottom] 1` has unimodular entries; by unitarity the modulus of the cell sum is at most
`2^(r/2) 2^((j'+1)/2)`. Hence on each cell `q(h2) mod 2^(j-1)` is equidistributed up to `O_j(2^(r/2))`, and so is
`x2 = y2 - 2 q(h2)` on the odd residues modulo `2^j`, with `y2 mod 2^j` fixed. Since `x1` is odd and
`2^s y1 y2 - [s = 0]` is even and constant modulo `2^j` on the cell, `N = x1 x2 + (2^s y1 y2 - [s = 0])` is
equidistributed on the odd classes modulo `2^j` up to `O_j(2^(r/2))`, on each of the at most `2^(2i-1)` cells. With
the cut-cell count of (B4) (at most `6 * 2^(n-i)` points; here `a = x1 <= b = 2^s y1` and `I = (h1, infinity)`) and
`i = ceil(n/4)`, the count on `B_n` is `2^(1-j)` of the total up to `O_j(2^(3n/4))` (trivially so for `n < 4j`).
Summing over the blocks `2^n <= Y/y1`, over the smaller states `y1 < Y^(1/2)` and over `s`, as in (C2), gives
`O_j(X^(7/8))`, because `sum_{y1 < Y^(1/2)} (Y/y1)^(3/4) = O(Y^(7/8))`.
*States by value.* The larger word is odd exactly when its bottom digit is 1, that is `v >= 1/2` in (B4)(i). The
proof of (D1) with this restriction replaces `F_s` by `F_s^odd(rho) = (1/rho) ln((3 rho/4 + 2^s)/(rho/2 + 2^s))`,
and `int_{1/2}^{1} F_s^odd = Li_2(-2^(-s-1)) - Li_2(-2^(-s-2)) + Li_2(-3 * 2^(-s-3)) - Li_2(-3 * 2^(-s-2))`. For
`s = 0` this is `0.52956 kappa_0`; the sum over `s >= 1` telescopes to `Li_2(-1/4) - Li_2(-3/8) = 0.51198 kappa_1`.
*Period 12.* Appending a digit `d` to a word with last digits `(h_0, h_1) = (a, b)` multiplies `(-1)^q` by
`(-1)^(d (a + b))` and gives the last digits `(d, a)`. In the order `(0,0), (0,1), (1,0), (1,1)` of the pairs the
signed transfer matrix is

    M = [[1, 0, 1, 0], [1, 0, -1, 0], [0, 1, 0, -1], [0, 1, 0, 1]]        (rows: from; columns: to),

and `S_n = e_(0,0)^T M^n 1`, `D_n = e_(0,0)^T M^n v` with `v(a, b) = (-1)^a` (because `x = 1 (mod 4)` iff
`h_0 + q(h)` is even). The rows of `M` are orthogonal of norm `sqrt 2`, so `M M^T = 2I`. The characteristic
polynomial is `lam^4 - 2 lam^3 + 2 lam^2 - 4 lam + 4 = (lam^2 - (1 + sqrt 3) lam + 2)(lam^2 + (sqrt 3 - 1) lam + 2)`,
with the roots `sqrt 2 e(+-1/24)` and `sqrt 2 e(+-7/24)`. `M / sqrt 2` is orthogonal, hence diagonalisable, and every
eigenvalue satisfies `lam^12 = 64 e(+-1/2) = -64`; so `M^12 = -64 I` (also by direct multiplication: `M^6 = 8 J` with
`J^2 = -I`), and `S_(n+12) = -64 S_n`, `D_(n+12) = -64 D_n`. The twelve initial values come from the recursion.
*Corank of `Q_n`.* A vector `v` is in the kernel iff, extended by `v_(-1) = v_0 = v_(n+1) = v_(n+2) = 0`, it
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
   `2 A_p / (kappa ln X)` to first order: far above zero when a p-divisible state is very small compared with `p^2`
   (`3, 5, 17, 11, 43, 257, 65537`: prime factors of `2^m + 1` with small `m`, the "Fermat boost" of HYP-9506 (B)),
   and negative when the divisible states are late or sparse (`13`, `29`, `31`, `37`, `61`, `73`, `127`; of these
   `13, 29, 37, 61` also divide some `2^m + 1`, so that property does not decide the sign). In the second case the
   excess of the distinguished class is below its generic value `p/(p^2 - 1)`; the class itself is under-represented
   only when `2f < 1/p^2` (to first order): on `[2^24, 2^27)` at `73`, `127` and, on sheet 1, `31` (the S90 audit
   observed it at 73 and 127), while at `13, 29, 37, 61` the class is still over-represented (`+0.044`, `+0.005`,
   `+0.007`, `+0.001` on sheet 0). The other classes are flat up to the power-saving term, which is visible only at
   slowly mixing primes (`7`: rate `0.932`).
3. *What the law does not say.* (A) is a statement about one prime, and (C), (D) are statements about a fixed odd
   prime or a fixed squarefree odd modulus ((E) about the powers of 2): the power saving is uniform in the region
   but not in the modulus. Uniformity in the modulus is OPEN (it needs at least the uniform gap of HYP-9660 (B), with a norm bound
   polynomial in the modulus); odd prime-power moduli are numerical
   only ((D3)). That no congruence excludes any odd number, 7 and 21 included, was KNOWN before (the one-apex
   family alone attains every residue modulo every odd modulus and every odd residue modulo every even
   modulus: Session 25; Session193 modular_sync Theorem 1;
   Session206 exact-modular-image); what this theorem adds is the size of the local densities of F1 at fixed odd
   primes, squarefree odd moduli and powers of 2. A positive local density does not give `r(N) > 0`. The
   representation count is over-dispersed: index of dispersion `Var(r)/mean(r)` about `2.8` on the sheet `s = 0`
   and `3.1` on the sheet `s >= 1`, stable from `2^18` to `2^27`; relative variance `Var/mean^2` from `0.56 / 0.51`
   on `[2^18, 2^21)` down to `0.45 / 0.39` on `[2^24, 2^27)` (audit, recompute-5128/a7_misses.out). And F1 has
   misses at every scale tested (`0.35 %` of the odd `N` in `[2^24, 2^27)`). A second-moment (minor-arc) bound
   would give at most a density-one statement for F1; completeness needs a pointwise lower bound, which F1 alone
   cannot satisfy.
