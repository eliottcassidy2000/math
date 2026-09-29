---
id: THM-4515
title: "Fair consecutive splits of a cycle necklace are circulant: for a cycle of y -> y/2, (qy+d)/2 of shape (K,X) with a j-fold fair split into arcs of shape (m,x) = (K/j, X/j), the cut elements n_i and share carries c_i satisfy N_k (2^m zeta^k - q^x) = d C_k in Z[zeta_j] and the clock factors as prod_(e|j) Phi_e(2^m, q^x); a fair 2-split always exists (discrete IVT) and splits the integrality congruence by CRT into (2^m - q^x) | d(c_0 + c_1), (2^m + q^x) | d(c_0 - c_1); a fair j-split exists iff the discrepancy K F(t) - tX is constant on a coset of mZ/KZ, for one bead per thief iff every run of t cyclic gaps lies in ((t-1)m, (t+1)m), and exactly m(m-1) of the 3 C(m,2) primitive three-bead necklaces of length 3m split fairly; fair antipodal cut elements satisfy n'/n = sqrt(P_0/P_1)"
status: >
  PROVED (elementary) + FINITE-EXACT. NOT independently audited (single
  session; self-checks below). Setting as in THM-4484: T_{q,d}(y) = y/2 (y even),
  (qy+d)/2 (y odd), q >= 3 odd, d odd, gcd(q,d) = 1; a cycle with primitive parity
  necklace w of shape (K,X) has y_0 = d c_w/(2^K - q^X), c_w = sum_{w_j=1} q^{ones after j} 2^j.
  (1) Circulant: if j | gcd(K,X) and the necklace has a j-fold fair consecutive
  split at rotation r (windows W_i = [r+im, r+(i+1)m) each with x = X/j odd letters),
  then with n_i = y_{r+im}, c_i = c_{w|W_i}: 2^m n_{i+1} = q^x n_i + d c_i (mod j), and
  for zeta = exp(2 pi i/j), N_k = sum zeta^{-ik} n_i, C_k = sum zeta^{-ik} c_i:
  N_k (2^m zeta^k - q^x) = d C_k. In particular (2^m - q^x) sum n_i = d sum c_i and,
  for j even, (2^m + q^x) sum (-1)^i n_i = -d sum (-1)^i c_i.
  (2) 2^K - q^X = prod_k (2^m - zeta^k q^x) = prod_{e|j} Phi_e(2^m, q^x). For j = 2 the
  factors are coprime and (2^K - q^X) | d c_w iff (2^m - q^x) | d(c_0 + c_1) and
  (2^m + q^x) | d(c_0 - c_1) (c_w = q^x c_0 + 2^m c_1): the CRT form of THM-4484(1).
  (3) Existence: j = 2 always (h(r) = F(r+m) - F(r) - X/2 has h(r+m) = -h(r) and
  unit steps); general j iff G(t) = K F(t) - tX is constant on {r + im}; X = j iff every
  run of t consecutive cyclic gaps (1 <= t <= j-1) sums to a value in ((t-1)m, (t+1)m);
  X = j = 3, K = 3m iff every gap <= 2m - 1, giving m(m-1) primitive good necklaces out
  of 3 C(m,2) (exactly one third fail, for every m).
  (5) Expected fair cuts (uniformly random word of shape (jm,jx), N_j = number of fair
  cut positions in [0,m)): E[N_j] = m C(m,x)^j / C(jm,jx) exactly; Stirling gives
  E[N_j] ~ sqrt(j) m^((3-j)/2) (2 pi rho(1-rho))^(-(j-1)/2), rho = x/m. Hence: j = 2
  always (IVT) with E[N_2] ~ sqrt(m/(pi rho(1-rho))); j = 3: E[N_3] -> sqrt(3)/(2 pi rho(1-rho))
  (1.1027 at rho = 1/2) while P(N_3 >= 1) -> 0 (EMPIRICAL law ~ 2.2/log m to m = 10^5;
  PROVED lower bound E[N_3]^2/E[N_3^2] with the exact second moment
  E[N_3^2] = (m/C(3m,3x)) sum_delta sum_s C(delta,s)^3 C(m-delta,x-s)^3 ~ 0.93 log m: the
  imbalance pair is a planar lattice walk visiting the origin); j >= 4: P(N_j >= 1) <= E[N_j] -> 0
  like m^(-(j-3)/2) (PROVED). A Pólya trichotomy for the one-arc-per-thief necklace problem.
  (4) Fair Eliahou identity: for a fair 2-split with cut elements n, n' and share
  products P_i = prod_{odd y in share i} (1 + d/(qy)): n'/n = sqrt(P_0/P_1), so
  |log(n'/n)| <= X|d| / (2(q n_min - |d|)); fairness cancels the drift q^x/2^m exactly.
  FINITE-EXACT: (1)-(2) verified in Z[zeta_j] on 49 primitive 3x+d cycles (|d| <= 100,
  least element <= 100|d|, gcd(K,X) > 1), the nine 7x+169 cycles of shape (9,3), and
  the Belaga-Mignotte 3x+17021 cycle (least element 5, shape (2140,1088), gcd 4: 79 fair
  antipodal cuts of 1070, no fair 4-split); (3) cross-checked against brute force on
  all 23,201 (primitive necklace, j) pairs with K <= 18, zero mismatches; m(m-1) for m <= 7.
  Every known 3x+-1 cycle has gcd(K,X) = 1, so none admits a fair split; (2) is
  CRT-equivalent to THM-4484 and yields no new necessary condition on Collatz cycles.
  The mechanism (fairness = Fourier diagonalisation of the cycle equation) refines
  the typing of necklace splitting in collatz_generic_price_topology_20260927.md section 2.
source: collatz-necklace-20260929 session (mac-mini), 2026-09-29; owner seed: K-bead dual multi-edges, two-colour orderings, the stolen-necklace problem
depends_on:
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md (periodic point d c_w/(2^K - q^X), shift criterion)
  - 01-canon/theorems/THM-4471-kawasaki-fixed-point-collatz-proof-refuted.md (Banach periodic points)
related:
  - 05-knowledge/results/collatz_generic_price_topology_20260927.md (necklace splitting typed sheet-blind; refined here)
  - 01-canon/theorems/THM-4516-perfect-power-clocks-are-fermat-catalan-identities.md (the Eisenstein-square instance 7x+169)
  - 01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md (rotation averaging on necklaces)
note: 05-knowledge/results/collatz_necklace_20260929_fair_splits_power_clocks_basins.md
scripts:
  - 04-computation/experiments/collatz_necklace_20260929_circulant_check.py
  - 04-computation/experiments/collatz_necklace_20260929_fairsplit_census.py
  - 04-computation/experiments/collatz_necklace_20260929_fairsplit_asymptotics.py
  - 04-computation/experiments/collatz_necklace_20260929_fairsplit_loglaw.py
outputs:
  - 05-knowledge/results/collatz_necklace_20260929_circulant_check_d100.out
  - 05-knowledge/results/collatz_necklace_20260929_fairsplit_criteria_K18.out
  - 05-knowledge/results/collatz_necklace_20260929_fairsplit_asymptotics.out
  - 05-knowledge/results/collatz_necklace_20260929_fairsplit_loglaw.out
extra_sha256: df66480e5af561fa99477d85edf03d57b914899eb1be1df0848bc65aef2c0f29 (fairsplit_asymptotics.py), 6b02b02a4b8281dbd0773edb9cb000ad47679e12b6638ea268ac23883f413ae9 (fairsplit_loglaw.py), 27c68009696d3efc9080532fd01b298cadb0d9a58cf2784e0f3758673d77f11d (asymptotics.out), ae58295f8b1ec309c675ff0091e57e330bfc0bae17a5c08ed1d6c273abb0fa47 (loglaw.out)
script_sha256: 8e8272b3b2c5027127185721dc99d17ded7f25531eee1aba0542912666108fe0 (circulant_check), 116cefe9c4708cc866f7306e2bd57e28567ebf98b67106aabc3d0153887fe7cd (fairsplit_census)
output_sha256: b8885a9ad4819611a59f23affbf954d23adc2450353ba684781e6d374bca1d01 (circulant_check_d100), c9f4ec136dad704b382bdba6ecdd6326f3177fb91086bdd761f2cd1e4cbdfe5f (fairsplit_criteria_K18)
hash_basis: raw LF bytes
audit: NOT independently audited; self-checks are exact identities in Z[zeta_j] via reduction modulo Phi_j (sympy), brute-force cross-checks of both existence criteria, and the closed count m(m-1).
---

# THM-4515 -- fair consecutive splits of cycle necklaces are circulant

**PROVED (elementary) + FINITE-EXACT; not independently audited.** Full note: [collatz_necklace_20260929_fair_splits_power_clocks_basins](../../05-knowledge/results/collatz_necklace_20260929_fair_splits_power_clocks_basins.md), section 1.

## 1. Setting

`T_{q,d}(y) = y/2` (`y` even), `(qy+d)/2` (`y` odd). A cycle `y_0 -> ... -> y_{K-1} -> y_0` has a primitive parity necklace `w` of shape `(K, X)`, carry `c_w = sum_{w_j = 1} q^{(ones after j)} 2^j`, clock `D = 2^K - q^X`, and `y_0 = d c_w/D` (THM-4484). Let `j | gcd(K,X)`, `m = K/j`, `x = X/j`. A **`j`-fold fair consecutive split** at rotation `r` cuts the cycle into the `j` windows `W_i = [r+im, r+(i+1)m)`, each containing `x` odd letters (the two-colour stolen necklace with one arc per thief). Cut elements `n_i = y_{r+im}`, share carries `c_i = c_{w|W_i}`.

## 2. The circulant system and its Fourier transform

Each share is one affine map of shape `(m,x)`: `2^m n_{i+1} = q^x n_i + d c_i`, indices mod `j`. Equal shapes make the system circulant, `(2^m P - q^x I) n = d c`, `P` the cyclic shift, and the DFT diagonalises it: with `zeta = e^{2 pi i/j}`,

`N_k (2^m zeta^k - q^x) = d C_k`, `N_k = sum_i zeta^{-ik} n_i`, `C_k = sum_i zeta^{-ik} c_i`, `k = 0, ..., j-1`.

*Proof.* Multiply the `i`-th share equation by `zeta^{-ik}` and sum; the left side is `2^m zeta^k N_k` after re-indexing `i+1 -> i`. ∎

`k = 0` gives `(2^m - q^x) sum n_i = d sum c_i`; for `j` even, `k = j/2` gives `(2^m + q^x) sum (-1)^i n_i = -d sum (-1)^i c_i`.

## 3. Clock factorisation and the CRT form

`2^K - q^X = prod_{k=0}^{j-1} (2^m - zeta^k q^x) = prod_{e | j} Phi_e(2^m, q^x)` (the eigenvalues of the circulant; homogenised cyclotomic polynomials). For `j = 2`: `c_w = q^x c_0 + 2^m c_1`, the factors `2^m -+ q^x` are coprime (odd; a common divisor divides `2 q^x`), `q^x = 2^m mod (2^m - q^x)` and `q^x = -2^m mod (2^m + q^x)`, so

`(2^K - q^X) | d c_w  <=>  (2^m - q^x) | d(c_0 + c_1)  and  (2^m + q^x) | d(c_0 - c_1)`.

The integrality congruence of THM-4484 splits along any fair antipodal cut into a congruence modulo the clock of the average share `(m,x)` and one modulo the clock of the multiplier `-q` (when `x` is odd). This is what the two thieves' split carries; it is CRT-equivalent to the single congruence and gives no new necessary condition.

## 4. Existence of fair splits

* **`j = 2`, always.** `h(r) = F(r+m) - F(r) - X/2` (`F` the prefix count, `K, X` even) satisfies `h(r+m) = -h(r)` and `|h(r+1) - h(r)| <= 1`, hence vanishes somewhere.
* **General `j`.** A fair `j`-split at `r` exists iff `G(t) = K F(t) - tX` takes one value on `{r, r+m, ..., r+(j-1)m}` (window counts are the increments of `F`).
* **One bead per thief, `X = j`.** With cyclic gaps `g_1, ..., g_j` and `P_i = g_1 + ... + g_{i-1}`, the cut `r` must lie in every `(P_i - im, P_i - (i-1)m]`; these `j` integer intervals meet iff each lower end is below each upper end, i.e. iff every run of `t` consecutive gaps (`1 <= t <= j-1`) has sum in `((t-1)m, (t+1)m)` (the complementary-run condition is the same family).
* **Three beads, `K = 3m`.** Iff every gap is `<= 2m-1`. Compositions of `3m` into three parts `<= 2m-1` number `binom(3m-1,2) - 3 binom(m,2) = 3m^2 - 3m + 1`; with the rotation-fixed `(m,m,m)` this gives `m^2 - m + 1` good necklaces, `m(m-1)` of them primitive, out of `(binom(3m-1,2) - 1)/3 = 3 binom(m,2)` primitive necklaces: exactly one third of the primitive three-bead necklaces of length `3m` have no fair 3-split, for every `m`.

## 5. How many thieves: expected fair cuts and a Pólya trichotomy

For a uniformly random word of shape `(jm, jx)` and `N_j` the number of fair cut positions `r in [0,m)`: each `r` is fair for exactly `binom(m,x)^j` words, so `E[N_j] = m binom(m,x)^j/binom(jm,jx)` exactly, and by Stirling `E[N_j] ~ sqrt(j) m^{(3-j)/2} (2 pi rho(1-rho))^{-(j-1)/2}`. Two thieves: always (section 4), `E[N_2] ~ sqrt(m/(pi rho(1-rho)))`. Four or more: `P(N_j >= 1) <= E[N_j] -> 0` like `m^{-(j-3)/2}`. Three: `E[N_3] -> sqrt(3)/(2 pi rho(1-rho))` (`1.1027` at `rho = 1/2`), but the fair cuts cluster — the imbalance pair `(F(r+m) - F(r) - x, F(r+2m) - F(r+m) - x)` is a planar lattice walk of `m` steps started at distance `~ sqrt(m)`, and a fair cut is a visit to the origin. Two fair 3-cuts at `0` and `delta` force the six window counts `s, x-s, s, x-s, s, x-s`, so `E[N_3^2] = (m/binom(3m,3x)) sum_{delta} sum_s binom(delta,s)^3 binom(m-delta,x-s)^3`, computed exactly to `m = 10^4` (`~ 0.93 log m` at `rho = 1/2`), giving the proved bound `P(N_3 >= 1) >= E[N_3]^2/E[N_3^2]` (`0.249, 0.180, 0.141` at `m = 10^2, 10^3, 10^4`). Monte Carlo to `m = 10^5`: `P(N_3 >= 1) log m = 2.13 .. 2.25` (flat) and `E[N_3 | N_3 >= 1] ~ 0.5 log m` — the law `P(N_3 >= 1) ~ 2.2/log m` is EMPIRICAL, its order `1/log m` from below is PROVED only at the computed `m` (`E[N_3]^2/E[N_3^2] = 0.249, 0.180, 0.141, 0.128` at `m = 10^2, 10^3, 10^4, 3·10^4`); that `E[N_3^2] = O(log m)` is EMPIRICAL [corrected by the independent audit `collatz_necklace_20260929_audit.md`, opus S23, 2026-09-29]. (`..._fairsplit_asymptotics.out`, `..._fairsplit_loglaw.out`.)

## 6. The fair Eliahou identity

One odd step multiplies `y` by `(q/2)(1 + d/(qy))`, so a share of shape `(m,x)` multiplies by `q^x 2^{-m} P`, `P = prod_{odd y in share} (1 + d/(qy))`. For a fair 2-split with cut elements `n, n'`: `n'/n = q^x 2^{-m} P_0` and `n/n' = q^x 2^{-m} P_1`, hence `n'/n = sqrt(P_0/P_1)` and `|log(n'/n)| <= X|d|/(2(q n_min - |d|))`. Fairness cancels the drift; the two thieves' cut values differ only by the carry ratio. For a hypothetical positive `3x+1` cycle with even `K, X` (`n_min >= 2^68`), antipodal fair cut elements agree to relative precision `X/(12 n_min)`.

## 7. Checks and boundary

Exact verification (`Z[zeta_j]` arithmetic modulo `Phi_j`) on 49 small `3x+d` cycles with `gcd(K,X) > 1`, the nine `7x+169` cycles of shape `(9,3)` (six with fair 3-splits; THM-4516), and the `3x+17021` cycle of shape `(2140,1088)` (79 fair antipodal cuts, no fair 4-split, so its factor `2^1070 + 3^544` has no share reading). Both existence criteria agree with brute force on all 23,201 `(necklace, j)` pairs with `K <= 18`. All known `3x+-1` cycles have `gcd(K,X) = 1` and admit no fair split. Collatz is OPEN; nothing here restricts a hypothetical cycle beyond THM-4484.
