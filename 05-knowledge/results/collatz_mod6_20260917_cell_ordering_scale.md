# Sandwich cell ordering across scales: E. Landau, Selberg-Delange, and the non-tournament

**Status:** PROVED scoped (no 4-vertex tournament has score profile max > mid = mid > min, by the parity identity a+2m+b=6 and by enumeration; the Selberg-Delange coefficient identities f_2 = M and g_2 = M - 5/6; the mirror involution k -> -k transposes the cell matrix and orients nothing, inherited SW2) + FINITE-EXACT (4x4 Omega matrices at K = 10^3..10^7, stabilization scales and their exact failure sets, ranking-change census, zero-run counts, crossovers, independence ratios, 3x3 conditionals) + PROVED-relative-to-cited-SD (the class-marginal asymptotics) + HEURISTIC (product model, predicted asymptotic order, crossover extrapolation, Hardy-Littlewood identification of R_PP) + REFUTED (the naive Landau ordering; its L = 2 crossover K = e^{e^2}/6 ~ 270; minimal witness N33(270) = 1, N22(270) = 42) + SCOPE (no intrinsic pairwise relation between cells beyond the mirror was found; no map from THM-4139/THM-4146/THM-3341/THM-3333 to the cell census). OPEN: every cell asymptotic, the mirror ties, and the predicted N33 > N22 crossover near x ~ 10^10-10^11. No twin-prime, Chen-refinement, Collatz, or Goldbach claim is made. The non-tournament verdict is classical and already inherited (no novelty claim); the score-sequence obstruction and the SD-F6 explanation of the observed order are this lane's contribution. Research result note, not a reserved canon ID.

**Session:** collatz-mod6-20260917 (mac-mini), lane `cell_ordering_scale`; script recovered from the agent transcript on 2026-09-21 and finalized after two verifier audits (recompute, proof-audit). Everything numeric below is quoted from the regenerated `collatz_mod6_20260917_cell_ordering_scale.out`.

## Inheritance and concept board

Read and not re-derived: [arithmetic_braids_20260917_divisors.md](arithmetic_braids_20260917_divisors.md), section SW1 (definition (SW1) of `N_ij(K) = #{1 <= k <= K : Omega(6k-1) = i, Omega(6k+1) = j}`; the exact CRT identity (SW2) `sum x^{h-} y^{h+} = prod (p-2+x+y)` and its covariance (SW3) `Cov(h_-,h_+) = -sum 1/p^2`; and the sentence "no intrinsic pairwise winner relation, so they are not tournaments", which is the **origin of the non-tournament verdict** restated here), section SW2 (the mirror `k -> -k` with `|6(-k)-1| = 6k+1`, swapping `(i,j) <-> (j,i)`; the class-parity law (SW4) `b(6k-1)` odd, `b(6k+1)` even), and section SW3 (the `K = 10^6` census table, used as a positive control, and the Chen 1973 / Pintz arXiv:1004.1065 citation audit). The sibling lane [collatz_mod6_20260917_sandwich_bias.out](collatz_mod6_20260917_sandwich_bias.out) owns the antisymmetric part `N_ij - N_ji`; it is cited for its independence prediction and the prime-race count, not re-derived. [arithmetic_braids_20260917_collatz.md](arithmetic_braids_20260917_collatz.md), [arithmetic_braids_20260917_summand.md](arithmetic_braids_20260917_summand.md) and the braids2 notes were read; the only overlap is the parallel guardrail in [arithmetic_braids2_20260917_floor_reciprocity.md](arithmetic_braids2_20260917_floor_reciprocity.md) (a symmetric relation carries "no tournament orientation"), the same principle applied here to the mirror involution.

Canon read: [THM-4057](../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md) (the repo's tournament-gauge guardrail: an ordering is a gauge unless an intrinsic pairwise observable exists; applied here, not extended), [THM-4139](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md) / [THM-4146](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md) (rational three-cycle of `x^2 - 29/16`, forced `3:4:5` class), [THM-3341](../../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md) / [THM-3333](../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md) (Gaussian squaring of triples, Pell hypotenuses). SCOPE: none of these shares a predicate with the cell census; THM-4146's related-list entry for THM-4057 is the same guardrail, not a mathematical link; **no map found** from cell counts or their ordering to a preperiodic set, a Pell orbit, or a Gaussian square. Nothing is imported from them.

Concept board. *Closest proved mechanism:* the Selberg-Delange secondary term, `g_2 = M - 5/6 < 0` on the classes coprime to 6, which is what makes the observed order differ from the E. Landau leading-term order at every computed scale (section 1). *Canonical hostile:* the random re-pairing of right endpoints at `K = 10^7`, which flattens every independence ratio to `1 +- 1.33/sqrt(N_ij)` and so proves the `R_ij` pattern is a property of the actual pairing (section 3). *Corrected near miss:* the draft's "SS first and CC fourth of six at every scale" and "the merged Landau order puts CC first", both false; in the six-class ranking CC is sixth or fifth at every decade and the merged Landau order has CC second for `2 < L < 4` (sections 2 and 5). *Least-used sidecar:* the exact failure sets of the 2x2 order, `{1..761}` and `{1..160} minus {128..136}`, which are intervals up to a nine-value gap and are the cheapest thing a future Lean or proof-of-record check could assert.

Conventions (SW1): rows describe `W-1 = 6k-1`, columns `W+1 = 6k+1`; tail `Omega >= 4` collapsed; `P, S, C` = `Omega` 1, 2, 3; `L = log log 6K`. Two different Landaus appear: **E. Landau** (the `pi_k(x)` asymptotic, 1900/1909) in sections 1-3 and **H. G. Landau** (the 1953 score-sequence theorem) in section 4.

## 1. The product heuristic, its predicted order, and why the leading term fails (HEURISTIC; f_2 = M PROVED)

E. Landau: `pi_k(x) ~ (x/log x) L^{k-1}/(k-1)!`. With class marginals `N_i.(K) ~ (3K/log 6K) L^{i-1}/(i-1)!` and a Hardy-Littlewood-type product `N_ij ~ c_ij N_i. N_.j / K` with a common singular series, the nine cells carry the weights

| cell | Landau weight |
|---|---|
| (1,1) | 1 |
| (1,2) = (2,1) | L |
| (2,2) | L^2 |
| (1,3) = (3,1) | L^2/2 |
| (2,3) = (3,2) | L^3/2 |
| (3,3) | L^4/4 |

Ties `(i,j) ~ (j,i)` are exact (symmetric weight); `(2,2)` vs `(1,3)` is the factor `1` vs `1/2`. The thresholds `(3,3) = (2,3) = (2,2)` all sit at `L = 2`, i.e. `K = e^{e^2}/6 = 269.7 ~ 270`. Predicted **nine-cell** order for `L > 2`: `CC > CS = SC > SS > CP = PC > PS = SP > PP`. Predicted **merged six-class** order (mixed weights doubled) for `2 < L < 4`: `SC+CS > CC > SS = PC+CP > PS+SP > PP`, with CC second; CC is first in the merged order only for `L > 4`, `x = exp(exp(4)) ~ 10^23.7` (both evaluated and asserted in the script). The user's "max, min, two equal middles" is the asymptotic 2x2 picture `SS > PS = SP > PP` of this heuristic; the tie is heuristic and asymptotic only.

**The leading term is wrong at every computed scale** (sections 2, 5). The mechanism is the Selberg-Delange secondary term. CITED: Sathe (1953) and Selberg (1954, J. Indian Math. Soc. 18; page numbers UNCITED-RECOLLECTION); Tenenbaum, *Introduction to analytic and probabilistic number theory*, Ch. II.5 (method) and II.6 (application to `pi_k`), chapter and theorem numbers UNCITED-RECOLLECTION: `sum_{n<=x} z^{Omega(n)} = x (log x)^{z-1} (F(z) + O(1/log x))` uniformly for `|z| <= 2 - delta`, with `F(z) = H(1,z)/Gamma(z)`, `H(1,z) = prod_p (1-1/p)^z (1-z/p)^{-1}`, so that `pi_k(x) ~ (x/log x) sum_{m<=k} f_m L^{k-m}/(k-m)!` where `F(z) = sum f_m z^m`.

**PROVED (f_2 = M, g_2 = M - 5/6).** `1/Gamma(z) = z + gamma z^2 + O(z^3)`, and `log H(1,z) = sum_p [z log(1-1/p) - log(1-z/p)] = z sum_p [log(1-1/p) + 1/p] + O(z^2) = (M - gamma) z + O(z^2)` by Mertens' second theorem in the form `M = gamma + sum_p [log(1-1/p) + 1/p]`. Multiplying, `f_0 = 0`, `f_1 = 1`, `f_2 = gamma + (M - gamma) = M`. Removing the primes 2 and 3 multiplies the Euler product by `(1 - z/2)(1 - z/3) = 1 - 5z/6 + z^2/6`, so for integers coprime to 6, `F_6(z) = F(z)(1 - z/2)(1 - z/3)` and `g_2 = f_2 - 5/6 = M - 5/6`, `g_3 = f_3 - 5 f_2/6 + f_1/6`. Numerically (mpmath, prime zeta, 30 digits; `f_2` checked against `mp.mertens` to `1e-18`): `f_2 = 0.261497212847643`, `f_3 = -0.562152927240038`, `f_4 = 0.305977961616984`; `g_2 = -0.571836120485691`, `g_3 = -0.613400604613073`. The negative `g_2` is the mechanism: on the classes `+-1 mod 6`, semiprimes and 3-almost-primes are far rarer than the leading term says.

**Class marginals (PROVED-relative-to-cited-SD).** The marginal `N_i.(K)` counts `n = 6k-1 <= 6K` with `Omega(n) = i` in the class `5 mod 6`; the SD count for integers coprime to 6 must be halved. This uses that exact-`Omega` counts equidistribute between the classes `1` and `5 mod 6`: the chi-twisted series `sum chi(n) z^{Omega(n)} n^{-s} = L(s,chi)^z x (holomorphic)` has no pole at `s = 1`, so the class difference is `O(x/log^A x)` for every `A`, below the secondary terms used here. The data agree: at `K = 10^7` the row marginals are `[1781238, 3763472, 2912338, 1542952]` and the column marginals `[1780875, 3763996, 2911912, 1543217]` (differences `363`, `-524`, `426`, `-265`, all below `0.03%`); `N_1. - N_.1 = 363` is the class-5 versus class-1 prime race (sibling lane, `pi_chi(60000001) = -363`).

**Leading-term accuracy (FINITE-EXACT).** Global `pi_k(x)` at `x = 10^6, 10^7, 6*10^7` (`pi(10^6) = 78498`, `pi(10^7) = 664579`, `pi_2(10^6) = 210035` asserted against known values): `actual/Landau` for `pi_2` is 1.1051, 1.1041, 1.1024; for `pi_3` 1.0053, 1.0196, 1.0279 (an accidental cancellation, `f_2 > 0 > f_3`); for `pi_4` 0.9069, 0.9231, 0.9339. The `k`-term SD polynomial gives 1.0050, 1.0092, 1.0108 for `pi_2` and 0.9703, 0.9779, 0.9825 for `pi_3`. On the sandwich marginals:

| K | N_2./Landau | N_3./Landau | N_2./SD-F6 | N_3./SD-F6 |
|---|---|---|---|---|
| 10^3 | 0.5858 | 0.1797 | 0.7963 | 0.8591 |
| 10^4 | 0.6746 | 0.2613 | 0.8858 | 0.8437 |
| 10^5 | 0.7268 | 0.3254 | 0.9329 | 0.8679 |
| 10^6 | 0.7579 | 0.3768 | 0.9571 | 0.8944 |
| 10^7 | 0.7787 | 0.4177 | 0.9712 | 0.9154 |

The leading term overestimates the class 3-almost-primes by a factor `1/0.4177 ~ 2.4` at `K = 10^7`; SD-F6 is within `0.84-0.92` for `i = 3` (non-monotone: 0.8437 at `10^4`) and within `0.80-0.97` for `i = 2`.

## 2. Exact matrices, rankings, stabilization (FINITE-EXACT)

`Omega` was sieved on `[1, 6*10^7+1]` (uint8; 3562115 primes up to 60000001). Audits, all with explicit `raise`: trial division on `2..30000`; `sympy.factorint` on 300 random `n`; an independent trial-division matrix at `K <= 1000`; equality with the inherited SW3 table at `K = 10^6`; the chunked cumulative totals against the direct count. All pass.

`K = 10^7` (rows `W-1`, columns `W+1`):

| L\\R | P | S | C | 4+ | row |
|---|---:|---:|---:|---:|---:|
| P | 280557 | 632905 | 539140 | 328636 | 1781238 |
| S | 632766 | 1380659 | 1118019 | 632028 | 3763472 |
| C | 539061 | 1118744 | 837851 | 416682 | 2912338 |
| 4+ | 328491 | 631688 | 416902 | 165871 | 1542952 |
| col | 1780875 | 3763996 | 2911912 | 1543217 | 10000000 |

`K = 10^5`: rows P (5330, 10091, 6537, 2614), S (10050, 17794, 10719, 3851), C (6561, 10751, 5600, 1664), 4+ (2583, 3844, 1666, 345). `K = 10^4`: P (810, 1338, 687, 206), S (1323, 1942, 915, 231), C (685, 935, 367, 62), 4+ (196, 233, 61, 9). `K = 10^3`: P (142, 182, 60, 13), S (167, 192, 71, 7), C (65, 71, 7, 2), 4+ (10, 9, 2, 0). `K = 10^6` equals the inherited SW3 table (positive control).

Nine-cell rankings: `10^3` `SS > PS > SP > PP > SC = CS > CP > PC > CC`; `10^4` `SS > PS > SP > CS > SC > PP > PC > CP > CC`; `10^5` `SS > CS > SC > PS > SP > CP > PC > CC > PP`; `10^6` `SS > SC > CS > PS > SP > CC > CP > PC > PP`; `10^7` `SS > CS > SC > CC > PS > SP > PC > CP > PP`.

Mirror-merged six classes: `10^3` `PS+SP > SS > PP = SC+CS > PC+CP > CC`; `10^4` `PS+SP > SS > SC+CS > PC+CP > PP > CC`; `10^5` `SC+CS > PS+SP > SS > PC+CP > CC > PP`; `10^6` and `10^7` `SC+CS > SS > PS+SP > PC+CP > CC > PP`.

Ranks at `K = 10^3..10^7` (1 = most frequent; asserted): nine-cell `SS` = `[1, 1, 1, 1, 1]`, nine-cell `CC` = `[9, 9, 8, 6, 4]`; six-class `SS` = `[2, 2, 3, 2, 2]`, six-class `CC` = `[6, 6, 5, 5, 5]`. So `SS` is first only in the nine-cell ranking, and `CC` is fourth only at `10^7` in the nine-cell ranking; in the merged ranking `CC` is sixth or fifth at every decade.

2x2 block: `PS > PP > SP > SS` (`K = 100`), `PS > SP > PP > SS` (270), `PS > SP > SS > PP` (500), `SS > PS > SP > PP` (every decade `>= 10^3`).

**Stabilization (running check over every `K <= 10^7`; failure sets asserted).** `SS > max(PS,SP)` fails for exactly 761 values of `K`, namely every `K` in `{1,...,761}`; `min(PS,SP) > PP` fails for exactly 151 values, namely `{1,...,160} minus {128,...,136}`. Hence `SS > {PS,SP} > PP` for every `K` in `[762, 10^7]`.

**Mirror differences (census only; mechanism in the sibling lane).** `D_PS-SP(K)`: 992 strict sign flips; `D = 0` at 14947 values of `K`, forming 1928 maximal zero runs (returns to zero, counting the initial run `K = 1..3` before the first mixed cell at `k = 4`); `max|D| = 650` at `K = 5512997`; final `+139` against `sqrt(N_PS + N_SP) = 1125`. `D_PC-CP`: 639 flips, 12330 zeros in 1368 runs (initial run `K = 1..20`), `max|D| = 570` at `K = 7331193`, final `+79` vs 1038. `D_SC-CS`: 470 flips, 4654 zeros in 1004 runs (initial run `K = 1..40`), `max|D| = 1323` at `K = 9539838`, final `-725` vs 1496. The final `+139` is not "pure fluctuation": the sibling lane's exact independence prediction at `K = 10^7` is `+229.95` (semiprime-class part 93.33, prime-class part 136.62) with residual `-90.95`, so the data are consistent with fluctuation around a small marginal-driven term. Consequently the nine-cell weak order changes 8918 times (last at `K = 8714402`, one past the sibling lane's last `K'` with `N_PS <= N_SP`, 8714401) and the 2x2 order 3885 times (same last `K`); the six mirror-merged classes change 259 times, last at `K = 920808` (the 2x2 mirror-merged order 33 times, same last `K`), and are stable through `10^7` afterwards.

**Crossovers for good (holding for all `K` in `(K_last, 10^7]`):** `SS > PS` after 761; `SS > SP` after 672; `SS > PS+SP` after 920807; `SC+CS > SS` after 16815; `SC+CS > PS+SP` after 66959; `CC > PS` after 1591506; `CC > SP` after 1559904; `PC+CP > PP` after 1445; `CC > PP` after 83683. `CC` never exceeds `SS`, `PS+SP`, or `PC+CP`, and neither `SC` nor `CS` ever exceeds `SS`, for any `K <= 10^7`.

**Model comparison (six-class order, symmetric product models with `c_ij := 1`).** The Landau-weight merged order `SC+CS > CC > SS = PC+CP > PS+SP > PP` (CC second, `SS` tied with `PC+CP`) differs from the actual order at all five decades: REFUTED. The SD-F6 product order matches at `10^4, 10^6, 10^7` (3 of 5); the measured-marginal product matches at `10^4, 10^5, 10^7` (3 of 5). The mismatches are the near-tie `SS` vs `PS+SP` and the tiny-count regime at `10^3`. The two `c = 1` models bracket the actual `SS > PS+SP` crossing (`K = 920807`, `x ~ 10^6.74`): SD-F6 marginals place it at `x ~ 10^5.68` (early) while the measured-marginal product still has `PS+SP > SS` at `10^6` (late). The observed ordering is therefore explained by the marginals; local correlation only nudges near-ties.

## 3. (3,3) versus (2,2), independence ratios, crossover (FINITE-EXACT counts, HEURISTIC extrapolation)

| K | L | N33 | N22 | N33/N22 | r- = N_3./N_2. | r+ = N_.3/N_.2 | c33/c22 | Landau L^2/4 |
|---|---|---:|---:|---|---|---|---|---|
| 10^3 | 2.1633 | 7 | 192 | 0.03646 | 0.33181 | 0.30837 | 0.35632 | 1.1699 |
| 10^4 | 2.3981 | 367 | 1942 | 0.18898 | 0.46452 | 0.45638 | 0.89142 | 1.4377 |
| 10^5 | 2.5881 | 5600 | 17794 | 0.31471 | 0.57943 | 0.57726 | 0.94089 | 1.6746 |
| 10^6 | 2.7477 | 72363 | 157420 | 0.45968 | 0.68303 | 0.68213 | 0.98661 | 1.8875 |
| 10^7 | 2.8854 | 837851 | 1380659 | 0.60685 | 0.77384 | 0.77362 | 1.01367 | 2.0813 |

(The intermediate scales `2, 5 x 10^j` are in the output.) The naive leading-term ratio `L^2/4` puts `N33 > N22` at `L = 2`, `K = e^{e^2}/6 ~ 270`: **REFUTED**, minimal witness `N33(270) = 1`, `N22(270) = 42` (also `N33(272) = 1`, `N22(272) = 43`; both asserted).

Crossover models for `N33 = N22` (HEURISTIC): B (SD-F6 marginal ratio `(L^2/2 + g_2 L + g_3)/(L + g_2)`, `c = c(10^7) = 1.0137`): `L* = 3.1589`, `K* = 2.800e9` (`log10 K* = 9.45`); C (linear fit `sqrt(r- r+) = 0.657 L - 1.123` on `K >= 10^5`, `c` fixed): `L* = 3.2197`, `K* = 1.224e10` (10.09); D (direct fit `sqrt(N33/N22) = 0.734 L - 1.338`): `L* = 3.1867`, `K* = 5.434e9` (9.74); E (`c` drifting linearly in `1/L` to `c_inf = 1.5541`): `L* = 3.1268`, `K* = 1.331e9` (9.12). Spread `log10 K* in [9.12, 10.09]`, i.e. `x = 6K` with `log10 x in [9.90, 10.87]`. The SD ratio overshoots the measured `r+-` at `10^7` (0.8209 vs 0.7738 / 0.7736), so B is probably early. Under SD-F6 with `c = 1`: `CC > SS` at `L = 3.1699`, `x ~ 10^10.3`; `CC > SC+CS` (the naive CC-first merged order) only at `L = 4.9285`, `x ~ 10^60.0`; for comparison the pure-Landau merged threshold is `L = 4`, `x ~ 10^23.7`. All extrapolations go 2-3 orders of magnitude beyond the data and assume `c33/c22` converges; none is a theorem.

**Independence ratios** `R_ij = N_ij K/(N_i. N_.j)` at `K = 10^7`:

| L\\R | P | S | C | 4+ |
|---|---|---|---|---|
| P | 0.8844 | 0.9440 | 1.0394 | 1.1955 |
| S | 0.9441 | 0.9746 | 1.0202 | 1.0882 |
| C | 1.0394 | 1.0206 | 0.9880 | 0.9271 |
| 4+ | 1.1955 | 1.0877 | 0.9279 | 0.6966 |

`R_PP = 0.9315, 0.8837, 0.8845, 0.8898, 0.8844` across the decades, against the Hardy-Littlewood value `4C_2/3 = 0.88022` (`C_2 = 0.6601619` from primes `<= 10^6`). Derivation (HEURISTIC identification): `2C_2 x/log^2 x` twin pairs up to `x`, all with center in `6Z` except `(3,5)`, divided by the product of marginals `(3K/log 6K)^2/K` with `x = 6K` gives `4C_2/3 = prod_{p>=5} (1 - 1/(p-1)^2)`; both sides use `li`-type main terms whose `1/log` corrections cancel to first order, which is why agreement to `0.5%` is plausible. The sign agrees with the inherited exact CRT covariance (SW1, eq. (SW3)) `-sum 1/p^2`, but that identity concerns hit counts of a finite prime set, not `Omega` cells, and does not derive `R_ij`. The matrix is nearly symmetric (mirror heuristic). **Hostile control:** randomly re-pairing the right endpoints (seed 20260917) gives every shuffled `R_ij` within `1.33/sqrt(N_ij)` of 1 (shuffled `R_PP = 0.9997`), so the pattern is a property of the actual pairing.

## 4. Tournament audit (repo guardrail)

**PROVED.** The score sequences of the 64 labelled tournaments on 4 vertices are exactly `(0,1,2,3), (0,2,2,2), (1,1,1,3), (1,1,2,2)` (enumeration; H. G. Landau 1953). One-line proof that no profile `max > mid = mid > min` occurs: the scores sum to `C(4,2) = 6`, so `(b, m, m, a)` with `0 <= a < m < b <= 3` would need `a + 2m + b = 6`, which has no integer solution in that range (the script lists the solution set: empty). The user's "6 arcs compressed because the 4 can be ordered into 3: a max, a min, and two equal middle values" is therefore not the score profile of any tournament; it is a weak order with one tie, a transitive comparability digraph with 5 arcs, not 6. Breaking the tie by any gauge yields the transitive tournament `(0,1,2,3)`, which carries no information beyond the ranking.

**Mirror involution (PROVED, inherited SW2) and SCOPE.** `iota: k -> -k` satisfies `|6(-k)-1| = 6k+1`, so it transposes the cell matrix: three fixed points `PP, SS, CC`, three 2-cycles. It is symmetric, so it orients nothing; and it leaves the positive prefix `1 <= k <= K`, so on the census it relates the positive and the signed census rather than members of two cells. Any orientation `(i,j) -> (j,i)` for `i < j` is the gauge "left endpoint < right endpoint". SCOPE: no other intrinsic pairwise relation between cells was found; this is a non-discovery, not a theorem of absence. The verdict "outcome table, not a tournament" is inherited from section SW1 of the divisors note; this lane adds the quantitative reason (score sequences) and the finite-scale evidence that the claimed tie is violated at almost every `K` (section 2).

Guardrail fields (THM-4057 style). Vertices: 4 (or 9) cells. Pairwise observable: none intrinsic; the only candidate, "`N_a(K) > N_b(K)` at this `K`", is scale-dependent and is a preorder on cells, not a relation between members. Orientation gauge: `6k-1 < 6k+1` or count order. Ties: the mirror pairs `(i,j) ~ (j,i)`, heuristic and violated at almost every finite `K`. Preserved target: nothing beyond the counts. Lost data: residue-factor parity (`b` odd/even, SW4), cofactors, everything arithmetic. Sidecar: the counts themselves. **Verdict:** the 4-tournament reading is a cosmetic ranking (a judgement supported by the above, not itself a theorem).

Typed analogies. (a) Cell frequencies -> tournament: map "orient `a -> b` iff `N_a > N_b`"; preserved: the ranking; lost: the tie; sidecar: none; test: the score sequences above. Not an analogy, a relabeling. (b) E. Landau `pi_k` -> cell weight: map = product with a common singular series; preserved: leading power of `L` and the tie structure; lost: the SD secondary terms (which set the order up to `~10^60`) and the local correlation `c_ij`; sidecar: `g_m` and the measured `R_ij`; test: the class-marginal table (leading term off by 2.4x). (c) Mirror `k -> -k`: inherited SW2; preserved: transposition; lost: the positive prefix; sidecar: the sign. (d) THM-3341/THM-3333 Gaussian squaring and Pell hypotenuses -> cell census: **no map found** (SCOPE).

## 5. The 3x3 (P,S,C) block (FINITE-EXACT)

Conditional on both endpoints having `Omega <= 3` (block totals 957, 9002, 83433, 769083, 7079702 of `K = 10^3..10^7` centers; fractions 0.9570, 0.9002, 0.8343, 0.7691, 0.7080):

| K = 10^7 | P | S | C | row |
|---|---|---|---|---|
| P | 0.03963 | 0.08940 | 0.07615 | 0.20518 |
| S | 0.08938 | 0.19502 | 0.15792 | 0.44231 |
| C | 0.07614 | 0.15802 | 0.11835 | 0.35251 |
| col | 0.20515 | 0.44244 | 0.35242 | |

At `K = 10^6`: P (0.04930, 0.10232, 0.07763), S (0.10178, 0.20469, 0.14617), C (0.07800, 0.14602, 0.09409). Exact fractions (e.g. `SS = 1380659/7079702`, `SC+CS = 2236763/7079702`, `CC = 837851/7079702` at `10^7`) are in the output for every decade and sum to 1 exactly (asserted). Adding `C` turns the 2x2 "max, min, two middles" into six mirror classes: three symmetric pairs `PS~SP`, `PC~CP`, `SC~CS` (asymptotically tied by the mirror heuristic) and three singletons `PP, SS, CC` (ordered). In the nine-cell ranking `SS` is first at every decade and `CC` rises from ninth (`10^3, 10^4`) to eighth, sixth and fourth (`10^7`); in the merged six-class ranking `SS` is second (third at `10^5`) and `CC` sixth (`10^3, 10^4`) then fifth (`10^5..10^7`). Only at `10^7` does the nine-cell order read `SS > CS > SC > CC > ...`. Under SD-F6 the CC-first merged order is a `~10^60` phenomenon.

## 6. Reproduction

```bash
python3 04-computation/experiments/collatz_mod6_20260917_cell_ordering_scale.py > 05-knowledge/results/collatz_mod6_20260917_cell_ordering_scale.out
```

About 8-9 s, peak RSS 0.86 GB (measured under `python3 -O`); `python3 -O` produces identical output modulo the `[t=...]` stamps (checked by diff). Final source sha256 `0cdebab552a502b92f8905b536958b78799125a41daa0928294cc5b92b1545ea`; output sha256 `c77be555253ac8a797666bf5157123c139264408b5fe1dad285936198ad67ceb` (395 lines). Universe: all centers `6k`, `1 <= k <= 10^7`, no filter. Controls: section 2 audits, the SW3 positive control, the shuffle control of section 3, and the audit-established assertions added at finalization (failure sets, zero-run counts, `N33/N22` at `K = 270` and `272`, the `a+2m+b=6` parity proof, the per-decade ranks, the merged Landau order at `L = 3`, `N_1. - N_.1 = 363`); all with explicit `raise`.

Audit provenance. The two verifier scripts `04-computation/experiments/collatz_mod6_20260917_cell_ordering_scale_audit_recompute.py` (smallest-prime-factor sieve, explicit SD series with A&S 1/Gamma coefficients, own tournament enumeration plus H. G. Landau's criterion) and `..._audit_proof-audit.py` (spf sieve, full-array cumulative census, hand-expanded SD coefficients, parity proof) reproduced every matrix, ratio, stabilization scale, crossover, census count, mirror statistic, model `L*`, block fraction and score sequence of the explorer's run. Their copies in this worktree are transcript-recovered drafts whose hashes differ from those the audits reported, and their outputs were not recovered; they are cited as provenance, not as regenerated evidence. The quantities they established beyond the explorer's output (failure sets, zero-run counts, `K = 270` witness, ranks, merged Landau order) were recomputed independently at finalization and are now asserted by the lane script itself.

## 7. Stopping boundary / next question

The ordering question is closed at the level available: exact through `10^7`, explained by the Selberg-Delange marginals on the classes coprime to 6, and not a tournament. Every cell asymptotic remains OPEN: `(1,1)` is the twin prime conjecture; Chen's theorem (inherited SW3 audit) gives only prime plus `P_2`, i.e. the union `PP u PS` (resp. `PP u SP` for `h = -2`), includes the twin cell, does not split the exact-`Omega` outcomes (parity obstruction, Pintz), and in its bare form does not restrict centers to `6Z`. The mirror ties `N_ij ~ N_ji` are conjectural for every pair. The first falsifiable prediction of the product model is `N33 > N22` near `x = 6K ~ 10^10-10^11` (`log10 K* in [9.12, 10.09]`); a segmented `Omega` sieve to `6*10^10` would test it, and a wrong decade would falsify the assumption that `c33/c22` converges. The next substantive step is to derive the heuristic constants `c_ij` for exact-`Omega` pairs `(n, n+2)` from a Hardy-Littlewood-type conjecture for almost-prime tuples and compare with `R_SS = 0.9746`, `R_SC = 1.0202`, `R_CC = 0.9880`, as `R_PP -> 4C_2/3` already matches to `0.5%`. Not attempted here: any Lean statement; the natural candidate is the finite fact "`SS > {PS,SP} > PP` for all `K` in `[762, 10^7]`" with its explicit failure sets, which is a decidable-by-computation claim.
