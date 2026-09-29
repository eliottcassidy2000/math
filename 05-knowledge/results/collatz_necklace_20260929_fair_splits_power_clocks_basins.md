# Fair splits of cycle necklaces are circulant, perfect-power clocks are Fermat–Catalan identities, and root-uniform basin bounds cannot exist: the corrected Krasikov–Lagarias argument

**Status: PROVED (elementary) for the three theorems THM-4515 / THM-4516 / THM-4517; FINITE-EXACT for every census and every density table; CITED for Darmon–Granville, the ten Fermat–Catalan solutions, Krasikov–Lagarias `X^0.84`, Barina's `2^68`, and the April 2025 survey of solved Fermat–Catalan signatures; UNVERIFIED for the owner's statement about `x^4 + y^3 = z^17`. Collatz OPEN. Single session, not independently audited.**

Session `collatz-necklace-20260929` (machine `mac-mini`, 2026-09-29). Scripts `04-computation/experiments/collatz_necklace_20260929_{basins_minus.c, entries_plus.c, fairsplit_census.py, power_clocks.py, circulant_check.py}`; outputs `05-knowledge/results/collatz_necklace_20260929_*.out` (hashes in the theorem files).

**Owner's seed.** A cycle of a Collatz-like graph with `K` total steps and `X` odd steps; in the dual of the planar unicyclic component the cycle face is joined to the unbounded face by `K` parallel edges; "differentiate and familiarize the possible orderings of lengths of the two colourings (even step / odd step)", "the chains extending from integers", the stolen-necklace problem, the Krasikov–Lagarias predecessor bound applied to any root and the `3x-1` sheet with its three cycles, and the statement "`x^4 + y^3 = z^17` has no solution with `xyz != 0`, `gcd(x,y) = 1`".

## Inheritance and concept board

* Closest proved mechanism: [THM-4484](../../01-canon/theorems/THM-4484-free-and-sporadic-cycles.md) — the periodic point of a parity word `w` of shape `(K,X)` under `T_{q,d}` is `d c_w/(2^K - q^X)`; a mixed shape is *free* iff the clock `2^K - q^X` divides `d`; `|2^K - q^X| = 1` only at Gersonides' four cases, hence the free cycles `{0}, {1,2}, {-1}, {-5,-7,-10}` and the sporadic `-17`.
* Canonical hostile example: the `3x-1` sheet (the negative `3x+1` sheet) with its three positive cycles through `1`, `5`, `17` ([counterexample portrait](collatz_mod6_20260922_counterexample_portrait.md)): every sheet-blind necessary condition is satisfied by them.
* Corrected near miss: the [topological reframe note](collatz_generic_price_topology_20260927.md) §2 typed necklace splitting as sheet-blind and found that Alon's `k(q-1)`-cut split of the doubled `-17` word gives shares that are unions of intervals and "carry no dynamics". Refined here: with two bead types the two thieves need only two cuts, each share is one orbit arc, and fairness is exactly what makes the cycle system circulant (§1).
* Least-used sidecar: the rotation action on carries `c_{rot w} = (q^{w_0} c_w + w_0 D)/2` and the run-length (`m`-cycle) form of the carry, `c_w = sum_r 2^{p_r} q^{X_{>r}} (q^{k_r} - 2^{k_r})` over the odd runs.
* Board: (C1) cycle necklace and its clock; (C2) fair consecutive split = circulant = DFT; (C3) trunk-entry spectrum and basins (the loops of the dual); (C4) run-length orderings and their gap discrepancy; (C5) perfect-power clocks = Fermat–Catalan identities.

## 0. Objects and the dual-graph reading

`T_{q,d}(y) = y/2` (`y` even), `(qy+d)/2` (`y` odd), `q >= 3` odd, `d` odd, `gcd(q,d) = 1`; `q = 3`, `d = +-1` are the two sheets of `3x+1`. A cycle `y_0 -> y_1 -> ... -> y_{K-1} -> y_0` has parity word `w = (y_i mod 2)`, a primitive binary necklace of length `K` with `X = sum w` ones; the carry is `c_w = sum_{j : w_j = 1} q^{(ones after j)} 2^j`, the clock is `D = 2^K - q^X`, and `y_0 = d c_w/D` (THM-4471/THM-4484). Rotating the word rotates the cycle.

Planar unicyclic component (cycle plus in-trees), dual graph: two vertices (the cycle face `f_c`, the unbounded face `f_oo`), `K` parallel edges `f_c f_oo` (one per cycle edge, i.e. one per bead of the necklace, coloured odd/even), and one loop at `f_oo` per tree edge. So the dual retains the necklace as a `K`-bond and remembers the trees only as a number of loops; "the chains extending from integers" (the doubling rays and their entry points) are those loops, and they are the subject of §3, invisible to the necklace lane of §1.

## 1. Fair consecutive splits are circulant (THM-4515)

Let `j | gcd(K, X)`, `m = K/j`, `x = X/j`. A **`j`-fold fair consecutive split** of the necklace at rotation `r` is a partition of the cycle into the `j` consecutive windows `W_i = [r + im, r + (i+1)m)` (indices mod `K`), each containing exactly `x` odd letters: `j` thieves, each taking one arc, each receiving the same number of beads of each colour. Write `n_i = y_{r+im}` (the cut elements) and `c_i = c_{w|W_i}` (the share carries, each of a word of shape `(m,x)`).

**1.1 Two thieves always split fairly (discrete IVT).** For `K, X` even and `m = K/2`, let `F` be the prefix count and `h(r) = F(r+m) - F(r) - X/2` on `Z/K`. Then `h(r+m) = -h(r)` and `|h(r+1) - h(r)| <= 1`, so `h` has a zero: a fair 2-split exists for every necklace with even `K` and `X`. (Verified exhaustively for all primitive necklaces with `K <= 18`.)

**1.2 Circulant system and its DFT.** Each share is one affine map: `2^m n_{i+1} = q^x n_i + d c_i` (indices mod `j`). Because all shares have the same shape the matrix is circulant, `(2^m P - q^x I) n = d c` with `P` the cyclic shift, and the discrete Fourier transform diagonalises it. With `zeta = e^{2 pi i/j}`, `N_k = sum_i zeta^{-ik} n_i`, `C_k = sum_i zeta^{-ik} c_i`:

`N_k (2^m zeta^k - q^x) = d C_k` in `Z[zeta]`, for `k = 0, ..., j-1`.

`k = 0`: `(2^m - q^x)(n_0 + ... + n_{j-1}) = d (c_0 + ... + c_{j-1})`. `j` even, `k = j/2`: `(2^m + q^x) sum (-1)^i n_i = -d sum (-1)^i c_i`. The `j` "twisted clocks" `2^m - zeta^{-k} q^x` are the eigenvalues; fairness is exactly the condition under which the cycle equation has this Fourier structure (an unfair split gives a non-circulant product of affine maps).

**1.3 Clock factorisation and the CRT form.** `2^K - q^X = prod_{k=0}^{j-1} (2^m - zeta^k q^x) = prod_{e | j} Phi_e(2^m, q^x)` (homogenised cyclotomic polynomials). For `j = 2` the two factors `2^m -+ q^x` are coprime (both odd, any common divisor divides `2 q^x`), and since `c_w = q^x c_0 + 2^m c_1`,

`(2^K - q^X) | d c_w  <=>  (2^m - q^x) | d (c_0 + c_1)  and  (2^m + q^x) | d (c_0 - c_1)`

(mod `2^m - q^x` one has `q^x = 2^m`, so `c_w = 2^m (c_0 + c_1)`; mod `2^m + q^x` one has `q^x = -2^m`, so `c_w = 2^m (c_1 - c_0)`; `2^m` is a unit modulo both). The integrality congruence of a cycle whose shape has even `K, X` splits along any fair antipodal cut into two congruences modulo the two algebraic factors of the clock: the "sum" factor `2^m - q^x` is the clock of the average share `(m,x)`, the "difference" factor `2^m + q^x` is the clock of the map with multiplier `-q` when `x` is odd (`2^m - (-q)^x`). This is the CRT decomposition of THM-4484's criterion and is what the split *carries*.

**1.4 When do `j` thieves split fairly?** (a) General criterion: a `j`-fold fair split at `r` exists iff the integer discrepancy `G(t) = K F(t) - tX` takes one common value at `t = r, r+m, ..., r+(j-1)m` (the window counts are the successive differences of `F`, and they all equal `x` iff `G` is constant on the coset). (b) One bead per thief, `X = j`: with cyclic gaps `g_1, ..., g_j` between consecutive odd letters, a fair split exists iff every run of `t` consecutive gaps (`1 <= t <= j-1`) has sum strictly between `(t-1)m` and `(t+1)m` (intersect the `j` admissible intervals `(P_i - im, P_i - (i-1)m]` for the cut `r`; they meet iff each lower end is below each upper end, which is exactly the run condition). (c) Three beads, `K = 3m`: a fair 3-split exists iff every gap is `<= 2m - 1`; compositions of `3m` into three parts `<= 2m-1` number `3m^2 - 3m + 1`, the rotation-fixed one is `(m,m,m)`, so `m^2 - m + 1` necklaces are good, of which `m(m-1)` are primitive out of the `3 binom(m,2)` primitive necklaces: **exactly two thirds of the primitive three-bead necklaces of length `3m` admit a fair 3-split, for every `m`.** Both criteria were cross-checked against brute force on all 23,201 `(necklace, j)` pairs with `K <= 18`, zero mismatches, and the count `m(m-1)` on `m <= 7`.

Census (primitive necklaces, `K <= 18`; file `..._fairsplit_criteria_K18.out`): `j = 2` never fails; `j = 3`: `(12,6)` 47 of 75, `(15,6)` 204 of 333, `(18,6)` 621 of 1026, `(18,9)` 1592 of 2700; `j = 4`: `(16,4)` 42 of 112, `(16,8)` 255 of 800; `j = 6`: `(18,6)` 107 of 1026; `j = 8`: `(16,8)` 30 of 800; `j = 9`: `(18,9)` 56 of 2700.

**1.5 The fair Eliahou identity.** Over one odd step `y -> (qy+d)/2 = (q y/2)(1 + d/(qy))`, so over a share of shape `(m,x)` the value is multiplied by `q^x 2^{-m} prod_{odd y in share} (1 + d/(qy))`. For a fair 2-split with cut elements `n, n'` and share products `P_0, P_1`: `n'/n = q^x 2^{-m} P_0` and `n/n' = q^x 2^{-m} P_1`, hence

`n'/n = sqrt(P_0/P_1)`, and `|log(n'/n)| <= X |d| / (2 (q n_min - |d|))`.

Fairness cancels the drift `q^x/2^m` exactly; the two thieves' cut elements differ only by the carry ratio. For a hypothetical positive `3x+1` cycle (`n_min >= 2^68`) with even `K, X`, antipodal fair cut elements agree to relative precision `X/(12 n_min)`.

**1.6 FINITE-EXACT checks (`..._circulant_check_d100.out`).** All identities of 1.2–1.3 verified exactly (in `Z[zeta_j]` via reduction mod `Phi_j`) on: 49 primitive cycles of `3x+d`, `|d| <= 100`, least element `<= 100|d|`, with `gcd(K,X) > 1` (first ones: `3x+7` `(4,2)`, `3x+11` `(6,2)` and `(14,8)`, `3x+13` `(24,15)` with `j = 3` and no fair 3-split, `3x-17` two cycles `(6,4)`); the nine primitive `7x+169` cycles of shape `(9,3)` (six admit fair 3-splits, the two-thirds law of 1.4c); and the Belaga–Mignotte cycle of `3x+17021` with least element `5`, shape `(2140, 1088)`, `gcd = 4`: 79 of its 1070 antipodal cuts are fair (the first at `r = 0`: cut elements `5` and `187226`, alternating sum `-187221 = -d (c_0 - c_1)/(2^1070 + 3^544)` exactly) and it has **no** fair 4-split, so the factor `Phi_4 = 2^1070 + 3^544` of its clock has no share reading. Every known `3x+-1` cycle has `gcd(K,X) = 1` (their densities are convergents or intermediates of `log_3 2`, THM-4484), so none of them admits any fair split; a hypothetical large cycle with `gcd(K,X) = j > 1` would.

**1.7 What this refines.** The [topological reframe note](collatz_generic_price_topology_20260927.md) is right that Borsuk–Ulam-type arguments are sheet-blind (its Proposition 2) and that Alon's split of the doubled `-17` word carries no dynamics. The two-colour, two-thief instance is different: it needs only the discrete IVT, its shares are orbit arcs, and fairness is the circulant symmetry that factors the clock. No new necessary condition on Collatz cycles follows (1.3 is CRT-equivalent to THM-4484's single congruence), but the mechanism — *fairness = Fourier diagonalisation of the cycle equation* — is exact and new to the repo.

## 2. Perfect-power clocks are Fermat–Catalan identities (THM-4516)

**2.1 Dictionary.** Let `q = p^s` be an odd prime power and `X >= 2`, `K >= 1`, `r >= 2`. An identity `2^K - q^X = +- m^r` is a primitive solution of `x^a + y^b = z^c` with exponent set `{K, sX, r}` (`m` is automatically coprime to `2p`); if `1/K + 1/(sX) + 1/r < 1` it is one of the finitely many primitive solutions of that signature (Darmon–Granville), conjecturally one of the ten known. By THM-4484(1) the mixed shape `(K, X)` is then *free* for `T_{q, +-m^r}`: every necklace of that shape is an integer cycle of `qy +- m^r`. Conversely each known Fermat–Catalan solution with a pure power of two as one term is such a clock.

**2.2 Census (`..._power_clocks_q201_K400.out`, `..._q3_K3000.out`).** For odd `q <= 201`, `1 <= X <= K <= 400`, `X >= 2`, the clocks with `|2^K - q^X| in {1} cup {m^r : r >= 2}` are exactly:

| identity | `q` | shape `(K,X)` | clock | free cycles of |
|---|---|---|---|---|
| Catalan `3^2 - 2^3 = 1` | 3 | `(3,2)` | `-1` | `3x+1`: `{-5,-7,-10}` (THM-4484) |
| Pythagorean family `2^K + (2^{K-2}-1)^2 = (2^{K-2}+1)^2` | `2^{K-2}+1` | `(K,2)` | `-(2^{K-2}-1)^2` | `(2^{K-2}+1) x - (2^{K-2}-1)^2`; `q = 5, 9, 17, 33, 65, 129` in range |
| `2^5 + 7^2 = 3^4` | 3 | `(5,4)` | `-7^2` | `3x-49`: `{65, 73, 85, 103, 130}` |
| same, read as `q = 9` | 9 | `(5,2)` | `-7^2` | `9x-49`: least elements `11, 13` |
| `7^3 + 13^2 = 2^9` | 7 | `(9,3)` | `13^2` | `7x+169`: nine primitive 9-cycles (least `67,71,79,85,93,95,109,121,137`) |
| same, read as `q = 13` | 13 | `(9,2)` | `7^3` | `13x+343`: least `15, 17, 29` primitive, `21 = 7 x (3x+49 ...)` scaled |
| `2^7 + 17^3 = 71^2` | 71 | `(7,2)` | `-17^3` | `71x-4913`: least `73, 75, 79` |

Nothing else appears; for `q = 3` alone, `K <= 3000` adds nothing with `X >= 2`. The family row is the trivial identity `(a+1)^2 - (a-1)^2 = 4a` at `a = 2^{K-2}`; `2^5 + 7^2 = 3^4` is its `K = 5` member, where the odd base `9` happens to be a square, which is why it is hyperbolic (`1/5 + 1/4 + 1/2 < 1`) and appears in the Fermat–Catalan list. Among the ten known Fermat–Catalan solutions exactly four have a pure power of two as a term (`1 + 2^3 = 3^2`, `2^5 + 7^2 = 3^4`, `7^3 + 13^2 = 2^9`, `2^7 + 17^3 = 71^2`), and the census finds exactly these (twice each where the odd base is a prime power in two ways). The `X = 1` rows (`2^K = q + m^r`, e.g. `2^7 = 3 + 5^3`, giving the cycle `{1, 64, 32, 16, 8, 4, 2}` of `3x+125`) are one-odd-step shapes and are listed in the output but not counted.

**2.3 Free cycles, verified (`..._circulant_check_d100.out` §4).** Each necklace of the listed shapes is an integer cycle; the non-primitive necklaces give the expected scalings (`001001001` of `7x+169` is `169 x {4,2,1}` of `7x+1`; `0101` of `5x-9` is `9 x {1,2}` of `5x-1`; `000001001` of `13x+343` has `gcd 7`; two of the three `17x-225` necklaces are scalings by `3` and `25`).

**2.4 The Eisenstein square behind `7^3 + 13^2 = 2^9`.** The shape `(9,3)` has `gcd = 3`, so the clock factors cyclotomically: `2^9 - 7^3 = (2^3 - 7) Phi_3(2^3, 7) = 1 x (64 + 56 + 49) = 13^2`, and the twisted clocks of the fair 3-split are `2^3 - 7 zeta_3^{+-1} = (3 - zeta_3^{+-1})^2` — squares of the Eisenstein primes above `13` (verified: `8 - 7z == (3 - z)^2 mod Phi_3(z)`). The third factor `2^3 - 7 = 1` is Gersonides' unit for `q = 7`, the trivial cycle `{1,4,2}` of `7x+1`. Six of the nine primitive `7x+169` cycles admit a fair 3-split, so for them the DFT identities of §1 read `N_k (8 zeta^k - 7) = 169 C_k`, `k = 0, 1, 2`, with `8 zeta^k - 7` an Eisenstein square for `k = 1, 2`.

**2.5 Beal's shadow and the owner's `(4,3,17)` statement.** A perfect-power clock `2^K - p^{sX} = +- m^r` with `K, sX, r >= 3` would be a Beal counterexample; so Beal's conjecture implies that no map `py +- m^r` with `r >= 3` has a free mixed shape `(K,X)` with `K >= 3`, `sX >= 3`, `X >= 2`. The `13x+343` and `71x-4913` examples have a square term (`13^2`, `71^2`) and are Fermat–Catalan, not Beal. The statement "`x^4 + y^3 = z^17` has no solution with `xyz != 0`, `gcd(x,y) = 1`" is **not** among the solved signatures of the April 2025 survey (arXiv:2412.11933v2, Table 1.1: among `(3,4,n)` only `n = 4, 5` are solved; the smallest open Beal signature is `(3,5,7)`, the smallest open Fermat–Catalan one `(2,5,7)`); a web search on 2026-09-29 found no later paper. It is recorded as **UNVERIFIED** (possibly a 2026 result unknown to this session). Its clock shadow is empty in any case: no term of `x^4 + y^3 = z^17` can be a pure power of two while the other two are powers of one odd base.

## 3. Root-uniform basin bounds cannot exist; the sheet cap; the loops of the dual (THM-4517)

**3.1 The pasted argument.** "Krasikov–Lagarias show at least `X^0.84` numbers up to `X` return to 1, and the same for any root, e.g. 1729. If for any root `r` the set of `x <= X` whose orbit hits `r` had density `c > 1/2`, Collatz would follow, because two roots cannot each eat more than half. But the arguments apply to `3x-1`, which has three cycles, so such techniques cannot prove more than `1/3`." Two errors and one correct kernel:

* *No root-uniform positive proportion can exist at all* (unconditionally, §3.2): the basins of the trunk-entry points `(4^i - 1)/3` are pairwise disjoint, so a uniform lower density `c > 0` over all roots `a != 0 (mod 3)` is impossible. Krasikov–Lagarias's sublinear `X^0.84` is the only kind of root-uniform statement that can be true. The premise is void before any sheet is mentioned.
* *Even `c > 1/2` for the basin of the trivial cycle alone would not prove Collatz*: it shows only that every other cycle basin and every divergent class has upper density `< 1/2`, not that they are empty. Only a root-uniform bound over cycle roots yields "at most one positive-density cycle", and nothing about density-zero cycles or divergent orbits.
* *Correct kernel*: a cycle-uniform lower bound proved by a sheet-blind method must hold for the three disjoint positive cycle basins of `3x-1`, so it is at most `1/3`, and in fact at most the smallest of the three densities, `0.3248` (the `{5,7,10}` basin) if the measured values converge. This is the repo's SHEET control with a number attached.

**3.2 Proposition (PROVED).** For `3x+1` in `T`-form let `B(a) = {n >= 1 : T^k(n) = a for some k >= 0}` and `a_i = (4^i - 1)/3` (`i >= 2`). (i) The sets `B(a_i)` are pairwise disjoint: after `a_i` the orbit is `2^{2i-1}, ..., 2, 1, 2, ...`, which contains no `a_{i'} >= 5`. (ii) `a_i != 0 (mod 3)` iff `3` does not divide `i` (`ord_9(4) = 3`); when `3 | i`, `B(a_i)` is the doubling ray of `a_i` (a multiple of three has no odd preimage) and has density zero. (iii) Hence there is no `c > 0` with `liminf |B(a) cap [1,X]|/X >= c` for all roots `a != 0 (mod 3)`: `M > 1/c` disjoint basins of lower density `>= c` would give a union of lower density `> 1` (lower density is superadditive on disjoint sets). (iv) `B(5) and B(32)` partition the positive integers minus the powers of two; more generally `B(a) = {a} cup B(2a) cup B((2a-1)/3)` (the last term when `a = 2 (mod 3)`), so along the trunk `B(4^i) = union_{i' >= i} B(a_{i'})` and the trunk-entry densities `e_i := dens B(a_i)` (when they exist) sum to `1` if and only if almost every orbit reaches `1`.

**3.3 FINITE-EXACT densities.**

`3x-1` on `[1, 2^29]` (`..._basins_minus_2e29.out`): basin of `{1}`: `0.3268569`; of `{5,7,10}`: `0.3247598`; of `{17,...,91}`: `0.3483833`. Top dyadic range `[2^28, 2^29)`: `0.326787 / 0.324876 / 0.348337`; the first basin drifts down (`0.3277` at `2^20`), the second up (`0.3243`), the third is flat to four decimals from `2^20` on. The three are not equal: the seven-odd-step cycle owns the largest basin. First ray element hit (odd cycle element `c`, exponent `k`, i.e. the entering odd number `(c 2^{k+1} + 1)/3`): `7*2^2 = 28` (from 19) `0.1984`; `61*2^2 = 244` (from 163) `0.1914`; `2^6` (from 43) `0.1878`; `2^4` (from 11) `0.1336`; `7*2^8` (from 1195) `0.0526`; `17*2` (from 23) `0.0460`; `25*2^2` (from 67) `0.0433`; `37*2^4` (from 395) `0.0335`; `5*2^5` (from 107) `0.0282`; `5*2^7` (from 427) `0.0264`; the ray of `1` is entered only at `2^k` with `k = 4, 6, 10, 12, 16, ...`, never at `k = 2, 8, 14` where `(2^{k+1}+1)/3` is a multiple of three.

`3x+1` on `[1, 2^30]` (`..._entries_plus_2e30.out`), first trunk element `2^{2i-1}` entered from `a_i`: `e_2 = 0.93796` (through `5`), `e_4 = 0.023647` (`85`), `e_5 = 0.037789` (`341`), `e_7 = 8.0e-5` (`5461`), `e_8 = 4.85e-4` (`21845`), `e_10 = 3.2e-5`, `e_11 = 2.1e-6`, `e_13 = 4e-7`, `e_14 = 1e-7`; `e_3 = e_6 = e_9 = e_12 = 0` exactly (3.2(ii)). Top range `[2^29, 2^30)`: `0.93794 / 0.02366 / 0.03780`, so `dens B(5) ~ 0.938` and `dens B(32) ~ 0.062`; `341` owns a larger basin than `85`. The tail `sum_{i >= 11} e_i` is below `3e-6` at this scale but is not converged (the basin of `a_i` is invisible below `~ 4^i`).

**3.4 The dual.** These entry statistics live entirely on the loops of the dual graph; the `K`-bond (the necklace) sees none of them, and conversely the necklace lane of §1 sees no basin. The two lanes are the two halves of the owner's picture.

## 4. Hypotheses filed

* [HYP-9165](../hypotheses/HYP-9165-basin-densities-exist-sheet-cap.md): the natural densities of the three `3x-1` basins and of `B(5)` exist, with the values above; in particular `dens B({5,7,10}) < dens B({1}) < dens B({17,...})`, so the sheet cap for cycle-uniform sheet-blind lower bounds is `0.3248`, not `1/3`.

## 5. Reproduction

```
gcc -O3 -o basins_minus  04-computation/experiments/collatz_necklace_20260929_basins_minus.c  && ./basins_minus 536870912
gcc -O3 -o entries_plus  04-computation/experiments/collatz_necklace_20260929_entries_plus.c  && ./entries_plus 1073741824
python3 04-computation/experiments/collatz_necklace_20260929_fairsplit_census.py 18 criteria
python3 04-computation/experiments/collatz_necklace_20260929_power_clocks.py 201 400
python3 04-computation/experiments/collatz_necklace_20260929_power_clocks.py 3 3000
python3 04-computation/experiments/collatz_necklace_20260929_circulant_check.py 100 100
```
(C sieves: `unsigned __int128` orbits, memoisation on the first value below the start; 1 GB for `2^30`; about a minute each. Python: `sympy` for cyclotomic reductions, `gmpy2` for perfect-power tests.)

## 6. Stopping boundary and next questions

* Does the CRT split of 1.3 combine with size (`n_min >= 2^68`, `|c_0 - c_1| < 2^m q^x`) or with the fair Eliahou identity to give a condition not already implied by Baker-type bounds on `|2^m - 3^x|`? Not attempted; the two factors have very different sizes (`2^m - 3^x` tiny at convergents, `2^m + 3^x ~ 2^{m+1}`), which is the natural next probe.
* Is the fraction of primitive necklaces of shape `(jm, jx)` admitting a fair `j`-split asymptotically a function of `j` and `x/m` only? The census suggests decay in `j` (`0.68`, `0.90`, `0.96`, `0.98` failing for `j = 4, 6, 8, 9` at `K = 16, 18`); a random-walk computation of the level criterion 1.4(a) would give the law.
* Perfect-power clocks for `q` composite with several prime factors are Fermat–Catalan only after rewriting; a census over all odd `q <= 10^4`, `K <= 60` (`X >= 2`) would test whether anything beyond the family and the four solutions appears — a cheap hostile probe of the Fermat–Catalan conjecture in this window.
* The trunk-entry spectrum `e_i` and the ray-entry spectra on the minus sheet are natural targets for the harmonic-mass formalism of the [Mazur digest](mazur_positive_density_20260928.md): `e_i` should be the harmonic mass of the Syracuse tree of `a_i`, and the observed ordering `e_5 > e_4` is a 3-adic statement about `341` versus `85`.

## 7. Cards used

Inheritance pass (THM-4484 as the mechanism, the minus sheet as the hostile example); "objects x operations": the necklace under the *split* operation instead of rotation; "type the connection" (source: necklace splitting; target: cycle equation; map: one arc per thief; preserved: shape; destroyed: everything about the trees; sidecar: the DFT; decisive test: exact identities on real cycles); "cheap hostile probe" (the perfect-power census, the criteria cross-check, the missing fair 4-split of the 17021 cycle). Candidate card, not promoted: *fairness = circulant symmetry* — whenever a cyclic object is cut into equal-shape consecutive pieces, the composed dynamics is circulant and the DFT factors the global invariant; counterindication: unequal shapes (no diagonalisation), non-abelian symmetry.
