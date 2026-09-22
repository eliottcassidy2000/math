# The 3n-1 sheet as a positive control: which descent statements survive the sign flip

**Status: PROVED (scoped, elementary): the parity-word map mod `2^J` is a bijection for both signs and the prefix-descent counts of the two sheets coincide for every `J` (S2); coefficient descent implies actual descent on the minus sheet and the converse holds on the plus sheet (S3); the residue reduction lemma for strict early descents (S3); the `B3` sheet-swapping bijection between `k=1` edges of one sheet and `k=2` edges of the other (S6). FINITE-EXACT: the `T_-` basin census on `[1, 10^7]` with no escape and no new cycle (S1); the prefix-descent table to `J = 20` (S2); no strict early descent for any positive `n` whose first descent time is `<= 20`, and none for `n <= 10^6` (S3); greedy `G_-` certification of `666665/666666` starts to `10^6` with the single hard start `m = 4` and BFS-shortest `E_-` paths from `1` to all nine odd cycle members (S4); the Berggren child/parent readings of all eleven cycle edges (S5); exact budget sums `sum q_i = 3 n_0` on fifteen starts (S7). HEURISTIC: the three basin fractions near `0.327`, `0.325`, `0.348`. OPEN: whether `sigma = sigma_c` holds on the minus sheet for all `n` (the mirror of Terras' coefficient stopping-time conjecture); `Q2_-` beyond `10^6`. SCOPE: the lane recalls no published basin fractions for `3n-1` with enough confidence to quote (UNCITED-RECOLLECTION declined); no Collatz convergence statement is claimed. REFUTED (as a proof mechanism): every density-one, ergodic, spectral or graph-reachability step of the pasted blueprint, because each survives the sign flip and the `3n-1` sheet has three cycles (S8).** Session collatz-mod6-20260922, lane `minus_sheet_positive_control`, script `04-computation/experiments/collatz_mod6_20260922_minus_sheet_positive_control.py`, output `collatz_mod6_20260922_minus_sheet_positive_control.out`.

## Inheritance and concept board

The control object is the map `T_-(n) = (3n-1)/2^v` on positive odd `n` (the `T`-form `n -> (3n-1)/2` on odd, `n/2` on even, is used for stopping times), whose three cycles with odd skeletons `(1)`, `(5,7)`, `(17,25,37,55,41,61,91)` are [arithmetic_braids2_20260917_signed_cycles](arithmetic_braids2_20260917_signed_cycles.md) section 5, with the gate `n_0 = B/(3^L - 2^K)` of [arithmetic_braids_20260917_collatz](arithmetic_braids_20260917_collatz.md) (B7). The closest proved mechanism is the sign-blind parity-vector law: Terras' residue determination of the halving word (CITED as recollection in [extended_collatz_scc](collatz_mod6_20260917_extended_collatz_scc.md) section 4, where the 3-adic mirror Theorem 4.4 and its conjugate Theorem 6.2 for `E_-` live), which S2 restates with a sign-independent proof. The canonical hostile is the minus sheet itself: every statement in the [synthesis](collatz_mod6_20260917_synthesis.md) sections 2-3 and 10 that is quantified over densities holds there verbatim, and the sheet does not converge to `1`. The corrected near miss is the lane's own first reading of the stopping-time comparison: on the plus sheet the carry `B > 0` gives `sigma >= sigma_c` and the equality is Terras' *conjecture* (UNCITED-RECOLLECTION: Terras 1976, surveyed by Lagarias 1985), not a theorem; on the minus sheet `B < 0` gives `sigma <= sigma_c` and the equality is the mirror conjecture, verified here for all first-descent times `<= 20` by a reduction to residues (S3). The least-used sidecar is the additive budget `sum q_i <= 3 n_0` of [collatz_guards_20260921_discrepancy](collatz_guards_20260921_discrepancy.md) (D8) and [glued_xor_20260921_blueprint](glued_xor_20260921_blueprint.md) (B9)-(B10): it is the only inherited object that is genuinely sign-specific, and S7 shows the three cycles exhaust it exactly, so even it does not exclude cycles. Also inherited by path: the Berggren transport table and its five sporadic readings ([berggren_edge_transport](collatz_mod6_20260921_berggren_edge_transport.md) Theorems 2.1 and 3.1), the `G_-` 2-cycle `{4,11}` and the `Q1_-`/`Q2_-` verification to `10^6` ([extended_collatz_scc](collatz_mod6_20260917_extended_collatz_scc.md) section 6), the `G` negatives lane ([g_negatives_joint_carry](collatz_mod6_20260921_g_negatives_joint_carry.md), whose S9-S10 are the `G`-side mirror of S7), and the blueprint audit ([collatz_blueprint_20260921_synthesis](collatz_blueprint_20260921_synthesis.md), quantifier gap). Canon THM-1858, THM-1370, THM-1745, THM-3756, THM-3341 were read for the paste's tournament vocabulary; no map to this lane's content was found (SCOPE).

Session lead probe (verified, S0): `8/pi^2 = 0.810569` has binary entropy `0.70028` bits, not `0.704`, and a positive entropy is not "zero entropy"; nothing else in the paste's section 3 is testable.

## 1. Basin census of `T_-` on `[1, 10^7]` (S1)

**S1 (FINITE-EXACT).** With the `T`-form map on all `n <= 10^7` (vectorised, first value below `n` recorded, then pointer-jumping to the cycle minimum), every `n` lands in the basin of `1`, `5` or `17`; no orbit exceeds the step cap `5000` (the largest first-descent time is `445`), no orbit value exceeds `30541433029400`, and the only `n <= 10^7` whose orbit never drops below `n` are exactly the three cycle minima `1, 5, 17`. No escape and no fourth cycle with minimum `<= 10^7`. The plus-sheet control run gives basin `{1}` for every `n <= 10^7`, largest first-descent time `401`, largest value `60342610919632`, and `1` as the only non-descending start.

```text
X          | basin{1}  all n   odd n | basin{5,7} all n  odd n | basin{17..91} all n  odd n
10000      |      3244 0.32440     1605 0.32100 |      3213 0.32130     1623 0.32460 |      3543 0.35430     1772 0.35440
100000     |     33030 0.33030    16553 0.33106 |     32104 0.32104    16026 0.32052 |     34866 0.34866    17421 0.34842
1000000    |    327679 0.32768   163486 0.32697 |    323351 0.32335   162122 0.32424 |    348970 0.34897   174392 0.34878
10000000   |   3273791 0.32738  1636054 0.32721 |   3244985 0.32450  1623149 0.32463 |   3481224 0.34812  1740797 0.34816
```

HEURISTIC: the fractions drift slowly (`0.324 -> 0.327` for `{1}`, `0.354 -> 0.348` for the 17-cycle over three decades) and the 17-cycle basin stays the largest; no limit is claimed. The all-`n` and odd-`n` columns agree to about `2e-4` because the basin of `2n` is that of `n` (trivial). Literature comparison: SCOPE, the lane declines to quote basin fractions it cannot cite.

## 2. Terras-type prefix-descent densities are sign-blind (S2)

Let `T_b(n) = (3n+b)/2` for odd `n`, `n/2` for even `n`, `b = +-1`. For a residue `r mod 2^J` let `w = (p_1, ..., p_J)` be the parity word of the first `J` steps and `a_i` the number of odd steps among the first `i`. Say the word has a *prefix descent* if `3^(a_i) < 2^i` for some `i <= J` (the coefficient stopping time is `<= J`).

**S2 (PROVED; no novelty claim, the argument is Terras' with the sign carried along).** For both signs the word of length `J` is a function of `n mod 2^J`, and the induced map `Z/2^J -> {0,1}^J` is a bijection. Hence, for every property of words, the number of residues mod `2^J` whose word has the property is the same on the two sheets; in particular the prefix-descent counts coincide for every `J`.

*Proof.* `T_b(n) mod 2^(J-1)` is determined by `n mod 2^J` on either branch (the branch is chosen by `n mod 2`, and `(3n+b)/2`, `n/2` are determined mod `2^(J-1)` by `n mod 2^J`); induction gives the word. If `n, n'` have the same word of length `J` with `a` odd steps, then `T_b^J(n) - T_b^J(n') = 3^a (n-n')/2^J`, the carries cancelling because they depend only on the word; the left side is an integer and `3^a` is odd, so `2^J | n - n'`. Injective on `Z/2^J`, hence bijective, and the count of residues with a word property equals the number of words with it, independent of `b`. QED

FINITE-EXACT confirmation to `J = 20` (both sheets computed separately, plus a DP count of words with no prefix descent; all three agree at every `J`):

```text
J  | #descend(+) #descend(-) | words(+) words(-) = 2^J | no-descent words (DP) | fraction
1  |           1           1 |        2        2        2 |          1 | 0.500000
2  |           3           3 |        4        4        4 |          1 | 0.750000
3  |           6           6 |        8        8        8 |          2 | 0.750000
4  |          13          13 |       16       16       16 |          3 | 0.812500
5  |          28          28 |       32       32       32 |          4 | 0.875000
6  |          56          56 |       64       64       64 |          8 | 0.875000
7  |         115         115 |      128      128      128 |         13 | 0.898438
8  |         237         237 |      256      256      256 |         19 | 0.925781
9  |         474         474 |      512      512      512 |         38 | 0.925781
10 |         960         960 |     1024     1024     1024 |         64 | 0.937500
11 |        1920        1920 |     2048     2048     2048 |        128 | 0.937500
12 |        3870        3870 |     4096     4096     4096 |        226 | 0.944824
13 |        7825        7825 |     8192     8192     8192 |        367 | 0.955200
14 |       15650       15650 |    16384    16384    16384 |        734 | 0.955200
15 |       31473       31473 |    32768    32768    32768 |       1295 | 0.960480
16 |       63422       63422 |    65536    65536    65536 |       2114 | 0.967743
17 |      126844      126844 |   131072   131072   131072 |       4228 | 0.967743
18 |      254649      254649 |   262144   262144   262144 |       7495 | 0.971409
19 |      509298      509298 |   524288   524288   524288 |      14990 | 0.971409
20 |     1021248     1021248 |  1048576  1048576  1048576 |      27328 | 0.973938
```

Consequence (the decisive hostile). Any argument whose input is the density of residues with a descending prefix, or the ergodic frequency of a residue event along orbits, reads the same numbers on the `3n-1` sheet. Applied there it would "prove" that all orbits converge to `1`, which S1 refutes at `n = 5`. The paste's section 4(3) ("ergodically forced to intersect the `0 mod 3` clock boundaries") is of this kind: it is sign-blind, so it cannot be a proof of anything sign-specific. This is not a new refutation of the blueprint (see the blueprint lanes cited above); it is the shortest one.

## 3. Actual versus coefficient stopping time: the two sheets face opposite ways (S3)

Write `sigma(n) = min{k : T_b^k(n) < n}` and `sigma_c(n) = min{k : 3^(a_k) < 2^k}`, with `T_b^k(n) = (3^(a_k) n + B_k)/2^k`, `B_k` the word carry, `B_k > 0` for `b = +1` and `B_k < 0` for `b = -1`.

**S3 (PROVED, elementary).** (i) Plus sheet: `sigma >= sigma_c`, since `3^(a_k) >= 2^k` and `B_k > 0` give `T^k(n) > n`. (ii) Minus sheet: `sigma <= sigma_c`, since `3^(a_k) < 2^k` and `B_k < 0` give `T^k(n) < n`. (iii) Reduction lemma (minus sheet): if some `n >= 1` has a *strict early descent* at time `k`, i.e. first descent at `k` with `3^(a_k) >= 2^k`, then some residue `r` in `[1, 2^k)` has a strict early descent at some time `<= k`.

*Proof of (iii).* The word `e^k` has coefficient `2^(-k) < 1`, so `n` is not `0 mod 2^k`; write `n = r + 2^k t`, `r` in `[1, 2^k)`, `t >= 0`; by S2 `r` has the same word of length `k`. Since `T^i` is affine on the class with slope `3^(a_i)/2^i`, `T^i(n) - n = [T^i(r) - r] + 2^(k-i) t (3^(a_i) - 2^i)` for `i <= k`. At `i = k` the left side is negative and the second term is `>= 0`, so `T^k(r) < r`. If `r` first descends at some `i < k`, then `T^i(n) - n >= 0` forces `3^(a_i) > 2^i`, a strict early descent of `r` at time `i`; otherwise `r` first descends at `k` with `3^(a_k) >= 2^k`. QED

**FINITE-EXACT.** For `n <= 10^6` on both sheets (`400`-step cap): the only starts with infinite `sigma` are `1` (plus) and `1, 5, 17` (minus); the count of `n` with `sigma < sigma_c` is `0` on both sheets, and no `n` has `sigma > sigma_c` on either. Residue level, minus sheet, `K = 20`: of the `1021247` residues `r` in `[1, 2^20)` that descend within `20` steps, all `1021247` have `sigma = sigma_c`, and the number of strict early descents is `0`. By (iii) this proves: **every positive `n` whose first `T_-`-descent occurs at time `k <= 20` has `3^(a_k) < 2^k`** (a statement about all `n`, not only `n <= 10^6`).

Reading. On the plus sheet the equality `sigma = sigma_c` for all `n > 1` is Terras' coefficient stopping-time conjecture (UNCITED-RECOLLECTION; it is known to imply the absence of nontrivial cycles), verified here to `10^6`. On the minus sheet the mirror equality is compatible with three cycles (their minima have `sigma = sigma_c = infinity`, the other members descend with `3^(a_k) < 2^k`, e.g. `7 -> 10 -> 5`), so the mirror conjecture carries no cycle-exclusion content: OPEN, and not a route to convergence. What the finite result does say is that on the minus sheet no orbit ever "beats the coefficient" using the negative carry before step `20`, i.e. the `B_k` term never decides a descent early; the same statement on the plus sheet is (i), a triviality of sign.

## 4. `E_-` reachability from `1` (S4)

`E_-` has arrows `n -> n/2` (even `n`) and `n -> 3n-1` (all `n`); `Q2_-` asks whether `1` reaches every `m` with `3 !| m` (multiples of `3` are transient singletons, [extended_collatz_scc](collatz_mod6_20260917_extended_collatz_scc.md) Theorem 1.2 and 6.1). The greedy inverse map is `G_-(m) = (2^k m + 1)/3` with `k = k_-(m mod 9)` from the table `(1,1),(2,0),(4,3),(5,0),(7,1),(8,2)` (checked minimal and landing in `{2,5} mod 9`).

**S4 (FINITE-EXACT).** Over the `666666` starts `2 <= m <= 10^6`, `3 !| m`: `666665` are certified by greedy chaining (the `G_-`-orbit drops below `m`, then induction), fraction `0.999998`; `0` hit `1` directly; `1` start cycles: `m = 4` on the 2-cycle `{4, 11}` (`G_-(4) = 11`, `G_-(11) = 4`); `0` hit the `2000`-step cap. Steps-to-certify histogram: `1: 444444`, `2: 148148`, `3: 24692`, `4: 24691`, ..., `35: 1` (full list in the `.out`; `444444/666666 = 2/3` certify in one step, which is the fraction of units mod 9 with `k_- <= 1` and a drop, HEURISTIC reading). Largest greedy peak before certification: `50331647 = 3*2^24 - 1` at `m = 797161`, the exact mirror of the plus-side compound peak reported in the cited SCC note at `m = 797162`. BFS fallback over single inverse arrows `x -> 2x`, `x -> (x+1)/3` rescues `m = 4` in `14` arrows: `4 -> 8 -> 16 -> 32 -> 64 -> 128 -> 43 -> 86 -> 29 -> 10 -> 20 -> 7 -> 14 -> 5 -> 2` (the cited note's compound rescue `4 -> 11 -> 59 -> 20 -> 7 -> 5 -> 2` is another `14`-arrow path). So `Q2_-` holds to `10^6` (reconfirming the cited note) and remains OPEN globally; the greedy-certified fraction is `1 - 1/666666`.

BFS-shortest forward `E_-` paths from `1` to every odd member of the three cycles (values capped at `5*10^6`):

```text
1 -> 5:   2 arrows: 1 -> 2 -> 5
1 -> 7:   4 arrows: 1 -> 2 -> 5 -> 14 -> 7
1 -> 41:  4 arrows: 1 -> 2 -> 5 -> 14 -> 41
1 -> 61:  6 arrows: 1 -> 2 -> 5 -> 14 -> 41 -> 122 -> 61
1 -> 91:  8 arrows: 1 -> 2 -> 5 -> 14 -> 41 -> 122 -> 61 -> 182 -> 91
1 -> 17: 13 arrows: 1 -> 2 -> 5 -> 14 -> 41 -> 122 -> 61 -> 182 -> 91 -> 272 -> 136 -> 68 -> 34 -> 17
1 -> 25: 15 arrows (via 17 -> 50 -> 25);  1 -> 37: 17 arrows;  1 -> 55: 19 arrows (along the cycle)
```

The paths enter the 5-cycle through the new arrow `2 -> 5` and the 17-cycle through `14 -> 41`: the even-to-`3n-1` arrows are exactly what makes the cycles reachable from `1` (and, by the cited `Q1_-`, what makes `1` reachable from them). This is the graph-theoretic form of the control: `E_-` has the same one-giant-SCC picture as `E` while `T_-` has three cycles, so the SCC picture does not see cycles at all.

## 5. Berggren `B3` sheet coupling of the cycle edges (S5, S6)

Edges are odd-skeleton steps `x -> y` with `3x + b = 2^k y`; children and readings follow the inherited Theorem 2.1 of the Berggren lane (case A `y > x`: `B1 = (y+2x, x)`, `B2 = (2y+x, y)`, `B3 = (2y-x, y)`; case B the mirror), a pair `(p,q)` being read as a `(3, b', j)` edge whenever `3s + b' = 2^j t` for `(s,t) = (p,q)` or `(q,p)`.

**S5 (FINITE-EXACT; every reading below is an exact identity).**

```text
edge (sheet, cycle-min): x->y k | children B1,B2,B3 with (3,+-1) readings | Berggren parent pair and its readings
(plus,1):   1->1  k=2 | B1=(3,1):3->1(minus,k=3) B2=(3,1):3->1(minus,k=3) B3=(1,1):1->1(plus,k=2),1->1(minus,k=1) | parent via B3-fixed: (1,1)
(minus,1):  1->1  k=1 | same pair (1,1): the two trivial cycles are one Berggren pair, its own B3 child
(minus,5):  5->7  k=1 | B1=(17,5):none B2=(19,7):19->7(minus,k=3) B3=(9,7):9->7(plus,k=2)     | parent via B3: (5,3):3->5(plus,k=1)
(minus,5):  7->5  k=2 | same pair (7,5): the sporadic reading E1 of the Berggren lane            | parent via B3: (5,3):3->5(plus,k=1)
(minus,17): 17->25 k=1 | B1=(59,17):none B2=(67,25):67->25(minus,k=3) B3=(33,25):33->25(plus,k=2) | parent via B3: (17,9):none
(minus,17): 25->37 k=1 | B1=(87,25):none B2=(99,37):99->37(minus,k=3) B3=(49,37):49->37(plus,k=2) | parent via B3: (25,13):none
(minus,17): 37->55 k=1 | B1=(129,37):none B2=(147,55):147->55(minus,k=3) B3=(73,55):73->55(plus,k=2) | parent via B3: (37,19):none
(minus,17): 55->41 k=2 | B1=(137,41):none B2=(151,55):none B3=(69,55):none                        | parent via B3: (41,27):27->41(plus,k=1)
(minus,17): 41->61 k=1 | B1=(143,41):none B2=(163,61):163->61(minus,k=3) B3=(81,61):81->61(plus,k=2) | parent via B3: (41,21):none
(minus,17): 61->91 k=1 | B1=(213,61):none B2=(243,91):243->91(minus,k=3) B3=(121,91):121->91(plus,k=2) | parent via B3: (61,31):none
(minus,17): 91->17 k=4 | B1=(125,17):none B2=(199,91):none B3=(165,91):none                       | parent via B1: (57,17):none
```

Counts: the minus-sheet cycles have `10` edges on `9` distinct unordered pairs (`k=1: 7`, `k=2: 2`, `k=4: 1`); the plus-sheet cycle has `1` edge (`k=2`). Edges with a child on the other sheet: `9` of `11` (all seven `k=1` minus edges via `B3`, plus both trivial edges through the fixed pair `(1,1)`). Edges that are Berggren children of an other-sheet edge: `5` of `11` (`7->5` and `55->41`, the two `k=2` minus edges, are `B3` children of the plus-sheet `k=1` edges `3->5` and `27->41`; `5->7` shares the pair of `7->5`; the two trivial edges are the fixed pair). The only uncoupled edge is `91 -> 17` (`k=4`, parent `(57,17)` reads as nothing). The `B2` children of the `k=1` minus edges are same-sheet `k=3` edges (`19->7`, `67->25`, ...), as Theorem 3.1 of the Berggren lane predicts.

**S6 (PROVED; the general law behind the table).** For `b = +-1`, `B3` restricted to the `k=1` edges of the sheet `b` is a bijection onto the `k=2` edges of the sheet `-b`, given by `x -> y  |->  (2x+b) -> y`. *Proof.* A `k=1` edge has `y > x`, so its `B3` child is `(2y-x, y) = (2x+b, y)`, and `3(2x+b) - b = 2(3x+b) = 4y`: a `k=2` edge of sheet `-b`. Conversely a `k=2` edge `u -> y` of sheet `-b` has `3u - b = 0 mod 4`, i.e. `u = b mod 4`... explicitly `u = 3 mod 4` for `b = +1` and `u = 1 mod 4` for `b = -1`, so `x = (u-b)/2` is odd and `3x + b = (3u - b)/2 = 2y` is a `k=1` edge of sheet `b` whose `B3` child is `(u, y)`; the two constructions are inverse. QED. FINITE-EXACT check on all `k=1` edges with `x <= 10^5` of both sheets (`25000` each, `x = 3 mod 4` on the plus sheet, `x = 1 mod 4` on the minus sheet): the `B3` images are exactly the `k=2` edges of the other sheet with source `<= 2*10^5 - 1`. Exact sheet coupling through the tree is therefore: every `k=1` edge of either sheet has one `k=2` child on the other sheet, every `k=2` edge is such a child, and `k >= 3` edges are coupled only through the five sporadic readings of the inherited Theorem 3.1 (none of which touches `91 -> 17`).

## 6. The additive budget is exhausted exactly by the cycles (S7)

**S7 (FINITE-EXACT, exact rationals; the identity is inherited (B10) and not re-proved).** With `q_i = 2^(K_i)/3^i` along the odd skeleton, `sum_(i>=0) q_i = 3 n_0` for `n_0 = 1, 3, 5, 7, 9, 11, 17, 25, 37, 41, 55, 61, 91, 27, 1001` (sums `3, 9, 15, 21, 27, 33, 51, 75, 111, 123, 165, 183, 273, 81, 3003`), the period ratios being `2/3`, `8/9`, `2048/2187 < 1`; on the plus sheet the cycle `{1}` has ratio `4/3 > 1` and a divergent sum (`D_L -> +infinity`). Gates: `n_0 = B/(3^L - 2^K)` with `(L, K, B) = (1, 1, 1)`, `(2, 3, 5)`, `(7, 11, 2363)`, i.e. `1/1`, `5/1`, `2363/139`. So the sign-specific budget of (D8)/(B9) holds on the sheet with cycles and is exhausted by them; "budget exhaustion" is a periodicity criterion there, not a convergence criterion.

## 7. Descent statements versus the minus sheet (S8)

**S8 (each row PROVED, FINITE-EXACT or CITED as marked).**

```text
statement (source by path)                                   | on the 3n-1 sheet                          | sign-specific? | usable alone for convergence?
density-1 prefix-descent / Terras (SCC note sec. 4, S2 here)  | identical counts, every J (PROVED S2)       | no             | no: the sheet has cycles
3-adic greedy density-1 stopping (SCC Thm 4.4, 6.2)           | holds verbatim by conjugation (CITED)       | no             | no; G_- even has the 2-cycle {4,11}
no bounded strip for q_j (guards sec. 2a)                     | q_j -> 0 on every orbit incl. cycles (D8)   | conclusion no  | no
dichotomy: periodic iff D_L -> +infinity (glued sec. 3)       | false as dichotomy: cycles have D_L -> -inf | yes (PROVED)   | no: an equivalence, both branches open
budget sum q_i <= 3 n_0 (guards D8, glued B9-B10)             | holds; exhausted exactly by the cycles (S7) | yes            | no: it is a periodicity test
E-graph one-SCC picture / Q1, Q2 (SCC note secs. 3, 6)        | same picture: Q1_-, Q2_- hold to 10^6 (S4)  | no             | no: E_- erases the cycles
G / G_- greedy certificates (SCC sec. 6, g_negatives)         | 666665/666666 certified (S4)                | no             | no
cycle gate n_0 = B/(3^L - 2^K) (braids B7; S7)                | same formula; sign law 2^K < 3^L            | sign law yes   | no: a cycle test, not a convergence test
Berggren transport (berggren lane; S5-S6)                     | B3 swaps sheets; k=1 <-> k=2 bijection      | no (it mixes) | no
coefficient = actual stopping time                            | mirror conjecture; exact for sigma <= 20    | direction yes  | no: compatible with three cycles
paste sec. 1-2 (spectra, Fano/Paley), sec. 4 chain            | same matrices and graphs exist for 3n-1     | no             | no (already REFUTED in the blueprint lanes)
```

The target inequality of the paste's section 4, "every odd `n > 1` has `L` with `3^L n + B_L < 2^(K_L) n`", is exactly the stopping-time formulation; on the minus sheet it fails for `n = 5` and `n = 17` only (among `n <= 10^7`, S1) and the density rows above cannot see those two integers. Any correct proof of Collatz must use an ingredient in the "yes" column, and the two proved sign-specific ingredients (the direction of `D_L` on cycles, the budget) are periodicity criteria that hold on the sheet with cycles. The control therefore isolates the missing object precisely: a sign-specific statement that is *not* a periodicity criterion, i.e. one that fails on the minus sheet at `5` or `17` while holding for every plus-sheet integer.

Typed control. Source: the `3n+1` sheet (`T_+`, `E`, `G`). Target: the `3n-1` sheet (`T_-`, `E_-`, `G_-`). Map: `n -> -n`, i.e. `T_-(n) = -T_+(-n)`, equivalently `b -> -b`. Preserved: the parity-word bijection and every residue count (S2), the 3-adic drift table up to `r -> 9-r` (CITED), the cycle gate formula, the Berggren tree with `B3` exchanging the sheets (S6), the `E`-reachability picture (S4). Lost: the sign of the carry `B` and with it the direction of `sigma` versus `sigma_c` (S3), the direction of `D_L` on cycles, and the convergence of `sum q_i` (S7). Sidecar: the sign of `B`, equivalently the sign law `2^K` versus `3^L` on cycles. Test: the three cycles; an argument that commutes with the map proves nothing sign-specific.

## 8. Reproduction block

```text
cd <worktree>
python3 04-computation/experiments/collatz_mod6_20260922_minus_sheet_positive_control.py \
    > 05-knowledge/results/collatz_mod6_20260922_minus_sheet_positive_control.out
python3 -O 04-computation/experiments/collatz_mod6_20260922_minus_sheet_positive_control.py 2>/dev/null | diff - \
    05-knowledge/results/collatz_mod6_20260922_minus_sheet_positive_control.out   # identical
```

Runtime about `8` s, peak RSS about `0.83` GB (numpy `int64` arrays of length `10^7`, freed between the two censuses); explicit `raise` only; timing lines go to stderr; the run ends with `ALL CHECKS PASSED`. Universes: `n <= 10^7` (S1), residues mod `2^J`, `J <= 20` (S2, S3 residue level), `n <= 10^6` with a `400`-step cap (S3), `m <= 10^6` with a `2000`-step greedy cap and a `10^6` BFS value cap (S4), forward BFS values `<= 5*10^6` (S4), `k=1` edges with `x <= 10^5` (S6). Source sha256 `455ba7149fcb3fd9d91317fc0c6cc7a4114b16e56b498d771dc9eac102df3165`, output sha256 `ccfd33c649c496dcd92b574e8cbe04df277d9c47aa5f7364481a10be9a7cae02` (computed after the run; not printed inside the `.out`).

## Stopping boundary / next question

The control is now exact: of everything the session has proved about descent, only the direction of `D_L` on cycles and the additive budget distinguish the sheets, and both are periodicity criteria that the three `3n-1` cycles satisfy. The lane stops at the mirror coefficient stopping-time statement (S3): it is exact for every first-descent time `<= 20` on the minus sheet by the residue reduction, and the same reduction bounds the search for a counterexample at any `k` to residues below `2^k`, so the next cheap question is whether the reduction can be closed for all `k` by a word-level inequality `|B_w| <= (3^a - 2^k) r_w` (with `r_w` the least positive minus-sheet representative of the word), which would prove `sigma = sigma_c` on the minus sheet outright; if it fails at some `k > 20`, the witness is the first orbit on either sheet whose carry decides a descent before its coefficient does. The other open item left untouched is `Q2_-` beyond `10^6`, where the only known obstruction is the greedy 2-cycle `{4,11}`.
