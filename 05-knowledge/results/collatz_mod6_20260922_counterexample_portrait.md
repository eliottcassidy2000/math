# Counterexample portrait: every proved necessary condition on a hypothetical 3n+1 cycle or divergent orbit, tested on the 3n-1 sheet

**Status: PROVED (elementary, from CITED inputs): the sheet criterion `T_+(-n) = -T_-(n)`; the odd/even-convergent refinement of the sign law; the length bound `L > 14878203146` for every nontrivial positive `3n+1` cycle from Barina `2^68` + Legendre + Hardy-Wright 171 alone; its minus-sheet mirror `L > 2738` for any fourth positive `3n-1` cycle from a new `10^7` census; the minus-sheet strip exclusion by positivity alone. FINITE-EXACT: the `10^7` three-basin census, gate censuses (`L <= 10`, 1-cycles `L <= 3000`, 2-cycles `L <= 24`), Terras densities in both conventions, squarefree prefix counts, all checks on the three minus cycles. CITED: Barina, Eliahou, Steiner, Simons-de Weger, Hercher (abstract read), Terras, Everett, Tao, Krasikov-Lagarias, Mihailescu, Legendre/Hardy-Wright. REFUTED (as sheet-transfers): Steiner's no-1-cycle (witness `{5,7}`), Simons-de Weger/Hercher no-`m`-cycle (witness the seven-cycle, `m=2`), "`D_L -> -infinity` separates cycles from escape" (all three minus cycles have `D_L -> -infinity`). Pairwise contradiction among the compiled conditions: NONE FOUND, on either sheet. The Collatz conjecture, completeness of the minus-sheet cycle list, and the existence question for divergent orbits remain OPEN.**

Session collatz-mod6 (machine `mac-mini`), lane `counterexample_portrait`, 2026-09-22. Script [collatz_mod6_20260922_counterexample_portrait.py](../../04-computation/experiments/collatz_mod6_20260922_counterexample_portrait.py), output [collatz_mod6_20260922_counterexample_portrait.out](collatz_mod6_20260922_counterexample_portrait.out). Claims are numbered S1-S15 to match the `.out` sections.

## Inheritance and concept board

Closest proved mechanism: the ordered cycle gate `n_0 = bB/Delta`, `q = |Delta|/gcd(B,|Delta|)`, integral iff `q | b`, with the sign law "every member has the sign of `b/Delta`" ([signed_cycles sec.3](arithmetic_braids2_20260917_signed_cycles.md), [catalan gate theorem](catalan_elliptic_20260921_catalan.md)), composed with the Eliahou product identity `2^K/3^L = prod(1 + b/(3n_i))` and the three-tier Legendre placement ([pillai S3-S4](collatz_mod6_20260921_pillai_convergents_cycle_gates.md)). Canonical hostile: the positive `3n-1` seven-cycle `{17,25,37,55,41,61,91}`, which satisfies every sheet-blind condition in the two tables below (gate `q=1`, product identity, necklace count, `D_L -> -infinity`, no bounded strip, Terras cylinder law, squarefree prefixes, Lyapunov obstruction, `G`-mirror) and therefore certifies that no conjunction of sheet-blind conditions can exclude a `3n+1` cycle. Corrected near miss (two): first, the glued_xor dichotomy "`D_L -> -infinity` iff the orbit escapes" is a plus-sheet statement only; on the minus sheet every orbit, periodic or not, has `D_L -> -infinity` and the discriminator is budget exhaustion `sum q_L = 3n_0` ([glued_xor sec.3](glued_xor_20260921_synthesis.md), verified exactly here on all three cycles: `3, 15, 51`); second, this lane's own first pass used a lazy factor `0.75` in the convergent-spacing inequality and left the clock `n=21` (`q_21 = 6586818670`) open below `L_A`; the sharp factor `e^{-L_A/(3N)}` closes it and gives the clean bound `L > 14878203146`. Least-used sidecar: the parity of the convergent index. Even-index convergents of `log_2 3` lie below it (`Delta < 0`), odd-index above (`Delta > 0`), so the sign law becomes "positive `3n+1` cycles live on odd-index convergents, positive `3n-1` cycles on even-index ones" wherever Legendre placement applies; the minus cycles sit on `p_0/q_0 = 1/1`, `p_2/q_2 = 3/2` and the extreme intermediate `11/7`, the plus cycle on `p_1/q_1 = 2/1`.

Typed analogy (the only map used). Source: the negative `3n+1` sheet. Target: the positive `3n-1` sheet. Map: `n -> -n`. Preserved: the accelerated map (`T_+(-n) = -T_-(n)`, S1), halving words, carries up to sign, the gate, `q`, `Delta`, `K/L`, cylinder densities, the greedy map `G` (S14). Lost: nothing (it is a conjugacy). Sidecar: the sign of the carry `bB` relative to `Delta`, i.e. which side of `log_2 3` the clock must lie on. Test: every condition below is checked on the three minus cycles; a condition is SHEET-BLIND if its proof never uses `n > 0` at `b = +1`, SIGN-SPECIFIC otherwise.

Concurrent notes cited by path, not re-proved: [synthesis](collatz_mod6_20260917_synthesis.md) sections 1-11, [pillai clocks](collatz_mod6_20260921_pillai_convergents_cycle_gates.md), [G on the negatives](collatz_mod6_20260921_g_negatives_joint_carry.md), [braids collatz](arithmetic_braids_20260917_collatz.md) (cylinder (B5), Lyapunov obstruction sec.5, gate (B7)), [braids2 squarefree](arithmetic_braids2_20260917_squarefree_symmetry.md) (SF3-SF5), [guards discrepancy 2a](collatz_guards_20260921_discrepancy.md), [glued_xor blueprint](glued_xor_20260921_blueprint.md) (B9-B10), [blueprint synthesis](collatz_blueprint_20260921_synthesis.md), and the Lean package [CollatzBlueprintAudit](../../04-computation/lean/CollatzBlueprintAudit/README.md) (which proves only the strong-induction reduction "reach 1 iff some positive-time iterate is below `n`", not any descent).

## 1. The sheet criterion (S1, PROVED)

**S1 (PROVED).** For every odd integer `n`, `T_+(-n) = -T_-(n)` where `T_b(n) = (3n+b)/2^{v_2(3n+b)}`: `3(-n)+1 = -(3n-1)` and `v_2` is sign-blind. Verified on all odd `n` in `[-20001, 20001]` (0 violations). Consequently the positive `3n-1` sheet is literally the negative `3n+1` sheet, and a condition proved for all odd integers under `3n+1` (both signs) is automatically sheet-blind; a condition is sign-specific exactly when its proof uses positivity of `n` (equivalently of the carry `bB`) at `b = +1`.

## 2. The three minus cycles through the gate, and the sign law (S2, FINITE-EXACT + PROVED)

| `b` | cycle | word | `K/L` | `B` | `Delta = 2^K - 3^L` | `q` | `n_0 = bB/Delta` |
|---|---|---|---|---|---|---|---|
| `+1` | `{1}` | `(2)` | `2/1` | `1` | `+1` | `1` | `1` |
| `-1` | `{1}` | `(1)` | `1/1` | `1` | `-1` | `1` | `1` |
| `-1` | `{5,7}` | `(1,2)` | `3/2` | `5` | `-1` | `1` | `5` |
| `-1` | `{17,25,37,55,41,61,91}` | `(1,1,1,2,1,1,4)` | `11/7` | `2363` | `-139` | `1` | `17` |

**S2 (PROVED, inherited).** The gate `q = 1` is SHEET-BLIND. The sign law is SIGN-SPECIFIC in its consequence: a positive `3n+1` cycle needs `Delta > 0` (`K/L > log_2 3`), a positive `3n-1` cycle needs `Delta < 0`. The three minus cycles have `Delta = -1, -1, -139`; the plus cycle has `Delta = +1`.

## 3. Convergents, intermediates, and where the known clocks sit (S3, FINITE-EXACT + CITED)

`log_2 3 = [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, 15, 1, 9, 2, 5, 7, 1, 1, 4, 8, ...]` (mpmath, 400 digits; determinant `p_{n-1}q_n - p_n q_{n-1} = +-1` checked for all 60 computed convergents). `Delta` at `p_n/q_n` alternates in sign with `n`: `-1, +1, -1, +13, -7153, +(18 digits), -(23 digits), +(144 digits), -(313 digits), +(7439 digits), -(15200 digits), +(37847 digits)` for `n = 0..11` (exact big integers), then by 400-digit logarithm through `n = 24`. Even index means `Delta < 0`, odd index `Delta > 0` (Hardy-Wright Thm 154 alternation, CITED classical).

`11/7 = (p_2 + p_3)/(q_2 + q_3) = (3+8)/(2+5)` with `a_4 = 2`: the unique intermediate fraction between `3/2` and `19/12`, not a convergent. Tier caps for the four known clocks:

| clock | `|Delta|` | tier-A cap `3^L/(4L)` | tier-B cap `3^L/(2L)` |
|---|---|---|---|
| `1/1` | `1` | `0.750` (OUT) | `1.500` (in) |
| `2/1` | `1` | `0.750` (OUT) | `1.500` (in) |
| `3/2` | `1` | `1.125` (in) | `2.250` (in) |
| `11/7` | `139` | `78.107` (OUT) | `156.214` (in) |

So the seven-cycle is exactly a tier-B, extreme-intermediate clock (pillai S3 tier B, PROVED), and the question in the task ("is 11/7 an intermediate fraction, not a convergent?") has the answer: yes, and consistently so, because the seven-cycle fails the tier-A hypothesis.

## 4. The Eliahou identity and the Legendre thresholds on both sheets (S4, PROVED)

**S4 (PROVED; SHEET-BLIND identity, SIGN-SPECIFIC hypothesis).** `2^K/3^L = prod_i (1 + b/(3n_i))` holds exactly on all four cycles (checked with rationals), with `|Delta|/3^L` = `0.3333, 0.3333, 0.1111, 0.0636` against the bounds `(1+1/(3N))^L - 1` (plus) or `1 - (1-1/(3N))^L` (minus) = `0.3333, 0.3333, 0.1289, 0.1294`. The smallest cycle minimum `N_min(L)` that forces tier A (hence a convergent clock):

| `L` | plus `N_min` | minus `N_min` | `4L^2/3` |
|---|---|---|---|
| 2 | 6 | 6 | 5.33 |
| 5 | 34 | 33 | 33.33 |
| 7 | 67 | 65 | 65.33 |
| 12 | 194 | 191 | 192.00 |
| 41 | 2248 | 2235 | 2241.33 |
| 53 | 3754 | 3737 | 3745.33 |

The seven-cycle has `N = 17 < 65`, so "the clock is a convergent" does not apply to it; `{5,7}` has `N = 5 < 6`, and `3/2` happens to be a convergent anyway. Verdict: the convergent theorem is SHEET-BLIND as a mechanism; what is SIGN-SPECIFIC is the numerical input `N >= 2^68` (Barina, plus sheet only). On the minus sheet the only available floor is this repo's census (next section).

## 5. Minus-sheet census to `10^7` (S5, FINITE-EXACT)

All `5000000` odd `n <= 10000000` under `3n-1` reach one of the three cycles; no fourth cycle and no escape below `10^7`. Basins: `{1}`: `1636054` (`0.327211`); `{5,7}`: `1623149` (`0.324630`); seven-cycle: `1740797` (`0.348159`). Longest segment before dropping below the start or hitting a cycle element: `172` accelerated steps at `n = 5960769`. Consequence (FINITE-EXACT): any fourth positive `3n-1` cycle has minimum `N >= 10000001`. Control: every odd `n <= 1000000` reaches `1` under `3n+1` (Barina 2020, CITED via pillai S4: all `n < 2^68`; Hercher's abstract, CITED below, uses `3 * 2^69`).

## 6. Cycle-length bounds from cited inputs alone (S6, PROVED given CITED inputs)

Chain, both sheets. (i) `|Delta|/3^L <= e^{L/(3N)} - 1` (plus) or `<= L/(3N)` (minus, Bernoulli). (ii) If this is `<= 1/(4L)`, tier A gives `|log_2 3 - K/L| < 1/(2L^2)`, so reduced `K/L = p_n/q_n` is a convergent (Legendre, CITED Hardy-Wright 184 via pillai) and `L = m q_n`. (iii) `|K - L log_2 3| = m|q_n log_2 3 - p_n| > m/(q_n + q_{n+1})` (Hardy-Wright Thm 171, CITED). (iv) Plus sheet, `Delta > 0`: `y <= e^y - 1 = Delta/3^L <= (L/(3N)) e^{L/(3N)}` with `y = (K - L log_2 3) ln 2`, hence `q_n(q_n + q_{n+1}) >= 3N ln 2 * e^{-L_A/(3N)}`. Minus sheet, `L >= 2`: `1 - e^{-y} <= 1/8` gives `y <= -ln(7/8)` and `1 - e^{-y} >= c_0 y` with `c_0 = 0.936109`, hence `q_n(q_n + q_{n+1}) >= c_0 * 3N ln 2`. `L = 1` is trivial on both sheets (`n_0 = b/(2^k - 3)` integral only for `k = 1, 2`).

| sheet | `N` | `L_A` (tier A holds for all `L <= L_A`) | spacing threshold | first admissible convergent (correct sign side) | conclusion |
|---|---|---|---|---|---|
| plus | `2^68` | `14878203146` | `6.137428e+20` (factor `1.000000`) | `n = 23`, `217976794617/137528045312` | no nontrivial positive `3n+1` cycle has `L <= 14878203146` |
| minus | `10^7 + 1` | `2738` | `1.946585e+07` (factor `0.936109`) | `n = 10`, `50508/31867` | no fourth positive `3n-1` cycle has `L <= 2738` |

On the plus sheet the clock `n = 21` (`q_21 = 6586818670 <= L_A`) has spacing `4.746e20 < 6.137e20` and is excluded; the first three surviving clocks are `n = 23, 25, 27` with `L = m * 137528045312, m * 5409303924479, m * 11571718688839`. On the minus sheet the survivors start at `n = 10, 12, 14` (`L = m * 31867, m * 111202, m * 10590737`), i.e. a fourth `3n-1` cycle with minimum `> 10^7` and `L <= 2738` is impossible, and with minimum `> 10^7` it must sit on an even-index convergent with `q_n >= 15601` (the first even-index one passing the threshold is `31867`; `665 * (665 + 15601)` fails it).

Literature (plus sheet): Simons-de Weger 2005 (CITED, Acta Arith. 117) no `m`-cycles `m <= 68`; **Hercher, "There are no Collatz m-cycles with m <= 91", J. Integer Seq. 26 (2023), art. 23.3.5 (CITED; abstract read 2026-09-22)**: `m >= 92`, verification range `3 * 2^69`, and at least `1.375e11` odd members in any nontrivial cycle. Hercher's bound is stronger than the `1.4878e10` derived here; the point of S6 is that the weaker bound follows from cited inputs by a chain every step of which is checked in the `.out`, and that the same chain on the minus sheet is limited only by the census depth.

## 7. Gate censuses on both sheets (S7-S8, FINITE-EXACT)

**S7.** All `250952` composition words with `L <= 10`, `L <= K <= 2L`, through the gate at `b = +-1`: positive cycles at `b = -1` are exactly `{1}, {5,7}, {17,25,37,41,55,61,91}` (as sets; `10, 10, 7` words including rotations and repetitions); at `b = +1` only `{1}`; the negative cycles at each sign are the negatives of the positive cycles at the other sign (S1). SHEET-BLIND mechanism.

**S8 (Steiner and `m`-cycles on the minus sheet).** 1-cycle words `(1^{L-1}, k)`, `L <= 3000`, `K` within `floor(L log_2 3) + {-1,0,1,2}`, positive and primitive: plus sheet only `(L,k) = (1,2)` giving `{1}`; minus sheet `(1,1)` giving `{1}` and `(2,2)` giving `{5,7}`. So `{5,7}` is a nontrivial 1-cycle: Steiner 1977 (CITED, no nontrivial `3n+1` 1-cycle) is SIGN-SPECIFIC and its transfer to `3n-1` is REFUTED by `{5,7}`. 2-cycle words `1^a k_1 1^b k_2` (`k_1, k_2 >= 2`), `L <= 24`, `9699` words: plus sheet reaches only `{1}` (non-primitive), minus sheet reaches `{5,7}` (repeated `(1,2)(1,2)`) and the seven-cycle `(1,1,1,2,1,1,4)`, a genuine `m = 2` cycle. Simons-de Weger's and Hercher's "no `m`-cycle, `m <= 68` / `91`" are SIGN-SPECIFIC; transfer REFUTED by the seven-cycle. Their Baker/continued-fraction machinery is sheet-blind; their conclusions are not, because the input verification (`2^68`, `3 * 2^69`) is plus-sheet only and the minus sheet genuinely has small cycles.

## 8. Unit gaps (Catalan) and the necklace structure (S9, FINITE-EXACT + CITED)

Clocks with `|Delta| = 1` for `L <= 400`: `(K,L,Delta) = (1,1,-1), (2,1,+1), (3,2,-1)`; clocks with `|Delta| <= 100`: `9` of them, `(1,1,-1), (2,1,1), (3,2,-1), (4,2,7), (4,3,-11), (5,3,5), (6,4,-17), (7,4,47), (8,5,13)` (pillai lists 26 in the larger box; the 9 here are the `K in {floor(L log_2 3), +1}` ones). Mihailescu 2004 (CITED; bibliographic details UNCITED-RECOLLECTION) makes `3/2` the last unit gap for all `L`. On a unit-gap clock `q = 1` for every word, so every composition is a cycle word: this is the entire "all-odd-word / Catalan" structure of the small cycles, and it is SHEET-BLIND (both signs of every unit gap are used: `(1) -> {1}` at `b=-1`, `{-1}` at `b=+1`; `(2) -> {1}` at `b=+1`; `(1,2),(2,1) -> {5,7}` / `{-5,-7}`). Necklace counts (rotation classes of compositions; `q` is rotation-invariant, PROVED in signed_cycles sec.3): `1/1: 1`, `2/1: 1`, `3/2: 1`, `11/7: 30` of `210` compositions (the seven-cycle is one necklace; its 7 rotations are the `7` `q=1` words on `11/7`, checked), `19/12: 2652` of `31824`, `65/41: 6113392816333320` of `250649105469666120`.

## 9. Divergent-orbit conditions on the minus sheet (S10, PROVED + FINITE-EXACT)

Identity (glued_xor, PROVED): `n_L q_L = n_0 + (b/3) sum_{i<L} q_i`, `q_i = 2^{K_i}/3^i`, verified exactly over three periods of each cycle.

| `b` | cycle | `D` per period `= K - L log_2 3` | `sum q_i` one period | ratio `2^K/3^L` | total `sum q` |
|---|---|---|---|---|---|
| `+1` | `{1}` | `+0.415037` | - | `4/3` | diverges |
| `-1` | `{1}` | `-0.584963` | `1` | `2/3` | `3 = 3n_0` |
| `-1` | `{5,7}` | `-0.169925` | `5/3` | `8/9` | `15 = 3n_0` |
| `-1` | seven-cycle | `-0.094738` | `2363/729` | `2048/2187` | `51 = 3n_0` |

Exact answer to the task's question: along the repeated seven-cycle `D_{7m} = -0.094738 m -> -infinity` linearly (`D_L` is NOT bounded; the periodic orbit is not a bounded-`D` object because `2^11 < 3^7`), and `q_{7m}` = `0.9364, 0.8769, 0.8212, 0.7690, 0.7201, 0.6744, 0.6315, 0.5914, 0.5538, 0.5186` for `m = 1..10`. What glued_xor's dichotomy says about periodic orbits on the minus sheet: `D_L -> -infinity` and `sum q_L <= 3n_0` for EVERY positive `3n-1` orbit (positivity alone, (B9)), with equality `sum q = 3n_0` iff eventually periodic ((B10)). So on the minus sheet "`D_L -> -infinity`" separates nothing; on the plus sheet it separates cycles (`D_L -> +infinity`) from escape (`D_L -> -infinity`, glued_xor with Garcia-Tal CITED). Verdict: SHEET-BLIND as a necessary condition on a divergent orbit, SIGN-SPECIFIC as a discriminator.

No bounded strip: on the plus sheet PROVED by source density ([guards 2a](collatz_guards_20260921_discrepancy.md)); on the minus sheet a strip `A <= q_j <= B` gives `n_j <= n_0/A`, a bounded hence eventually periodic orbit whose `D_L -> -infinity` (sign law `2^K < 3^L`) contradicts `q_j >= A`. SHEET-BLIND conclusion with different proofs (the same split as [g_negatives S9/S10](collatz_mod6_20260921_g_negatives_joint_carry.md) for `G`).

## 10. Terras/Everett stopping densities, both conventions (S11, FINITE-EXACT + CITED)

Odd-step convention (`J` accelerated steps, density among odd starts, exact rationals): `d_1 = 1/2, d_2 = 5/8, d_3 = 3/4, d_4 = 51/64, d_5 = 109/128, d_6 = 7/8, d_7 = 911/1024, d_8 = 3729/4096`, then `d_12 = 0.947876`, `d_20 = 0.980746`, `d_30 = 0.992925`, `d_40 = 0.997014`. Terras's own convention (`U(n) = (3n+1)/2` or `n/2` on all `n`, `J` total steps, density mod `2^J`): `F_1 = 1/2, F_2 = 3/4, F_4 = 13/16, F_8 = 237/256`, `F_12 = 0.944824`, `F_20 = 0.973938`, `F_30 = 0.988106`, `F_40 = 0.994177`; `F_12 = 0.9448` and `F_20 = 0.9739` reproduce the three_adic lane's quoted values exactly. Terras 1976 / Everett 1977 (CITED): `F_J -> 1`. SHEET-BLIND: the cylinder of a word has odd-source density `2^{-K_j}` on both sheets (the 2-adic prefix law, `2000` random tests at `L = 6`, `0` violations, including the residue identity `3^L n + bB = 2^{K_L} mod 2^{K_L+1}`), and a word with `2^{K_j} > 3^j` gives descent for ALL `n` at `b = -1` and for all `n > B/(2^{K_j} - 3^j)` at `b = +1`; the minus sheet is slightly easier.

## 11. Squarefree growth prefixes, Lyapunov obstruction, `G`-mirror (S12-S14, PROVED inherited + FINITE-EXACT)

**S12.** Growth prefixes `n_j = 3^j 2^{L+1-j} q - b` (all `k = 1`; on the minus sheet the start is `2^{L+1} q + 1`). `L = 5`, `q <= 100000`: all `6` nodes squarefree for `49167` starts (plus) and `49173` (minus), densities `0.4917, 0.4917`, against the sieve product `delta_L = 0.4918` (`p < 20000`, `nu_L(3) = 1`, `nu_L(p) = min(L+1, ord_{p^2}(3/2))` for `p >= 5`); sieve cross-checked against trial division for `q <= 300` (`148` starts on each sheet). SHEET-BLIND ([braids2 SF5](arithmetic_braids2_20260917_squarefree_symmetry.md), the local counts do not see the sign).

**S13.** `n_0 = 2^{L+1} M q - b`, `M = 6`, `L = 12`, `q = 1`: plus `n_0 = 49151`, minus `n_0 = 49153`, both with word `(1^12)`, `T^L(n_0) = 6377291` resp. `6377293 > n_0`, all nodes in one residue class mod 6 (`5` resp. `1`). SHEET-BLIND ([braids sec.5](arithmetic_braids_20260917_collatz.md)).

**S14.** `G(-m) = -G_-(m)` for all `m <= 3000` coprime to 3 (0 violations): the E-graph/`G` mirror, word law mod `3^{J+1}`, drift `log(2/3)` and tail `(7/9)^{J-1}` transfer verbatim ([g_negatives](collatz_mod6_20260921_g_negatives_joint_carry.md)). SHEET-BLIND.

## 12. The portrait: Table A (hypothetical nontrivial positive `3n+1` cycle)

| # | condition | status | source | 3n-1 check | verdict |
|---|---|---|---|---|---|
| A1 | ordered word `w`, `n_0 = B/Delta`, `q = 1` | PROVED | signed_cycles sec.3; catalan gate | all three cycles have `q = 1` (S2) | SHEET-BLIND |
| A2 | sign law: `Delta > 0`, i.e. `K/L > log_2 3` | PROVED | signed_cycles sec.3 | minus cycles have `Delta = -1, -1, -139 < 0` | SIGN-SPECIFIC (side of `log_2 3`) |
| A3 | `q` rotation-invariant; gate is a test per necklace | PROVED | signed_cycles sec.3; S9 | `11/7`: 30 necklaces, 7 words, one cycle | SHEET-BLIND |
| A4 | `|Delta| = 1` only at `2/1` (`3n+1`) | CITED Mihailescu | S9 | unit gaps `1/1`, `3/2` carry `{1}`, `{5,7}` | SHEET-BLIND (both signs used) |
| A5 | minimum `N >= 2^68` | CITED Barina 2020 (via pillai S4); Hercher abstract `3 * 2^69` | S5 control | minus floor is `N >= 10000001` from this census only | SIGN-SPECIFIC (input) |
| A6 | `K/L` a convergent for `L <= 14878203146` | PROVED from A5 + Legendre (pillai S3-S4) | S4, S6 | seven-cycle `11/7` is an intermediate; `N = 17 < 65` so the hypothesis fails; minus version needs `N >= 10^7`, `L <= 2738` | SHEET-BLIND mechanism, SIGN-SPECIFIC hypothesis |
| A7 | `K/L` an odd-index convergent | PROVED (A2 + alternation) | S3 | minus cycles on even index `0, 2` and an intermediate | SIGN-SPECIFIC refinement |
| A8 | `q_n(q_n + q_{n+1}) >= 3N ln 2 e^{-L_A/(3N)}`; hence `L > 14878203146` | PROVED from A5 + A6 + HW 171 | S6 | minus: `L > 2738`, survivors `31867, 111202, 10590737` | SHEET-BLIND mechanism |
| A9 | no 1-cycle | CITED Steiner 1977 | S8 | `{5,7}` is a 1-cycle | SIGN-SPECIFIC, transfer REFUTED |
| A10 | no `m`-cycle, `m <= 68` (`91`) | CITED Simons-de Weger 2005 (Hercher 2023) | S8 | seven-cycle is a 2-cycle | SIGN-SPECIFIC, transfer REFUTED |
| A11 | at least `1.375e11` odd members | CITED Hercher 2023 abstract | S6 | no minus analogue published | SIGN-SPECIFIC (input) |
| A12 | `D_L -> +infinity`, `sum q_L` diverges along the cycle | PROVED glued_xor | S10 | minus cycles: `D_L -> -infinity`, `sum q = 3n_0` | SIGN-SPECIFIC |
| A13 | Eliahou product identity | PROVED | pillai S4 | holds on all four cycles | SHEET-BLIND |
| A14 | word cylinder mod `2^{K+1}` (2-adic prefix law) | PROVED braids (B5) | S11 | 0 violations both signs | SHEET-BLIND |
| A15 | `G`-mirror gate `m_0(2^J - 3^L) = B'` | PROVED three_adic/g_negatives | S14 | conjugation exact | SHEET-BLIND |

## 13. Table B (hypothetical divergent positive `3n+1` orbit)

| # | condition | status | source | 3n-1 check | verdict |
|---|---|---|---|---|---|
| B1 | `liminf D_L = -infinity` (packing), `D_L -> -infinity` and `sum q_L < infinity` (with Garcia-Tal) | PROVED / CITED input | glued_xor sec.3 | every minus orbit has `D_L -> -infinity`, `sum q_L <= 3n_0`; strict `<` iff infinite | SHEET-BLIND as necessary condition; SIGN-SPECIFIC as discriminator |
| B2 | no bounded strip for `q_j` | PROVED guards 2a | S10 | positivity alone excludes it | SHEET-BLIND conclusion, different proof |
| B3 | start in the density-zero complement of Terras/Everett stopping set | CITED Terras 1976, Everett 1977 | S11 | same cylinder densities; minus is easier | SHEET-BLIND |
| B4 | start in the log-density-zero exceptional set of Tao: for every `f -> infinity`, `Col_min(N) < f(N)` for almost all `N` in logarithmic density | CITED Tao arXiv:1909.03562 Thm 1.3 (quantifiers as stated) | - | statement published for `3n+1` only; whether the proof is sign-blind: UNCITED-RECOLLECTION, left OPEN | SIGN-SPECIFIC as cited |
| B5 | the tree of `1` has `>= x^{0.84}` elements `<= x` | CITED Krasikov-Lagarias 2003 (details UNCITED-RECOLLECTION) | S5 | minus sheet has three trees with finite basin fractions `0.327, 0.325, 0.348` at `10^7`; no exponent computed | SIGN-SPECIFIC as cited |
| B6 | no fixed-modulus periodic weight or bounded correction of `log n` decreases every step | PROVED braids sec.5 | S13 | same construction with `+1` | SHEET-BLIND |
| B7 | all-squarefree growth prefixes of every length exist | PROVED braids2 SF5 | S12 | same counts | SHEET-BLIND |
| B8 | 2-adic prefix law and E-graph/`G` mirror | PROVED | S11, S14 | exact | SHEET-BLIND |
| B9 | finite-congruence Lyapunov obstruction (any `M`) | PROVED braids sec.5 | S13 | exact | SHEET-BLIND |
| B10 | `G`-strip exclusion `3^j/2^{K_j}` | PROVED g_negatives S9/S10 | cited | proved on both sheets there | SHEET-BLIND |

## 14. Contradictions and the minimal unresolved conjunction (S15)

Pairwise contradiction search among A1-A15: NONE FOUND. Witness of joint satisfiability at clock level: the surviving odd-index convergents `n = 23, 25, 27` violate none of A1-A15 a priori; the only untested predicate is A1's integrality `q = 1` on a clock `(m p_n, m q_n)` with all members `>= 2^68` (A5) and at least `1.375e11` members (A11). Among B1-B10: NONE FOUND; a density-zero divergent orbit with `D_L -> -infinity`, `sum q_L < infinity`, unbounded `n_L`, avoiding the Terras and Tao exceptional sets' complements, is consistent with all of them.

Minimal unresolved conjunction (cycle): exists an odd-index convergent `p_n/q_n` with `q_n > 14878203146`, `m >= 1`, and a composition of `m p_n` into `m q_n` parts with `q = 1` whose members are all `>= 2^68`. (Divergence): exists a positive odd orbit with `n_L` unbounded, `D_L -> -infinity`, `sum q_L < infinity`, whose start lies in the log-density-zero exceptional set.

**What a proof must add (consistent with the blueprint audit).** The seven-cycle satisfies every SHEET-BLIND row of both tables, so a proof cannot consist of sheet-blind conditions alone; it must use the sign of the carry at `b = +1` (A2/A7/A12: `Delta > 0`, `D_L -> +infinity` on cycles, `sum q_L` grows), or plus-sheet numerical inputs (A5, A11), or a genuinely new mechanism. The pasted blueprint's objects (adjunction-bipartite operator, `K_5`/`K_{3,3}` minors, `S_2 x S_3` monodromy, `B^3 = -I`, the `8/pi^2` entropy) are sign-blind constructions, hence by this portrait cannot exclude a `3n+1` cycle or escape even if they were Collatz objects, which [collatz_blueprint synthesis](collatz_blueprint_20260921_synthesis.md) already refutes; the Paley/Fano facts in the session lead's probe (21 arcs, 14 cyclic triples, two Fano planes) are correct as tournament facts and irrelevant to either table (no map found). The entropy of `8/pi^2` is `0.70028` bits, not zero and not `0.704` (session lead probe; not recomputed here).

## Reproduction

```bash
cd /tmp/math-wt-collatz-mod6-b
python3 04-computation/experiments/collatz_mod6_20260922_counterexample_portrait.py \
  > 05-knowledge/results/collatz_mod6_20260922_counterexample_portrait.out
python3 -O 04-computation/experiments/collatz_mod6_20260922_counterexample_portrait.py | \
  diff - 05-knowledge/results/collatz_mod6_20260922_counterexample_portrait.out   # only the two timing lines differ
```

About 64 s, under 200 MB; mpmath at 400 digits, exact `Fraction` arithmetic elsewhere; every `check(...)` raises on failure.

## Stopping boundary / next question

Stopped at: the compiled portrait, the sheet verdicts, the two derived length bounds, and the negative contradiction search. Not done: (i) pushing the minus-sheet census beyond `10^7` (each factor of `100` in `N` adds a factor `10` to the minus `L_A`; matching Hercher's method on the minus sheet would need a published-scale verification that does not exist); (ii) the residue-DP gate test on the first surviving minus clocks `50508/31867` and `176251/111202` (`|Delta|` has about 15200 and 53000 digits; the pillai lane's DP is feasible only for `|Delta| <= 30000`), which is the one place where a fourth `3n-1` cycle could be excluded by finite computation; (iii) whether Tao's argument is sign-blind (read the paper's dependence on the sign of the carry; UNCITED-RECOLLECTION here). Next question: is there ANY sign-specific condition beyond A2/A7/A12 that is not a numerical input, i.e. a proved statement about positive `3n+1` orbits whose proof fails for negative ones for a structural reason other than `Delta > 0`? The portrait says no such condition is currently on file.
