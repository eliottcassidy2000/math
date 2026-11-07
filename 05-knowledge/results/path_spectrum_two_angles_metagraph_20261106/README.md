# The Hamiltonian-path spectrum of tournaments seen from both ends — the exact set of counts at each order, the exact set of orders at each count — and how the two meet the levels of the tournament metagraph

**Session:** klein-2026-11-06 (branch `codex/session-spectrum-20261104`). Exact censuses of every tournament of orders 3–11, resolved by the number of cyclic triangles (the level `k = c3` of the metagraph `G_n`, [`METAGRAPH-ATLAS`](../../../00-navigation/METAGRAPH-ATLAS.md)), by strong connectivity, and (orders ≤ 10) by the number of Hamiltonian cycles. The order-10 census (9,733,056 classes) and the order-11 census (903,753,248 classes) are new; the repository previously had exhaustive spectra through order 9 (THM-4861 (b), THM-4744) and samples at orders 10 and 11 (`h_spectrum_n10.out`, `h_spectrum_n11.out`).

Notation. `h(T)` = number of Hamiltonian paths of a tournament `T` (odd, Rédei). `S(n)` = `{h(T) : |T| = n}`, `Ss(n)` = the same over strong `T`. `N(n) = S(n) − S(n−1)` = the values **first attained** at order `n`. `max H(n)` = `max S(n)` (OEIS A038375). `f(n)` = Busch's minimum over strong `n`-tournaments (3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 625, …; THM-1370, THM-4744). `n(h) = min{n : h ∈ S(n)}` (least order), `m(h) = min{n : max H(n) ≥ h}`, `B(h) = max{n : f(n) ≤ h}` (Busch bound), `W(h) = {n : h ∈ Ss(n)}` (strong window). `Λ(m)` = the largest number of cyclic triangles of an `m`-tournament: `(m³ − m)/24` for odd `m`, `(m³ − 4m)/24` for even `m`. `φ_n(k)` = least `h` at level `k` of order `n` (level floor), `ψ_n(k)` = greatest (level ceiling).

## 0. Status

- **KNOWN, used.** `S(n) ⊆ S(n+1)` (pad with a sink: same paths; THM-1370). `h` is the product over the strong components (THM-1370 (1)). Insertion lemma `h(T − v) ≤ h(T)` (THM-1370 (2)). A strong `m`-tournament has at least `m − 2` cyclic triangles (Moon's theorem that a strong tournament of order `m` has at least `m − l + 1` cycles of each length `l`, `3 ≤ l ≤ m`, at `l = 3`; J. W. Moon, *Topics on Tournaments*, 1968; sharp: Burzio–Demaria count the `2m − 4` strong tournaments with exactly `m − 2`; checked here at every order 3–11: the least level with a strong class is `m − 2`). Busch's floor `f(m)` and the chains attaining it (THM-4744). The spectrum omits 7 and 21 and nothing else is known to be omitted (THM-1370; completeness is the open Chan–Pak conjecture, MISTAKE-559).
- **FINITE-EXACT (this packet).** The spectra `S(n)`, `Ss(n)`, the level-resolved tables and the cycle-resolved tables for `n ≤ 10`; the level-resolved spectrum `S(11)`; the least-order function `n(h)` for every `h ≤ max H(11)`; the strong windows `W(h)` for every odd `h < 375` (complete by Busch, since `f(12) = 375`); the level floors `φ(k)` for `k ≤ 9` as minima over **all** orders (§3.3); the flip and deletion gradients of orders ≤ 9.
- **PROVED (elementary, §3.3, §2.2).** `φ(k) = min over all tournaments with exactly k cyclic triangles` is attained at order `≤ k + 2` (Moon's bound + multiplicativity), so the census of orders `≤ k + 2` determines it. The strong window of `h` ends at or before `B(h)`. The level support of `h` ends at or before `max{k : φ(k) ≤ h}`.
- **Pre-registered and tested at order 11** ([`PREREGISTERED.md`](PREREGISTERED.md), committed 14:22 before the census finished): P1–P6, P9, P11, P12 confirmed, P8 off by one level (the band reaches level 52), P7 and the second half of P10 **refuted** (the first gap at order 11 is 80,453 = 0.846 · 95,095, below HYP-9452's 0.85), P13 untested; §6.
- **OPEN (HYP-9481).** The laws of §2 and §3 for all `n`: `n(h) − m(h) ∈ {0, 1}`; `N(n) ⊆ Ss(n)`; `W(h)` an interval ending at `B(h)` or `B(h) − 1`; the zone law `φ_n(Λ(n−1) + 1) > max H(n−1)`; monotone level floors; top-level separation.
- **Sequences not in OEIS (checked 2026-11-06):** `|S(n)|` = 1, 1, 2, 3, 7, 19, 77, 320, 1520, 7663, 43911; `|N(n)|` = 1, 1, 4, 12, 58, 243, 1200, 6143, 36248; unattained below the maximum 0, 0, 0, 2, 16, 9, 157, 208, 3635 (orders 3–11); the level floor `φ(k)` = 1, 3, 5, 9, 11, 15, 23, 33, 41, 51, 69, 97, 107, 141, 171, 213, 227, 321 (exact to `k = 9`, bounds beyond); the level widths of order 7 (1, 5, 7, 8, 17, 23, 39, 40, 59, 61, 79, 52, 47, 15, 3) — none matched.

## 1. Angle one: the exact set of counts at each order

| `n` | classes | strong | distinct `h` | strong-distinct | `max H(n)` | new values `|N(n)|` | unattained odd `< max` (besides 7, 21) | solid prefix `s(n)` | `s/max` |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 2 | 1 | 2 | 1 | 3 | 1 | 0 | 3 | 1 |
| 4 | 4 | 1 | 3 | 1 | 5 | 1 | 0 | 5 | 1 |
| 5 | 12 | 6 | 7 | 4 | 15 | 4 | 0 | 15 | 1 |
| 6 | 56 | 35 | 19 | 13 | 45 | 12 | 2 | 33 | 0.733 |
| 7 | 456 | 353 | 77 | 67 | 189 | 58 | 16 | 61 | 0.323 |
| 8 | 6,880 | 6,008 | 320 | 297 | 661 | 243 | 9 | 609 | 0.921 |
| 9 | 191,536 | 178,133 | 1,520 | 1,482 | 3,357 | 1,200 | 157 | 2,881 | 0.858 |
| 10 | 9,733,056 | 9,355,949 | 7,663 | 7,590 | 15,745 | 6,143 | 208 | 14,683 | 0.933 |
| 11 | 903,753,248 | 884,464,590 | 43,911 | 43,775 | 95,095 | 36,248 | 3,635 | 80,451 | 0.846 |

(`s(n)` = the largest odd number such that every odd value `≤ s(n)` other than 7 and 21 is attained at order `n`.) The census reproduces A000568 (classes), A051337 (strong classes: 9,355,949 and 884,464,590 at orders 10, 11) and A038375 (maxima; `15,745` at order 10 confirms Kagey's value, `95,095` at order 11 Royle's), and at orders ≤ 9 the recorded spectra class by class (THM-4861 (b), THM-4744's `h_multiplicities.out`). The ratio `s/max` alternates with the parity of the order (0.733, 0.323, 0.921, 0.858, 0.933, 0.846 at orders 6–11): at odd orders the top level is the few regular tournaments (3, 15, 1,223 at orders 7, 9, 11), at even orders the many near-regular ones (85; 13,333).

- **The unattained values are 7, 21 and a top band.** Order 6: 35, 39. Order 7: 63, 107, 119, 149, 161–169, 173, 177–187 (the lowest-order exception to "solid up to a top band": 63 is isolated, 0.33 of the maximum). Order 8: 611, 615, 617, 619, 623, 625, 635, 647, 655. Order 9: 157 values in `[2883, 3355]` ([`data/unattained_9.txt`](data/unattained_9.txt)). Order 10: 208 values in `[14685, 15743]` ([`data/unattained_10.txt`](data/unattained_10.txt)); its largest block is 14735–14875 (71 values). Order 11: 3,635 values in `[80453, 95093]` ([`data/unattained_11.txt`](data/unattained_11.txt)). **HYP-9452 (1) (every odd `h ≤ 0.85 · max H(n)` attained for `n ≥ 8`) holds at order 10** (`s(10) = 14,683 = 0.933 · 15,745`) **and fails at order 11** (`s(11) = 80,451 = 0.846 · 95,095`; the first gap 80,453 is below 80,830): the constant is not 0.85, and the pattern of §3.5 says the ratio is governed by the sparsity of the top level, lower at odd orders.
- **Values attained only by non-strong tournaments** (`S(n) − Ss(n)`): all values below `f(n)` (padding), and a few above it whose strong window ended one order early (§2.2): at order 10 these are 127, 129, 131, 133, 137, 139, 141, 143, 149, 151, 157, 163, 173 (and 1–123 below `f(10) = 125`): 73 values; at order 11, 136 values, the largest 331, every attained value `≥ 333` being strongly attained (so THM-4104's strong interval `[429, 80265]` is confirmed and is in fact `[333, 80451]`).
- **Every value first attained at order `n` is attained there by a strong tournament** (`N(n) ⊆ Ss(n)`, orders 3–11, 43,910 values): a non-strong tournament's count is a product of counts of smaller orders, and no such product was ever new. No new value lies below `f(n)`.

## 2. Angle two: the orders — and the strong orders — of each count

The set of all orders attaining `h` is `[n(h), ∞)` (padding), so this angle is the least-order function `n(h)` together with the strong window `W(h)`, which is finite because `f(n) → ∞`. The per-value tables are [`data/least_order_10.txt`](data/least_order_10.txt) and [`data/least_order_11.txt`](data/least_order_11.txt) (every odd `h ≤ 15,745`, resp. `≤ 95,095`: `n(h)`, `m(h)`, the strong window with class counts, and the level windows at each order).

### 2.1 Least orders: `n(h) = m(h)` or `m(h) + 1`

For every one of the 43,911 values with `n(h) ≤ 11`: `n(h) = m(h)` except for 392 values, and for those `n(h) = m(h) + 1`. The 392 are exactly the unattained-below-the-maximum values of orders 6–10 (2 + 16 + 9 + 157 + 208): **every value missing from the top band of order `n` is attained at order `n + 1`** (orders 6→7 through 10→11; the 3,635 band values of order 11 are the first open cases). Hence `n(h) = m(h) + [h is in the top band of order m(h)]`, which is HYP-9452 (2) with the band exception made exact. The least orders of the Fermat primes (THM-4861 (a)) are the instances 3, 4, 6, 8 and (by §6) 11 of this rule.

### 2.2 Strong windows: intervals ending at the Busch bound or one before it

For odd `h < 375 = f(12)` the window `W(h)` is complete from orders ≤ 11: all 184 such values (every odd number in `[3, 373]` except 7 and 21) have a strong realiser, every window is an interval of length 1–4 (14, 27, 84, 59 values), `n_s(h) = n(h)`, and `e(h) = B(h)` for 139 values, `B(h) − 1` for 45, never less. In detail below 225 (complete from orders ≤ 10): all 109 such values (every odd number in `[3, 223]` except 7 and 21) have a strong realiser; **every window is an interval** `[n(h), e(h)]` of length 1 (14 values), 2 (27), 3 (58) or 4 (10); **`n_s(h) = n(h)`** (the least strong order is the least order, for all 109, and for all 7,663 values of `S(10)`); and the last strong order satisfies `e(h) = B(h)` for 90 values and `e(h) = B(h) − 1` for the 19 values 47, 55, 67 (`B = 8`), 77, 79, 83 (`B = 9`), 127, 129, 131, 133, 137, 139, 141, 143, 149, 151, 157, 163, 173 (`B = 10`). Never `B(h) − 2`. (The upper bound `e(h) ≤ B(h)` is Busch's theorem; the lower bound `e(h) ≥ B(h) − 1` is the empirical law.) The 26 further early-ending windows in `[225, 373]` are the values above `f(11)` that order 11 attains only non-strongly (227–241, 247–253, 257, 259, 263, 269, 271, 275, 277, 283, 287, 293, 299, 311, 317, 331, …).

The strong class counts per order (the coefficients of THM-4744's polynomials `P_d`) agree with the recorded `strong_ears_375.out` for all 125 recorded `d = 125..373` at orders ≤ 11, with the nine predicted multiplicities `N(10, h)` of THM-4744 (45: 96, 63: 6, 75: 70, 81: 106, 105: 59, 125: 43, 135: 146, 189: 102, 223: 82) and `N(11, h)` (157, 8, 121, 170, 85, 67, 278, 148, 117), with its order-11 inventory (3,536 strong classes with `H ≤ 625`, least strong value 225), and with the polynomial laws `N(n, 3) = n − 2`, `N(n, 9) = (n − 4)(n − 1)/2`, `N(n, 15) = (n − 4)(n − 3)`, `N(n, 27)`, `N(n, 33)`, `N(n, 45)`, `N(n, 49)`, `N(n, 63)`, `N(n, 75)` at `n = 10, 11` (receipts `spectrum_tables_10/11`).

## 3. Where the two angles meet the metagraph: the levels

`G_n` is stratified by the level `k = c3(T) ∈ [0, Λ(n)]`. The census gives every level's class count (the level widths, a named gap of the atlas), strong count, value set, floor and ceiling ([`data/least_order_10.txt`](data/least_order_10.txt) tail, receipt `level_tables_10`).

### 3.1 Level widths (classes per level; strong classes; distinct values)

- Order 9: `1, 7, 16, 23, 46, 83, 149, 234, 405, 569, 1042, 1307, 2200, 2993, 4339, 5380, 8169, 8679, 12526, 13801, 16868, 17119, 21413, 17482, 19089, 15181, 11988, 6645, 3310, 457, 15` (levels 0–30); strong from level 7: `32, 84, 180, 309, 567, 1019, 1752, 2558, 3900, 6311, 7435, 11312, 13277, 16698, 17119, 21413, …` (all classes are strong from level 21 on); distinct values `1, 1, 2, 3, 6, 11, 18, 29, 49, 66, 97, 121, 179, 218, 268, 317, 374, 417, 472, 470, 473, 437, 408, 365, 336, 291, 234, 174, 116, 63, 13`.
- Order 10: `1, 8, 22, 36, 71, 140, 262, 440, 793, 1240, 2264, 3220, 5672, 8080, 13416, 18232, 30218, 38020, 61172, 76420, 116350, 137784, 206344, 230868, 333066, 370112, 490744, 520012, 690172, 665864, 837084, 785312, 874245, 757500, 794222, 582472, 507631, 306180, 190310, 63724, 13333` (levels 0–40; the widest level is 32, not the top); strong from level 8 (`64, 192, 453, 844, …`; all classes strong from level 31); distinct values `1, 1, 2, 3, 7, 11, 22, 35, 66, 94, 145, 202, 280, 367, 501, 623, 798, 965, 1194, 1392, 1617, 1798, 2068, 2231, 2378, 2511, 2605, 2579, 2467, 2341, 2276, 2100, 1981, 1853, 1743, 1423, 1185, 931, 756, 472, 314` (most distinct values at level 26).
- Order 11: `1, 9, 29, 54, 107, 224, 432, 793, 1432, 2447, 4461, 7111, 12275, 19275, 32413, 48504, 80867, 115386, 186844, 266621, 411844, 563429, 875223, 1155919, 1717900, 2260002, 3248210, 4122118, 5889316, 7187436, 9897018, 12083630, 16003679, 18643977, 24626425, 27753156, 35000261, 38819276, 47043786, 49897941, 59479028, 59249707, 66650749, 65078847, 67552905, 60951613, 59927541, 48244798, 41473941, 29324760, 20765018, 10704181, 4964716, 1305963, 98427, 1223` (levels 0–55; widest level 44; strong from level 9: `128, 432, 1110, 2239, …`, all classes strong from level 41); distinct values per level rise to 14,427 at level 36 and fall to 427 at level 55 (the regular tournaments, A096368).
- Moon's bound is exact: the least level with a strong class is `n − 2` at every order 3–11.

### 3.2 Level floors and ceilings

`φ_n(k)` (least `h` at level `k`):

| `n` | `φ_n(0), …, φ_n(Λ(n))` |
|---|---|
| 6 | 1, 3, 5, 9, 11, 15, 23, 33, 41 |
| 7 | 1, 3, 5, 9, 11, 15, 23, 33, 41, 51, 69, 99, 121, 147, 171 |
| 8 | 1, 3, 5, 9, 11, 15, 23, 33, 41, 51, 69, 97, 107, 141, 171, 263, 325, 395, 463, 549, 621 |
| 9 | 1, 3, 5, 9, 11, 15, 23, 33, 41, 51, 69, 97, 107, 141, 171, 213, 227, 321, 415, 549, 621, 891, 1105, 1275, 1491, 1743, 1953, 2223, 2621, 2931, 3159 |
| 10 | 1, 3, 5, 9, 11, 15, 23, 33, 41, 51, 69, 97, 107, 141, 171, 213, 227, 321, 411, 445, 467, 657, 847, 1143, 1357, 1605, 1953, 2223, 2621, 2931, 3159, 4923, 5657, 6351, 7209, 8601, 9771, 11083, 12261, 13623, 14877 |
| 11 | …, 321 (k ≤ 17 as at order 10), 411, 445, 467, 657, 831, 895, 933, 963, 1401, 1839, 2343, 2781, 3159, 4509, 5657, 6351, 7209, 8601, 9771, 11083, 12261, 13623, 14877, 20479, 23397, 26315, 30441, 34567, 40573, 44967, 49591, 55739, 60589, 66099, 72207, 79359, 85935, 91395 |

`ψ_n(k)` (greatest `h` at level `k`): order 9: `1, 3, 9, 27, 27, 51, 81, 129, 195, 297, 399, 465, 645, 831, 1017, 1233, 1479, 1635, 1683, 1833, 1995, 2205, 1983, 2151, 2283, 2535, 2595, 2721, 2919, 3109, 3357`; order 10: `1, 3, 9, 27, 45, 51, 99, 153, 257, 387, 585, 783, 1049, 1341, 1777, 2259, 2889, 3297, 4001, 4779, 5557, 5965, 6669, 6673, 7189, 7675, 8369, 8885, 9509, 8627, 9205, 9717, 10773, 10905, 11549, 11963, 12657, 13261, 14049, 14733, 15745`.

- **The level floor is nondecreasing in `k` at every order 5–11** (more cyclic triangles never lower the least path count). The ceiling is not monotone from order 8 on (order 10: `ψ(28) = 9509 > ψ(29) = 8627`; order 11: `ψ(39) = 53361 > ψ(40) = 49193`).
- **The floor falls with the order at fixed `k`** (padding gives `φ_{n+1}(k) ≤ φ_n(k)`), and it keeps falling late: level 20 has floor 621 at orders 8 and 9 and 467 at order 10; order 11 keeps every floor of order 10 for `k ≤ 21` and lowers levels 22–29 (847 → 831, 1143 → 895, 1357 → 933, 1605 → 963, …) and 31 on. So "stable for two orders" is not stabilisation.

### 3.3 The global level floor `φ(k) = min over all orders` — exact to `k = 9`

**Lemma.** A tournament with exactly `k` cyclic triangles and `h = φ(k)` has all its nontrivial strong components of order `≤ k_i + 2` where `k_i` is the component's number of cyclic triangles (Moon's bound), the `k_i` sum to `k` (cyclic triangles lie inside strong components) and `h` is the product of the components' counts. Hence `φ(k) = min(φ_s(k), min over compositions k = Σ k_i of Π φ_s(k_i))` where `φ_s(k)` is the least count of a **strong** tournament with `k` cyclic triangles, attained at order `≤ k + 2`. The census of orders `≤ 11` therefore gives `φ(k)` exactly for `k ≤ 9`:

`φ(0..9) = 1, 3, 5, 9, 11, 15, 23, 33, 41, 51` (every minimiser is strong; the best products `Π φ_s(k_i)` are 9, 15, 25, 33, 45, 69, 99, 123 at `k = 2..9`, never below the strong floor). The census values `φ(10..17) = 69, 97, 107, 141, 171, 213, 227, 321` are upper bounds (unchanged at orders 9, 10 and 11, but see the late fall above).

The minimisers at their least order ([`data/level_floor_witnesses.txt`](data/level_floor_witnesses.txt), gentourng strings and scores) are transitive tournaments with one vertex `v` added and at most one further arc reversed: `v` beats the first `a` vertices of the transitive order and loses to the last `b` (`k = ab` cyclic triangles): `(a, b)` = (1,1), (2,1), (2,2), (3,2), (3,3), (4,3), (4,4) give `φ` = 3, 5, 11, 23, 51, 107, 227 at `k` = 1, 2, 4, 6, 9, 12, 16; the others use one more reversal. The regular 7-tournament with 171 paths is the floor of level 14 (`φ(14) = 171` is attained at order 7 only, by the regular tournament whose count is the least of 171, 175, 189).

### 3.4 The zone law: inherited values stop exactly at the top level of the previous order

Split the values of order `n` into **inherited** (`h ∈ S(n−1)`) and **new** (`h ∈ N(n)`). Counting, per level, the inherited and new values (receipt `level_tables_10` (3)):

| `n` | `Λ(n−1)` | highest level with an inherited value | lowest level with a new value | `φ_n(Λ(n−1) + 1)` | `max H(n−1)` |
|---|---|---|---|---|---|
| 5 | 2 | 2 | 3 | 9 | 5 |
| 6 | 5 | 5 | 4 | 23 | 15 |
| 7 | 8 | 8 | 6 | 51 | 45 |
| 8 | 14 | 14 | 8 | 263 | 189 |
| 9 | 20 | 20 | 13 | 891 | 661 |
| 10 | 30 | 30 | 16 | 4,923 | 3,357 |
| 11 | 40 | 40 | 22 | 20,479 | 15,745 |

So at every order the inherited values occupy exactly the levels `0..Λ(n−1)` (level `Λ(n−1)` always carries one: the padded extremal `(n−1)`-tournament), and the levels above carry only new values. The content of the law is the inequality **`φ_n(Λ(n−1) + 1) > max H(n−1)`: every `n`-tournament with more cyclic triangles than any `(n−1)`-tournament has more Hamiltonian paths than any `(n−1)`-tournament** (margins 1.80, 1.53, 1.13, 1.39, 1.35, 1.47, 1.30 at orders 5–11). The lower boundary is the trivial one: a level carries a new value only if its ceiling reaches the first gap `u(n−1)` of order `n−1` (`ψ_10(16) = 2889 ≥ u(9) = 2883`, and 2889 is the third value of that top band; `ψ_11(21) = 13353 < u(10) = 14685 ≤ ψ_11(22) = 15819`).

The metagraph is thus cut into three zones by `k`: `k ≤ λ_new(n) − 1` purely inherited (levels 0–15 at order 10, 0–21 at order 11), a mixed middle, and `k > Λ(n−1)` purely new (levels 31–40 at order 10, 41–55 at order 11).

### 3.5 The top band's anatomy and top-level separation

| `n` | `ψ_n(Λ−1)` | `φ_n(Λ)` | top-level classes / distinct values | band values in the inter-level gap `(ψ(Λ−1), φ(Λ))` / inside the top level's range / below `ψ(Λ−1)` |
|---|---|---|---|---|
| 6 | 37 | 41 | 5 / 3 | 1 / 0 / 1 |
| 7 | 159 | 171 | 3 / 3 | 5 / 7 / 4 |
| 8 | 613 | 621 | 85 / 16 | 3 / 5 / 1 |
| 9 | 3,109 | 3,159 | 15 / 13 | 24 / 87 / 46 |
| 10 | 14,733 | 14,877 | 13,333 / 314 | 71 / 121 / 16 |
| 11 | 89,357 | 91,395 | 1,223 / 427 | 1,018 / 1,424 / 1,193 |

- **Top-level separation (orders 6–11): `φ_n(Λ(n)) > ψ_n(Λ(n) − 1)`** — every tournament with the maximum number of cyclic triangles (regular at odd `n`, near-regular at even `n`) has more Hamiltonian paths than every tournament with fewer. Equivalently the values of the top level occur at no other level (all 3, 3, 16, 13, 314, 427 of them). The margin is 1.6 % at order 9 (3109 → 3159) and 2.3 % at order 11 (89,357 → 91,395), where the 1,223 regular tournaments span `[91395, 95095]` — the pre-registered test P12.
- The band is made of the inter-level gap (its largest block: 14735–14875 at order 10; 1,018 values at order 11), of the sparsity of the top level (15 regular 9-tournaments carry 13 values over 100 odd numbers; 13,333 near-regular 10-tournaments carry 314 values over 435; 1,223 regular 11-tournaments carry 427 values over 1,851), and of values below the second level's ceiling (16 at order 10, all in `[14685, 14731]`, where only level 39 reaches; 1,193 at order 11, from 80,453, where levels 52 and 53 overlap thinly). The pigeonhole frontier (least `x` with fewer classes of count `≥ x` than odd numbers in `[x, max]`) is 143, 3,041, 15,551 at orders 7, 9, 10 — above the band's onset, so the band is a level phenomenon, not a counting necessity.
- The classes above the first gap lie at the top two levels (orders 8, 10), three (order 9) or four (order 11: 317 classes at level 52, then 1,305,828, 98,427, 1,223).

### 3.6 Level windows of a value are not intervals; the `(n, k)`-support of a value is a staircase

`K(n, h) = {k : a class of order n at level k has h paths}` is **not** always an interval: 34 exceptions at order 8, 130 at order 9, 896 at order 10, 3,866 at order 11 (e.g. 39 at levels {5, 7} and 41 at {6, 8} from order 7 on; 51 at {5, 6, 7, 9}). Across orders the support only grows (padding keeps `(h, k)`): the least level of `h` falls with `n` and its greatest level is fixed and bounded by `max{k : φ(k) ≤ h}` (attained for 45, 47, 63, 93, not for 35). Section (5) of receipt `level_tables_10` lists the staircase of every `h < 160`.

## 4. The gradients along the metagraph's two kinds of edge

**Flips** (within `G_n`; `hflip.c`, receipt `flip_stats`): all `(class, arc)` pairs of orders 6–9 (840; 9,576; 192,640; 6,895,296).

- **Sign law.** A flip that raises the level (`Δc3 > 0`) raises `h` in 98.7 %, 97.6 %, 98.0 %, 98.2 % of cases (orders 6–9), and symmetrically for lowering; the level change ranges over `[−(n−2), n−2]` with a symmetric distribution. This quantifies the atlas's "H-gradient mostly increasing".
- **Magnitude.** `max |Δh| / h` at level 0 is exactly `2^{n−2}` (the extreme-arc reversal of the transitive tournament, THM-4749 (4)) and decreases with the level to 0.44, 0.16, 0.27, 0.09 at the top level (orders 6–9). The fraction of `h`-preserving flips falls from 12.6 % (order 6) to 1.2 % (order 9) and, within an order, from level 0 (`1/(n−2)`-like: 0.333, 0.286, 0.250, 0.222) to about 1 % at the top. (Every attained value is joined to every other by a chain of flips — trivially, since `G_n` is connected.)

**Deletions** (between `G_n` and `G_{n−1}`; `hdel.c`, data `dels_<n>.txt`): for every class, the least and greatest `h(T − v)`.

- `1 ≤ h(T) / max_v h(T − v) ≤ 4.09, 4.27, 4.38, 5.35` (orders 6–9; the trivial bound is `n`, from `h(T) = Σ_v (paths ending at v) ≤ Σ_v h(T − v)`); the order-9 maximum is a strong class at level 27 with `h = 2721` and all nine deletions at 509.
- `h(T) / min_v h(T − v)` is unbounded (23, 51, 107, 227 = `h` itself when some deletion is transitive; deleting a vertex from a Hamiltonian path need not leave a path, so there is no upper insertion bound).
- **A value new at order `n` needs growth ≥ 1.64 over every deletion:** the least `h(T) / max_v h(T − v)` over classes whose value is new is 1.727 (`h = 19`), 1.690 (49), 1.646 (107), 1.681 (711) at orders 6–9.

## 5. Hamiltonian cycles alongside

Orders ≤ 10 record `(h, hc)` jointly (`C` lines of `census_<n>.txt`). `max hc` = 3, 5, 24, 54, 249, 971 (orders 5–10). The path range of each small cycle count at order 10: `hc = 1`: `[125, 467]` (91 values), `hc = 2`: `[205, 503]` (113), `hc = 3`: `[225, 681]`, `hc = 4`: `[295, 751]`, `hc = 5`: `[333, 895]`, `hc = 6`: `[345, 849]`; the floors 125, 205 are Busch's `f(10)` and the two-cycle floor `f₂(10)` of HYP-9472, as they must be. The `hc ≤ 2` values of each order (receipt `spectrum_tables_10` (6)) are the finite slices of HYP-9472's family.

## 6. Order 11: the pre-registered column

The census ran from 14:07 to 17:11 (11,011 s wall on eight parts, the machine at load 30 from other sessions; 47 µs per class at full speed) and read 903,753,248 classes, 884,464,590 of them strong, with 43,911 distinct values (43,775 strong), 36,248 new values, maximum 95,095, and 3,635 unattained odd values below the maximum, from 80,453 ([`data/census_11.txt.gz`](data/census_11.txt.gz), [`data/least_order_11.txt`](data/least_order_11.txt), [`data/unattained_11.txt`](data/unattained_11.txt), receipts `*_11`, run log `receipts/census_11_run.log`). Verdicts on [`PREREGISTERED.md`](PREREGISTERED.md):

| # | verdict |
|---|---|
| P1 | **confirmed**: `max H(11) = 95,095`, 903,753,248 classes (and 884,464,590 strong = A051337) |
| P2 | **confirmed**: all 208 band values of order 10 are attained at order 11; `n(h) − m(h) ∈ {0, 1}` for every odd `h ≤ 15,745` (and now for all `h ≤ 95,095` with `n(h) ≤ 11`: 392 exceptions, all with difference 1) |
| P3 | **confirmed**: all 36,248 new values are strongly attained |
| P4 | **confirmed**: `φ_11(41) = 20,479 > 15,745` (margin 1.30); the inherited values occupy exactly levels 0–40 |
| P5 | **confirmed**: `φ_11(k)` nondecreasing; `φ_11(k) = φ_10(k)` for all `k ≤ 21` (the weaker guess `k ≤ 17` included); the floors of levels 22–29 and 31–40 fell |
| P6 | **confirmed**: `N(11, h)` = 157, 8, 121, 170, 85, 67, 278, 148, 117 as predicted; 3,536 strong classes with `H ≤ 625`; least strong value 225; every `z^11` coefficient of the 125 recorded `P_d` matches |
| P7 | **refuted** in its HYP-9452 half: the first gap is 80,453 < 80,830 = 0.85 · 95,095 (solid prefix 80,451 = 0.846 · max); THM-4104's strong interval `[429, 80,265]` is confirmed (every attained value `≥ 333` is strongly attained) |
| P8 | **confirmed in substance, off by one level**: 3,635 > 1,000 unattained values; the classes above the first gap lie at levels 52–55, not 53–55 (level 52 contributes 317 classes) |
| P9 | **confirmed**: the 184 strong windows of odd `h < 375` are intervals (lengths 1–4: 14, 27, 84, 59) with `e(h) = B(h)` (139) or `B(h) − 1` (45); `n_s(h) = n(h)` throughout |
| P10 | distinct values 43,911 ∈ [40,000, 48,000] **confirmed**; `s(11) ≥ 0.85 · max` **refuted** (0.846) |
| P11 | **confirmed**: strong classes from level 9; the level widths sum to 903,753,248 |
| P12 | **confirmed**: `φ_11(55) = 91,395 > ψ_11(54) = 89,357` — every regular 11-tournament has more Hamiltonian paths than every non-regular one (margin 2.3 %); the 427 values of the regular tournaments occur at no lower level |
| P13 | untested (no deletion data at order 11) |

So the laws (a)–(f) of §7 hold one order further, and the one quantitative guess that failed — the 0.85 constant — fails for the reason §3.5 gives: at odd orders the top level is the regular tournaments, few and spread out (1,223 classes over 1,851 odd numbers), so the band begins lower.

## 7. Proof status and next obligations

- **PROVED here:** the lemma of §3.3 (exactness of `φ(k)` from orders `≤ k + 2`); `e(h) ≤ B(h)` (Busch); the level-support bound; the counting identity `h(T) ≤ Σ_v h(T − v)`.
- **FINITE-EXACT:** everything in §§1–6 as stated, from the censuses (`hspec.c` reads `gentourng`'s default format; receipts byte-identical for `python3` / `python3 -O` where Python is involved).
- **OPEN — HYP-9481.** (a) `n(h) − m(h) ∈ {0, 1}` for all odd `h ∉ {7, 21}` (the top band of order `n` fills at order `n + 1`; exact through 10→11); (b) `N(n) ⊆ Ss(n)` for all `n` (exact to 11); (c) `W(h)` is an interval with `e(h) ∈ {B(h) − 1, B(h)}` (exact for `h < 375`); (d) the zone law `φ_n(Λ(n−1) + 1) > max H(n−1)` (exact to 11); (e) `φ_n(k)` nondecreasing in `k` (exact to 11); (f) top-level separation `φ_n(Λ(n)) > ψ_n(Λ(n) − 1)` (exact to 11). (a) with Chan–Pak would make `n(h)` a function of A038375 alone. A proof of (d) or (f) would be a lower bound on Hamiltonian paths in terms of cyclic triangles — the natural next target, since `h` and `c3` are both determined by the same local data (`c3` by the scores, `h` by Rédei's recursion).
- **Next computations.** The deletion gradient at orders 10 and 11 (P13); the exact `φ(10..13)` need the strong census of orders 12–15 restricted to `c3 ≤ 13` (feasible by a generator by level); the order-12 top levels (`gentourng -d5 -D6 12`) for (d)–(f) one order further; the lower part of `S(12)` by ears of the order-11 band for (a).

## 8. Files

- [`scripts/hspec.c`](scripts/hspec.c) — the census (`gentourng -q n | hspec n [hc]`; `hspec n witness h k`); [`scripts/merge_census.py`](scripts/merge_census.py) merges `res/mod` parts; [`scripts/spectrum_tables.py`](scripts/spectrum_tables.py) (§1, §2, §5), [`scripts/level_tables.py`](scripts/level_tables.py) (§3), [`scripts/value_tables.py`](scripts/value_tables.py) (the per-value table and the level widths), [`scripts/hflip.c`](scripts/hflip.c) + [`scripts/flip_stats.py`](scripts/flip_stats.py) (§4 flips), [`scripts/hdel.c`](scripts/hdel.c) (§4 deletions).
- [`data/census_<n>.txt`](data/) (`n = 3..11`): `T h c3 strong count` lines (one per realised triple), `W h witness` / `S h witness` (the first tournament, and the first strong tournament, with `h` paths, in gentourng's upper-triangle format), `C h hc count` (orders ≤ 10). [`data/least_order_10.txt`](data/least_order_10.txt), [`data/least_order_11.txt`](data/least_order_11.txt), [`data/unattained_9.txt`](data/unattained_9.txt), [`data/unattained_10.txt`](data/unattained_10.txt), [`data/level_floor_witnesses.txt`](data/level_floor_witnesses.txt), [`data/dels_<n>.txt`](data/) (orders 6–9), [`data/flips_<n>.txt`](data/) (orders 5–8; order 9's 60 MB per-class file is regenerated by `hflip 9` in two minutes, its summary is in the receipt).
- [`receipts/`](receipts/): `spectrum_tables_10/11`, `level_tables_10/11`, `value_tables_10/11`, `independent_recount_10` (normal and `-O`, byte-identical), `flip_stats`, `deletion_gradient_summary`, `flips_9_summary`, `census_summaries`, and the order-11 run log `census_11_run.log`. [`scripts/independent_recount.py`](scripts/independent_recount.py) recounts 1,628 witnesses (all of orders ≤ 8, 200 per order at 9 and 10 plus the maxima and least new values, and the level-floor minimisers) with an independent pure-Python DP: all agree.
