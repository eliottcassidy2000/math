# Independent audit of `collatz_coalescence_20261004` (Probes A, B, E, F; Proposition E; the typed claims of the note, HYP-9174 and HYP-9175)

**VERDICT: SOUND WITH CORRECTIONS.** Every number was re-derived blind and reproduces (0 of 33 checks
fail; the share dynamic programme to `k = 6000` is now confirmed by exact big-integer counts, not
floats). Proposition E, the three-place fixed-point-group statement and the sheet-cancellation
arithmetic are correct. Six corrections are needed: one substantive (HYP-9174's title clause "of a
modulus q prime to 6" is false -- THM-4522's own `1/10`, `1/14` and the values `7/104`, `1/15`
contradict it -- and P3's ordered list holds only on moduli prime to 6), one wrong constant
(HYP-9175 (iii), low by a factor 10-22), two wrong modulus labels, one mis-identified pair of
periodic points, one typing label.

**Auditor:** independent subagent, 2026-10-04, worktree `collatz-synthesis-20261004`.
**Code (written before opening the session's script):**
`04-computation/experiments/collatz_coalescence_20261004_audit.py`, output
`collatz_coalescence_20261004_audit.out` beside it (20 s, `python` + numpy). Exact integers and
`Fraction`s throughout; the only floating-point numbers are the shares at `k = 12000, 24000` (the
float DP agrees with the exact DP to `1e-12` at `k = 6000`).
**Method for `I`:** `I(m/q) = g_q(m)/q`, where `g_q(m)` is the minimum of `D(x) = min(x, q-x)` over
the forward orbit of `m` under `x -> 2x, 3x (mod q)`; computed for all `q` (including `q` divisible
by 2 or 3, where the orbit contains non-units and `0`) by an `O(q)` reverse propagation in increasing
`D`, cross-checked against a per-`m` forward BFS and against the coset formula on `q <= 300`.
**Read:** the coalescence note, HYP-9174, HYP-9175, THM-4522 (statement), THM-4479, the Spitzer
note's `W_k` table, S23's `collatz_three_mirrors_character_spectrum_20260929.out` (it lives in
`05-knowledge/results/`, not beside the S23 script), S19's level-18 table.

## 1. Check table (my number vs the claimed number)

| # | check | mine | claimed | result |
|---|---|---|---|---|
| 1a | reduced `m/q`, `q <= 600`, `I >= 1/14` | 88 points | 88 | PASS |
| 1a | denominators | `5,7,10,11,13,14,26,28,33,52,56` | same | PASS |
| 1a | values x multiplicity | `1/5 x4, 1/7 x6, 1/10 x4, 1/11 x20 (q=11,33), 1/13 x30 (13,26,52), 1/14 x24 (14,28,56)` | same | PASS |
| 1b | coset formula = per-`m` orbit BFS = table algorithm, `q <= 300` prime to 6 | 99 moduli, 13686 fractions, 0 mismatches | identity (PROVED) | PASS |
| 1c | top 12 below `1/14`, `q <= 4000` prime to 6 | `5/73, 1/17, 1/19, 5/97, 13/259, 7/145, 29/601, 11/247, 31/697, 19/431, 1/23, 1/25` | same | PASS |
| 1c | note section 6's top 20 and quoted indices | values 13-20: `1/29, 13/385, 1/31, 7/241, 1/35, 1/37, 13/485, 23/865`; indices 73:2, 97:2, 259:3, 145:2, 601:8, 247:3, 697:8, 431:10, 23:2, 385:2, 241:2, 485:8, 865:4 | same | PASS |
| 1c | same census over ALL `q <= 4000` | `5/73, 7/104, 1/15, 1/17, 1/19, 5/97, 5/99, 11/219, 13/259, 77/1539, 1/20, 7/145, 25/518, 29/601` | not claimed; see C1 | -- |
| 1d | primes `5 <= q <= 20000` | 2260 | 2260 | PASS |
| 1d | index histogram of `<2,3>` | `{1:1598, 2:465, 3:86, 4:42, 5:16, 6:19, 7:5, 8:13, 9:3, 10:4, 12:3, 14:1, 15:1, 16:1, 17:1, 24:1, 56:1}` (lcm(ord 2, ord 3) formula = explicit subgroup on primes `<= 3000`) | same | PASS |
| 1d | fraction with `<2,3> = (Z/q)^x` | `1598/2260 = 0.7071`; naive product `prod_l (1 - 1/(l^2(l-1))) = 0.69750`; difference `0.0096 = 0.99` SE | `0.7071`, `0.6975`, "within one SE" | PASS |
| 1d | primes `q = 1 mod 24`, index exactly 2 | 186; `q I_max(q)` = least QNR in all 186 (QNR histogram `5:96, 7:57, 11:21, 13:8, 17:4`) | 186, no exception | PASS |
| 1d | largest `I_max` over primes with proper `<2,3>` | `5/73(2), 5/97(2), 29/601(8), 19/431(10), 1/23(2), 197/6563(17), 7/241(2), 5/193(2), 11/439(3), 109/4513(12), 29/1201(4), 157/6553(56)` | same values; indices 17, 12, 56 as quoted | PASS |
| 2a | `|Bad_k|`, `k = 4..24`, word DP = direct residue enumeration (numpy, `2^23` residues at `k = 24`) | `3,4,8,13,19,38,64,128,226,367,734,1295,2114,4228,7495,14990,27328,46611,93222,168807,286581` | same; `= W_4..W_20` of the Spitzer table | PASS |
| 2b | `|Bad_k cap {1,3 mod 8}|` | `1,1,2,3,4,8,13,26,45,71,142,247,395,790,1386,2772,5019,8468,16936,30485,51280` | same | PASS |
| 2b | `Bad_k subset {3, 7 mod 8}` | residues mod 8 present: `{3, 7}` for every `k = 4..24` | claimed | PASS |
| 2c | `<3> = {1, 3 mod 8}` in `(Z/2^k)^x`, `-1 notin <3>`, `k = 3..20`; `2^(L/2) = -1 mod 3^n`, `n = 2..20` (Probe F) | holds | claimed | PASS |
| 2d | shares `s_k`, EXACT big integers, `k = 50, 100, 200, 400, 800, 1600, 3200, 6000` | `0.169832, 0.165041, 0.162544, 0.161188, 0.160470, 0.160098, 0.159908, 0.159819` | `0.16983, 0.16504, 0.16254, 0.16119, 0.16047, 0.16010, 0.15991, 0.15982` | PASS |
| 2d | float DP `k = 12000, 24000` | `0.159767, 0.159741` | -- | -- |
| 2d | `1/k` extrapolation | `0.159717` (1600, 6000), `0.159716` (6000, 24000); successive differences `5.2e-5`, `2.6e-5` confirm the `1/k` law (`c = 0.62`) | `0.15972`; "0.1597 to four digits" | PASS |
| 2d | `C_k = W_k k^(3/2) 2^(-hk)`, `h = H_2(log_3 2) = 0.94996` | `4.62 (k=24), 7.66 (100), 9.84 (400), 10.76 (1600), 11.05 (3200), 10.98 (6000)`; `W_6000/W_5999 = 1.9251`, `2^h = 1.9318` | Spitzer window `[9.66, 11.05]` | consistent |
| 2d | HYP-9175 (iii): fraction of resisting exponent classes `|Bad_k cap <3>| / 2^(k-2)` | exact `51280/2^22 = 1.22e-2` at `k = 24`; the HYP's `0.16 x 2 x 2^(-(1-h)k) k^(-3/2)` gives `1.18e-3`; ratio `10.3` at `k = 24`, `21.9` at `k = 6000` | "`~ 0.16 x 2 x ...`" | see C2 |
| 2e | unique `a mod 2^(k-2)` with `3^a = -5 mod 2^k`, `k = 3..14` | `1,3,3,11,11,11,11,11,267,267,1291,3339` | same | PASS |
| 2e | `a mod 2^37` | `68011179275` | same | PASS |
| 2e | the labels of those digits | `267 mod 512` and `mod 1024`, `1291 mod 2048`, `3339 mod 4096` | note: "`1291 mod 1024, 3339 mod 2048`" | see C3 |
| 2e | `{3^a mod 2^k} = {1, 3 mod 8}`, `k <= 20`; cycle points | `-5,-7,-37,-55,-61,1` are `1` or `3 mod 8`; `-1,-17,-25,-41` are `7`, `-91` is `5`; all eleven lie on `T`-cycles | same | PASS |
| 2f | `|det(A^j - I)| = |L_j - 1 - (-1)^j|`, `tr A^j = L_j`, `j <= 60` | holds; `tr - |det| = 1 + (-1)^j` | same | PASS |
| 2g | carry census `A <= 18` (Probe B) | integral carries exactly at the 18 points of the known `T`-cycles (`0`; `1, 2`; `-1`; `-5, -7, -10`; the 11 points of the `-17` cycle); clock `242461 = 2^18 - 3^9` hit only by `x_0 = 1, 2` | same (odd points listed) | PASS |
| 4.5 | table inputs from S23's output, `n = 8, 10, 12, 14, 16` | `sum_prim S = -0.057189, -0.039532, -0.036036, -0.017939, -0.015797`; `parity sum = 8.324134, 18.73390, 42.16941, 94.86095, 213.43393` | `-0.0572, ..., -0.0158`; `8.324, ..., 213.434` | PASS |
| 4.5 | `E_n = ((E+O) + (E-O))/2`, `O_n`, `(0.975/6)(3/2)^n`, ratio | `4.134, 9.347, 21.067, 47.42, 106.71`; `-4.191, ..., -106.72`; `4.165, 9.371, 21.084, 47.44, 106.74`; `1.4e-2, 4.2e-3, 1.7e-3, 3.8e-4, 1.5e-4` | same | PASS |
| 4.5 | `E_18` | `(1/2)[Delta rho_18(1) + Delta rho_18(-1)] = (1/2)[~ -0.01 + 0.9748 (1/3)(3/2)^18] = 240.1` (S19: `rho_18(1) = 0.2736`, `rho_18(-1) = 0.9748 (3/2)^18`) | `~ 240` | PASS |

## 2. Part 3 findings

**(i) Proposition E -- correctly proved.** `Z_2^x = {+-1} x (1 + 4Z_2)`; `1 + 4Z_2` is pro-cyclic and
any element `5 mod 8` (so `-3`) generates it topologically, because `(1 + 4Z_2)^2 = 1 + 8Z_2`. Writing
`3^a = (-1)^a (-3)^a`, the closure of `3^Z` is `{(-1)^a (-3)^a : a in Z_2}`: `a` even gives
`(-3)^a in 1 + 8Z_2`, `a` odd gives `-(5 + 8Z_2) = 3 + 8Z_2`, and both classes are filled. So the
closure is `(1 + 8Z_2) cup (3 + 8Z_2)`, the "diagonal" index-2 subgroup of `Z/2 x Z_2`, and
`a -> 3^a` is a topological isomorphism `Z_2 -> closure` (injective since `3^a = 1` forces `a` even
and `(-3)^a = 1`), which is why `log_3 x` is unique, as the note says. The word of `-5` is `(110)^inf`
with `3^(o_j) > 2^j` for all `j` (`9^t > 8^t`), so the class `log_3(-5) mod 2^(k-2)` lies in
`Bad_k`: proved. The finite instances (`k <= 20`; the digits to `2^37`) reproduce.

**(ii) The three-place solenoid statement -- correct.** `Sigma = (R x Q_2 x Q_3)/Z[1/6]` is the
Pontryagin dual of the discrete group `Z[1/6]`; for a dual endomorphism `u-hat`, `Fix(u-hat)` is the
annihilator of `(u-1)Z[1/6]`, i.e. the dual of `Z[1/6]/(u-1)Z[1/6]`; since `2^A` is a unit of
`Z[1/6]`, `(u-1)Z[1/6] = (3^p - 2^A)Z[1/6]`, and `Z[1/6]/N Z[1/6] = Z/N` exactly because `N` is prime
to 6, so `Fix(u-hat) = dual(Z/N) = Z/|2^A - 3^p|`. (On the 2-adic solenoid `dual(Z[1/2])` the same
`u` is not a unit of `Z[1/2]`, so `u-hat` is an endomorphism, not an automorphism; the fixed-point
group is the same. At three places `u` is a unit, as the note says.) The cat-map side is the
classical `Fix(A^j) = Z^2/(A^j - I)Z^2 = Z[phi]/(phi^j - 1)` of order `|N(phi^j - 1)| =
|L_j - 1 - (-1)^j| = |F(2,j)^ab|`; checked to `j = 60`. One mis-identification: C4.

**(iii) Sheet cancellation -- arithmetic correct, inputs verified.** `E_n = (1/2)[Delta rho_n(1) +
Delta rho_n(-1)]`; with `rho_n(-1) = 0.975 (3/2)^n` at consecutive levels,
`Delta rho_n(-1) = 0.975 (1 - 2/3)(3/2)^n = (0.975/3)(3/2)^n`, hence `E_n = (0.975/6)(3/2)^n + O(0.03)`
and `E_18 = 0.1625 x 1477.9 = 240.1`. The five table rows are exactly the `sum_prim S` and `parity
sum` lines of S23's output (table above); S19 gives `rho_18(1) = 0.2736` and `0.9748` at level 18
(increasing in `n`, which shifts `E_18` by well under 1).

**(iv) PROVED / FINITE-EXACT inventory.** PROVED: Proposition E (proved), `I(m/q) = n_q(m)/q` on
moduli prime to 6 (proved; the orbit of a unit is a coset because 2, 3 have finite order), the
three-place `Fix` statement (proved), the cat-map identities (classical; reproduced), S23's Theorem
1(ii) (cited, numbers verified). FINITE-EXACT: the 88-point census, the value spectrum to 4000 (on
moduli prime to 6), the index census and the 186-prime law to 20000, the `I_max` list, the carry
census to `A <= 18`, `|Bad_k cap <3>|` to `k = 24` -- all reproduced by own code; the DP to 6000 is
float in the session but now exact here (C6). Typing issues: C1 (a PROVED label attached to a
false universal sentence), C2 (a wrong constant inside a `~`), C6, and remarks R1-R3.

## 3. Corrections

**C1 (substantive overclaim; HYP-9174 title, and the note's P3).** Title: "every value of I(t) =
inf ||2^j 3^k t|| is n/q with n the least absolute residue of a <2,3>-coset of a modulus q prime to
6 (I(m/q) = n/q exactly; PROVED)". The PROVED identity is for `q` prime to 6 only. Values on mixed
moduli are not of this form: THM-4522's own `1/10` and `1/14` (a value `n/q'` with `q'` prime to 6
cannot equal `1/10`), and below `1/14` on `q <= 4000` the values `7/104, 1/15, 5/99, 11/219,
77/1539, 1/20, 25/518` -- e.g. `I(1/15) = 1/15` (the orbit of `1 mod 15` contains `1`) and
`I(7/104) = 7/104` lie strictly between `5/73` and `1/17`. The same unqualified sentence is P3
("first 5/73, then 1/17, 1/19, 5/97, 13/259, ...") and HYP-9174's test (i) ("no other mechanism");
section 6 of the note and the HYP body are correctly qualified ("moduli prime to 6"). Fix: title ->
"on moduli q prime to 6 every value is n/q with n the coset minimum (PROVED); on q = 2^e 3^f q' the
value is the minimum over the intermediate denominators and can be new (1/10, 1/14, 1/15, 7/104)";
P3 -> insert "on moduli prime to 6"; test (i) -> "with the mixed-denominator mechanism". The first
value below `1/14` over all rationals is still `5/73`; the all-`q` list begins `5/73, 7/104, 1/15,
1/17, 1/19, 5/97, 5/99, 11/219, 13/259, 77/1539, 1/20, 7/145, 25/518, 29/601`.

**C2 (wrong constant; HYP-9175 status, consequence (iii)).** "the resisting exponent classes at
depth k number |Bad_k cap <3>| = s_k W_k among the 2^(k-2) classes, i.e. a fraction ~ 0.16 x 2 x
2^(-(1-h)k) k^(-3/2) of all exponents". The exact fraction is `4 s_k W_k / 2^k = 4 s_k C_k
2^(-(1-h)k) k^(-3/2)` with `C_k = W_k k^(3/2) 2^(-hk) -> ~ 11` (THM-4495's observed window
`[9.66, 11.05]`; `4.6` at `k = 24`, `10.98` at `k = 6000`), so the prefactor is `~ 4 x 0.16 x 11 = 7`,
not `0.32`: the stated formula is low by a factor `10` at `k = 24` (`1.18e-3` against the exact
`51280/2^22 = 1.22e-2`) and `22` asymptotically. Fix: "a fraction `4 s_k W_k / 2^k ~ 7 x
2^(-(1-h)k) k^(-3/2)` (exactly `51280/2^22 = 1.2 x 10^-2` at `k = 24`)".

**C3 (wrong modulus labels; note section 6, Probe E).** "267 mod 512, 1291 mod 1024, 3339 mod
2048". `1291 > 1024` and `3339 > 2048`. The residues are `267 mod 512` (and `mod 1024`), `1291 mod
2048`, `3339 mod 4096` (`k = 11, 12, 13, 14`). Fix the two moduli (the value list `1, 3, 3, 11, 11,
11, 11, 11, 267, 267, 1291, 3339` in HYP-9175 is right, and `68011179275 mod 2^37` is right).

**C4 (mis-identified points; note 4.1).** "the cat map loses exactly the shift's two points 0^inf
and (01)^inf when j is even". The count `tr(A^j) - |det(A^j - I)| = 1 + (-1)^j` is right, but `0^inf`
is the shift's unique fixed point and codes the origin, which is fixed by every `A^j`; the two points
lost for even `j` are the period-two codings `(01)^inf` and `(10)^inf` (`j = 2`: the shift has
`L_2 = 3` points, the torus `|det(A^2 - I)| = 1`). Fix: "`(01)^inf` and `(10)^inf`".

**C5 (mis-stated fact; note 4.4 table, row THM-4522 + HYP-9174).** "top {1/5, ..., 1/14} where 2 is
a primitive root". `2` is not a primitive root mod 7 (order 3; `3` generates), and the top also lives
on `10, 14, 26, 28, 33, 52, 56`. Fix: "where `<2,3> = (Z/q')^x` for the 6-free part `q' in
{5, 7, 11, 13}`".

**C6 (typing; note status line).** "FINITE-EXACT (... |Bad_k cap <3>| exactly to k = 24 and by a
dynamic programme to k = 6000)". The session's DP is floating point (section 6 says so). The audit's
big-integer DP (`04-computation/experiments/collatz_coalescence_20261004_audit.py`, 5 s) gives the
exact shares at `k = 50..6000`; they round to every digit the note prints. Fix: "exactly to
`k = 6000` (audit's big-integer DP; the session's float DP agrees to `1e-12`)".

**Remarks (no change required).** R1: HYP-9174's 186-prime law is a one-line theorem, not a census:
for a prime `q = 1 mod 24`, `2` and `3` are squares, index 2 forces `<2,3>` = the squares, `-1` is a
square, so the extremal coset is the non-residues and its least absolute residue is the least QNR;
its typing can be PROVED (the run is a software check). R2: "the seven points of the -17 cycle
((11,7)) with its seven rotations" -- all 11 rotations of the word hit the zero class of the clock
139 (7 odd points and the 4 even points `-82, -136, -68, -34`); "seven" counts odd points only, say
so. R3: HYP-9175's "a rational cycle point x" should read "odd rational cycle point" (even points are
not 2-adic units). R4: the harmonic-function formula for `s_inf` is consistent with the measured
value (`h(y_0)/h(y_1) = (s/(1-s)) (r/(1-r)) = 0.3249`) and is correctly typed CONJECTURED; the
derivation needs the factor `e^(-theta*) = (1-r)/r` from the Cramer change of measure, which is what
the `(1-r)` versus `r` weights supply.

## 4. Summary

Everything computable in the note and the two hypothesis files reproduces from independent code:
the 88-point census with its denominators and multiplicities, the coset identity on 13686 fractions,
the twenty largest values below `1/14` on moduli prime to 6 with their indices, the prime index
histogram, the two-generator Artin fraction `1598/2260` (0.99 standard errors above the naive
product `0.69750`), the 186 primes with the least-non-residue law, the `I_max` list to 20000,
`|Bad_k|` and `|Bad_k cap <3>|` to `k = 24` by two methods, the exact shares to `k = 6000`
(`s_6000 = 0.159819`, `1/k` limit `0.15972`), the digits of `log_3(-5)` to `2^37`, the Lucas identity
to `j = 60`, the carry census to `A <= 18`, and the sheet-cancellation table from S23's own output.
Proposition E and the three-place `Fix` statement are proved as stated. What needs changing is
language, not mathematics: HYP-9174's title and P3 state the coset mechanism as if it covered all
rationals (it covers moduli prime to 6; the values `1/10, 1/14, 1/15, 7/104` need the
intermediate-denominator mechanism the HYP body already names), HYP-9175's consequence (iii) carries
a prefactor that is too small by an order of magnitude, two modulus labels and one pair of periodic
points are mislabeled, and the float DP should not be called FINITE-EXACT until (now) an exact
count stands behind it.
