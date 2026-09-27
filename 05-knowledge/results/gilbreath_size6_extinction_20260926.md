# The size-6 extinction law of Gilbreath's automaton: exact Markov reduction, exact rationals, and the dyadic anatomy of the re-emission excess

**Session:** opus, `gilbreath6-collatz-precision-20260926`, 2026-09-26.
**Owner's directive:** "attack the size 6 extinction law next".
**Inherits:** [THM-4511](../../01-canon/theorems/THM-4511-gilbreath-size-four-wall-theorem.md)
(wall theorem: a `4` never crosses the first sea `2` on its left; lone size-4
extinction exactly `2^(1-F)`), the S10 note
([`gilbreath_fermat_platonic_20260926.md`](gilbreath_fermat_platonic_20260926.md):
finite-context table `p_6(F)` exceeding the front-only tail `F 2^(1-F)` by
`6`-`10` per cent), the S9 automaton reading (fronts at speed `1`, the `0/2`
sea is XOR, the Lucas kernels).

**Status: PROVED (the exact Markov reduction, Theorem A; the front-only law,
recalled) + FINITE-EXACT (exact rationals `p_6(F)` for `F <= 11` and `p_8(F)`
for `F <= 11`, floating values to `F = 17`, the exact decomposition of the
size-6 excess by the number of zeros before the wall for `F <= 11` and to
`F = 16` in floating point) + CONJECTURE (HYP-9163: the dyadic law of the
excess). No claim on Gilbreath's conjecture.** Script:
`04-computation/experiments/gilbreath_size6_exact_20260926.py` -> `.out`.

## 0. Result in one paragraph

For a lone defect of size `d` at distance `F` in a uniform `0/2` sea, the
extinction probability is an exact rational computable from a finite
absorbing Markov chain (Theorem A). For `d = 4` it returns `2^(1-F)` exactly
(THM-4511, a consistency check). For `d = 6` the excess over the front-only
tail `F 2^(1-F)` is `2^(-F) c_6(F)` with `c_6(F)` an exact dyadic rational
that decomposes by the number `z` of zeros between the defect and the first
sea `2` on its left:

```text
c_6(F) = sum_z c_z(F),   c_0(F) = (2/3)(1 - 4^-(2^k - 1)),  c_1(F) = c_2(F) = (2/15)(1 - 16^-(2^(k-1) - 1))
         for 2^k < F <= 2^(k+1)  (k >= 1; c_1, c_2 start at F = 5),
         c_4 = 1/2 (F = 7, 8), 261/512 (F >= 9);  c_3 = c_5 = c_6 = 1/128 (F >= 9);  c_7 = 0;
         c_(z+8)(F) = c_z(F - 8) for the computed range  (z = 0, 1, 2, 4;  F <= 16).
```

The terms are truncated geometric series whose lengths jump at powers of
two (the Frobenius tower of the sea), new terms switch on at `z = 2^k`, and
on the computed range the total is `c_6(F) = 1.4 log_2 F + O(1)` (`0.5,
0.91, 1.41, 1.47, 1.97, 2.37, 2.87` at `F = 3, 5, ..., 15`; the one doubling
`9 -> 17` adds `1.41`), so the excess ratio `p_6(F)/(F 2^(1-F)) - 1 = c_6(F)/(2F)`
(`0.083` at `F = 3`, `0.100` at `7`, `0.096` at `15`, `0.085` at `17`) is flat at
`6`-`10` per cent for `F <= 17`; that it eventually returns to zero like
`0.7 log_2 F / F` is HYP-9163's conjecture, not a result (the audit removed a
first draft's "settled"). For `d = 8` the
exact denominators are `2^m` times products of Fermat numbers
(`255 = F_0 F_1 F_2`, `65535 = F_0 F_1 F_2 F_3`), the same dyadic geometric
series in a different dress. A full proof of the size-6 law is not given;
its mechanism (the phases of the defect column and the toggling of the wall
at the Lucas hit times) is identified in section 4.

## 1. Theorem A: the exact Markov reduction (PROVED)

Columns `1..F-1` carry the left sea word, column `F` the defect `d`, and
columns `> F` the right sea. The rule looks only to the right, so the
values at columns `1..F` at row `s + 1` are determined by their values at
row `s` and the value `u_s` of column `F + 1` at row `s`. The right sea is a
pure `0/2` word, so `u_s = 2 XOR_(j subset s) b(F + 1 + j)` (sea kernel
theorem), which contains the fresh bit `b(F + 1 + s)` with coefficient `1`:
conditionally on everything at rows `< s`, `u_s` is uniform on `{0, 2}`. Hence
the row vector `(a(1), ..., a(F))` is a Markov chain with two equally likely
transitions `new[i] = |a[i] - a[i+1]|` (`i < F`), `new[F] = |a[F] - u|`,
`u in {0, 2}`; the states with `a(1) >= 4` are absorbing with extinction
probability `1` (the leading `1` becomes `|1 - a(1)| != 1` at the next row),
the states with all `a(i) <= 2` are absorbing with probability `0`, and
`p_d(F)` is the mean over the `2^(F-1)` initial states of the absorption
probability, a rational number (finite chain). ∎

The chain is small (`30`, `97`, `253`, `669`, `1639`, `3952`, `8795`, `20084`,
`44559`, `101200`, `220506`, `484666`, `1047839`, `2271338`, `4766630` reachable
states for `d = 6`, `F = 3..17`); exact rationals by sparse elimination over
`Q` for `F <= 11`, floating value iteration (to `10^-15`) beyond, the two
agreeing to nine digits where both exist. All values reproduce the
finite-context numbers of S10 (`R = 10`, `R = 9` at `F = 12`, where S10's
`0.006338` is the truncated value and the chain gives `0.0063395`).

## 2. Exact values (FINITE-EXACT)

Size `4`, all `F = 3..12`: `p_4(F) = 2^(1-F)` (floating chain values to `10^-15`,
agreeing with THM-4511's exact law).

Size `6` (`F`, exact `p_6(F)`, excess over `F 2^(1-F)`, `2^F *` excess):

```text
 3  13/16            1/16          1/2
 4  17/32            1/32          1/2
 5  349/1024         29/1024       29/32  = 0.90625
 6  413/2048         29/2048       29/32
 7  493/4096         45/4096       45/32  = 1.40625
 8  557/8192         45/8192       45/32
 9  159469/4194304   12013/2^22    12013/8192 = 1.466431
10  175853/8388608   12013/2^23    1.466431
11  196333/16777216  16109/2^24    16109/8192 = 1.966431
12  (float)          4.800856e-4   1.966431
13  (float)          2.896339e-4   2.372697
14  (float)          1.448169e-4   2.372697
15  (float)          8.766726e-5   2.872697
16  (float)          4.383363e-5   2.872697
17  (float)          2.192064e-5   2.873238
```

Two exact regularities: the excess **halves exactly from every odd `F` to
`F + 1`** (`2^F *` excess is constant on the pairs `(3,4), (5,6), ..., (15,16)`),
and it changes only at odd `F`, by the switching-on of new `z`-terms
(section 3). Ratios `p_6/(F 2^(1-F))`: `1.083, 1.0625, 1.091, 1.076, 1.100,
1.088, 1.081, 1.073, 1.089, 1.082, 1.091, 1.085, 1.096, 1.090, 1.085`
(`F = 3..17`).

Size `8` (exact, `F = 3..11`): `1`, `29/32`, `703/960`, `16243/30720`,
`388363/1044480`, `509743/2088960`, `2576519/16711680`, `202620131/2139095040`,
`64449019483/1099494850560`; excess ratios `1.036`-`1.082` (`F = 4..17`),
denominators `2^m * 15`, `2^m * 255`, `2^m * 65535` (products of Fermat
numbers): the denominators suggest that the size-8 law is a sum of geometric
series with ratios `2^-(2^i)` (not shown).

## 3. The dyadic anatomy of the size-6 excess (FINITE-EXACT, exact rationals for `F <= 11`)

Group the left words by `z` = number of zeros between the defect and the
first `2` on its left (the wall, at column `c = F - 1 - z`). Exact
`2^F *` excess by `z`:

```text
 F :   z=0        z=1        z=2       z=3     z=4        z=5     z=6     z=7   z=8     z=9   z=10  z=12
 3,4   1/2        0          0
 5-8   21/32      1/8        1/8       0       1/2 (F>=7) 0       0       0
 9-11  5461/8192  273/2048   273/2048  1/128   261/512    1/128   1/128   0     1/2 (F>=11)
 12    (float)    same                                                                 0.5
 13,14 same       same       same      same    same       same    same    0     21/32   1/8   1/8   0
 15,16 same       ...                                                                                 1/2
```

Read as series:

* `z = 0`: `1/2 = (2/3)(1 - 4^-1)`, `21/32 = (2/3)(1 - 4^-3)`,
  `5461/8192 = (2/3)(1 - 4^-7)`: the partial sums of `(1/2) sum_i 4^-i` with
  `1, 3, 7` terms for `F` in `(2, 4], (4, 8], (8, 16]`. The number of terms is
  `2^k - 1` for `2^k < F <= 2^(k+1)`: a Mersenne count, the Frobenius
  tower again (the left region of width `F - 2` is refreshed by the sea's
  `2^k`-step kernels).
* `z = 1, 2`: `1/8 = (2/15)(1 - 16^-1)`, `273/2048 = (2/15)(1 - 16^-3)`: the
  partial sums of `(1/8) sum_i 16^-i` with `1, 3` terms for `F` in `(4, 8]`,
  `(8, 16]`.
* `z = 4`: `1/2`, then `261/512 = 1/2 + 1/128 + 1/512`; `z = 3, 5, 6`: `1/128`
  each from `F = 9`; `z = 7`: nothing to `F = 16`.
* Self-similarity: `c_8(F) = c_0(F - 8)`, `c_9 = c_10 = c_1(F - 8)`,
  `c_12(F) = c_4(F - 8)` for all computed `F <= 16` (the eight zeros cost
  exactly `2^-8`, and the hit pattern of `z + 8` repeats that of `z` for the
  first eight rows after arrival); the analogous `c_4(F) = c_0(F - 4)` holds
  for `F = 7, 8` and fails at `F = 9` (`261/512` against `21/32`), where more
  than four hit rows matter.

So `c_6(F) = 2/3 + 4/15 + (1/2 + 3/128 + ...)` plus one new block of about
`1.4` at each `F = 2^k + 3, 2^k + 5, 2^k + 7` (`z = 2^k, 2^k + 1, 2^k + 2,
2^k + 4`), i.e. `c_6(F) = 1.4 log_2 F + O(1)` on the computed range, which is
the observed growth
(`0.5, 0.91, 1.41, 1.47, 1.97, 2.37, 2.87` at `F = 3, 5, 7, 9, 11, 13, 15`).

## 4. Mechanism (PROVED pieces, marked)

* The defect column `F` is `6` while its right neighbour is `0`, becomes
  `4` at the first right `2`, stays `4` while the neighbour is `0`, becomes
  `2` at the next right `2`, and is then in `{0, 2}` for ever: phases of
  geometric lengths (PROVED, as in THM-4511).
* During the `6`-phase the region between the wall and the defect is the
  `{0, 6}` Pascal cone, and the wall column `c` sees a `6` exactly at the
  rows `s` with `z subset s` (Lucas); each such hit toggles it `2 <-> 4`
  (`|2 - 6| = 4`, `|4 - 6| = 2`) (PROVED for the pure-cone rows, i.e. while
  the cone's light cone stays left of `F + 1`).
* The first front: the `6` converts to a `4` at the wall and dies at the
  next `2` on its path; the path cells are a unit-triangular image of the
  left word, so exactly `F` words are destroyed by the first front alone
  (S10, PROVED), giving the tail `F 2^(1-F)`.
* Re-emissions: a wall at `4` with a `0` to its left releases a `4`; the
  column left of the wall cycles `0, 2, 2, 0, 4` (period `4`) under an
  alternating wall, and in the `4`-phase of the defect column the wall
  column itself alternates `0, 4` (`|4 - 4| = 0`, `|0 - 4| = 4`), releasing
  `4`s every two rows; each released `4` needs a zero path to column `1`
  through the evolved sea. The geometric ratios `1/4` (`z = 0`) and `1/16`
  (`z = 1, 2`) are the costs of keeping the defect column alive for two,
  respectively four, more rows per re-emission; the truncations at
  `2^k - 1` terms are where the released `4`'s path leaves the region
  refreshed by the sea's `2^k`-kernels. These identifications are read off
  the exact data and the small cases (`F = 3` by hand: the single weight-2
  word `(0, 2)` is destroyed iff the first two right bits are `0`,
  probability `1/4`, giving `13/16`); they are not a proof of the closed
  forms.

## 5. What this settles and what it does not

Settled: the size-6 extinction probability is an exact, computable rational
(Theorem A), and its excess over the front-only tail is a dyadic quantity
organised by the Frobenius tower; on the computed range its relative size
is `c_6(F)/(2F)`, flat at `6`-`10` per cent for `F <= 17`. Whether it decays
like `0.7 log_2 F / F`, which would make the front-only law `F 2^(1-F)`
asymptotically exact for size `6`, is HYP-9163's conjecture, supported by
one doubling of `F` and by the dyadic anatomy, not proved. Open: the
closed forms of section 3 beyond the computed range (HYP-9163), and the same
programme for sizes `>= 8` (whose denominators already show the Fermat
products). For Gilbreath's conjecture itself nothing changes: all sizes
`>= 6` in a real row interact through their `4`-phases, and the certified
statement remains THM-4511's (`0/2/4` rows).

Hypothesis file: [HYP-9163](../hypotheses/HYP-9163-gilbreath-size6-excess-dyadic-law.md).

## 6. Independent audit (2026-09-26)

Auditor subagent (own chain builder, exact values by solving the linear
system modulo two 61/62-bit primes and by Fraction elimination for `F <= 8`,
exhaustive full-automaton test of Theorem A over all `2^(F-1+14)` contexts
for `F <= 7`): `04-computation/experiments/gilbreath_size6_collatz_precision_20260926_audit.py`
(sha256 `4fb9ad9389589a6a1197bdb817ed651162fa2ff6018bcf04ea175f437c4837e2`) ->
`05-knowledge/results/gilbreath_size6_collatz_precision_20260926_audit.out`
(sha256 `52c57be11ff0fcdf7a5eea8251d7a027f0fb59e97542fa3da7d664e744de1db6`). CONFIRMED: Theorem A
(exact; absorption probability one on every transient state), every exact
rational `p_6(F)`, `p_8(F)` for `F <= 11`, the `2^F`-excess values, the per-`z`
table and its closed forms for `F <= 15`, the shift law, the halving, the
state counts, the size-8 denominators. CORRECTED (applied above): the
asymptotic statements were presented as settled (the decay of the excess
ratio is HYP-9163, item 4; the printed formula `0.7 log_2 F / F` evaluates to
`0.37, 0.28, 0.18, 0.17` at `F = 3, 7, 15, 17`, while the quoted numbers were
`c_6(F)/(2F)`); `1.063` for `1.0625`; the S10 comparison at `F = 12`; the
size-`4` chain values are floating; the size-8 geometric-series reading is
inferred from denominators only. Verdict for this note: HAS GAPS -> repaired.
