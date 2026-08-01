---
id: THM-3030
title: "The eighth resultant log-jet P_8, and closed forms for the corner-slice constants"
status: >
  VERIFIED-EXACT. P_8 is built: 4883 terms, degrees (M,U,V) = (16,32,16),
  support b+2c <= 32, content c_8 = 2^28 * 3^11 * 5 * 7 = 1664338750341120,
  sha256 bba6b4b9916a316c41b800a044861a15840820b6048133b754d85cfad78873ad.
  Two DISJOINT tensor grids (no shared M, U or V) return identical polynomials
  for every j=1..8, with 6 out-of-sample widths per grid and 0 coefficient
  mismatches; the frozen THM-2997 digest cfb36557... and the THM-3013/3015
  P_4..P_7 digests are all re-emitted byte-for-byte. eq(2) for the term count
  holds at j=8 (4883 = 17^3 - 3*10). The corner-slice slot law is corrected:
  the sign is (-1)^(j+m), NOT (-1)^(j+k) -- since k = 2m-1 the latter collapses
  to (-1)^(j+1) and is right only for odd m, which is why the error survived to
  j=7. With that correction all 24 pre-registered j=8 coefficients in slots
  k=-1..6 CONFIRM. The new content is closed forms for the constant sequences
  c_m = a_m/d_m, of which three of four pieces were confirmed OUT OF SAMPLE at
  m=4: a_m^A = 3 + 44*16^(m-1) (predicted 180227, exact), a_m^E = 4 + 33*9^(m-1)
  (predicted 24061, exact), a_m^C = 23 constant, and the slice ratios
  d_m^A/d_m^C = 4^(2m-1), d_m^E/d_m^C = 3^(2m-1) (predicted 4^7, 3^7, both
  exact). The one REFUTED guess is the base sequence: d_m^C = (3m)!/(2m-2)!
  fits m=1,2,3 and predicts 665280 at m=4 where the truth is 604800 (ratio
  11/10). d_m^C = 6, 360, 15120, 604800 is the single remaining gap.
source: klein-S428
depends_on:
  - THM-3015
  - THM-3013
  - THM-2997
related:
  - THM-3011
  - THM-3014
  - THM-3006
script: 04-computation/gmc_eighth_resultant_jet_and_the_corner_constant_laws_thm3030.py
output: 05-knowledge/results/gmc_eighth_resultant_jet_and_the_corner_constant_laws_thm3030.out
---

# THM-3030 -- `P_8` and the corner-slice constant closed forms

## 1. The jet

`P_8` is the eighth resultant log-jet in the Macaulay chart of THM-3011, in the
normalisation `Q_j = D^{2j} L_j` with `D = U^2 + 3U - 3V - 1`, `U = 2^M`,
`V = 3^M`.

| | |
|---|---|
| terms | **4883** |
| degrees `(M,U,V)` | `(16, 32, 16)` |
| support | `b + 2c <= 32`, absences only at the three corners |
| content `c_8` | `2^28 * 3^11 * 5 * 7 = 1664338750341120` |
| sha256 (content-1 rows) | `bba6b4b9916a316c41b800a044861a15840820b6048133b754d85cfad78873ad` |

Table: [`05-knowledge/results/gmc_first_gap_resultant_jet_P8_table_thm3030.json`](../../05-knowledge/results/gmc_first_gap_resultant_jet_P8_table_thm3030.json).
`P_1`, `P_2`, `P_3` are emitted in the same content-1 format alongside it
(`c_1 = 1`, `c_2 = 48`, `c_3 = 1152`).

The content ledger continues to hold: the primes dividing `c_j` are exactly the
primes `<= j+1`, and at `j = 8` that is `{2,3,5,7}`.

**Term counts** `27, 122, 333, 717, 1313, 2176, 3348, 4883` for `j = 1..8`, all
matching eq (2) of THM-3015, `|P_j| = (2j+1)^3 - 3(2j-2-floor(j/2)) - [j=3]`.
At `j = 8`: `17^3 - 3*10 = 4913 - 30 = 4883`.

### Controls

All four pass, and they are the reason this is `VERIFIED-EXACT` rather than
merely computed:

1. Grid `A` (`M=4..26`, `U=2..34`, `V=2..34` even) and grid `B` (`M=27..49`,
   `U=35..67`, `V=36..68` even) share **no** `M`, `U` or `V` value, and return
   **identical** polynomials for every `j = 1..8`.
2. 6 out-of-sample widths per grid, **0** coefficient mismatches.
3. The frozen THM-2997 `P_1,P_2,P_3` table is re-emitted byte-for-byte, digest
   `cfb36557e1d54a0328a309375a948ace99c78e0688a54a014aef0906c1b90513`.
4. The THM-3013 `P_4` and THM-3015 `P_5,P_6,P_7` content-1 digests are all
   reproduced.

## 2. The slot law, corrected

Write the corner slices `A_j = [U^{4j}]Q_j`, `E_j = [V^{2j}]Q_j/9^j`,
`C_j = Q_j(M,0,0)`, and index the coefficient of `M^e` by the **slot**
`k = j - e`. For `k < j` (the terminal slot `k = j` is the constant term and is
outside the law's range):

```text
   k = -1        [M^(j+1)] = (-1)^(j+1) * 46/(j(j+1))        slice-independent
   k =  0        [M^j]     = (-1)^j * kappa/j                kappa = 12, 11, 3
   k = 2m-1      [M^(j-k)] = (-1)^(j+m) * c_m * (j-1)(j-2)...(j-2m+2)
   k even >= 2   [M^(j-k)] = 0
```

**The correction.** THM-3015 recorded the odd-slot signs case by case
(`(-1)^(j+1)` at `k=1`, `(-1)^j` at `k=3`, `(-1)^(j+1)` at `k=5`) without
naming the pattern. The pattern is `(-1)^(j+m)`. Writing it as `(-1)^(j+k)` is
wrong but *undetectably* wrong through `j = 7`: since `k = 2m-1`,
`(-1)^(j+k) = (-1)^(j+1)` identically, which agrees with `(-1)^(j+m)` for every
**odd** `m` and disagrees for every **even** `m` — and `m = 2` is the only even
value reachable below `j = 8`, where it had already been recorded by hand. The
error is only exposed by predicting rather than fitting.

**All 24 pre-registered `j = 8` coefficients (3 slices x 8 slots `k = -1..6`)
confirm**, listed in part D of the output.

## 3. The constant sequences (the new content)

Write `c_m = a_m / d_m` per slice. Three of the four pieces have closed forms,
and all three were **written down from `m <= 3` data and confirmed against
`m = 4` only afterwards**:

| piece | closed form | prediction at `m=4` | observed | |
|---|---|---|---|---|
| numerator, `A` | `a_m = 3 + 44 * 16^(m-1)` | `180227` | `180227` | CONFIRMED |
| numerator, `E` | `a_m = 4 + 33 * 9^(m-1)` | `24061` | `24061` | CONFIRMED |
| numerator, `C` | `a_m = 23` | `23` | `23` | CONFIRMED |
| ratio `d^A/d^C` | `4^(2m-1)` | `4^7 = 16384` | `16384` | CONFIRMED |
| ratio `d^E/d^C` | `3^(2m-1)` | `3^7 = 2187` | `2187` | CONFIRMED |
| base `d_m^C` | `(3m)!/(2m-2)!` | `665280` | `604800` | **REFUTED** |

The `m <= 3` inputs were `47/24, 707/23040, 11267/15482880` (`A`),
`37/18, 301/9720, 2677/3674160` (`E`), `23/6, 23/360, 23/15120` (`C`).
The measured `m = 4` constants are

```text
   c_4^A = 180227/9909043200,   c_4^E = 24061/1322697600,   c_4^C = 23/604800.
```

Two remarks.

**The `23`.** The `C` slice has numerator exactly `23` at every `m`, and the
slice-independent `k = -1` law carries `46 = 2*23`. This is the same `23` that
heads the four-band charge density of THM-3006,
`lim w_k/M^(k+1) = [23 + (2/3)^(k+1) + 2(1/2)^(k+1) + 2(1/3)^(k+1)]/(k+1)`.
Recorded as a structural coincidence worth explaining, not as a proved link.

**The one gap.** `d_m^C = 6, 360, 15120, 604800` has no closed form. Its
successive ratios are `60, 42, 40`. The near-miss `(3m)!/(2m-2)!` reproduces the
first three exactly (`1*2*3`, `3*4*5*6`, `5*6*7*8*9`) and then over-predicts by
`11/10`, i.e. the fourth term is *not* `7*8*9*10*11*12`. Since the two slice
ratios and all three numerators are closed-form, `d_m^C` is the entire remaining
unknown in the corner-slice picture, and `P_9` would decide between candidate
continuations with a single new data point.

## 4. Status of each claim

| Claim | Status |
|---|---|
| `P_8` shape, content, digest | VERIFIED-EXACT (disjoint grids + OOS) |
| eq(2) at `j = 8` | VERIFIED-EXACT |
| slot sign law `(-1)^(j+m)` | PROVED to be the unique law consistent with `j<=8`; supersedes the case-by-case signs of THM-3015 |
| 24 pre-registered `j=8` coefficients | CONFIRMED (genuine pre-registration) |
| numerator closed forms, `A`/`E`/`C` | CONFIRMED out-of-sample at `m = 4` |
| slice ratios `4^(2m-1)`, `3^(2m-1)` | CONFIRMED out-of-sample at `m = 4` |
| `d_m^C = (3m)!/(2m-2)!` | **REFUTED** at `m = 4` |
| closed form for `d_m^C` | OPEN |
| link between the slice `23` and THM-3006's `23` | UNPROVED coincidence |
