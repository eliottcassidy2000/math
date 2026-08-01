---
id: THM-3030
title: "The eighth resultant log-jet P_8, and closed forms for the corner-slice constants"
status: >
  FINITE-EXACT TABLE-INTERNAL (j<=8) + INDEPENDENTLY HOSTILE-AUDITED;
  all-order continuation OPEN.  The frozen P_8 table has 4883 terms, degrees
  (M,U,V)=(16,32,16), support b+2c<=32, content
  2^28*3^11*5*7, and LF hash bba6b4b9....  Its historically reported two-grid
  interpolation and out-of-sample build are NOT replayed by the canonical
  companion and no build artifact is stored, so that provenance is reported,
  not independently verified.  Direct table arithmetic verifies all visible
  slot laws, repairs the odd-slot sign to (-1)^(j+m), and proves on j<=8 the
  exact Faulhaber identity p_j=46 sum_(s<M)s^j+20M^j+K_j.  Consequently the
  measured C-slice constants are 46|B_(2m)|/(2m)! for m=1..4, closing the
  observed denominator sequence 6,360,15120,604800.  Extension to all m is a
  conjecture.  P_10, not P_9, is its first new nonterminal test; P_12 would see
  the first Bernoulli-numerator break through B_12=-691/2730.
source: klein-S428
audit: >
  Independent hostile audit ACCEPTS the frozen-table claims after scope repair.
  A dependency-free referee checks all eight table hashes and shapes, 48 visible
  odd slots, 36 even-zero slots, the j=8 sign hostile, eight exact width
  recurrences, and the Bernoulli-Faulhaber identity.  It identified three truth
  boundaries recorded in MISTAKE-342: finite pattern uniqueness is not an
  all-order proof; the advertised interpolation-grid/OOS controls are absent
  from the executable companion; and P_9 cannot test c_5 because k=9 is its
  exceptional terminal slot.  Primary and independent normal/-O/stored
  transcripts byte-match.
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
script_sha256: c730b9d6c49ff2909abd274b7413bc7f5957b9348bffb83b340ab1d1a93c70b6
output_sha256: e32ee36076e26be2fde8eb0dc28b0c8c21a06f1ac8dca2601e8315eec640220c
independent_script: 04-computation/gmc_eighth_resultant_jet_bernoulli_faulhaber_referee_thm3030.py
independent_output: 05-knowledge/results/gmc_eighth_resultant_jet_bernoulli_faulhaber_referee_thm3030.out
independent_script_sha256: 045282d09a255a55ccdb6a42af86948c88ed1e03115b926c3f0d72f83f807913
independent_output_sha256: 1edb2bca1c15c686f8f4f86200ea36c2fdffb93d8966b44137f817df1229a96a
hash_basis: LF-normalized bytes
---

# THM-3030 -- `P_8` and the corner-slice constant closed forms

**FINITE-EXACT TABLE-INTERNAL (`j<=8`) + INDEPENDENTLY HOSTILE-AUDITED;
ALL-ORDER CONTINUATION OPEN.**

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

### Controls and evidence boundary

The canonical primary companion reproducibly checks the frozen tables: it
re-emits THM-2997's `P_1,P_2,P_3` rows byte-for-byte, reproduces the
THM-3013/3015 `P_4,...,P_7` digests, checks `P_8`'s shape and support, and
replays every displayed corner coefficient.  The independent referee reads no
primary code; it checks all eight table hashes and shapes directly.

The historical build report says that grid `A` (`M=4..26`, `U=2..34`,
`V=2..34` even) and disjoint grid `B` (`M=27..49`, `U=35..67`, `V=36..68`
even) returned identical jets and that six out-of-sample widths per grid had
zero mismatches.  **Those grid computations are not present in either
canonical script, and the referenced interpolation pickle is not stored.**
They are therefore `REPORTED BUILD PROVENANCE`, not a reproducible control of
this promotion.  The theorem's audited status is deliberately table-internal.

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

**The finite correction.** THM-3015 recorded the odd-slot signs case by case
(`(-1)^(j+1)` at `k=1`, `(-1)^j` at `k=3`, `(-1)^(j+1)` at `k=5`) without
naming the pattern. The pattern is `(-1)^(j+m)`. Writing it as `(-1)^(j+k)` is
wrong but *undetectably* wrong through `j = 7`: since `k = 2m-1`,
`(-1)^(j+k) = (-1)^(j+1)` identically, which agrees with `(-1)^(j+m)` for every
**odd** `m` and disagrees for every **even** `m` — and `m = 2` is the only even
value reachable below `j = 8`, where it had already been recorded by hand. The
error is only exposed by predicting rather than fitting.

**All 24 pre-registered `j = 8` coefficients (3 slices x 8 slots `k = -1..6`)
confirm**, listed in part D of the output.  Across all frozen tables the
independent referee checks `16` visible odd slots per slice and `36` visible
even-zero slots.  Thus `(-1)^(j+m)` is **FINITE-EXACT through `j=8`**.  No
finite list makes it a proved all-order law; that continuation remains open.

## 3. The Bernoulli--Faulhaber mechanism

Put

    p_j^C(M)=(-1)^(j-1) j C_j(M).

The independent referee finds the following **exact polynomial identity for
every frozen order `1<=j<=8`**:

    p_j^C(M)=46 sum_(s=1)^(M-1) s^j + 20 M^j + K_j,      (3)

where `K_j` is independent of `M`.  Equivalently,

    p_j^C(M+1)-p_j^C(M)=26M^j+20(M+1)^j.                (4)

All eight recurrences are checked coefficientwise over `Q`; this is not a
floating fit.  Faulhaber's formula now gives, whenever `2m-1<j<=8`,

    [M^(j-2m+1)] C_j
      =(-1)^(j+1) 46 B_(2m)/(2m)! (j-1)...(j-2m+2).    (5)

Since `sign B_(2m)=(-1)^(m+1)`, (5) is exactly the corrected sign law with

    c_m^C = 46 |B_(2m)|/(2m)!,        m=1,2,3,4.       (6)

Thus the formerly unexplained measured sequence is

    c_m^C = 23/6, 23/360, 23/15120, 23/604800,

and its reduced denominators are `6,360,15120,604800`.  The failed guess
`(3m)!/(2m-2)!` was approximating the first Bernoulli values; it is replaced,
on the measured range, by the exact formula (6).

The separately pre-registered slice formulas also survive at `m=4`:

| finite formula | predicted at `m=4` | observed | status |
|---|---:|---:|---|
| `a_m^A=3+44*16^(m-1)` | `180227` | `180227` | CONFIRMED |
| `a_m^E=4+33*9^(m-1)` | `24061` | `24061` | CONFIRMED |
| reduced `a_m^C=23` | `23` | `23` | CONFIRMED through `m=4` |
| `d_m^A/d_m^C=4^(2m-1)` | `16384` | `16384` | CONFIRMED |
| `d_m^E/d_m^C=3^(2m-1)` | `2187` | `2187` | CONFIRMED |

These are finite confirmations, not all-order closed forms.  In particular,
the Bernoulli continuation predicts that the reduced `C` numerator will not
remain `23`: `B_12=-691/2730` gives a factor `691` at `m=6`.

**What the `23` now means.**  Within the frozen `C` slice, `23` is no longer an
unexplained numerator: it is half of the exact bulk density `46` in (3), after
the Bernoulli coefficient and log-jet normalization.  Its equality with the
terminal band density `23` in THM-3006 remains an unproved connection between
two different objects.

## 4. Sharp future tests

Formula (6) beyond `m=4` is `CONJECTURAL`.  Its first new nonterminal slot is
`k=9` at **`j=10`**, not `j=9`; at `j=9`, `k=j` is the exceptional constant
term and lies outside the law.  The exact preregistered predictions are

    P_10: c_5^C = 46|B_10|/10! = 23/23950080,
    P_12: c_6^C = 46|B_12|/12! = 15893/653837184000.

The `P_12` value is the first discriminator between the Bernoulli mechanism
and the observed reduced-numerator-`23` extrapolation because
`15893=23*691`.  A more structural proof target is the width quotient suggested
by (4): establish the `C`-corner log-jet recurrence at the determinant level,
which would prove (3)--(6) simultaneously at every order.

## 5. Status of each claim

| Claim | Status |
|---|---|
| `P_8` frozen shape, content, digest | FINITE-EXACT table-internal |
| historical disjoint grids + OOS | REPORTED; build/pickle absent from canon |
| eq(2) at `j = 8` | FINITE-EXACT |
| slot sign law `(-1)^(j+m)` | FINITE-EXACT on every visible slot through `j=8`; all-order OPEN |
| 24 pre-registered `j=8` coefficients | CONFIRMED (genuine pre-registration) |
| numerator/slice-ratio formulas | CONFIRMED at the new `m=4` datum only |
| `d_m^C = (3m)!/(2m-2)!` | **REFUTED** at `m = 4` |
| `c_m^C=46|B_(2m)|/(2m)!` | PROVED from frozen tables for `m<=4`; all-order OPEN |
| first new test | `P_10`; `P_9` is terminal and cannot test it |
| link between the slice `23` and THM-3006's `23` | UNPROVED coincidence |

## 6. Reproduction

    python3 04-computation/gmc_eighth_resultant_jet_and_the_corner_constant_laws_thm3030.py
    python3 -O 04-computation/gmc_eighth_resultant_jet_and_the_corner_constant_laws_thm3030.py
    python3 04-computation/gmc_eighth_resultant_jet_bernoulli_faulhaber_referee_thm3030.py
    python3 -O 04-computation/gmc_eighth_resultant_jet_bernoulli_faulhaber_referee_thm3030.py

The primary transcript is `7,363` LF bytes and the independent transcript is
`1,264` LF bytes.  Their hashes and both script hashes are frozen in
frontmatter.  Normal, optimized, and stored output agree byte-for-byte.

**QED for the finite table statements.**
