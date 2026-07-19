---
id: THM-1258
title: Scale-thirty-six complementary four-nine fibre obstruction
status: PROVED STRUCTURAL + FINITE-EXACT — all 206,725,596 labelled support/order contexts are exhausted. Scalar capacity leaves 82,332 rows; the sound Z/4 relaxation leaves 7,824 rows live at all six owners and the sound Z/9 relaxation leaves 12,222, but their all-owner intersection is empty. Normal and optimized exact runs agree byte-for-byte; the two generic logical consumers are sorry-free Lean.
source: codex-2026-07-19-S78 scale-36 sporadic-frontier audit
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-1249]
related: [HYP-6820]
verification:
  - 04-computation/lrc13_scale_thirty_six_hamming_six_complementary_fibre_obstruction_codex_c36.py
  - 05-knowledge/results/lrc13_scale_thirty_six_hamming_six_complementary_fibre_obstruction_codex_c36.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCScaleThirtySixComplementaryFibre.lean
---

# THM-1258 — complementary four- and nine-fibres close scale thirty-six

> **The primitive proper AP-centred common-scale-36 Hamming-six face is
> empty.**

The proof visits every labelled support and hereditary effective-order word.
Unit-independent scalar owner capacity leaves 82,332 rows.  Two complementary
upper relaxations then finish the bank: one retains all order-dividing-four
masks exactly, while the other retains all order-dividing-nine masks exactly.
A literal cover must survive both relaxations at all six owners, but no scalar
row does.

This is a finite face theorem in the AP-centred shallow Hamming-six chart.  It
does not prove AP extraction, arbitrary-height ballot rigidity, or uniform
emptiness of the twelve-runner sporadic branch.

## 1. Hereditary divisor grammar

THM-860 associates to the six replaced labels an effective-order word

```text
(D_1,...,D_6),                 D_i | 36,
lcm(D_i : i != j)=36           for every j.                 (1)
```

The divisor and exact-unit ledgers are

```text
D          1  2  3  4  6  9  12 18 36
# units    1  1  2  2  2  6   4  6 12.                    (2)
```

Because `36=4*9`, (1) is equivalent to having at least two coordinates
divisible by four and at least two divisible by nine.  Split the nine divisors
according to whether they carry the full factors four and nine.  The four
category sizes are

```text
neither=4,       four-only=2,       nine-only=2,       both=1.
```

Inclusion-exclusion therefore gives

```text
9^6 - 2(6^6+6*3*6^5) + 65,536 = 223,729             (3)
```

hereditary order words.  Weighting each coordinate by `phi(D)` gives

```text
1,904,124,672                                                (4)
```

literal unit words per labelled support.  Across all `binom(12,6)=924`
supports, the exact banks thus contain

```text
206,725,596 labelled support/order contexts,
1,759,411,196,928 raw labelled unit states.                  (5)
```

The program checks the carrier description against all six literal
leave-one-out lcm equations on every one of the `9^6` order words.  It does
not quotient supports or order words by symmetry.

## 2. Literal masks and the scalar residual

For provider label `a`, owner label `b`, effective order `D`, and exact unit
`u`, let `B` be the unique CRT class

```text
B = D*a (mod 13),                   B = u (mod D).            (6)
```

The owner-sheet mask on `Z/36Z` is

```text
M(a,D,u;b)
 = {t : <B*(b^(-1)+13t)>_(13D) lies in (-D,D]}.              (7)
```

The script constructs (6) both algebraically and by bounded literal search,
then checks every bit of (7), its exact cardinality, and its period `D`.
The unit-independent cardinality vectors, indexed by the ratio `r=a/b` and
with columns ordered as in (2), are

```text
r= 1: 36 18 12  9  6 8 6 6 6      r= 7:  0 18  0  9 6 4 6 6 6
r= 2:  0  0  0  0  6 8 6 6 6      r= 8:  0  0 12  0 6 8 6 6 5
r= 3:  0  0  0  9  6 4 6 6 6      r= 9:  0  0 12  9 6 4 6 4 5
r= 4:  0  0 12  9  6 4 6 4 5      r=10:  0  0  0  9 6 4 6 6 6
r= 5:  0  0 12  0  6 8 6 6 5      r=11:  0  0  0  0 6 8 6 6 6
r= 6:  0 18  0  9  6 4 6 6 6      r=12:  0  0  0  0 0 4 3 4 5. (8)
```

Summing these exact cardinalities over the six providers is a necessary
capacity test at each owner.  The number of contexts having exactly `k`
owners with scalar capacity at least 36 is

```text
k=0:  2,400,018       k=1: 49,776,708       k=2: 77,652,252
k=3: 56,649,096       k=4: 17,848,230       k=5:  2,316,960
k=6:     82,332.                                             (9)
```

Those 82,332 rows occur on 860 supports, have 312 order-multiplicity profiles
and 37,247 distinct ordered capacity vectors, and represent 698,516,352 exact
unit words.  Every other context already fails a literal owner because union
cardinality is at most the scalar sum.

## 3. Two sound thick-fibre relaxations

For `q in {4,9}`, define the anchor orders

```text
A_4={1,2,4},                       A_9={1,3,9}.              (10)
```

If `D|q`, the period-`D` law checked above implies period `q`; hence every
anchor mask is a union of full fibres of `Z/36Z -> Z/qZ`.  This periodicity
makes the exact anchor bank small, but soundness uses only the following
elementary union inequality.

Fix a surviving row, an owner, and an exact tuple of anchor units, and put

```text
Q_q = union_(D_i in A_q) M_i(u_i).
```

For every simultaneous choice of the remaining units,

```text
|union_i M_i(u_i)|
 <= |Q_q| + sum_(D_i notin A_q) |M_i(u_i) minus Q_q|
 <= |Q_q| + sum_(D_i notin A_q)
                  max_e |M_i(e) minus Q_q|.                  (11)
```

Finally maximize the right side over the finite exact anchor bank and call
the result `U_q`.  This deliberately forgets shared nonanchor units and all
mutual nonanchor overlaps.  Both losses can only increase (11), so

```text
literal cover at owner b  ==>  U_4(b)>=36 and U_9(b)>=36.    (12)
```

The complete live-owner histograms on the 82,332 scalar rows are

```text
# live owners     0      1      2     3     4     5     6
Z/4           44,946 11,208 10,722 5,172 2,280   180  7,824
Z/9           45,600 10,752  4,428 5,436 2,550 1,344 12,222. (13)
```

Neither flag alone closes the face.  Their intersection does:

```text
# rows with six Z/4-live owners and six Z/9-live owners = 0. (14)
```

More explicitly, the 7,824 all-Z/4-live rows split as
`(7,338; 468; 18)` by `0,1,2` Z/9-live owners.  The 12,222 all-Z/9-live rows
split as `(10,278; 1,392; 480; 72)` by `0,1,2,3` Z/4-live owners.  There is
no hidden `(6,6)` cell.

If a literal packet covered all 36 sheets at every owner, its context would
survive (9), and (12) would put its row in the missing cell (14).  This
contradiction proves the theorem.

## 4. Tournament and carrier audit

The useful vertices are the six owner obligations, not the six runners or
providers.  In each quotient gauge compare the lexicographic key

```text
(1[U_q>=36], U_q, scalar capacity, anchor-bank size),        (15)
```

with owner-coordinate order as the tie Hamiltonian path.  Every completed
tournament is transitive: score word `(0,1,2,3,4,5)`, no directed triangle,
six singleton SCCs, and one Hamiltonian path.  There are 449 joint
tie/flip fingerprints; changing from the four-fibre to the nine-fibre gauge
flips between zero and fourteen edges.

This telemetry does not prove (14).  Runner vertices retain speed order but
lose the owner threshold; provider vertices retain ratio data but lose the
common sheet union; divisor vertices retain ramification but lose the
embedded support.  The paired owner-coloured quotient masks are the smallest
carrier here that preserves the absolute noncoverage predicate.  The
challenged assumption is that a single prime-power quotient must suffice:
at scale 36, neither does, while their fibre product is terminal.

## 5. Replay, formal consumer, and scope

Normal and `python -O` executions reproduce the stored output byte-for-byte.
Frozen SHA-256 values are

```text
Python source    2fe0755fa880c22f94c30ceec88dca38eb1cbc00c1801a51ff8d4f04aeebfd76
Python output    d8a19d72873957c3a16429356cacfa88a75730190491938c1263f3ff768d84c1
Lean consumer    fee02fe0ec5d3bdbf9a2e87cb72bc763373321d2afd18f8db67a60125fc4f5ea
```

`LRCScaleThirtySixComplementaryFibre.lean` is sorry-free and checks two
generic facts: an anchor/deviation ceiling below 36 cannot equal the full
36-sheet set, and one owner failing either of two necessary upper-bound gauges
prevents all-owner cover.  The 206,725,596-context enumeration remains the
external exact certificate.

This advances the primitive proper AP-centred H6 common-scale frontier through
36.  Scale 37 is already excluded by THM-983's uniform prime-scale theorem;
the next untreated numerical common scale on this line is 38.  Smooth H5
ramification, non-AP-centred shallow packets, deep two-sheet packets, higher
sheets, and uniform n=12 sporadic emptiness remain open.  ∎
