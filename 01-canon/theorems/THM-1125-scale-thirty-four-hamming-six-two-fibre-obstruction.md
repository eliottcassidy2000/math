---
id: THM-1125
title: Scale-thirty-four Hamming-six two-fibre obstruction
status: PROVED STRUCTURAL + INDEPENDENT FINITE-EXACT — independent algebraic/literal-CRT Python and bounded-literal-CRT standard-library C++20 implementations exhaust all 3,249 hereditary words and 3,002,076 labelled contexts representing 1,271,272,375,296 raw unit states. Scalar capacity leaves 552 rows on 36 supports. The sound Z/2 thick-fibre relaxation bounds all 3,312 surviving owner obligations by 29<34, so every row is impossible at every owner. Normal/Python-optimized and C++ O3/O0/sanitized runs replay the frozen outputs byte-for-byte.
source: codex-2026-07-18 c34 frontier continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-994, THM-1072, THM-1090, THM-1096, THM-1124]
related: [HYP-6820]
verification:
  - 04-computation/lrc13_scale_thirty_four_hamming_six_two_fibre_obstruction_codex_c34.py
  - 05-knowledge/results/lrc13_scale_thirty_four_hamming_six_two_fibre_obstruction_codex_c34.out
  - 04-computation/lrc13_scale_thirty_four_hamming_six_two_fibre_referee_codex_c34.cpp
  - 05-knowledge/results/lrc13_scale_thirty_four_hamming_six_two_fibre_referee_codex_c34.out
---

# THM-1125 — scale thirty-four has a terminal two-fibre deficit

The theorem is:

> **The primitive proper AP-centred common-scale-34 Hamming-six face is
> empty.**

After scalar capacity, retain every order-one and order-two mask exactly as a
union of the two thick fibres of `Z/34 -> Z/2`.  For a fixed exact anchor union,
allow every order-seventeen and order-thirty-four provider to choose its unit
independently.  This forgets global unit compatibility and nonanchor overlaps,
so it is an upper relaxation.  Its largest value is `29<34` on every one of the
3,312 scalar-surviving owner obligations.  Thus every scalar row fails at every
owner, a stronger terminality statement than the hereditary bridge requires.

## 1. Squarefree hereditary grammar

The effective orders and exact-unit counts are

```text
D                 1   2   17   34
# units            1   1   16   16.
```

For an order word `(D_1,...,D_6)`, hereditary primitivity requires

```text
lcm(D_i : i != j)=34                    for every j.          (1)
```

Since `34=2*17` is squarefree, (1) is equivalent to

```text
at least two coordinates carry 2,
at least two coordinates carry 17.                              (2)
```

Each prime carrier can be any subset of six coordinates except the empty set
and the six singletons.  Hence the number of order words is

```text
(2^6-1-6)^2 = 57^2 = 3249.                                   (3)
```

The unit-weighted count per labelled support factors across the two carriers:

```text
(2^6-1-6) * (17^6-1-6*16) = 1,375,835,904.                   (4)
```

Across all `binom(12,6)=924` supports, this gives

```text
3,002,076 labelled support/order contexts,
1,271,272,375,296 raw labelled unit states.                    (5)
```

Both implementations check (2) against the literal leave-one-out lcm on all
`4^6` words and check the two factorizations independently.

## 2. Literal masks and scalar capacity

For provider label `a`, owner label `b`, effective order `D`, and exact unit
`u`, let `B` be the unique CRT representative

```text
B = D*a (mod 13),                    B = u (mod D).             (6)
```

The owner-sheet mask is

```text
M(a,D,u;b)
 = {t in Z/34 : <B*(b^(-1)+13t)>_(13D) lies in (-D,D]}.        (7)
```

Writing `r=a/b` in `F_13^*`, literal interval counting gives the
unit-independent cardinalities below, with columns `D=1,2,17,34`:

```text
r= 1: 34 17 6 6          r= 7: 0 17 6 5
r= 2:  0  0 4 5          r= 8: 0  0 4 5
r= 3:  0  0 6 5          r= 9: 0  0 6 6
r= 4:  0  0 6 6          r=10: 0  0 6 5
r= 5:  0  0 4 5          r=11: 0  0 4 5
r= 6:  0 17 6 5          r=12: 0  0 4 5.                     (8)
```

Summing (8) over the six providers is a necessary owner-capacity test.  Over
all contexts in (5), the number of owners reaching capacity 34 has histogram

```text
0:33112, 1:387180, 2:1283004, 3:1099712,
4:180156, 5:18360, 6:552.                                    (9)
```

Exactly 552 rows on 36 supports survive at all six owners.  The support
histogram is

```text
0:888, 15:24, 16:12.                                        (10)
```

All scalar survivors have precisely two order-two coordinates and no order-one
coordinate.  Their five multiplicity profiles are

```text
(0,2,0,4):36, (0,2,1,3):144, (0,2,2,2):216,
(0,2,3,1):144, (0,2,4,0):12.                                (11)
```

They represent 36,175,872 literal unit words and have 356 distinct scalar
capacity vectors.  Every row omitted here is already impossible at some owner
by the scalar union bound.

## 3. The sound `Z/2` relaxation

Replacing `t` by `t+D` in (7) changes the CRT argument by a multiple of `13D`.
Every order-`D` mask is therefore pulled back from `Z/DZ`.  When `D|2`, it is
in particular an exact union of the two fibres of `Z/34 -> Z/2`.

Fix an owner and an exact anchor-unit tuple for the providers with
`D in {1,2}`.  Write

```text
Q = union_(D_i|2) M_i(u_i).
```

For arbitrary choices of all nonanchors,

```text
|union_i M_i(u_i)|
 <= |Q| + sum_(D_i in {17,34}) max_e |M_i(e) minus Q|.        (12)
```

Maximizing (12) over the exact anchor bank gives a sound owner-local upper
bound `U_2`.  The independent nonanchor maxima discard their shared global unit
word and their mutual intersections, but both losses can only increase the
right-hand side.

For the 552 scalar rows, the exact anchor bank has size one at all 3,312 owner
obligations.  The complete bound histogram is

```text
25:48, 26:1056, 27:408, 28:1728, 29:72.                     (13)
```

In particular

```text
U_2 <= 29 < 34                                               (14)
```

at **every** scalar-surviving owner.  Split by the owner's own order, the
stronger carrier data are

```text
own order 2:   25:48, 26:288, 27:408, 28:288, 29:72
own order 17:  26:384, 28:672
own order 34:  26:384, 28:768.                              (15)
```

Thus all 552 rows have zero live owners after the relaxation.  A literal unit
word capable of global coverage would have to cover all 34 sheets at all six
owners, contradicting (14).  Equations (1), (9), (12), and (14) prove the
stated emptiness.

## 4. Complementary flags and the loss ledger

The primary evaluates two nonessential competitors on the same complete scalar
bank.  The transverse `Z/17` flag retains orders `1,17` exactly and relaxes
orders `2,34`; it leaves 180 rows live at all six owners and has bounds as high
as 40.  A mixed flag retaining all proper orders `1,2,17` has the same threshold
histogram as `Z/2`, but exact anchor banks as large as 11,292 rather than one.

So the terminal quotient is selected by both proof strength and compression:
`Z/2` loses more fine unit information yet preserves the absolute predicate
needed by (12).  It destroys shared nonanchor units and mutual needle overlaps.
Those losses are sound only because the surviving number is an upper bound;
the same quotient would not preserve exact reachability or a lower bound.

## 5. Tournament Analysis and challenged vertices

Runner tournaments are not the proof carrier.  Use the six owner obligations
as vertices and compare the lexicographic cost key

```text
(1[U_q>=34], U_q, scalar capacity, anchor-bank size)          (16)
```

in the `q=2` and `q=17` gauges; the harder obligation wins and coordinate order
breaks ties.  All 552 completed tournaments are transitive, with score word
`(0,1,2,3,4,5)`, no directed triangle, six singleton SCCs, and one Hamiltonian
path.  Switching gauges flips as many as eight edges.  This ordering telemetry
does not imply (14); the absolute owner threshold does.

An alternate tournament uses the three quotient flags `Z/2`, `Z/17`, and
mixed anchors as vertices, orienting toward the smaller upper bound and then
the smaller anchor bank.  `Z/2` is first on all 3,312 obligations.  This
challenges several vertex assumptions:

- provider vertices keep label ratios but lose the owner-cover threshold;
- sheet vertices keep incidence but obscure the thick-fibre compression;
- quotient flags compare proof cost but forget which owner is impossible;
- owner-coloured `Z/2` fibre incidence is the faithful finite carrier.

The tournaments are deliberately reported as lossy sidecars rather than
mistaken for the proof object.

## 6. Independent replay, integrity, and scope

The Python primary checks CRT bases both algebraically and by bounded literal
search, verifies every mask cardinality and period, exhausts the complete
labelled scalar bank, and evaluates the `Z/2`, `Z/17`, and mixed flags plus the
tournament sidecars.  Normal and `python -O` runs reproduce its frozen output.
Its SHA-256 values are

```text
Python source   6a91cb4f190a0c5ff4bc836a31ff20b417f64296f1efbca8bb06f70fab88a632
Python output   45c62249efefb185732e427562275b46276e1b70c8c41b63e91d254681344350
```

The independently developed C++20 referee starts from bounded literal CRT
search, enumerates supports as twelve-bit words and order words in base four,
and intentionally omits the primary's algebraic CRT, NumPy enumeration,
competing flags, and tournament telemetry.  It independently reproduces every
proof-facing count in (3)--(15).  Warning-clean O3/O0 and ASan+UBSan builds
replay its frozen output byte-for-byte.  Its SHA-256 values are

```text
C++ source      320e9452791824133b1ec12236bd2337a71cbd7412a7048f7c06082266df76a3
C++ output      f1c52da15c52a01a46b95a3d922efa28544a64a1b3ad72ac11afd264eaf6b9f8
```

This theorem closes only the primitive proper AP-centred common-scale-34
Hamming-six face.  It does not close the smooth Hamming-five bank,
non-AP/deep sheets, the global sporadic branch, or LRC(14).  Scale 35 is the
next untreated common-scale face on this line.  ∎
