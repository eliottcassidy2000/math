---
id: THM-1124
title: Scale-thirty-three Hamming-six three-fibre obstruction
status: PROVED STRUCTURAL + INDEPENDENT FINITE-EXACT — a deterministic all-labelled algebraic/literal-CRT Python primary and a separately developed standard-library C++20 bounded-literal-CRT referee agree on all 3,249 hereditary words, 3,002,076 labelled support/order contexts representing 1,171,996,056,000 literal unit states, the complete 1,344-row scalar bank, every Z/3 anchor-bank and owner-bound histogram, and the terminal own-order-eleven ceiling 29. Hereditary lcm forces at least two such owners. Normal/Python-optimized and C++ O3/O0 runs replay both frozen outputs byte-for-byte.
source: codex-2026-07-18 c33 frontier continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-994, THM-1072, THM-1090, THM-1096]
related: [HYP-6820]
verification:
  - 04-computation/lrc13_scale_thirty_three_hamming_six_frontier_scout_codex_c33.py
  - 05-knowledge/results/lrc13_scale_thirty_three_hamming_six_frontier_scout_codex_c33.out
  - 04-computation/lrc13_scale_thirty_three_hamming_six_three_fibre_referee_codex_c33.cpp
  - 05-knowledge/results/lrc13_scale_thirty_three_hamming_six_three_fibre_referee_codex_c33.out
---

# THM-1124 — scale thirty-three has a terminal three-fibre deficit

The theorem is:

> **The primitive proper AP-centred common-scale-33 Hamming-six face is
> empty.**

After the complete scalar filter, retain every order-one and order-three mask
exactly as a union of the three thick fibres of

```text
Z/33 -> Z/3.
```

Let every order-eleven or order-thirty-three provider choose its unit
independently outside that exact anchor union.  This is a sound upper
relaxation.  At any owner whose own coordinate has order `11` or `33`, the
relaxed union has size at most `29<33`.  Hereditary leave-one-out lcm forces
at least two coordinates to have one of those two orders.  Those matching
owners are therefore terminal, and no global unit word exists.

## 1. Hereditary product grammar

The effective orders and unit counts are

```text
D                 1   3   11   33
# units            1   2   10   20.
```

For an order word `(D_1,...,D_6)`, hereditary primitivity supplies

```text
lcm(D_i : i != j)=33                    for every j.          (1)
```

Because `33=3*11` is squarefree, (1) is equivalent to the pair of carrier
conditions

```text
at least two coordinates are divisible by 3,
at least two coordinates are divisible by 11.                (2)
```

Indeed, deleting any one coordinate preserves a prime factor precisely when
at least two coordinates carry it.  For one prime, its carrier set can be
any subset of six coordinates except the empty set and six singletons.  The
exact number of order words is consequently

```text
(2^6-1-6)^2 = 57^2 = 3249.                                (3)
```

The unit-weighted count factors independently across the two prime carriers.
For a prime `p`, its contribution is

```text
sum_(|S|>=2) binom(6,|S|)(p-1)^|S|
 = p^6-1-6(p-1).
```

Thus the exact number of literal unit words per support is

```text
(3^6-1-6*2)(11^6-1-6*10)
 = 716*1771500
 = 1268394000.                                             (4)
```

Across all `binom(12,6)=924` labelled supports, (3)--(4) give

```text
3002076 labelled support/order contexts,
1171996056000 raw labelled unit states.                       (5)
```

The primary checks the literal leave-one-out lcm predicate against (2) on all
`4^6` order words and checks both factorizations.

## 2. Literal masks and scalar capacities

For provider label `a`, owner label `b`, effective order `D`, and unit `u`,
let `B` be the unique CRT representative

```text
B = D*a (mod 13),                   B = u (mod D).             (6)
```

The owner-sheet mask is

```text
M(a,D,u;b)
 = {t in Z/33 : <B*(b^(-1)+13t)>_(13D) lies in (-D,D]}.       (7)
```

Writing `r=a/b` in `F_13^*`, direct interval counting gives the
unit-independent cardinality

```text
|M(a,D,u;b)|
 = (33/D) * #{x in (-D,D] : x = D*r (mod 13)}.                (8)
```

The complete scalar table, with columns `D=1,3,11,33`, is

```text
r= 1: 33 11 6 6          r= 7: 0  0 3 5
r= 2:  0  0 6 5          r= 8: 0 11 6 5
r= 3:  0  0 6 5          r= 9: 0 11 6 5
r= 4:  0 11 6 5          r=10: 0  0 6 5
r= 5:  0 11 6 5          r=11: 0  0 6 5
r= 6:  0  0 3 5          r=12: 0  0 3 5.                     (9)
```

Summing (9) over providers is a necessary owner-capacity test.  Over all
contexts in (5), the number of owners reaching scalar capacity `33` has
histogram

```text
0:16086, 1:446568, 2:1446252, 3:891684,
4:171594, 5:28548, 6:1344.                                  (10)
```

Exactly `1,344` rows on `150` supports survive at all six owners.  They have
eleven order-multiplicity profiles, `600` capacity vectors, and represent
`40,022,400` literal unit words.  The support histogram by survivor count is

```text
0:774, 2:48, 3:24, 4:24, 9:24, 18:12, 36:18.                (11)
```

Every context omitted from this bank is already impossible at some owner by
the scalar union bound.

## 3. The sound `Z/3` relaxation

Replacing `t` by `t+D` in (7) changes its CRT argument by a multiple of
`13D`.  Hence every order-`D` mask is pulled back from `Z/DZ`.  When `D|3`,
the mask is in particular periodic under `t -> t+3` and is an exact union of
the three thick fibres for `Z/33 -> Z/3`.

Fix an owner and retain every provider with `D in {1,3}`.  For one exact
anchor-unit tuple write

```text
Q = union_(D_i|3) M_i(u_i).
```

For arbitrary choices of all nonanchors,

```text
|union_i M_i(u_i)|
 <= |Q| + sum_(D_i in {11,33}) |M_i(u_i) minus Q|
 <= |Q| + sum_(D_i in {11,33}) max_e |M_i(e) minus Q|.        (12)
```

Therefore

```text
U_3 = max_Q [|Q| + sum_(D_i in {11,33})
                         max_e |M_i(e) minus Q|]              (13)
```

is a sound upper bound on every literal owner union.  The independent
nonanchor maxima forget their shared unit word and their mutual overlaps;
both losses only enlarge attainable coverage.  Thus `U_3<33` is terminal.

Across all `1344*6=8064` owner obligations, the exact anchor-union bank sizes
are exceptionally small:

```text
1:444, 2:3336, 3:4284.                                      (14)
```

The complete relaxed-bound histogram is

```text
25:996, 26:4248, 27:384, 28:288,
29:228, 30:96, 33:1824.                                     (15)
```

The number of owners per scalar row at which `U_3` still reaches 33 is

```text
0:306, 1:360, 2:624, 4:54.                                  (16)
```

In particular no row remains live at five or six owners.  The sharper bridge
between the hereditary and fibre layers is the owner-order refinement:

| own order | exact `U_3` distribution over owner obligations |
|---:|---|
| `1` | `33:168` |
| `3` | `25:216, 26:2376, 27:12, 28:192, 29:96, 30:96, 33:1656` |
| `11` | `25:624, 26:912, 27:228, 28:48, 29:84` |
| `33` | `25:156, 26:960, 27:144, 28:48, 29:48` |

Every own-order `11` or `33` obligation therefore has

```text
U_3 <= 29 < 33.                                             (17)
```

Condition (2) forces at least two coordinates to carry eleven, so their two
matching owners satisfy (17).  A global unit word would have to cover all 33
sheets at every owner.  Equations (10), (12), and (17) prove the stated
emptiness.

## 4. Complementary flags and the loss ledger

The primary also evaluates two nonessential competitors on the same complete
scalar bank.  The independent referee deliberately omits them and reconstructs
only the theorem-bearing `Z/3` carrier.

The transverse `Z/11` flag retains orders `1,11` exactly and relaxes orders
`3,33`.  Its live-owner histogram is

```text
0:30, 1:192, 2:912, 3:72, 4:72, 6:66.                       (18)
```

Thus it leaves 66 all-owner rows and cannot prove the theorem alone.  This is
the same complementary-prime phenomenon seen at scale thirty, but here the
smaller `Z/3` quotient is already terminal and no second pass is required.

A mixed anchor bank retains every proper order `1,3,11` exactly and relaxes
only order `33`.  It sharpens some numerical bounds, but its threshold live
histogram is again exactly (16).  Its largest exact anchor bank has `2,025`
unions, versus only three for the proof-facing `Z/3` carrier.  This identifies
the useful quotient by compression as well as by terminality: the mixed bank
adds information without improving the binary theorem predicate.

Multiplication by `F_13^*` gives scalar-context orbit sizes

```text
6:4, 12:110,
```

and support orbit sizes

```text
6:1, 12:12.
```

These are telemetry only.  No orbit quotient skips a labelled row.

## 5. Tournament Analysis and challenged vertices

Runner tournaments are not the proof object.  Use the six owner obligations
as vertices and compare the lexicographic cost key

```text
(1[U_q>=33], U_q, scalar capacity, anchor-bank size)          (19)
```

in either the `q=3` or `q=11` gauge; the harder obligation wins and support
coordinate order is the tie Hamiltonian path.  Every one of the 1,344
completed tournaments is transitive, with score word `(0,1,2,3,4,5)`, no
directed triangle, six singleton SCCs, and one Hamiltonian path.  Switching
from the `Z/3` to the `Z/11` gauge flips between zero and nine of the fifteen
edges.  The stable transitivity is therefore an ordering artefact, while the
absolute inequality (17) is the obstruction.

There is a second alternate-vertex tournament on the three quotient flags

```text
Z/3, Z/11, mixed proper-order anchors.
```

At each owner, the smaller sound upper bound wins; flag order resolves ties.
All 8,064 flag tournaments are again transitive with score word `(0,1,2)`,
no directed triangle, three singleton SCCs, and one Hamiltonian path.  The
`Z/3` flag defeats `Z/11` on every obligation.  The mixed flag defeats `Z/3`
only 192 times and never changes its threshold verdict.

This explicitly challenges several vertex assumptions:

- runner or provider vertices retain label ratios but lose the absolute
  owner-cover threshold;
- bare sheet vertices retain point incidence but obscure the hereditary
  prime-carrier bridge;
- quotient flags compare proof strength but lose which owner order forces the
  deficit;
- the faithful finite carrier is the owner-coloured `Z/3` anchor incidence,
  decorated by the own-order eleven-carrier bit and the numerical bound.

The relaxation destroys shared nonanchor units and all mutual needle
overlaps.  It remains sound precisely because (12) is an upper bound.  A
completed tournament destroys even the numerical threshold and is telemetry
only.

## 6. Exact replay, integrity, and scope

The deterministic Python primary:

1. checks every CRT base algebraically and by bounded literal search;
2. checks ratio covariance, mask cardinality (8), and order periodicity;
3. checks the hereditary grammar against all literal leave-one-out lcms;
4. visits all 3,002,076 labelled scalar contexts without quotienting;
5. evaluates all three exact anchor systems on every one of the 8,064
   scalar-survivor owner obligations; and
6. checks closure and sizes of every multiplication orbit and both tournament
   sidecars.

Normal and `python -O` executions reproduce the frozen 62-line primary output
byte-for-byte.  The primary SHA-256 values are

```text
Python primary source  356318648ffbef93702d7fbcd040c07df15ffcc8a36ee8ead7bc730fad46de65
primary output         893cf4d9b0ff65b10fcf13746ac91b45104e2933ba223065a32061c53b64334a
```

The internal theorem-bearing stream digests are

```text
grammar       2900299ec2c09d25e41e3281575ca5de5fa5077388a34b8694bc0eda40d067c3
CRT bases     dc9e41ee276778a8b43d911c2e823109050c771285998307131d8fd39c148859
literal masks bc85dd8ecfc3d5a178b7da29be292b203472d2a5d1dfbf1b3cbefda8f42e0554
scalar bank   a1d26583a73929cd8f3bd66563741e09b550100dfe414d8334d71989cef280a4
Z/3 rows      af6ad015c829d6b8dd6ffad44752297f51ac6f2535548f1b3dea3b7b3b900564
```

Replay with

```bash
python3 04-computation/lrc13_scale_thirty_three_hamming_six_frontier_scout_codex_c33.py |
  cmp - 05-knowledge/results/lrc13_scale_thirty_three_hamming_six_frontier_scout_codex_c33.out
python3 -O 04-computation/lrc13_scale_thirty_three_hamming_six_frontier_scout_codex_c33.py |
  cmp - 05-knowledge/results/lrc13_scale_thirty_three_hamming_six_frontier_scout_codex_c33.out
```

The independent standard-library C++20 referee constructs CRT bases only by
bounded literal search, verifies ratio covariance against every literal mask,
enumerates the `4^6` grammar in base-four order, and reproduces precisely the
small theorem-facing surface:

```text
hereditary words / scalar survivors / supports     3249 / 1344 / 150
Z/3 anchor-bank sizes                              1:444,2:3336,3:4284
own-order 11 bound maximum                         29
own-order 33 bound maximum                         29.
```

It does not import or reproduce the nonessential `Z/11`, mixed-anchor, orbit,
or tournament telemetry.  It shares no code, NumPy tables, digests, or
assertion stream with the primary.  Warning-clean `-O3` and `-O0` builds both
reproduce its frozen 14-line output byte-for-byte.  Its SHA-256 values are

```text
C++ referee source  b54507c6a7244ff8db7bffecab59024c934fa3d7c6e1e24a5dd5bffcd791a7f0
referee output      1a4c62f9736b52e561228a81951dea7fc84fc1c50acc5b3ff56c4fc326d8f664
```

Reproduce with

```bash
clang++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc13_scale_thirty_three_hamming_six_three_fibre_referee_codex_c33.cpp \
  -o /tmp/c33-referee
/tmp/c33-referee |
  cmp - 05-knowledge/results/lrc13_scale_thirty_three_hamming_six_three_fibre_referee_codex_c33.out
```

The two implementations therefore agree independently on every finite step
used by (17), while inequality (12) supplies the language-independent
soundness bridge.

This theorem closes only the primitive proper AP-centred common-scale-33
Hamming-six face.  It does not close the Hamming-five ramified bank,
non-AP/deep sheets, the global sporadic branch, or LRC(14).  Scale 34 is the
next untreated common-scale Hamming-six face on this line.  ∎
