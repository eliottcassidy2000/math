---
id: THM-1090
title: Scale-thirty Hamming-six two-flag obstruction
status: PROVED STRUCTURAL + MULTIPLY INDEPENDENT FINITE-EXACT — two separately developed Python primaries and two standard-library C++ referees agree on the 185,193-word hereditary grammar, all 171,118,332 labelled scalar contexts representing 588,280,492,800 raw states, the 54,050-row scalar bank, the 120-row Z/6 residual, and all 720 terminal Z/10 owner deficits. The referees additionally certify exact reachable masks or all 18,874,368 literal residual owner/unit words; all frozen outputs replay byte-for-byte across their recorded build matrices.
source: codex-2026-07-18-S67 scale-thirty continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-994, THM-1072]
related: [THM-983, HYP-6820]
verification:
  - 04-computation/lrc13_scale_thirty_hamming_six_two_flag_obstruction_codex_c30.py
  - 05-knowledge/results/lrc13_scale_thirty_hamming_six_two_flag_obstruction_codex_c30.out
  - 04-computation/lrc13_scale_thirty_hamming_six_two_flag_referee_codex_c30.cpp
  - 05-knowledge/results/lrc13_scale_thirty_hamming_six_two_flag_referee_codex_c30.out
  - 04-computation/lrc13_scale_thirty_hamming_six_complementary_fibre_obstruction_codex_c30.py
  - 05-knowledge/results/lrc13_scale_thirty_hamming_six_complementary_fibre_obstruction_codex_c30.out
  - 04-computation/lrc13_scale_thirty_complementary_fibre_referee_codex_c30.cpp
  - 05-knowledge/results/lrc13_scale_thirty_complementary_fibre_referee_codex_c30.out
---

# THM-1090 — scale thirty has a terminal complementary-fibre deficit

The theorem is:

> **The primitive proper AP-centred common-scale-30 Hamming-six face is
> empty.**

The proof-facing carrier changes once.  The first pass retains the masks whose
orders divide six, hence the thick fibres of

```text
Z/30 -> Z/6.
```

That pass kills all but 120 scalar rows.  Those rows have no divisor-of-six
anchor at all, so the transverse quotient

```text
Z/30 -> Z/10
```

is the correct second carrier.  Retaining all divisor-of-ten masks makes every
remaining owner obligation miss at least two sheets, even after every
nonanchor provider is allowed to choose its unit independently.  This is a
sound upper relaxation and is therefore terminal.

The Python primary contains multiple internal reconstructions, but those do
not count as independence.  The separate C++ referee was written without
reading the primary source or output, uses a different scalar scan and owner
gauge, and independently agrees on every theorem-bearing count.  Section 9
records the exact cross-check.

## 1. Squarefree hereditary grammar

At `c=30=2*3*5`, the effective orders and literal unit counts are

```text
D                 1  2  3  5  6  10  15  30
# units            1  1  2  4  2   4   8   8.
```

For an order word `(D_1,...,D_6)`, hereditary leave-one-out lcm means

```text
lcm(D_i : i != j) = 30                  for every j.       (1)
```

Since 30 is squarefree, (1) is equivalent to the three independent conditions

```text
at least two coordinates are divisible by 2,
at least two coordinates are divisible by 3,
at least two coordinates are divisible by 5.                (2)
```

Indeed deletion preserves the exponent one of a prime exactly when at least
two coordinates carry it.  Conversely (2) preserves all three prime factors
after every deletion.

For each prime, a carrier set can be any subset of the six coordinates except
the empty set and the six singletons.  Thus the exact unweighted count factors:

```text
(2^6-1-6)^3 = 57^3 = 185193.                               (3)
```

The unit-weighted count factors just as cleanly.  For a squarefree order,
`phi(D)=prod_(p|D)(p-1)`, so the contribution of one prime is

```text
sum_(|S|>=2) binom(6,|S|)(p-1)^|S|
 = p^6-1-6(p-1).
```

At `p=2,3,5` these factors are `(57,716,15600)`, giving

```text
636667200 literal state words/support,
924*636667200 = 588280492800 raw labelled states.             (4)
```

The primary checks (1) against (2) on all `8^6=262144` order words and checks
both factorizations (3)-(4).

## 2. Literal masks and the scalar table

For provider label `a`, owner label `b`, order `D`, and unit `u`, let `B` be
the CRT representative

```text
B = D*a (mod 13),             B = u (mod D).
```

The local mask is

```text
M(a,D,u;b)
 = {t in Z/30 : <B*(b^(-1)+13t)>_(13D) lies in (-D,D]}.       (5)
```

Writing `r=a/b in F_13^*`, literal interval counting gives the unit-independent
formula

```text
|M(a,D,u;b)|
 = (30/D) * #{x in (-D,D] : x = D*r (mod 13)}.               (6)
```

The resulting exact table is

```text
D=1 : 30 on {1}, otherwise 0;
D=2 : 15 on {1,6,7}, otherwise 0;
D=3 : 10 on {1,4,5,8,9}, otherwise 0;
D=5 :  0 on {4,9,12}, otherwise 6;
D=6 :  0 on {12}, otherwise 5;
D=10:  6 on {1,2,3,6,7,10,11}, otherwise 3;
D=15:  6 on {1,6,7}, otherwise 4;
D=30:  5 on {1,3,4,6,7,9,10}, otherwise 4.                  (7)
```

The primary constructs every CRT base algebraically and by bounded literal
search, constructs all `12*30*12=4320` label/state/owner masks from (5), and
checks every cardinality against both (6) and (7).

Summing (7) over the six providers is a necessary scalar capacity test.  The
complete

```text
924*185193 = 171118332
```

labelled support/order census has feasible-owner histogram

```text
0:1401966, 1:36143640, 2:66874158, 3:49478260,
4:15326622, 5:1839636, 6:54050.                              (8)
```

Thus exactly `54050` rows, on `772` supports, reach scalar capacity 30 at all
six owners.  They have 244 order-multiplicity profiles and 28965 distinct
capacity vectors.  Their unit fibres contain only `64678912` literal unit
words.  Every other labelled context is already impossible at some owner by
the union bound.

## 3. The sound flag relaxation

For every effective order `D`, replacing `t` by `t+D` in (5) changes the CRT
argument by a multiple of `13D`.  Hence an order-`D` mask is the pullback of a
subset of `Z/DZ`.  In particular, if `D|q|30`, the mask is periodic under
`t -> t+q` and is a valid anchor for `Z/30 -> Z/q`.

Fix an owner and a flag `q`.  Retain every provider with `D_i|q`.  For an
anchor-unit tuple let

```text
Q = union_(D_i|q) M_i(u_i).
```

For arbitrary choices of the remaining units,

```text
|union_i M_i(u_i)|
 <= |Q| + sum_(D_i not| q) |M_i(u_i) minus Q|
 <= |Q| + sum_(D_i not| q) max_e |M_i(e) minus Q|.             (9)
```

Therefore

```text
U_q = max_Q [|Q| + sum_(D_i not| q) max_e |M_i(e) minus Q|]   (10)
```

is a sound upper bound on every literal union.  The independent nonanchor
maxima in (9)-(10) forget their shared unit word and their mutual overlaps;
both losses can only enlarge attainable coverage.  Thus `U_q<30` is terminal.

## 4. The complete `Z/6` pass

On all 54050 scalar rows, retain exactly the orders

```text
D in {1,2,3,6}.
```

The 324300 owner obligations have saturated bound histogram

```text
16:144, 17:144, 18:384, 19:840,
20:7632, 21:9264, 22:26496, 23:60972,
24:104388, 25:31368, 26:42480, 27:22464,
28:6600, 29:300, 30:10824.                                  (11)
```

Their exact anchor-union banks have sizes only `1,2,3,4,5,7`, with histogram

```text
1:114588, 2:123024, 3:42996, 4:35760, 5:5244, 7:2688.        (12)
```

The number of owners per context at which the relaxed bound still reaches 30
is

```text
0:45110, 1:7536, 2:1284, 6:120.                              (13)
```

Consequently 53930 of the 54050 scalar rows already have at least one
impossible owner.  The 120 remaining rows have no divisor-of-six anchor at
all.  In order coordinates

```text
(n1,n2,n3,n5,n6,n10,n15,n30)
```

their exact multiplicity census is

```text
(0,0,0,0,0,2,3,1):24,
(0,0,0,0,0,4,0,2):48,
(0,0,0,1,0,3,0,2):48.                                      (14)
```

These rows represent only `3145728` literal unit words.

## 5. The terminal `Z/10` pass

On the complete 120-row residual, retain exactly the orders

```text
D in {1,2,5,10}.
```

The 720 owner obligations have anchor-union bank-size histogram

```text
10:144, 34:48, 52:96, 70:96, 86:48, 94:288.                 (15)
```

Their sound upper bounds are

```text
U_10: 23:48, 24:120, 25:192, 26:240, 27:72, 28:48.          (16)
```

Every entry in (16) is strictly below 30.  More structurally, the three
profiles in (14) have uniform sharp relaxation ceilings

```text
D10^2 D15^3 D30       : U_10 <= 28,
D10^4 D30^2           : U_10 <= 26,
D5 D10^3 D30^2        : U_10 <= 27.                          (17)
```

Thus all 120 residual contexts have zero surviving owners under the second
relaxation.  A global unit word would have to cover all thirty sheets at every
owner, while (16)-(17) make even each owner-local projection impossible.
Together with (8) and (13), this proves the emptiness implication.  The
independent referee reproduces all three gates.

The competing `Z/15` flag is weaker:

```text
U_15: 27:48, 28:144, 29:144, 30:384;
live owners/context: 0:24, 2:48, 6:48.                        (18)
```

On every residual owner obligation `U_10<=U_15`; the inequality is strict on
672 and an equality on 48.  The two complementary CRT pictures are

```text
Z/30 = Z/6 x Z/5,              Z/30 = Z/10 x Z/3.
```

The first sees the broad bank; the second detects the anchor-free residual.

## 6. Exact residual sharpness

An immutable-union DP separately retains every provider mask at all 720
residual owner obligations.  It constructs `3401088` reachable masks, with
largest bank `12936`, and finds

```text
exact maximum union:
23:72, 24:144, 25:192, 26:264, 27:48;

exact feasible owners/context:
0:120.                                                       (19)
```

The `Z/10` upper bound equals the exact maximum on 576 owner rows, exceeds it
by one on 48, and exceeds it by two on 96.  Every discrepancy remains below
threshold.  Hence the deliberately lossy second flag is threshold-exact on
the residual bank.

This DP is a sharpness sidecar, not an independent certificate.  It shares
the primary's CRT and mask construction.

## 7. Covariance and orbit audit

No orbit quotient is used in the scalar or flag scans.  As telemetry, the
primary nevertheless checks all twelve multipliers in `F_13^*`, including
ownerwise capacity and flag/exact-summary covariance.  The full scalar bank
has context orbit histogram

```text
4:2, 6:9, 12:4499
```

and support orbit histogram

```text
4:1, 6:2, 12:63.
```

The residual consists of ten context orbits of size twelve.  After normalizing
an owner to label one, there are only 60 owner keys.  Representatives of the
ten context orbits, written `label:order`, are

```text
1:5  2:10 4:10 5:30 9:30 11:10
1:5  2:10 4:10 8:30 9:30 11:10
1:5  2:10 4:30 5:30 9:10 11:10
1:5  2:10 4:30 8:30 9:10 11:10
1:10 2:10 4:10 5:30 9:30 11:10
1:10 2:10 4:10 8:30 9:30 11:10
1:10 2:10 4:30 5:30 9:10 11:10
1:10 2:10 4:30 8:30 9:10 11:10
1:10 2:15 4:30 6:15 7:15 12:10
1:10 2:15 6:15 7:15 9:30 12:10.                     (20)
```

Equation (20) is a compact audit view, not the proof enumeration.

## 8. Cayley and tournament telemetry

The scalar high-cardinality switches already distinguish the relevant label
geometry.  Directing `a->b` when `b/a` is an off-diagonal high ratio gives

```text
       arcs  symmetric edges  reciprocal  triangles  SCCs
D=2      24       24              0          0       one
D=3      48       36             12         16       one
D=10     72       48             24         64       one
D=30     72       48             24         64       one.       (21)
```

These are Cayley relations, not tournaments; missing and reciprocal pairs are
part of the carrier.

For an owner-obligation tournament, compare the exact key

```text
(1[U_10>=30], U_10, exact maximum, scalar capacity,
 U_15, exact-bank size),                                      (22)
```

let the harder key win, and use support-coordinate order as the tie
Hamiltonian path.  All 120 completed tournaments are transitive: score word
`0,1,2,3,4,5`, no directed cycle, six singleton SCCs, and one Hamiltonian
path.  Their tie histogram is `1:96,2:24`.

A second tournament uses quotient flags `{6,10,15}` as vertices.  The smaller
absolute upper bound wins, with flag order as the tie path.  All 720 ownerwise
flag tournaments are again transitive with score word `0,1,2`; `U_10` defeats
`U_6` everywhere and defeats or ties `U_15` everywhere.

Neither tournament proves the theorem: orientation discards the absolute
threshold 30.  The proof-faithful vertices are owner obligations or quotient
fibres carrying their numerical bounds.  Runners/providers lose fibre
overlap; isolated residues and wall events lose mask incidence; Fourier modes
lose monotone union; a completed tournament loses both missing edges and the
terminal threshold.  This challenges the assumption that tournament vertices
must be runners or arcs: here quotient flags are useful diagnostic vertices,
while owner-coloured coset fibres are the actual proof object.

## 9. Independent replay and reproducibility

The deterministic primary:

1. reconstructs every CRT base algebraically and by literal search;
2. checks the squarefree grammar against all leave-one-out lcms;
3. checks every mask against (6)-(7) and proves all flag periodicities;
4. scans all 171118332 labelled scalar contexts without quotienting;
5. evaluates `U_6` on all 324300 scalar-survivor owner obligations;
6. evaluates `U_10`, competing `U_15`, and exact DP on all 720 residual owner
   obligations; and
7. checks every multiplication covariance, orbit, Cayley, and tournament
   fingerprint.

Frozen artifact hashes are

```text
Python primary source  c1c25e32b3c762287e787b31ceb9952c4a21877a0cd9eec84e57745fe4980901
primary output         10729d3d29f377ef44c9b8e04c439772ff29f405bc9533d7078fbed641f8e6e4
```

Internal stream digests are

```text
grammar          d5ec0f6ac4c8f74ed399fbce0e8b50f232789f82619f556f7606b2e41fb6d55f
CRT bases        c94c62de38c0238f686e34cfc77293359e8638418485e59e49bed509f992dc3e
literal masks    dcff92077f5e320927a98c7f5a86df40e9c6baed2075e7cad46542d40d3e8da9
scalar bank      6be8cfd5291e955df7cffd5112471bc0d2e232e816cc68d25de3031034e4dd93
U6 rows          bf5be40ba8f998c194319593abbddfbe277a15df3805c2ebff0f24d3b81e3342
residual rows    e8c698b8d62f804c67925977c4dbc397b26af00a5908b553fb5fb4f0eb10d59d
U10/U15 flags    c1607cc6959e5bd1aea3bbcf63353dbee6eadce5b30709ee98864c108ca3ea32
exact banks      955669e80f11fc3239ea076416966fc4f27ad6fb5e9dc9453f1e692979c959d8
```

The independent C++ referee uses only the standard library and was developed
without opening the Python source or output.  Its construction is deliberately
different:

1. it builds CRT representatives only by bounded literal search;
2. it proves the ratio/rotation gauge against every literal
   `(label,owner,D,unit)` mask, then evaluates normalized owner keys;
3. it enumerates the base-eight grammar and separately checks its weighted and
   unweighted inclusion-exclusion counts;
4. it expands scalar capacities in packed six-byte tables rather than NumPy;
5. it reconstructs every anchor from its `Z/6` or `Z/10` fibre signature; and
6. it computes each exact residual bank in both provider orders and requires
   the sorted forward and reverse banks to agree.

It independently obtains

```text
hereditary words / weighted states       185193 / 636667200
labelled contexts / raw states           171118332 / 588280492800
scalar survivors / supports              54050 / 772
U6 live owners per context               0:45110 1:7536 2:1284 6:120
all-owner U6 residual                    120 rows
U10 terminal bounds                      23:48 24:120 25:192
                                         26:240 27:72 28:48
exact maxima                             23:72 24:144 25:192
                                         26:264 27:48
U10 minus exact                          0:576 1:48 2:96
exact reachable masks / largest bank     3401088 / 12936.
```

The primary reports the `U6` histogram after saturation at 30.  The referee
retains the raw tail

```text
30:10608, 31:72, 32:96, 34:48;
```

and `10608+72+96+48=10824`, exactly the primary's saturated `30:10824`.
This is a representation difference, not a count discrepancy.

The independent frozen hashes are

```text
C++ referee source  8ca947d10371153bf69f13faf5c5eed331f465817efd8610d98e84f4cde1942d
referee output      af3ddc65b9ba927062e65f6a08880e99380590b2b23b4e521c2fb3f20b34263f
```

Apple clang 17 builds at `-O3`, `-O0`, ASan+UBSan, and libc++ debug hardening
all produce this byte-identical output; the clang static analyzer emits no
diagnostics.  `/usr/bin/g++` is the same Apple clang on this host, so no claim
of a genuine GCC replay is made.

Reproduction uses

```bash
python3 04-computation/lrc13_scale_thirty_hamming_six_two_flag_obstruction_codex_c30.py
python3 -O 04-computation/lrc13_scale_thirty_hamming_six_two_flag_obstruction_codex_c30.py
clang++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc13_scale_thirty_hamming_six_two_flag_referee_codex_c30.cpp \
  -o /tmp/c30-referee
/tmp/c30-referee
```

Both Python modes and the independently compiled referee have been replayed
against their stored outputs byte-for-byte.

## 10. A second independent certificate pair

A separately written NumPy/literal-CRT primary and standard-library C++
referee independently reconstruct the same grammar, scalar census, `U6`
live-owner histogram

```text
0:45110, 1:7536, 2:1284, 6:120,
```

and terminal `U10` histogram

```text
23:48, 24:120, 25:192, 26:240, 27:72, 28:48.
```

That referee goes beyond the relaxed bank: it checks all `18,874,368`
literal residual unit words on the 720 owner obligations pointwise, finds
zero covers, and obtains exact literal maxima
`23:72,24:144,25:192,26:264,27:48`.  Its optimized, unoptimized, and
ASan+UBSan outputs are byte-identical.  Frozen hashes are

```text
second Python source  42acbec7a6d5131d4ad34f7d51e4fc1b41d6f87c6c33f6ae4ccf569ff31ad253
second Python output  942ae6038530378a32d43ff2fa3bd91bb73bf4720a2a32e4175aa81fd064fc6d
second C++ source     16d272548d305ba376f7ef8780f0056e1058753b3651aa0023c8b2aed7e5c83e
second C++ output     83b49a729b595652fba0edb98adba20f176f9a439a36903eeefeeb31caa7a33b
```

This theorem concerns only the AP-centred Hamming-six common-scale-30 face.  It
does not close Hamming five, non-AP/deep sheets, the global sporadic branch, or
the full LRC(14) theorem.  Scale 31 is already prime-excluded by THM-983, so
THM-1096 now closes scale 32 separately; the next untreated composite on this
particular common-scale H6 line is 33.
