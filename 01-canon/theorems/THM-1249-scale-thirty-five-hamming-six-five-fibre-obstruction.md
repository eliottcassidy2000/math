---
id: THM-1249
title: Scale-thirty-five Hamming-six five-fibre obstruction
status: PROVED STRUCTURAL + FINITE-EXACT — the complete 3,249-word/3,002,076-context scalar bank leaves 216 rows on 24 supports; the sound Z/5 thick-fibre relaxation is at most 31<35 at all 1,296 surviving owner obligations, so the primitive proper AP-centred common-scale-35 Hamming-six face is empty. Normal and optimized exact runs agree; the generic ceiling consumer is sorry-free Lean. This closes one AP-centred slice, not the uniform n=12 sporadic branch.
source: codex-2026-07-19-S78 sporadic-frontier audit
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-994, THM-1072, THM-1090, THM-1096, THM-1124, THM-1125]
related: [HYP-6820]
verification:
  - 04-computation/lrc13_scale_thirty_five_hamming_six_five_fibre_obstruction_codex_c35.py
  - 05-knowledge/results/lrc13_scale_thirty_five_hamming_six_five_fibre_obstruction_codex_c35.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCScaleThirtyFiveFibre.lean
---

# THM-1249 — scale thirty-five has a terminal five-fibre deficit

> **The primitive proper AP-centred common-scale-35 Hamming-six face is
> empty.**

The proof retains every effective-order-one and effective-order-five mask
exactly as a union of thick fibres of `Z/35Z -> Z/5Z`.  For each resulting
anchor union, it allows every order-seven and order-thirty-five provider to
choose its unit independently outside the anchor.  This is an upper
relaxation: it forgets shared unit compatibility and nonanchor intersections.
Its maximum is only `31<35` at every scalar-surviving owner obligation.

This is a finite face theorem inside the AP-centred shallow Hamming-six chart.
It is not AP extraction, shallow all-height rigidity, the deep content law, or
uniform emptiness of the twelve-runner sporadic branch.

## 1. Squarefree hereditary grammar

At common scale `35=5*7`, the effective orders and exact-unit counts are

```text
D                 1   5   7   35
# units            1   4   6   24.
```

For an order word `(D_1,...,D_6)`, hereditary primitivity requires

```text
lcm(D_i : i != j)=35                    for every j.          (1)
```

Because the scale is squarefree, (1) is equivalent to at least two
coordinates carrying the factor five and at least two carrying the factor
seven.  Each carrier subset is therefore any subset of six coordinates except
the empty set and the six singletons.  Hence

```text
# hereditary words = (2^6-1-6)^2 = 57^2 = 3,249.             (2)
```

The unit-weighted count per labelled support factors independently:

```text
(5^6-1-6*4)(7^6-1-6*6) = 1,834,747,200.                     (3)
```

Across the `binom(12,6)=924` supports this is

```text
3,002,076 labelled support/order contexts,
1,695,306,412,800 raw labelled unit states.                  (4)
```

The exact program checks the carrier condition against all six literal
leave-one-out lcm equations on all `4^6` words.

## 2. Literal owner masks and scalar capacity

For provider label `a`, owner label `b`, order `D`, and exact unit `u`, let
`B` be the unique CRT representative

```text
B = D*a (mod 13),                    B = u (mod D).            (5)
```

The owner-sheet mask is

```text
M(a,D,u;b)
 = {t in Z/35 : <B*(b^(-1)+13t)>_(13D) lies in (-D,D]}.       (6)
```

The script constructs (5) both algebraically and by literal bounded search,
then checks every bit of (6), its cardinality, and its order-periodicity.  If
`r=a/b` in `F_13^*`, the unit-independent cardinality vectors, in columns
`D=1,5,7,35`, are

```text
r= 1: 35  7 10  6          r= 7:  0  7  5  5
r= 2:  0  7  5  6          r= 8:  0  7  5  6
r= 3:  0  7  5  5          r= 9:  0  0  5  5
r= 4:  0  0  5  5          r=10:  0  7  5  5
r= 5:  0  7  5  6          r=11:  0  7  5  6
r= 6:  0  7  5  5          r=12:  0  0  5  5.                (7)
```

Summing (7) over the six providers gives a necessary scalar capacity at each
of the six owner labels.  Across all contexts in (4), the number of owners
reaching capacity 35 has histogram

```text
0:8,522, 1:1,135,668, 2:1,246,950, 3:485,440,
4:114,720, 5:10,560, 6:216.                                  (8)
```

Exactly 216 rows on 24 supports survive all six scalar tests; every such
support has nine rows.  The surviving order-multiplicity profiles
`(#1,#5,#7,#35)` are

```text
(0,2,0,4):24,       (0,2,1,3):96,       (0,2,2,2):96.        (9)
```

Thus every survivor has exactly two order-five anchors and no order-one
provider.  The rows have 130 distinct capacity vectors and represent
286,654,464 literal unit words.

## 3. The sound `Z/5` relaxation

Every order-`D` mask in (6) is periodic by `D`.  For `D in {1,5}` it is
therefore pulled back from `Z/5Z`.  Fix an owner and an exact anchor-unit tuple
for those providers, and put

```text
Q = union_(D_i in {1,5}) M_i(u_i).                            (10)
```

For arbitrary choices of the remaining units,

```text
|union_i M_i(u_i)|
 <= |Q| + sum_(D_i in {7,35}) max_e |M_i(e) minus Q|.         (11)
```

This is just the union bound after splitting off `Q`; independent maximization
can only increase the right side.  Maximizing (11) over the exact anchor bank
gives the complete histogram over the `216*6=1,296` surviving owner
obligations:

```text
U_5 = 27:576, 28:336, 29:240, 30:96, 31:48.                  (12)
```

The exact anchor bank has size four on 432 obligations and ten on 864.  In
particular,

```text
U_5 <= 31 < 35                                                (13)
```

at every obligation.  Split by the owner's own effective order, (12) is

```text
own 5:   27:288, 28:144
own 7:   29:144, 30:96, 31:48
own 35:  27:288, 28:192, 29:96.                              (14)
```

A literal cover would have to cover all 35 sheets at every owner, contradicting
(13).  Equations (1), (8), (11), and (13) prove the theorem.

The transverse `Z/7` relaxation is strictly weaker: it leaves 24 rows live at
all six owners and reaches 37.  Retaining all proper-order anchors reproduces
the `Z/5` bounds but enlarges the anchor bank to as many as 210 states.  Thus
the five-fibre quotient is both the stronger and the smaller proof carrier.

## 4. Tournament and assumption audit

The six owner obligations, not the runners, are the useful tournament
vertices.  In either the `Z/5` or `Z/7` gauge, compare the lexicographic key

```text
(1[U_q>=35], U_q, scalar capacity, anchor-bank size),         (15)
```

using coordinate order as the tie Hamiltonian path.  All 216 tournaments are
transitive, with score word `(0,1,2,3,4,5)`, zero directed triangles, six
singleton SCCs, and one Hamiltonian path.  Gauge changes flip up to six edges.
This telemetry does not prove (13); the absolute owner threshold does.

Challenging the vertex choice explains the loss precisely:

- runner vertices retain speed order but lose the owner threshold;
- provider vertices retain label ratios but lose the common sheet union;
- quotient-flag vertices compare proof cost but forget which owner fails;
- owner-coloured five-fibre incidence retains exactly the noncoverage
  predicate used in (11).

The relaxation deliberately destroys shared nonanchor units and mutual
needle overlaps.  That loss is sound only because (11) is an upper bound.

## 5. Exact replay, Lean consumer, and scope

Normal and `python -O` runs reproduce the stored output byte-for-byte.  Frozen
SHA-256 values are

```text
Python source   2c87857e3b296677d8809db2551a50714c52e5342ca6425f0a93a344db6f082a
Python output   75ba6e5a64e4ccd2b168ab894099a9549bad5d6c7f98eb46801caa475a7890fd
Lean consumer   833a00acde6707c126be64ee250bd64124685988ad436beb842d82d2b6814647
```

`LRCScaleThirtyFiveFibre.lean` imports the generic, kernel-pure
`LRCNestedFibreRelaxation` theorem and proves that any certified relaxation
ceiling at most 31 prevents equality with `Finset.univ : Finset (Fin 35)`.
The large row bank remains the external exact certificate; the Lean file has
no `sorry` and no `native_decide`.

This advances the AP-centred Hamming-six common-scale frontier from 35 to 36.
It does **not** close the finite smooth Hamming-five bank, non-AP-centred
shallow packets, two-sheet or higher deep packets, arbitrary-height ballot
rigidity, AP extraction, or uniform n=12 sporadic emptiness.  THM-860 still
bounds this particular AP-centred H6 scale line by `c<=840`; the next untreated
common scale on that line is `c=36`.  ∎
