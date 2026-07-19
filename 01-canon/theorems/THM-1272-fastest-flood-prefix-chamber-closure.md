---
id: THM-1272
title: FASTEST-FLOOD PREFIX-CHAMBER CLOSURE -- the top-dominated branch has d6/c<798
status: PROVED / COMPUTER-EXACT (exact Pbar-by-tooth-count suprema; exhaustive THM-1233 prefix-chamber dynamic program; five basic e-branch cuts; functional private-count closure of all e=5 chambers from d6/c=798; top-dominated projective bound d6/c<798; optimization-safe referee; sorry-free Lean arithmetic consumer).  This closes a projective tail only under d6>=(7/2)d5; it does not prove six-comb noncoverage or LRC(14)
source: codex-2026-07-19 tail-tax branch closure audit
depends_on: [THM-1198, THM-1233, THM-1267, THM-1275]
related: [THM-1199, THM-1241, THM-1266]
script: 04-computation/lrc14_fastest_flood_prefix_chamber_closure_thm1272.py
output: 05-knowledge/results/lrc14_fastest_flood_prefix_chamber_closure_thm1272.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCFastestFloodPrefixChamberClosure.lean
script_sha256: 691edce90f5071f1ecb199e8129cfa829ed97f38914ad9d4715070542a9bfff3
output_sha256: ecb5d0cf24e746f1101b947df892689959ca527a35797a695293a221cc0f98e2
formalization_sha256: df2e9c543c3e4da2394db00f7135c32393ef773d19d66477a12c652b291b5d55
---

# THM-1272 -- fastest-flood prefix-chamber closure

## 1. Statement

Let a complete `c`-safe gap be covered by the strict danger combs of

```text
c<d1<d2<d3<d4<d5<h,                       h=d6.       (1)
```

Put

```text
x_i=d_i/c,                 x=h/c,
e=#{j<=5:h>6d_j}.                                      (2)
```

Assume the top-dominated branch

```text
h>=(7/2)d5.                                           (3)
```

THM-1275's fastest-owner flood/turn theorem gives, for the number `K` of
selected fastest teeth,

```text
F6>=1/(8x)[ceil(K/(e+1))-1],
K>=ceil(x/7),                                         (4)
K>=ceil(36x eta_h/7),                                 (5)

F6=sum_(i=1)^6 Pbar(6d_i/(7c))-1,
eta_h=1-sum_(i=1)^5 Pbar(6d_i/(7c)).                  (6)
```

Only the dominated corollary (4)--(5) is used below; the result is therefore
independent of whether THM-1275 presents its stronger invoice as disjoint
flood/selected-turn regions or as the layered whole-seam multiplicity law.

The basic count in (4), coupled to THM-1233's prefix chambers, gives the exact
certified branch cuts

```text
e=0  => 45h<563c,
e=1  => h<=42c,
e=2  => 20h<1407c,
e=3  => 5h<1078c,
e=4  => 2h<1323c,
e=5  => h<959c.                                      (7)
```

The functional private count (5) sharpens the last row to

```text
e=5  => h<798c.                                      (8)
```

Consequently every cover in the top-dominated branch satisfies

> **Top-dominated projective closure.**
>
> ```text
> d6/c<798.                                          (9)
> ```

This improves THM-1233's universal `d6/c<2345` by almost a factor three on
the branch (3).  It does not treat `d6<(7/2)d5`.

## 2. Tooth-count chambers retain the missing one-body information

Write

```text
T(y)=ceil((6y+1)/7).                                 (10)
```

For an integer `m>=2`, the chamber `T(y)=m` is exactly

```text
(7m-8)/6<y<=(7m-1)/6.                               (11)
```

Let

```text
q_m=sup{Pbar(6y/7):T(y)=m}.                          (12)
```

Here and below `Pbar(6y/7)` is abbreviated to `Pbar(y)` in the ratio
coordinate.  Intersecting (11) with THM-1198's twelve exact pieces and its
BV ray gives

```text
q_2=7/36,
q_3=8/45,
q_4=4/21,
q_m=1/7+1/(7m-8),                    m>=5.           (13)
```

The anomalously large `q_4` is essential.  The exact compact envelope has
value `1/7` at `y=7/2`, but THM-1198's BV majorant on the open ray `y>7/2`
has right limit `4/21`.  Dropping that upward jump would manufacture a false
gain.

Every putative cover has `x_1>1`.  THM-1267 gives

```text
x_1<563/270<13/6,                                   (14)
```

so `T(x_1)=2`, and its actual load is **strictly** below `q_2`: the value
`7/36` occurs only at the excluded endpoint `x_1=1`.  Therefore, for every
prefix chamber word,

```text
sum_(i=1)^r Pbar(x_i)<sum_(i=1)^r q_(m_i),
m_i=T(x_i).                                         (15)
```

THM-1233's individual ratio caps give

```text
m_1<=2, m_2<=9, m_3<=15, m_4<=76, m_5<=232.         (16)
```

The counts are nondecreasing.  They also retain the prefix-component
geometry.  If `a_m=(7m-8)/6`, every feasible word must satisfy

```text
a_(m2)<(91/29)(1+m1),
a_(m3)<(91/22)(1+m1+m2),
a_(m4)<(49/15)(1+m1+m2+m3),
a_(m5)<(21/8)(1+m1+m2+m3+m4).                       (17)
```

These are necessary conditions, not a relaxation in the wrong direction:
`x_i>a_(m_i)`, while THM-1233 bounds `x_i` strictly above by the corresponding
right side.

For depth `r` and count floor `S`, define the finite exact maximum

```text
Q_r(S)=max sum_(i=1)^r q_(m_i),                      (18)
```

over monotone words obeying (16)--(17) and `sum m_i>=S`.  The exact dynamic
program keeps `(sum m_i,last m_i)`; those are precisely the coordinates on
which the next feasibility test depends.  Its compressed state counts are

```text
1, 8, 84, 872, 9040                                 (19)
```

at depths one through five.

## 3. Converting `e` into a prefix count floor

When `e<5`, definition (2) gives

```text
h<=6d_(e+1).                                         (20)
```

Combining (20) with the depth-`e` row of THM-1233 yields

```text
e=1: x<(546/29)(1+m1),
e=2: x<(273/11)(1+m1+m2),
e=3: x<(98/5)(1+m1+m2+m3),
e=4: x<(63/4)(1+m1+m2+m3+m4).                       (21)
```

For `e=5`, THM-1233's last-tooth row is

```text
x<7(1+m1+...+m5).                                   (22)
```

Thus a lower bound on `x` forces an exact integer lower bound on the prefix
count sum.  The important chamber jumps are

```text
x>=1078/5  and e=3 => sum_(i<=3)m_i>=11,
x>=1323/2  and e=4 => sum_(i<=4)m_i>=42,
x>=959     and e=5 => sum_(i<=5)m_i>=137.            (23)
```

For example, `(1078/5)/(98/5)=11`, so (21) says the integer count sum is
strictly greater than ten.  The strict signs in (21)--(22) are exactly what
turn equality in (23) into the next integer chamber.

The exact maxima needed at the jumps are

```text
Q_1(2)   =7/36                    at (2),
Q_2(4)   =7/18                    at (2,2),
Q_3(11)  =61/108                  at (2,4,5),
Q_4(42)  =112571/158004           at (2,4,5,31),
Q_5(137) =9882995/11534292        at (2,4,5,31,95).  (24)
```

The immediately preceding `e=3` and `e=5` chambers are respectively

```text
Q_3(10)=145/252                 at (2,4,4),
Q_5(136)=476431/549252          at (2,4,4,31,95).    (25)
```

The load drops in (24)--(25) are not a rounding convenience.  They are the
discrete prefix geometry that a naked continuous ratio box misses.

## 4. The basic five branch contradictions

Every threshold in (7), except the already smaller `e=0` row, is above
`x=21`.  If `j>e`, then (20) gives `x_j>=x/6>7/2`; hence THM-1198's decreasing
BV ray gives

```text
Pbar(x_j)<=1/7+1/x.                                  (26)
```

The fastest load is `1/7+1/(6x)`.  If `Q` is the appropriate prefix maximum,
equations (15) and (26) give the strict upper bound

```text
F6<Q-(e+1)/7+(31-6e)/(6x).                          (27)
```

Nested ceilings in (4) simplify exactly:

```text
ceil(ceil(x/7)/(e+1))-1
 =ceil(x/[7(e+1)])-1=:B_e(x).                        (28)
```

Thus a cover requires

```text
B_e(x)/(8x)<=F6.                                     (29)
```

At the five boundary banks, the cross-multiplied margin

```text
B/8-{x_0[Q-(e+1)/7]+(31-6e)/6}                      (30)
```

is exactly

| `e` | claimed cut | count floor | `Q` | `B` | margin (30) |
|---:|---:|---:|---:|---:|---:|
| 1 | `x>42` | 2 | `7/36` | 3 just right of 42 | `1/24` |
| 2 | `x>=1407/20` | 4 | `7/18` | 3 | `0` |
| 3 | `x>=1078/5` | 11 | `61/108` | 7 | `29/216` |
| 4 | `x>=1323/2` | 42 | `112571/158004` | 18 | `11503/5016` |
| 5 | `x>=959` | 137 | `9882995/11534292` | 22 | `1185455/411939` |

The zero in the `e=2` row still contradicts a cover because (27) is strict.
At `e=1,x=42`, the basic tax has only `B=2`; this is why the conclusion is
the strict `x<=42`, not `x<42`.

For each row, `Q-(e+1)/7<0`.  Beyond the displayed cut, the count floor can
only rise, `Q_e` can only fall, and `B_e` can only rise.  Between jumps the
cross-multiplied left side has negative slope.  Hence checking (30) at the
boundary proves the whole corresponding ray, establishing the five nonzero
rows of (7).

For completeness, `e=0` means `h<=6d1`.  THM-1267's integer form gives

```text
45h<=270d1<=563c-1<563c,                             (31)
```

which is the first row of (7) with no envelope relaxation.

## 5. Functional private mass moves the last cut from 959 to 798

Assume `e=5`.  Put

```text
s=floor(x/7).                                        (32)
```

THM-1233 gives `x<2345`, so `x>=798` implies

```text
114<=s<=334.                                         (33)
```

On the ratio chamber

```text
7s<=x<7(s+1),                                        (34)
```

equation (22) forces `sum_i m_i>=s`.  Write `Q_s=Q_5(s)`.  The strict prefix
load bound (15) gives

```text
eta_h>1-Q_s.                                         (35)
```

Using `x>=7s` in THM-1275's functional private count (5), define

```text
k_s=floor(36s(1-Q_s))+1,
b_s=ceil(k_s/6)-1.                                   (36)
```

Then every cover in the chamber has `K>=k_s`, so (4) gives

```text
F6>=b_s/(8x).                                        (37)
```

On the other hand,

```text
F6<Q_s-6/7+1/(6x).                                  (38)
```

After multiplication by `x`, the worst endpoint of (38) is `7(s+1)` when
`Q_s>6/7`, and `7s` when `Q_s<=6/7`.  The exact `221`-chamber census
`s=114,...,334` proves

```text
max_endpoint x(Q_s-6/7)+1/6 < b_s/8                 (39)
```

in every row.  The smallest margin occurs already at `s=114`, where

```text
Q_114=2847097/3276540       at (2,4,4,25,79),
k_114=538,
b_114=89,
margin=1921973/1310616>0.                            (40)
```

The preceding chamber is an exact guardrail.  At `s=113`, the maximizing
word is `(2,2,4,26,79)` and the same certificate has negative margin

```text
-113823/63220.                                       (41)
```

Thus `798=7*114` is the genuine onset of this functional prefix-chamber
certificate, not a decimal rounding of the basic `959` cut.  Equations
(37)--(40) prove (8).  Combining `e=0,...,5` now proves (9).

## 6. Kakeya, tournament, and information-loss audit

The protected slow gap remains the Kakeya needle, but the useful vertices in
this proof are **prefix tooth-count chambers**.  Give two such obligations the
observable `q(m_i)-q(m_j)`, switch by its sign, and break ties first by count
and then by chronological speed order.  The resulting tournament is
transitive: its score histogram is `0,1,...,r-1`, it has no directed cycles,
all SCCs are singletons, and the tie-resolved Hamiltonian path is unique.

This quotient preserves the prefix component feasibility (17), the one-body
load price (13), and the integer jumps (23).  It destroys the five phases,
fastest-tooth addresses, locations of floods and turns, and every seam clock.
THM-1275 supplies those destroyed geometric facts before the quotient; this
theorem consumes only their scalar dominated tax.

We challenged runners, gaps, fixed circle sections, boundaries, teeth,
wall-crossing events, residues, cover arcs, Fourier modes, return packets,
prefix components, and proof obligations as vertices.  The challenged
assumption is that a runner or arc tournament must be the faithful final
carrier.  Here the decisive object is a monotone word of **load-bearing proof
obligations**, whose count jumps expose structure invisible in real ratios.

## 7. Verification, formalization, and scope

The dependency-free referee reconstructs all `231` chamber suprema
`m=2,...,232`, including the `m=4` BV jump.  Its exact dynamic program checks
compressed state banks of sizes (19), the maximizers (24)--(25), every
strict/equality ceiling in (20)--(31), all five basic margins, and all `221`
functional `e=5` chambers including the negative `s=113` guardrail.  It has
zero Python `assert` nodes, and ordinary and optimized outputs are
byte-identical.

The sorry-free Lean module kernel-checks the chamber endpoint algebra, the
`e=0` integer chain, the strict count jumps, the nested-ceiling arithmetic
consumers, all five basic boundary margins, the functional `s=113/114`
guardrail, and the abstract contradiction/monotonicity steps.  The exact
Pbar envelope, finite dynamic-program exhaustiveness, and THM-1275's
multiplicity theorem remain explicit paper/referee providers.

Frozen hashes are

```text
source         691edce90f5071f1ecb199e8129cfa829ed97f38914ad9d4715070542a9bfff3
output         ecb5d0cf24e746f1101b947df892689959ca527a35797a695293a221cc0f98e2
formalization  df2e9c543c3e4da2394db00f7135c32393ef773d19d66477a12c652b291b5d55
```

THM-1272 closes the formerly enormous projective tail of the
`d6>=(7/2)d5` branch and identifies `798` as the first chamber where the
current functional certificate becomes uniform.  It does **not** exclude the
remaining compact band `d6/c<798`, treat `d6<(7/2)d5`, prove six-comb
noncoverage, empty the sporadic branch, or prove LRC(14).  The highest-leverage
next use is to attack the complementary near-top cluster, where `d5>(2/7)d6`
and the THM-1241 macroscopic cut must occur inside the first five speeds.
