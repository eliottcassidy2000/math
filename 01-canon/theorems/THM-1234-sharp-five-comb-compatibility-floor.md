---
id: THM-1234
title: The sharp compatibility-sensitive five-comb pair floor
status: PROVED by an analytic all-range reduction and an exact-integer projective spanning-tree bank.  Every five distinct radius-1/14 danger combs have total pair overlap at least 44/273, and equality is attained.  Consequently R_6>=22/91 and R_7>=22/65, with sharpened slow-gap gcd, common-period, quadratic-hole, and Fano-period consequences.  This does not close the phase-coupled six-comb cover or LRC(14)
source: codex-2026-07-18 sharp five-speed compatibility session
depends_on: [THM-1166, THM-1179]
related: [THM-1191, THM-1196, THM-1199, THM-1205, THM-1231]
script: 04-computation/lrc14_five_comb_sharp_compatibility_floor_thm1234.py; 04-computation/lrc14_five_comb_sharp_compatibility_floor_thm1234_bank.cpp
output: 05-knowledge/results/lrc14_five_comb_sharp_compatibility_floor_thm1234.out
---

# THM-1234 -- the sharp compatibility-sensitive five-comb pair floor

For a positive integer `s`, put

```text
D_s={t in R/Z: ||st||<1/14},
rho(a,b)=measure(D_a intersect D_b),
R(S)=sum_{{a,b} subset S} rho(a,b).                    (1)
```

THM-1231 proved `R(S)>=13/81` for five distinct speeds and found the
telemetry row `(1,12,13,27,156)` at `44/273`.  The rounding and heavy-ratio
threshold in that certificate intentionally stopped at `13/81`.  Enlarging
the exact ratio alphabet, retaining the rows lost to fixed-grid rounding,
and resolving those rows rationally proves that the telemetry value was the
sharp constant after all.

> **Theorem (sharp five-comb floor).**  Every set `S` of five distinct
> positive integer speeds satisfies
>
> ```text
> R(S)>=44/273.                                         (2)
> ```
>
> Equality is attained, for example, by
>
> ```text
> S=(1,12,13,27,156).                                  (3)
> ```

Thus `44/273` is the optimal universal five-speed constant.

## 1. The sharp heavy-graph threshold

THM-1166 supplies the universal pair and triple bounds

```text
rho(a,b)>=1/91,
rho12+rho13+rho23>=1/24.                               (4)
```

Summing the triple inequality over the four triples of a four-speed set
counts each pair twice and gives

```text
R_4>=1/12.                                             (5)
```

Put

```text
T=44/273,
c=(T-1/12)/4=85/4368,
theta=1/49-c=29/30576,                                (6)
```

and call a pair **heavy** when `rho(a,b)<c`.  A nonheavy pair contributes at
least `c`.  If the heavy graph on the five labels is disconnected, its
component-size partition gives the following ledger.  Inside components of
sizes `1,2,3,4` use respectively `0,1/91,1/24,1/12`; every cross-component
edge is nonheavy.

```text
components       cross edges       lower bound          excess over T
4+1                   4               44/273                  0
3+2                   6              185/1092                3/364
3+1+1                 7                37/208               73/4368
2+2+1                 8                97/546                3/182
2+1+1+1               9               271/1456             109/4368
1+1+1+1+1            10               425/2184              73/2184. (7)
```

The first row is the exact threshold identity

```text
1/12+4(85/4368)=44/273.                               (8)
```

Therefore every hypothetical row with `R<T` has a **connected** heavy
graph and hence a heavy spanning tree.  Notice that this argument is for a
strict counterexample.  It does not classify equality rows whose heavy graph
is disconnected; no such classification is needed for (2).

## 2. The finite all-range ratio alphabet

For a reduced coprime pair `1<=u<v`, write `F(r)=r(14-r)` for the
representative `0<=r<14`.  The exact folded formula and its defect bound are

```text
rho(u,v)
 =[4uv+F((u+v) mod 14)-F((v-u) mod 14)]/(196uv),
rho(u,v)>=1/49-1/(7v).                                (9)
```

If the pair is heavy, (6) and (9) imply

```text
theta<1/(7v),
v<1/(7theta)=4368/29,
v<=150.                                               (10)
```

Exact enumeration of all coprime pairs in (10) leaves exactly

```text
72 heavy unoriented ratio types.                      (11)
```

Their actual largest upper coordinate is `125`.  The ordered `u:v`
serialization has SHA-256

```text
785e3b41577df477c6c04a284efaeb99a319382cea130ae128dc4d6ce87520c3. (12)
```

The referee independently compares (9) with exact interval geometry on all
`6,857` coprime rows `1<=u<v<=150`.  Thus (10), not an experimental speed
ceiling, makes the bank exhaustive.

## 3. Exact projective spanning-tree enumeration

Take a heavy spanning tree, root it, and divide all speeds by the root
speed.  Beginning at the singleton coordinate `1`, attach each new vertex to
an existing coordinate `x` by

```text
x -> x(v/u) or x(u/v),                                (13)
```

where `u:v` is one of the 72 types.  After each attachment, clear
denominators, divide the full row by its gcd, and sort.  Induction along a
rooted tree enumerates the projective class of every connected-heavy
five-tuple.  Common scaling is harmless because (9) is dilation invariant.

The pruning calculation contains no floating point.  At scale `Q=10^9`, set

```text
L_Q(A)=Q^(-1) sum_{{a,b} subset A} floor(Q rho(a,b)). (14)
```

For a partial row of size `n`, each of the
`10-binomial(n,2)` unformed pairs is at least `1/91`.  Hence a row is
discarded only if the rigorous lower bound

```text
L_Q(A)+[10-binomial(n,2)]/91>=44/273                  (15)
```

already proves that every completion is safe.  Every comparison in (15) is
unsigned-integer cross multiplication.

The surviving projective banks are

```text
n=2:             72 rows       (144 attempts, 0 pruned),
n=3:         10,048 rows    (20,592 attempts, 180 pruned),
n=4:      1,189,651 rows (4,300,386 attempts, 1,524,044 pruned). (16)
```

Attaching the fifth vertex checks

```text
678,030,510 terminal extensions.                      (17)
```

Of these, `678,030,462` are certified directly by (14)--(15).  There are
only `48` occurrences, or four unique projective rows, for which fixed-grid
lower rounding lies below the target:

```text
(1,12,13,27,156),
(1,12,13,156,351),
(4,9,108,117,1404),
(9,52,108,117,1404).                                  (18)
```

Every row in (18) has the exact ten-pair word

```text
2*(1/63+17/819+1/84+1/91+23/1092)=44/273.            (19)
```

The least fixed-grid numerator is `161172154/10^9`; its apparent deficit is
exactly

```text
44/273-161172154/10^9=1958/(273*10^9).                (20)
```

Equation (19) shows that this is only the sum of the ten downward-rounding
errors, not a counterexample.  Thus every terminal row is at least `T`,
contradicting `R<T` and proving (2).  Row (3) and (19) prove sharpness.

The integer-width audit is also explicit.  Since the actual alphabet has
upper coordinate at most `125`, four attachments give coordinates at most
`125^4`; a reduced pair denominator is bounded by

```text
196*125^8=11682510375976562500<2^64.                  (21)
```

The implementation additionally checks coordinate and denominator overflow
before multiplication and uses `unsigned __int128` for scaled products.

## 4. Six-speed and seven-speed consequences

For `m>=5`, sum (2) over all five-subsets of an `m`-speed set.  Each pair is
counted `binomial(m-2,3)` times, so

```text
R_m>=T*binomial(m,5)/binomial(m-2,3)
   =11m(m-1)/1365.                                    (22)
```

In particular,

```text
R_6>=22/91,                 R_7>=22/65.               (23)
```

For six faster combs covering an `a`-slow gap, retain THM-1179's notation

```text
delta=a sum_i 1/b_i-1,
g_ij=gcd(b_i,b_j),
eta_ij=rho_ij(1-rho_ij).                              (24)
```

Its complete-graph debt inequality and (23) give

```text
a sum_(i<j) eta_ij/g_ij +(18/49)delta
 >=(6/7)R_6>=132/637.                                 (25)
```

Since `eta_ij<=13/196`, this implies the sharper reciprocal-gcd law

```text
a sum_(i<j)1/g_ij>=528/169-(72/13)delta.              (26)
```

The optimal six-comb quadratic majorant gives

```text
measure(union_i D_bi)
 <=6/7-(1/3)R_6<=212/273.                             (27)
```

If all six killers have common gcd `G_0`, the periodic-component argument
therefore sharpens to

```text
G_0/a< (7/6)(212/273)=106/117.                        (28)
```

The inequality is strict for an actual open-danger cover of the closed slow
gap.

For seven combs, THM-1166's quadratic certificate yields

```text
measure(uncovered)>=2R_7/7>=44/455.                   (29)
```

On its common-dilate protected-needle branch this gives

```text
G/m<=7-2R_7<=411/65.                                  (30)
```

On the covered lower-LRC/protected-needle branch, the same scalar improvement
propagates through the metric-Fano clock.  From THM-1166,

```text
(2LR_7)/7<=sum_ell e_ell/G_ell,
L>=1/(7m),                  e_ell<=11/128.            (31)
```

Consequently every labelled Fano plane obeys

```text
sum_ell m/G_ell>=256R_7/539>=512/3185,                (32)
```

and at least one of its seven lines satisfies

```text
G_ell/m<=22295/512.                                   (33)
```

## 5. What this does and does not change at the frontier

The improvement is real but deliberately scalar.

* **Six-comb `H`-drift.**  Equations (25)--(28) tighten the complete-graph
  and common-period faces of the remaining `r=6` cone.  They do not couple
  the six phases, select the nerve spanning tree, or remove the `1/c`
  toothpick drift in THM-1199.
* **Fano/`chi_7`.**  Equations (32)--(33) give the still-open Fano probe a
  stronger global line-period budget.  The proof does not say which Fano
  line carries the debt or correlate it with `chi_7` colour, endpoint owner,
  or beat chronology.
* **The `j=4` flood tail.**  The larger `R_7` increases the global quadratic
  hole in the chamber whose pointwise slack peaks at multiplicity four.
  THM-1196's phase-local component-span law `d_5/d_4<189/8`, however, comes
  from the four-prefix dual survivor and two-tooth chronology, not from a
  global pair floor.  Its constant is unchanged.

Thus (2) closes the optimal **five-speed pair functional**, not the six-comb
coverage problem.  The missing object remains a phase-coupled transport
certificate combining loads, seam owners, beat meshes, and line periods.

THM-1205's concurrent active-pair formula is complementary rather than
redundant: its determinant `D` retains the located maximizer and phase-side
integer that the projective overlap carrier deliberately quotients out.  The
present theorem controls all ten magnitudes `rho_ij` but forgets `D`; THM-1205
retains one active `D` but does not control the six-comb cover.  A useful
synthesis must decorate the slow-gap transport stalk with an active-pair
determinant instead of expecting either scalar quotient to recover the other.

## 6. Tournament and carrier audit

The pairwise observable is the signed bulk defect

```text
w_ij=1/49-rho(s_i,s_j).                               (34)
```

On the sharp witness (3), orient an edge from the lower-speed vertex to the
higher when `w_ij>0` and reverse it when `w_ij<0`.  Zero edges, if present,
are switched by the declared tie Hamiltonian path `(0,1,3,2,4)`.  The witness
has no ties and fingerprint

```text
scores                         (3,2,2,1,2),
score multiset                  {1,2,2,2,3},
directed triangles              014,024,123,234,
strong components               one of size 5,
Hamiltonian paths               13,
flips from increasing order     4.                    (35)
```

This tournament preserves the sign topology of pair defects and correctly
exhibits a strongly connected compatibility packet.  It destroys the defect
magnitudes and exact reduced ratios, so it cannot distinguish the target
from the weaker THM-1231 constant.

Runner vertices were therefore challenged.  Alternatives considered were
runners, gaps, fixed circle sections, section boundaries, wall events,
residues, cover arcs, Fourier modes, Fano lines, matroid circuits, and proof
obligations.  The faithful carrier for (2) is

```text
projective speed coordinates + a rooted heavy spanning-tree ratio word.  (36)
```

It preserves every ratio needed to reconstruct all ten pair densities and
quotients out only common scaling.  It loses phase, absolute gcd, slow-gap
location, endpoint ownership, and cover chronology.  Those are precisely
the coordinates that must be restored before (2) can contribute to a full
six-comb contradiction.

## 7. Exact replay and scope

Run

```bash
python3 04-computation/lrc14_five_comb_sharp_compatibility_floor_thm1234.py
```

The Python referee checks (9) against exact interval geometry on all `6,857`
coprime cutoff rows, recomputes the heavy alphabet and hash, all six
disconnected partitions, the sharp and fallback rows, (22)--(33), and the
tournament fingerprint.  It then compiles and runs the integer-only C++
spanning-tree bank and requires its complete output byte for byte.  The stored
output is

```text
05-knowledge/results/lrc14_five_comb_sharp_compatibility_floor_thm1234.out.
```

Normal and `PYTHONOPTIMIZE=1` replays are required to be byte-identical.
Frozen SHA-256 hashes are

```text
Python referee  1300338cedfeff1df17e2beb1f5674d0d79eae605d4ed7420b4a1780129636fd
C++ exact bank  5a16e4ab322b30dbbbf5e1bb79f9cd14955117cd3d0b38539c88b16a3f0c87d9
stored output   248e7356dbe531042acf0b52284957cb47f6cb74fb3d6173cc5acef3a23a2dcf
```

This theorem sharpens every scalar consumer of the five-subset average, but
does not prove universal six-comb noncoverage, empty the full mixed-period
slow-gap cone, or prove LRC(14).
