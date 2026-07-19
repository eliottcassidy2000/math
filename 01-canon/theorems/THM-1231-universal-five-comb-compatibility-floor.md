---
id: THM-1231
title: The compatibility-sensitive five-comb pair floor
status: PROVED by an analytic all-range reduction and an exact-integer finite bank.  Every five distinct radius-1/14 danger combs have total pair overlap at least 13/81.  Consequently every six have pair mass at least 13/54, restoring the stronger THM-1179 slow-gap gcd inequality by a valid five-way compatibility route; every seven have pair mass at least 91/270.  No claim that 13/81 is the optimal five-comb constant is made
source: codex-2026-07-18 five-speed pair-compatibility session
depends_on: [THM-1166, THM-1179, THM-1191]
related: [THM-1176]
script: 04-computation/lrc14_five_comb_compatibility_floor_codex_20260718.py; 04-computation/lrc14_five_comb_compatibility_floor_bank_codex_20260718.cpp
output: 05-knowledge/results/lrc14_five_comb_compatibility_floor_codex_20260718.out
---

# THM-1231 -- compatibility-sensitive five-comb pair floor

For a positive integer `s`, put

```text
D_s={t in R/Z: ||st||<1/14},
rho(a,b)=measure(D_a intersect D_b).                    (1)
```

If `S` is a finite set of distinct speeds, write

```text
R(S)=sum_{{a,b} subset S} rho(a,b).                     (2)
```

The four-comb route proposed before THM-1191 was false: a `13`-adic quartet
has pair mass below `13/135`.  Five-way compatibility nevertheless supplies
exactly the constant that route had sought at six speeds.

> **Theorem (universal five-comb floor).**  Every five distinct positive
> integer speeds satisfy
>
> ```text
> R(S)>=13/81.                                          (3)
> ```

The proof has an analytic all-range reduction followed by a finite exact
bank.  Floating point is not used in the certificate.

## 1. The heavy graph must be connected

THM-1166 proves, for distinct speeds,

```text
rho(a,b)>=1/91                                         (4)
```

and, for every triple,

```text
rho12+rho13+rho23>=1/24.                               (5)
```

Summing (5) over the four triples of a four-speed set counts every pair
twice, hence

```text
R_4>=1/12.                                             (6)
```

Set

```text
theta=71/63504,
c=1/49-theta,                                          (7)
```

and call a pair **heavy** when

```text
rho(a,b)<c.                                            (8)
```

Thus a nonheavy pair contributes at least `c`.  Form the simple heavy graph
on the five speed labels.  If it is disconnected, let the component sizes
give a partition of five.  Cross-component pairs are nonheavy.  Inside a
component of size `1,2,3,4`, use respectively the lower bounds

```text
0, 1/91, 1/24, 1/12.                                  (9)
```

Every disconnected partition then gives the following exact ledger.

```text
component sizes    cross edges    lower bound for R      excess over 13/81
4+1                    4               13/81                   0
3+2                    6             1655/9828              233/29484
3+1+1                  7              229/1296                7/432
2+2+1                  8             2599/14742             233/14742
2+1+1+1                9             2419/13104            2843/117936
1+1+1+1+1             10              125/648                 7/216. (10)
```

For example, the first line is the identity

```text
1/12+4(1/49-71/63504)=13/81.                          (11)
```

Consequently a counterexample with `R(S)<13/81` would have a **connected**
heavy graph.  This is the compatibility step missing from a black-box
average of four-subset constants.

## 2. Every heavy edge belongs to a finite ratio alphabet

Reduce a pair to coprime integers `u<v`.  THM-1166's folded formula and
defect bound are

```text
rho(u,v)
 =[4uv+F((u+v) mod 14)-F((v-u) mod 14)]/(196uv),
F(r)=r(14-r),                                          (12)

rho(u,v)>=1/49-1/(7v).                                 (13)
```

If the pair is heavy, (7)--(8) and (13) imply

```text
theta<1/49-rho(u,v)<=1/(7v),
v<1/(7theta)=9072/71<128.                              (14)
```

Thus `v<=127`.  Exact evaluation of (12) over all coprime
`1<=u<v<=127` leaves exactly

```text
61 heavy unoriented ratio types.                       (15)
```

Their actual largest upper coordinate is `97`.  Serialize the rows in
increasing `v`, then increasing `u`, as ASCII words `u:v` separated by one
space.  This ordered serialization has SHA-256

```text
b18b3cac5db13715b062292e288e9d061c21e805f348eacb2d9e6e04c549b151. (16)
```

This is a genuine all-range reduction: (14), rather than an experimental
speed ceiling, proves that no ratio beyond the bank can be heavy.

## 3. The projective spanning-tree bank

A connected heavy graph has a heavy spanning tree.  Root such a tree and
divide all speeds by the root speed.  Starting from the singleton coordinate
`1`, attach each new tree vertex to an existing coordinate `x` by

```text
x -> x(v/u) or x(u/v),                                 (17)
```

where `u:v` is one of the 61 types.  After each attachment, clear
denominators, divide the whole row by its gcd, and sort.  Induction on a
rooted-tree order proves that this procedure enumerates the projective class
of every possible counterexample.  Common scaling is harmless because every
pair overlap in (12) is dilation invariant.

The bank prunes only with a rigorous lower bound.  At integer scale

```text
Q=10^9,
L_Q(A)=Q^(-1) sum_{{a,b} subset A} floor(Q rho(a,b)).  (18)
```

For a partial row of size `n`, every one of the
`10-binomial(n,2)` unformed pairs contributes at least `1/91`.  Therefore
the row can be discarded only when

```text
L_Q(A)+[10-binomial(n,2)]/91>=13/81.                   (19)
```

Each floor in (18) is performed by exact unsigned `128`-bit integer division;
the comparison in (19) is unsigned `128`-bit cross multiplication.  No
approximate comparison can remove a row.  Along a four-edge tree, clearing
denominators introduces at most one factor at most `127` per edge, so

```text
max(A)<=127^(n-1)<=127^4=260,144,641.                 (19a)
```

Consequently `196ab<1.33*10^19<2^64` at the terminal level; every larger
fixed-grid/cross-product accumulator is explicitly promoted to unsigned
`128`-bit arithmetic.  The program exits nonzero on a coordinate or
numerator overflow.

The exact bank sizes are

```text
n=2:         61 surviving projective rows,
n=3:      7,179 surviving projective rows,
n=4:    720,678 surviving projective rows.             (20)
```

Attaching the fifth vertex scans

```text
347,319,154 terminal extensions.                       (21)
```

Every extension passes (18) at the terminal level.  The least fixed-grid
lower numerator is

```text
161172154/10^9,                                        (22)
```

and its exact margin over the target is certified by

```text
161172154*81-13*10^9=54,944,474>0.                    (23)
```

Hence every enumerated terminal row has `R>13/81`.  By the connected-heavy
reduction, this contradicts the existence of a counterexample and proves
(3).

As telemetry, the lexicographically first fixed-grid minimizer is

```text
(1,12,13,27,156),                                     (24)
```

with exact pair word

```text
2*(1/63+17/819+1/84+1/91+23/1092)=44/273,
44/273-13/81=5/7371.                                  (25)
```

The bank proves (3), not that `44/273` is the optimal universal constant.
The rounding in (18) was deliberately chosen only to separate `13/81`.

## 4. Six-speed and all-cardinality averaging

Let `T` contain six speeds.  There are six five-subsets, and each pair lies
in exactly four of them.  Therefore

```text
sum_(A subset T, |A|=5) R(A)=4R(T).                   (26)
```

Applying (3) to all six terms gives the compatibility-sensitive consequence

```text
R(T)>=13/54.                                           (27)
```

This does **not** revive the false four-comb floor refuted by THM-1191.
It bypasses that guardrail with a separately proved five-comb theorem.

More generally, for `m>=5`, each pair belongs to `binomial(m-2,3)`
five-subsets, so

```text
R_m>=13m(m-1)/1620.                                   (28)
```

In particular

```text
R_6>=13/54,             R_7>=91/270.                  (29)
```

THM-1166's seven-comb quadratic certificate therefore improves from global
uncovered mass `1/12` to

```text
measure(uncovered by seven combs)>=(2/7)(91/270)=13/135. (30)
```

## 5. The restored slow-gap consequences

In THM-1179, let six faster combs cover an `a`-slow gap and put

```text
delta=a sum_i 1/b_i-1,
g_ij=gcd(b_i,b_j),
eta_ij=rho_ij(1-rho_ij).                              (31)
```

The complete-graph debt theorem says

```text
a sum_(i<j) eta_ij/g_ij +(18/49)delta
 >=(6/7)R_6.                                           (32)
```

Using (27) in (32) gives

```text
a sum_(i<j) eta_ij/g_ij +(18/49)delta>=13/63.          (33)
```

Since `eta_ij<=13/196`, this implies the restored reciprocal-gcd law

```text
a sum_(i<j) 1/g_ij>=28/9-(72/13)delta.                (34)
```

This is the stronger constant whose earlier quartet-based justification was
invalidated by THM-1191; (26)--(27) now provide a valid proof frame.

The optimal six-comb quadratic majorant also yields

```text
measure(union_i D_bi)
 <=6/7-(1/3)R_6<=881/1134.                            (35)
```

If the six killers share `G_0=gcd(b_1,...,b_6)`, the periodic-component
argument of THM-1179 sharpens to

```text
G_0/a< (7/6)(881/1134)=881/972.                       (36)
```

The inequality is strict in the actual open-danger convention because the
union must contain the closed slow gap.

## 6. Tournament and alternate-carrier audit

The pairwise observable is the signed bulk defect

```text
w_ij=1/49-rho(s_i,s_j).                               (37)
```

On the telemetry row (24), orient from the smaller-speed vertex to the
larger when `w_ij>0`, and reverse when `w_ij<0`.  Zero edges, if present,
are switched according to the declared tie Hamiltonian path

```text
(0,1,3,2,4).                                           (38)
```

There are no ties in (24).  Its tournament fingerprint is

```text
score sequence             (3,2,2,1,2),
score multiset              {1,2,2,2,3},
directed triangles          014, 024, 123, 234,
strong components           one of size 5,
Hamiltonian paths           13,
flips from increasing order 4.                         (39)
```

This tournament preserves the sign topology of cheap versus expensive
pairs and exposes the strongly connected compatibility packet.  It destroys
the defect magnitudes and the exact reduced ratios, so it cannot certify
(3) by itself.

Runner vertices were explicitly challenged.  The faithful finite carrier is

```text
projective speed coordinates + a heavy spanning-tree ratio word.          (40)
```

It preserves every ratio needed to reconstruct all ten pair densities and
quotients out only common scaling.  That quotient destroys absolute gcds,
slow-gap location, phase, and endpoint ownership; those data re-enter only
when (27) is inserted into THM-1179.  Alternate vertices considered were
runners, gaps, fixed circle sections, section boundaries, wall events,
residues, cover arcs, Fourier modes, Fano lines, matroid circuits, and proof
obligations.  Gaps and sections preserve local coverage but are irrelevant to
the global pair functional; residues alone lose ratio size; the rooted-tree
proof obligations give a useful enumeration order but (40) is the minimal
carrier that retains exact pair mass.

## 7. Exact replay and scope

Run

```bash
python3 04-computation/lrc14_five_comb_compatibility_floor_codex_20260718.py
```

The Python referee independently checks (12) against 780 exact interval-
geometry rows, recomputes (10), the 61-type bank and its hash, (25)--(36),
and (39).  It then compiles and runs the integer-only C++ spanning-tree bank.
The stored output is

```text
05-knowledge/results/lrc14_five_comb_compatibility_floor_codex_20260718.out.
```

Normal and `PYTHONOPTIMIZE=1` replays are byte-identical.  SHA-256 values are

```text
Python referee  b6a12198350e2578e09d9e320a9770bf316b3415690dfff61616618325874dcf
C++ exact bank  addf2fa8ba4959fd69148ad7e7d42914c81f31885ffcde151fbbd432b942af91
stored output   69882188ab6787f3fd833adfa56394037eee2d7531f073252a456819a8277b5b.
```

This theorem materially strengthens the global pair and slow-gap gcd
ledgers, but does not by itself close the remaining mixed-period slow-gap
cone or LRC(14).
