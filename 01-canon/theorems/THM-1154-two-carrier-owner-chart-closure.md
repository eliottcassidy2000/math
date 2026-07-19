---
id: THM-1154
title: The residue-owner chart strengthens the carrier reduction and closes all thirty two-carrier switch obstructions
status: PROVED — an exact nonuniform refinement of the translated thirteen-grid argument.  With rho carriers, a broad-cap shift cM<=13s loses exactly d=#{p in P:cp>s} distinct core-owner multipliers and at most 2(6-rho) further multipliers, so d<2rho suffices.  In particular rho>=4 is automatic once the carrier phases are compatible, and rho=3 needs only d<=5.  On every exceptional pair from THM-1151, the large-carrier shift c=15 is compatible and has d=1, leaving at least three exact lonely-time charts.  Consequently every r=6 row with exactly two 13-divisible killers is lonely, with no Covering hypothesis and no bound on the noncarriers
source: codex-2026-07-18-S67 two-carrier residual attack
depends_on: [THM-1136, THM-1151]
related: [THM-1135, THM-1137]
script: 04-computation/lrc14_r6_two_carrier_owner_chart_referee_codex_S67.py
output: 05-knowledge/results/lrc14_r6_two_carrier_owner_chart_referee_codex_S67.out
---

# THM-1154 — the residue-owner chart closes the two-carrier branch

Let `P` be a seven-element subset of `{1,...,12}` and put `M=max(P)`.  Let the six
killers be distinct integers greater than `13M`, exactly two of which are
divisible by `13`.  Write those two carriers as

```text
z=13u<h=13v,       M<u<v.                              (1)
```

THM-1151 proves that the integer-shift chart of THM-1136 supplies a lonely
time unless

```text
8<=M<=12,       u=M+1,       13M+14<=v<=15M-1.         (2)
```

The thirty triples in (2) were obstructions to the *uniform core cap* used
there, not obstructions to the underlying translated grid.  Keeping the
owners of boundary residues instead of replacing all of them by `M` gives a
general strengthening of THM-1136 and closes the thirty triples.

## The general residue-owner chart

Let `O` be the set of normalized carriers, so the carrier killers are
`{13o:o in O}`, and put `rho=|O|`.  There are `6-rho` noncarrier killers.
Choose `s in O` and a positive integer `c` satisfying the broad core cap

```text
cM <= 13s.                                             (G1)
```

For `a=1,...,12`, consider

```text
t_a=a/13+c/(14*13s)=a/13+c/(182s).                    (G2)
```

> **Residue-owner lemma.**  Suppose every normalized carrier `o in O`
> satisfies
>
> ```text
> dist(co/s,14Z)>=1.                                   (G3)
> ```
>
> Put
>
> ```text
> d=#{p in P: cp>s}.                                   (G4)
> ```
>
> Then the core forbids exactly `d` of the twelve charts (G2), every
> noncarrier forbids at most two, and a lonely chart exists whenever
>
> ```text
> d+2(6-rho)<12,
> equivalently d<2rho.                                 (G5)
> ```

### Exact core owners and the no-wrap sign

Fix `p in P` and write `b=ap mod 13` in `{1,...,12}`.  Its phase at (G2) is

```text
b/13 + cp/(182s).                                      (G6)
```

The cap (G1) bounds the positive shift by `1/14`.  If `b<=11`, then (G6)
lies between `1/13` and

```text
11/13+1/14=167/182<13/14,                              (G7)
```

so it is safe.  If `b=12`, then (G6) is at most
`12/13+1/14=181/182<1`, so there is no wrap through the integer.  On this
last residue it is safe exactly when

```text
12/13+cp/(182s) <= 13/14
  iff cp<=s.                                           (G8)
```

Thus `p` forbids one chart exactly when `cp>s`, namely the unique chart

```text
a=-p^(-1) mod 13.                                      (G9)
```

Different core speeds have different inverses, so these owner charts are
distinct.  Their union therefore has cardinality exactly `d`, not merely at
most `|P|`.

For a carrier `13o`, the base term in (G2) is integral and its remaining
phase is `co/(14s)`, so (G3) is precisely its safety condition.  For a
noncarrier `k`, multiplication by `k` permutes `F_13^*`; hence its twelve
phases form a translate of the punctured `1/13`-grid.  An open arc of length
`1/7` contains at most two such points because `2/13>1/7`.  This proves the
two-point noncarrier bound.  Taking the union of the `d` singleton owner
hyperedges and `6-rho` hyperedges of size at most two proves (G5).

The scalar cap in THM-1136 was the special case `cM<=s`, which forces
`d=0`.  The owner ledger permits the thirteen-times-broader cap (G1) as long
as the crossed boundary owners fit in the multiplier budget.  In particular:

```text
rho=1: d<=1 suffices;
rho=2: d<=3 suffices;
rho=3: d<=5 suffices;
rho>=4: every d<=|P|=7 suffices.                        (G10)
```

The last line is conditional only on finding `(s,c)` satisfying (G1)--(G3);
it does not assert that every multiple-carrier ratio set has such a switch.

## The owner-sensitive exceptional chart

> **Owner-chart theorem.**  Assume (2).  Among the twelve integers
> `a=1,...,12`, at least three make
>
> ```text
> t_a = a/13 + 15/(14h) = a/13 + 15/(182v)             (3)
> ```
>
> safe at level `1/14` for the core, both carriers, and all four
> noncarriers.

### Both carriers are safe

For the large carrier, (3) gives

```text
h t_a = av + 15/14,
||h t_a|| = 1/14.                                      (4)
```

For the small carrier,

```text
z t_a = au + 15u/(14v).                                (5)
```

The endpoints in (2), together with `u=M+1`, give

```text
13u < v < 15u.                                         (6)
```

Hence

```text
1 < 15u/v < 15/13 < 2,                                 (7)
```

so the fractional term in (5) lies strictly between `1/14` and `1/7`.
Thus both carriers are safe, with (4) recording the only endpoint equality.
There is no hidden wrap: (7) lies in the first safe segment after zero.

### The core loses exactly one multiplier

Fix `p in P` and put

```text
r = ap mod 13 in {1,...,12},
y = 15p/(182v)>0.                                      (8)
```

If `r<=11`, then

```text
1/13 < r/13+y
       <= 11/13+15/182
       = 13/14,                                        (9)
```

where `p<v` was weakened to `p/v<=1`.  Thus this phase is safe.  If `r=12`,
then `12/13+y<1`: indeed `15p<14v` follows already from `p<=M<v` and
`15M<14v`.  On this last residue the safety condition is therefore exactly

```text
12/13+y <= 13/14
  iff 15p <= v.                                        (10)
```

For `p<=M-1`, the lower endpoint in (2) gives

```text
v-15p >= 13M+14-15(M-1)=29-2M>=5.                     (11)
```

For `p=M`, the upper endpoint gives `v<=15M-1<15M`.
Consequently a core speed can be endangered only when

```text
p=M and aM=-1 mod 13.                                  (12)
```

There is exactly one such multiplier `a` because `M` is nonzero modulo
`13`.  The other eleven multipliers are safe for the entire core.  This is
the information discarded by THM-1151's scalar cap: the cap charged all core
speeds as though they simultaneously owned residue `-1`, although a
multiplier has only one such owner.

### Four arbitrary noncarriers cost at most eight multipliers

Let `k` be any of the four noncarriers, so `13` does not divide `k`.  As `a`
runs through `1,...,12`, the residues `ak mod 13` run bijectively through
`1,...,12`.  Its twelve phases at (3) are therefore the translated grid

```text
{b/13 + 15k/(182v): b=1,...,12} modulo 1.              (13)
```

The open danger arc has length `1/7`.  Three distinct points of the
`1/13`-grid have circular span at least `2/13>1/7`, so (13) puts at most two
points in that arc.  Thus one noncarrier forbids at most two multipliers,
independently of its size, residue lift, or position among the killers.

The four noncarriers forbid at most eight multipliers in total.  Adding the
single core-owner loss (12) gives at most

```text
1 + 4*2 = 9 < 12.                                      (14)
```

At least three multipliers survive, and every corresponding time (3) is
safe for all thirteen speeds.  This proves the owner-chart theorem.

## Uniform two-carrier closure

> **Corollary.**  Every seven-core/six-killer row with exactly two
> `13`-divisible killers is `1/14`-lonely.

Outside (2), THM-1151 constructs a carrier-safe switch and THM-1136 absorbs
the four noncarriers.  Inside (2), the owner-chart theorem above applies.
Neither branch uses Covering.  In particular, the thirty rows listed by
THM-1151 are now an empty mathematical residual rather than a finite search
obligation.

For a concrete boundary example, take

```text
P=(1,2,7,9,10,11,12),
K=(157,158,159,160,169,2210).
```

Here `(M,u,v)=(12,13,170)`, the first triple in the last exceptional block.
The owner chart has seven surviving multipliers; for example `a=2` gives

```text
t=955/6188,
min_{w in P union K} ||w t|| = 1/14,
```

with the large carrier `2210` as the unique equality owner.

## Tournament Analysis and assumption challenge

The speed-order tournament on the two carriers is again the unique transitive
two-vertex tournament: score multiset `{0,1}`, no directed cycle, singleton
SCCs, and one Hamiltonian path.  It cannot see the closure above.  The same is
true of a tournament on the four noncarriers ordered by magnitude: its
orientation contains no translated-grid incidence data.

The proof-bearing vertices are instead the twelve **multiplier charts**.  A
core boundary owner marks one vertex, while each noncarrier adds a hyperedge
of size at most two.  The surviving-chart predicate is the complement of the
union of these five hyperedges.  A pairwise orientation can record, for
example, which chart has the larger number of obstructing owners, but it
destroys the hyperedge union whose capacity inequality is exactly (14).

Candidate vertex sets challenged here are runners, carriers, carrier ratios,
danger arcs, residue owners, multiplier charts, and proof obligations.  The
faithful quotient uses multiplier charts as vertices and owner sets as
hyperedges.  It preserves existence of the explicit lonely time (3), while
discarding continuous component geometry and killer order, neither of which
the certificate needs.

## Exact replay

The dependency-free referee uses only integer and `Fraction` arithmetic.  It
checks 1,586,779 seven-core broad-cap charts and 1,649,232 individual phase
decisions, verifies the carrier condition (G3) on 1,087,262 candidate charts,
and exhausts the abstract owner/noncarrier hyperedge capacities for every
`rho=2,...,6`.  It separately checks all thirty endpoint triples, both carrier
inequalities and wrap orientations, the exact one-owner core mask, and every
noncarrier in a large lift box.  The displayed witness is recomputed exactly.
Ordinary and optimized executions must be byte-identical.

```text
04-computation/lrc14_r6_two_carrier_owner_chart_referee_codex_S67.py
SHA-256 1a876069e3b9f96d5f291524f6a53ce339a9f1e480c70b1c56fc3d0a9f0d4bbd

05-knowledge/results/lrc14_r6_two_carrier_owner_chart_referee_codex_S67.out
SHA-256 b9233ef2be521912eb854ea2f92cae7f002c964554a3b2399c91b85785e33197
```
