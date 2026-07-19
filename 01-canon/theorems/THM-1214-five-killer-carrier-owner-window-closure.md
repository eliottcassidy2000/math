---
id: THM-1214
title: Five-killer carrier-owner windows close the clustered r=5 stratum
status: PROVED uniformly.  The five-killer owner ledger has two more multiplier slots than the six-killer ledger.  Zero and one carrier close directly; two carriers close by an endpoint integer contradiction; three carriers reduce to 478 exact first-comb rows and three bounded two-comb candidates, with zero point-only residuals and zero covers; four and five carriers close by density and one nested first-safe window.  This proves the entire eight-core/five-killer clustered stratum, not general LRC(14)
source: codex-2026-07-18-S78
depends_on: [THM-1136, THM-1154, THM-1155, THM-1169, THM-1212]
related: [THM-1148, THM-1168, THM-1213]
script: 04-computation/lrc14_r5_carrier_owner_window_referee_codex_S78.py
output: 05-knowledge/results/lrc14_r5_carrier_owner_window_referee_codex_S78.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCFiveKillerCarrierWindow.lean
---

# THM-1214 -- five-killer carrier-owner windows

Let

```text
P={p1<...<p8} subset {1,...,12},       M=p8,
13M<k1<...<k5.                                           (1)
```

Call a killer a **carrier** when it is divisible by `13`.  If there are
`rho` carriers, write their normalized speeds as

```text
O={o1<...<o_rho},       x=o1,       carrier killers={13o:o in O}.  (2)
```

> **Theorem.** Every family (1) has a time at which all thirteen speeds in
> `P union {k1,...,k5}` have distance at least `1/14` from the integers.

Thus the complete clustered `r=5` stratum is uniformly closed.  No Covering
hypothesis and no upper bound on a killer are used.

## 1. The five-killer owner ledger

The owner argument must be restated because THM-1154 has seven core speeds
and six killers, whereas (1) has eight and five.  The proof mechanism is the
same, but its obstruction budget is genuinely different.

Suppose `rho>=1`.  In the normalized carrier clock `u`, choose a carrier-safe
boundary

```text
u=c/(14s),       s in O,       c a positive integer,                 (3)
```

such that

```text
cM<=13s.                                                           (4)
```

For each `a=1,...,12`, put

```text
t_a=a/13+u/13=a/13+c/(182s).                                      (5)
```

Let

```text
d=#{p in P:cp>s}.                                                   (6)
```

> **Five-killer residue-owner lemma.** If (3) is safe for every normalized
> carrier, then exactly `d` distinct charts in (5) are forbidden by the
> core, each noncarrier forbids at most two charts, and a chart survives
> whenever
>
> ```text
> d+2(5-rho)<12,
> equivalently d<2rho+2.                                           (7)
> ```

**Proof.**  Fix `p in P` and write `b=ap mod 13` in `{1,...,12}`.  Its phase
at (5) is

```text
b/13+cp/(182s).                                                     (8)
```

By (4), the positive shift is at most `1/14`.  If `b<=11`, (8) lies in the
closed safe band.  If `b=12`, there is still no wrap, and (8) is safe exactly
when `cp<=s`.  Thus a crossed owner `p` forbids precisely the chart
`a=-p^(-1) mod 13`.  Distinct core speeds have distinct inverses, so the
core loss is exactly `d`.

For a noncarrier `k`, the twelve phases in (5) form a translate of the
punctured `1/13` grid.  An open arc of length `1/7` contains at most two grid
points, since `2/13>1/7`.  Hence all `5-rho` noncarriers cost at most
`2(5-rho)` charts.  A carrier `13o` has phase `ou`, so it is safe by the
hypothesis.  The union bound on forbidden charts proves (7). ∎

The useful budgets are

```text
rho=1: d<=3 suffices;       rho=2: d<=5 suffices;
rho=3: d<=7 suffices;       rho=4,5: every d<=8 suffices.            (9)
```

These are two owner slots stronger than the corresponding six-killer
inequality `d<2rho`.

## 2. Reading a boundary from a first-safe window

Put

```text
D_o={u:||ou||<1/14},
I_x=[1/(14x),13/(14x)].                                             (10)
```

The least carrier `x` is safe throughout the closed interval `I_x`.  More
generally, for `x<=q<13x`, every carrier `o` with `x<=o<=q` is safe on

```text
J(x,q)=[1/(14x),13/(14q)].                                         (11)
```

Indeed `1/14<=ou<=13/14` there, with safe equality endpoints.

If a closed subinterval `J subset I_x` is not covered by the open danger
sets of the later carriers, its carrier-safe subset is compact and nonempty.
Take its least point `u`.  Either `u` is the left endpoint, or it is a danger
boundary of some carrier `s`; in both cases it has form (3).  Since
`u<=13/(14x)`,

```text
cx<=13s,       hence cM<13s.                                       (12)
```

If also `u<=1/(14p_j)`, then `cp_j<=s`; the first `j` core speeds are not
crossed and

```text
d<=8-j.                                                            (13)
```

Thus every case below consists of finding a carrier-safe point in a window
whose right endpoint retains the required number of core owners.

We use THM-1155 in the form: if `m` danger combs of speeds `W` cover a closed
interval of length `L`, then

```text
(7-m)L<=sum_{w in W}1/w.                                            (14)
```

## 3. Zero and one carrier

If `rho=0`, every speed in (1) is nonzero modulo `13`.  The single time

```text
t=1/13                                                             (15)
```

puts every phase at a nonzero thirteenth residue, hence at distance at least
`1/13>1/14`.  This row is immediate.

If `rho=1`, take `u=1/(14x)`, so `(s,c)=(x,1)`.  Since every core speed is
less than `x`, (6) has `d=0`.  The four noncarriers cost at most eight charts,
leaving at least four survivors.  This is also the five-killer specialization
of THM-1136's unique-carrier mechanism.

## 4. Two carriers: a one-tooth integer contradiction

Let `O={x<y}` and put `p=p3`.  An eight-subset of `{1,...,12}` satisfies

```text
3<=p<=7,       M>=p+5,       x>=p+6.                              (16)
```

Use the owner-truncated window

```text
J=[1/(14x), min(13/(14x),1/(14p))].                                (17)
```

If `x>=13p`, this is the full first-safe window, of length `6/(7x)`, which
is longer than one open `y`-tooth of length `1/(7y)`.

Suppose `x<13p`.  Then

```text
|J|=(x-p)/(14xp).                                                   (18)
```

If `x>=3p`, (18) is at least `1/(7x)>1/(7y)`, so again `D_y` cannot cover
the closed interval.  It remains to consider `x<3p`.  If one open `y`-tooth,
centred at `n/y`, contained `J`, endpoint strictness would give

```text
x(14n-1)<y<p(14n+1).                                               (19)
```

The right endpoint forces `n>=1`.  Eliminating `y` gives

```text
14n(x-p)<x+p.                                                       (20)
```

But the left side is at least `14*6=84` by (16), whereas
`x+p<4p<=28`.  This contradiction proves that (17) contains a point safe for
both carriers.

Its boundary satisfies `cp3<=s`, so `d<=5`.  The three noncarriers cost at
most six charts, and `5+6=11<12`.  Hence every two-carrier row closes.

The referee independently checks the twenty possible residual `(p,x)` rows,
117 bounded `y` rows, and the exact open-tooth endpoints; it finds no cover.

## 5. Three carriers: the tiny exact wall-cell lemma

Let `O={x<y<z}` and now put `p=p1`.  Then

```text
1<=p<=5,       M>=p+7,       x>=p+8.                              (21)
```

Use

```text
J_{p,x}=[1/(14x), min(13/(14x),1/(14p))].                          (22)
```

If the two later combs covered this interval, (14) would imply

```text
5|J_{p,x}|<=1/y+1/z<2/(x+1).                                      (23)
```

For `x>=13p`, the full-window identity

```text
6/(7x)-2/[5(x+1)]
  =(16x+30)/[35x(x+1)]>0                                          (24)
```

contradicts (23).  For `x<13p`, density closes exactly when

```text
5(x-p)(x+1)>28xp.                                                  (25)
```

The complementary integer rows are

| `p=p1` | cores with this `p` | first `x` closed by (25) | residual `x` |
|---:|---:|---:|---:|
| 1 | 330 | 9  | none |
| 2 | 120 | 13 | 10..12 |
| 3 | 36  | 19 | 11..18 |
| 4 | 8   | 26 | 12..25 |
| 5 | 1   | 33 | 13..32 |

> **Finite three-carrier lemma.** None of the residual windows (22) is
> covered by two distinct later carrier combs.

Here is the exhaustive reduction.  From (23) and `z>y`, a cover forces

```text
y<2/[5|J_{p,x}|].                                                   (26)
```

There are 478 resulting `(p,x,y)` rows, by stratum

```text
(0,6,51,141,280).                                                   (27)
```

Subtract the open `y`-comb exactly.  Every row retains a positive-length
closed cell; there are zero one-comb covers and zero point-only complements.
If `ell` is the longest retained cell, the `z`-comb could cover it only if

```text
z<1/(7ell).                                                        (28)
```

Only three integer candidates survive:

```text
(p,x,y,z)=(5,13,14,15),(5,13,14,16),(5,13,15,16),
ell=4/455.                                                         (29)
```

The three terminal rows do not need a black-box endpoint decision.  In every
one of them, `y,z in {14,15,16}` and the window is

```text
J_{5,13}=[1/182,1/70].
```

For every `q in {14,15,16}` and every `u in J_{5,13}` there is no wrap and

```text
1/13=14/182 <= qu <= 16/70=8/35 <13/14.                (29a)
```

Thus the **entire** terminal window is safe for both later carriers.  The two
independent exact tests provide a redundant endpoint audit: one checks every
rational wall and one midpoint of every open cell; the other merges strict
open-tooth overlaps without merging touching intervals.  They agree with
(29a) and find no cover.  This proves the lemma.

A least carrier-safe boundary in (22) has `cp1<=s`, so `d<=7`.  The two
noncarriers cost at most four charts; `7+4=11<12`.  Every three-carrier row
therefore closes.

The finite lemma is a complete exact endpoint theorem, not a sampled height
claim: (26) and (28) prove both integer cutoffs, and the zero point-only count
closes the only route by which `z` could escape the second cutoff.

## 6. Four and five carriers

For `rho=4`, suppose the three later carriers covered `I_x`, whose length is
`L=6/(7x)`.  Equation (14) would give

```text
4L<=sum_{i=2}^4 1/o_i<3/x,
```

but `4L-3/x=3/(7x)>0`.  Thus `I_x` has a common carrier-safe point.  The
boundary ledger has `d<=8`, while the one noncarrier costs at most two;
`8+2=10<12`.

For `rho=5`, write the carriers as `x<y<z<w<v`.  If the four later combs
covered `I_x`, (14) would first force

```text
y<14x/9.                                                           (30)
```

The prefix window `J(x,y)` from (11) is safe for `x,y`.  If the remaining
three combs covered it, (14) would give

```text
4|J(x,y)|<3/y.                                                      (31)
```

Yet

```text
4|J(x,y)|-3/y=(5x-2y)/(7xy)>0,                                    (32)
```

because (30) is stronger than `y<5x/2`.  This contradiction supplies a
carrier-safe boundary.  There are no noncarriers and `d<=8<12`.

This is the five-carrier step of THM-1212's nested-window recursion; no
finite shape bank is involved.

## 7. Consequence and proof-map correction

The six carrier cardinalities now form a complete dispatch:

```text
rho=0: mod-13 time 1/13;
rho=1: least-carrier boundary;
rho=2: third-owner window and the one-tooth contradiction;
rho=3: first-owner window and the 478-row/3-pair exact lemma;
rho=4: three-comb density;
rho=5: one nested prefix window.
```

> **Uniform clustered r=5 closure.** For every eight-element
> `P subset {1,...,12}` and all five ordered killers above `13 max(P)`, the
> resulting thirteen-speed family is `1/14`-lonely.

THM-1148, THM-1168, and THM-1213 remain valid independent metric/shape
theorems, but their proportional four-comb residual is no longer an open
branch of the LRC(14) proof map.  This theorem does not close the adjacent
nonclustered or global inverse/covering obligations and is not a proof of
LRC(14).

## 8. Tournament and alternate-vertex audit

Orient the six carrier-cardinality strata by proof-obligation size, with
`rho` as the tie gauge.  The telemetry tournament is transitive, with score
sequence `(5,4,3,2,1,0)`, no directed cycles, six singleton SCCs, one
Hamiltonian path, and fifteen flips when the gauge is reversed.

This tournament is only a dispatcher.  The predicate-preserving object is

```text
(twelve multiplier charts,
 distinct core-owner singleton hyperedges,
 noncarrier hyperedges of size at most two,
 labelled carrier wall cells with open/closed endpoint data).       (33)
```

Runner or carrier tournaments destroy the hyperedge union in (7).  A wall
order destroys the metric lengths used in (23), (28), and (32).  A residue
graph without boundary owners cannot recover `d`.  Candidate vertices
explicitly challenged here are runners, carriers, core gaps, wall events,
residues, multiplier charts, and cover obligations.  The decorated owner
hypergraph/wall-cell complex in (33) is the smallest carrier used by the
proof.

## 9. Exact replay and formal boundary

The dependency-free referee uses only integers and `Fraction`; every check
uses `require`, so ordinary and optimized Python execute the same
obligations.  It verifies both core-order histograms, every owner budget,
the rho-two endpoint contradiction, the complete rho-three cutoffs and two
independent cover tests, the rho-four/five identities, and the tournament
fingerprint.  Ordinary and optimized output are byte-identical to the frozen
result.

```text
04-computation/lrc14_r5_carrier_owner_window_referee_codex_S78.py
SHA-256 cd6f615d1e3806125ecdd5ff0bfc508fd67559b9b540344f318a1b75918c8490

05-knowledge/results/lrc14_r5_carrier_owner_window_referee_codex_S78.out
SHA-256 94447a6dcbf63905007df0d5fb3d63f430d78bc00a58feb286a988ce011313cf
```

The reusable formal boundary is small: the owner-capacity inequality, the
rho-two endpoint contradiction, the rho-three density cutoffs, and the
rho-four/five rational inequalities.  The terminal three-row geometry is no
longer part of that boundary: `LRCFiveKillerCarrierWindow.lean` kernel-checks
the uniform safe band (29a).  It also instantiates the multi-comb density
theorem on the actual closed interval `I_x` and proves the rho-four
three-comb cover impossible as a set inclusion.  What remains external is the
exact 478-row reduction to the terminal rows, the analogous full rho-five
set-level composition, and the assembly of the owner charts.
