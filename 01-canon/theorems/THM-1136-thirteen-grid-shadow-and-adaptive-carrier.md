---
id: THM-1136
title: The thirteen-grid shadow gives a full sharp component and absorbs every noncarrier at arbitrary scale
status: PROVED — exact F_13 grid geometry plus interval topology.  If the first five killers are nonzero modulo 13, their residual contains a component longer than 1/(7k5).  More generally, after a small integer shift from a thirteen-grid observer, each noncarrier forbids at most two of the twelve multipliers regardless of its magnitude; only compatibility among the 13-divisible carriers remains.  In particular every row with a unique 13-carrier is lonely.  Combining the least-carrier obstruction with THM-1135's s=5 ratio horn shows that the least carrier of any surviving row is one of k1,...,k4.  A weighted lift-aware residue quotient and the earlier s,d cones are retained as corollaries.  The complete translated-grid cells, residue universes, weighted incidence assignments, and endpoint inequalities are independently replayed by an exact Fraction referee
source: codex-2026-07-18-S67 sharp-horn uniformity continuation, with the adaptive-carrier inequality from the concurrent tail audit
depends_on:
  - THM-1132   # sharp final-comb target; this theorem proves two uniform pieces
  - THM-1135   # the s=5 proper-prefix ratio horn locates the least carrier
related:
  - THM-1121   # exact bounded r=6 atlas
  - THM-1128   # four-center thirteen-grid rectangle
  - THM-1134   # sequential scalar-tail refutation
script: 04-computation/lrc14_r6_q13_shadow_adaptive_referee_codex_S67.py
output: 05-knowledge/results/lrc14_r6_q13_shadow_adaptive_referee_codex_S67.out
---

# THM-1136 — the thirteen-grid shadow and carrier-ratio reduction

For an integer `k>0`, put

```text
D_k={t in R/Z: ||kt||<1/14}.
```

The strict inequality is important: `D_k` is a union of **open** teeth.  Let
`P` be a nonempty subset of `{1,...,12}`, put `M=max(P)`, and let

```text
13M<k1<k2<...<k6.                                      (1)
```

These are the seven-core/six-killer coordinates of the `r=6` clustered
stratum, although neither theorem below needs the Covering hypothesis.

## 1. Exact endpoint topology: equality already suffices

The sharp final-comb criterion has a useful non-strict form which is hidden
if endpoints are discarded.

> **Closed-component lemma.**  If `I` is a closed circular interval and
> `|I|>=1/(7k)`, then `I` is not contained in `D_k`.

Indeed, a connected subset of `D_k` lies in one open tooth, whose length is
exactly `1/(7k)`.  A closed interval of the same length cannot be contained
in that open interval: strict containment at both endpoints would make its
length strictly smaller.  Hence `I` contains a point with `||kt||>=1/14`.

Thus the exact closed-safe residual needs

```text
L>=1/(7k),                                             (2)
```

not `L>1/(7k)`.  The strict version in THM-1132 remains a convention-free
sufficient statement, but equality is already decisive for the actual LRC
predicate.  The referee records one tooth closure whose two endpoints both
have distance exactly `1/14`.

## 2. No carrier among the five removed killers gives a sharp component

> **Theorem A (five-comb thirteen-grid shadow).**  If
>
> ```text
> 13 does not divide ki,       i=1,...,5,              (3)
> ```
>
> then
>
> ```text
> G={t: ||pt||>=1/14 for p in P,
>       ||ki t||>=1/14 for i=1,...,5}
> ```
>
> contains a component of length strictly greater than `1/(7k5)`.
> Consequently `P union {k1,...,k6}` is `1/14`-lonely.

### Proof

Work in `F_13^*`, which has twelve multipliers.  A nonzero residue `r`
forbids only the two multipliers

```text
a=+r^(-1), -r^(-1),                                   (4)
```

if the goal is to make `ar` avoid `{+1,-1}`.  The five killer residues
therefore forbid at most ten multipliers.  Choose `a` outside their union.
Then

```text
ak_i mod 13 is not in {0,+1,-1},
||a k_i/13||>=2/13,       i=1,...,5.                  (5)
```

Put

```text
t0=a/13,                 h=1/(14k5),
J=[t0-h,t0+h].                                           (6)
```

For `p in P` and `|u|<=h`, reverse Lipschitz continuity gives

```text
||p(t0+u)||
 >= 1/13-p|u|
 >= 1/13-M/(14k5)
 > 1/13-1/182
 = 1/14,                                                 (7)
```

where the strict inequality is exactly `k5>13M`.  For every removed killer,

```text
||ki(t0+u)||
 >=2/13-ki/(14k5)
 >=2/13-1/14
 =15/182
 >1/14.                                                 (8)
```

Thus the whole closed interval `J`, of length exactly `1/(7k5)`, is
**strictly** safe for the core and the five removed killers.  Its minimum
clearance has positive slack; compactness and continuity therefore thicken
`J` slightly while preserving all twelve inequalities.  The component of
`G` containing it has length strictly greater than `1/(7k5)`.

Since `k6>k5`, this length is also greater than `1/(7k6)`.  One open tooth of
`D_k6` cannot cover the component, so a point survives all thirteen speeds.
This proves Theorem A.  Notice that the proof produces a metric interval,
not merely a rational witness.

The exact example replayed by the referee is

```text
P=(1,2,4,7,9,11,12),
(k1,...,k5)=(171,173,175,177,179),
a=1,
J=[2493/32578,2519/32578],
|J|=1/1253=1/(7*179).
```

The exact lower bound throughout `J` is `1175/16289>1/14`, owned by speed
`12`.  Appending `k6=182`, either endpoint is already `182`-safe with
distance `13/179`.

## 3. The translated grid absorbs every noncarrier, at every scale

Let

```text
C={ki:13 divides ki}                                    (9)
```

be the carrier set, and suppose it is nonempty.  For a real number `x`, write
`dist(x,14Z)=min{|x-14m|:m in Z}`.

> **Theorem B (integer-shift carrier reduction).**  Choose `z in C` and a
> positive integer `c` such that
>
> ```text
> c<=floor(z/(13M)).                                   (10)
> ```
>
> If every carrier `h in C` satisfies
>
> ```text
> dist(c h/z,14Z)>=1,                                  (11)
> ```
>
> then some `a in F_13^*` makes
>
> ```text
> t=a/13+c/(14z)                                       (12)
> ```
>
> a `1/14`-lonely time for all of `P union {k1,...,k6}`.

### Proof

For every core speed `p`, condition (10) and reverse Lipschitz continuity
give

```text
||pt|| >= ||ap/13||-pc/(14z)
       >= 1/13-Mc/(14z)
       >= 1/13-1/182
        = 1/14.                                        (13)
```

The possible equality is safe because the danger teeth are open.

Now fix one noncarrier `k`.  As `a` runs through `F_13^*`, the residue
`ak mod 13` runs bijectively through `1,...,12`.  Thus its twelve possible
phases are the translated grid

```text
{b/13+c k/(14z): b=1,...,12} modulo 1.                 (14)
```

One open danger tooth has length `1/7`.  No such tooth contains three points
of a `1/13`-spaced grid, because three points would span at least
`2/13>1/7`.  Therefore this noncarrier forbids at most two multipliers `a`,
**independently of the size of `k/z`**.  Since `z` is a carrier, there are at
most five noncarriers.  Their forbidden sets have union of size at most
`2*5=10<12`, so choose a multiplier outside that union.

Finally, if `h in C`, then `ah/13` is integral and

```text
||ht||=||c h/(14z)||=dist(c h/z,14Z)/14>=1/14           (15)
```

by (11).  This proves the theorem.  Notice the quantifier improvement: the
noncarrier ratio has disappeared completely.  The only surviving obstruction
is a finite compatibility problem among the carriers themselves.

For `c=1`, the chosen carrier `z` satisfies (11) at equality.  We immediately
obtain two useful corollaries.

> **Unique-carrier corollary.**  If exactly one killer is divisible by `13`,
> the full thirteen-speed family is lonely, with no upper bound on any killer.

> **First-carrier corollary.**  If `z=min(C)` and
> `dist(h/z,14Z)>=1` for every `h in C`, the family is lonely.  In particular,
> this holds if every carrier lies in `[z,13z]`.

Equivalently, the `c=1` test fails only when another carrier ratio lies in an
open band

```text
(14m-1,14m+1),                 m>=1.                    (16)
```

The endpoints are accepted.  An exact covering example far beyond the old
ratio cone is

```text
P=(1,2,4,7,9,11,12),
(k1,...,k6)=(157,160,169,196,1000,10001),
C={169}, z=169, c=1, a=1,
t=183/2366.
```

Here `k6/z=10001/169>59`, yet the minimum distance is exactly `1/14`, owned
by `169`.  The integer-shift freedom is genuine as well.  For

```text
P=(1,2,4,7,9,11,12),
(k1,...,k6)=(312,313,315,316,350,4212),
C={312,4212}, z=312,
```

the shift `c=1` fails because the second carrier has distance `1/28`, while
`c=2=floor(z/(13M))` works.  Taking `a=1` gives `t=13/168`; the minimum
distance is `1/14`, owned by speeds `12` and `4212`.

## 4. The weighted lift-aware quotient and the old cones

There is a useful one-sided quotient of Theorem B which exposes how the
earlier adaptive argument was losing wrap-around.  Put `z=min(C)`.  For each
distinct nonzero residue class `r` occupied by a killer, choose an integer
`d_r in {0,...,11}` such that

```text
k/z <= (14d_r+1)/13       for every k congruent to r mod 13.             (17)
```

If

```text
k/z<=13 for every k in C,              sum_r d_r<12,    (18)
```

then the time (12) with `c=1` is lonely.  Indeed, forbid the top `d_r`
residues `{13-d_r,...,12}` for class `r`.  At most `sum_r d_r<12`
multipliers are lost.  For an allowed multiplier and a noncarrier in class
`r`, its unreduced phase obeys

```text
1/13 < (ak mod 13)/13+k/(14z)
     <=(12-d_r)/13+(14d_r+1)/182
      =13/14.                                          (19)
```

The carrier phases lie in `[1/14,13/14]` by (18).  This proves the weighted
criterion, including `d_r=0` when a lift is small enough that no top residue
must be removed.

If `S` is the set of distinct nonzero classes, `s=|S|`, and one uses the same
`d` on every class, (17)--(18) recover the simpler condition

```text
1<=d<=11,       sd<12,       k6/z<=(14d+1)/13.          (20)
```

Thus the old clean cones remain valid:

```text
s<=5, d=2:   k6/z<=29/13;
s<=3, d=3:   k6/z<=43/13;
s<=2, d=5:   k6/z<=71/13;
s<=1, d=11:  k6/z<=155/13.                             (21)
```

The weighted quotient can leave more explicitly named multipliers available,
but Theorem B is stronger for existence: circular wrap-around caps each
individual noncarrier cost at two at arbitrary scale.

The previous boundary-bearing cone example remains useful for checking (19):

```text
P=(1,2,4,7,9,11,12),
(k1,...,k6)=(157,158,160,169,179,338),
z=169, S={1,2,4,10}, d=2, a=1,
t=183/2366.
```

Its minimum distance is exactly `1/14`, owned by `z=169`; speed `338` sits
at distance `1/7`.

## 5. What is now removed from the r=6 frontier

In a Covering family the core cannot supply a multiple of `13`, so `C` is
nonempty.

- If the first five killers are noncarriers, Theorem A gives the requested
  sharp component, and Theorem B also gives a direct witness.
- Every family with `|C|=1` is closed, at arbitrary scale.
- More generally, any `(z,c)` satisfying (10)--(11) closes the family.

Consequently a surviving counterexample must have at least two carriers and
must satisfy the following finite obstruction statement:

```text
for every z in C and 1<=c<=floor(z/(13M)),
there is h in C with dist(c h/z,14Z)<1.                 (22)
```

In particular, with `z0=min(C)` and `c=1`, some later carrier obeys

```text
h/z0 in (14m-1,14m+1) for an integer m>=1; hence h>13z0.                (23)
```

This strict ratio combines with THM-1135 to locate the first carrier.

> **First-four carrier corollary.**  In every surviving `r=6` row,
>
> ```text
> min(C) is one of k1,k2,k3,k4.                         (24)
> ```

Indeed, if `min(C)=k6`, then `C={k6}` and the unique-carrier corollary
closes the row.  If `min(C)=k5`, survival first forces a second carrier, which
can only be `k6`; then (23) gives `k6/k5>13`.  This is precisely THM-1135's
ratio horn with `s=5`, since `6(8+5)/(5+1)=13`, and again closes the row.
The endpoint issue is aligned on both sides: (23) is strict, while the horn
requires a strict ratio.

This is substantially smaller than the former residual
`z<=k5, 29z<13k6`: noncarriers no longer participate in the obstruction at
all.  THM-1135 independently bounds the mixed-scale integer box.  THM-1136
does not yet prove that the multiple-carrier ratio-avoidance problem (22) is
empty, and therefore does not by itself close all of `r=6`.

The sampled worst row of THM-1132 has five removed residues
`{2,4,6,8,10}`, all nonzero modulo `13`.  It lies in Theorem A and is now
closed uniformly by a metric component argument, not by extrapolating its
sampled component ratio.

## 6. Conflict graphs and Tournament Analysis

The useful vertices are not runners or teeth.  For Theorem A they are the
twelve observer multipliers `a in F_13^*` and the occupied nonzero residue
classes.  For Theorem B the left vertices must be the **lifted noncarrier
obligations** `(k,c k/(14z))`, because two speeds with the same residue can
translate their grids differently.  Join an obligation to the multipliers it
makes dangerous.  Every left degree is at most two, so the faithful finite
object is a sparse bipartite conflict graph

```text
lifted noncarrier obligations  <-->  observer multipliers.              (25)
```

The weighted quotient in Section 4 merges obligations with the same residue
and remembers the one-sided lift budget `d_r`; this preserves a sufficient
free-observer certificate but deliberately forgets wrap-around.  The raw
residue quotient without either a translate or a lift budget is not faithful.

The genuinely remaining object has the carriers as vertices.  For each
admissible switch `c`, put a directed obstruction edge

```text
z -[c]-> h     iff     dist(c h/z,14Z)<1.               (26)
```

A row of this labelled digraph with no outgoing obstruction edge is exactly
the carrier condition of Theorem B.  This relation is not generally a
tournament: a pair can obstruct in both directions or in neither, and the
available switch set depends on the source through `floor(z/(13M))`.
Forcing an orientation by ordinary speed order erases precisely the
Diophantine bands (16).

For the repository's Tournament Analysis audit, orient the twelve multiplier
vertices by `1<...<12`.  This gauge gives a transitive tournament: score
histogram `{0,...,11}`, no directed cycles, twelve singleton SCCs, and one
Hamiltonian path.  It is pure telemetry.  It forgets every neighborhood in
(25), every label in (26), the core cap (10), and endpoint topology.  Candidate
vertex sets explicitly challenged here were runners, combs, teeth, residues,
observer multipliers, carrier ratios, core sections, wall events, and proof
obligations.  The faithful carrier is the sparse conflict graph (25) followed
by the labelled ratio obstruction graph (26), not a naked runner tournament.

## Verification

The dependency-free referee uses `fractions.Fraction` for every metric
quantity.  It exhausts

```text
1,585 distinct nonzero residue subsets of sizes 1,...,5;
3,707 admissible (S,d) parameter rows with |S|<=5 and |S|d<12;
4,220,504 weighted residue/incidence assignments with sum d_r<12;
all 24 translated-grid boundary cells and their 24 interiors;
11 symbolic d-inequality rows.
```

The translated-grid census has maximum danger degree `2`, including every
wrap cell.  The referee checks the component example, the unique-carrier
arbitrary-scale example, the `c=2` repair example, the weighted boundary
example, the closed-tooth endpoint identity, and the conflict/tournament
fingerprints.  Ordinary and `PYTHONOPTIMIZE=1` executions are byte-identical.

```text
04-computation/lrc14_r6_q13_shadow_adaptive_referee_codex_S67.py
SHA-256 459279cb7f97ed27d1aa4390f0cd71196a2e66faeffe894656fbcc366e3c4980

05-knowledge/results/lrc14_r6_q13_shadow_adaptive_referee_codex_S67.out
SHA-256 3c936ffdcfa6bbf60bf56883501321ef478db96c239f33cddc2053c4a5a39b9e
```
