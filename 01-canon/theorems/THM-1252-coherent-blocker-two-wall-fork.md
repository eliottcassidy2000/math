---
id: THM-1252
title: A COHERENT IRREDUNDANT BLOCKER EXPORTS A TWO-WALL FORK — two chronological lcm seams, or a finite detuned toothpick rung
status: PROVED (common-irredundant-subcover blocker gauge; all-edge target-tooth containment in the slow gap; minimum-speed fallback; ascent-only source-safe protection; two adjacent disjoint full-invoice lcm handoffs on every edge; distinct-provider rank-two fork; same-provider multiplicity-two / exact detuned toothpick rung; slowest-spoke positioned-Hunter strengthening; same-edge unconditional binary-descent coherent-cycle composition; dependency-free exact referee; paper topology with sorry-free Lean arithmetic core)
source: codex-2026-07-19-S78 mixed-circuit / coherent-word synthesis
depends_on: [THM-1233, THM-1238, THM-1240, THM-1244, THM-1248, THM-1250, THM-1253, THM-1254]
related: [THM-1156, THM-1243, THM-1246, HYP-7870]
script: 04-computation/lrc14_coherent_blocker_two_wall_fork_thm1252.py
output: 05-knowledge/results/lrc14_coherent_blocker_two_wall_fork_thm1252.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCMinimalBlockerTwoWallFork.lean
script_sha256: 2227fbdc7fc7f101313ccf4fea6c78d0643b99e99fc97359376240d9f4793680
output_sha256: fb29e9e6d9f315303a0d333f08a825ca397c669b5fb5298e56530bc46e714d47
formalization_sha256: 3f4ca6f473f5345a2f1ff3012addf53625f1d36b6dd97f4577996b67afa3d88d
---

# THM-1252 — a coherent irredundant blocker exports a two-wall fork

## 1. Put every centered blocker on the chronological word

Let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)]
```

be a complete closed safe gap of the positive integer speed `c`, and suppose
the six strict danger combs of

```text
c<d1<d2<d3<d4<d5<d6                                  (1)
```

cover `G`.  Choose one deletion-minimal finite cover of `G` by individual
open fast teeth and order it by left endpoint:

```text
I_1,...,I_A.                                           (2)
```

At the THM-1240 centered spoke `t_i` of each `d_i`, both `c` and `d_i` have
depth greater than `1/4`.  Since (2) covers `t_i`, choose one selected tooth
through it and let `b(i)` be that tooth's owner.  Then `b(i)!=i`, so `b` has
a directed cycle.  This is THM-1254's coherent gauge: every marked blocker
tooth belongs to the same irredundant chronological word.  THM-1248 permits
arbitrary strict blocker choices, so its finite relative-address word,
contracting affine transport, and gcd sheet apply without alteration.

Fix any marked blocker edge

```text
i -> j,                                                  (3)
```

and let

```text
J=((14N-1)/(14d_j),(14N+1)/(14d_j))                  (4)
```

be the selected `d_j`-tooth containing `t_i`.  The whole tooth lies strictly
inside `G`, without any speed-order assumption.  Indeed, nearest-integer
selection gives

```text
dist(t_i,boundary(G))
 >=3/(7c)-1/[2(c+d_i)]>5/(28c),
diam(J)=1/(7d_j)<1/(7c)=4/(28c).                     (5)
```

Since `t_i` lies strictly in `J`, every endpoint of `J` is less than one full
tooth width from it.  Thus

```text
closure(J) subset int(G).                            (5a)
```

There is a stronger ascent-only statement.  Put `S_i` for the closed
`d_i`-safe component through `t_i` and `K_i=G intersect S_i`.  If
`d_i<d_j`, THM-1248 gives

```text
closure(J) subset int(K_i).                          (5b)
```

Every blocker cycle contains both a speed ascent and a speed descent.  The
use of one common irredundant subcover turns (5a) into a two-wall theorem on
**every** edge; (5b) adds third-support protection on ascents.

## 2. Both walls are consecutive chronological handoffs

Write the two walls of `J` as `x_-<x_+`.  At either wall, the carrier `c` is
strictly safe and `d_j` is equality-safe.  The strict six-cover therefore
supplies labels

```text
h_-,h_+ !=d_j,
x_- in D_(h_-),                 x_+ in D_(h_+).        (6)
```

On a speed ascent, (5b) says the source is strictly safe too, so then
`h_-,h_+ notin {d_i,d_j}`.  On a descent one of the wall providers is allowed
to be the source; this does not affect either seam invoice.

Choose the provider teeth from the same selected subcover (2).  Neither can
cover all of `J`: because `J` itself is selected and lies wholly in `G`, this
would make `J` deletable, contradicting irredundancy.

More is true.  A selected provider overlaps `J` on its inward side.  By
THM-1253's minimality separation, nonconsecutive selected teeth do not
overlap.  Hence `J=I_a` is an internal tooth of (2), the left provider is
exactly `I_(a-1)`, the right provider is exactly `I_(a+1)`, and the two wall
seams are the raw chronological handoffs `W_(a-1),W_a`.  THM-1253 also proves
that all raw handoffs are pairwise disjoint.

This makes both wall intersections *proper*.  If the left provider tooth has
address `M_-`, its component inside `J` adjacent to `x_-` has length

```text
omega_-=(14M_-+1)/(14h_-)-(14N-1)/(14d_j)>0.         (7)
```

If the right provider tooth has address `M_+`, the corresponding length is

```text
omega_+=(14N+1)/(14d_j)-(14M_+-1)/(14h_+)>0.         (8)
```

Both endpoint pairs are fast-tooth walls in `int(G)`.  Their positive
integer numerators are divisible by the relevant gcd, so

```text
omega_- >=gcd(d_j,h_-)/(14d_jh_-)
         =1/[14 lcm(d_j,h_-)],
omega_+ >=gcd(d_j,h_+)/(14d_jh_+)
         =1/[14 lcm(d_j,h_+)].                        (9)
```

Thus the two raw intersections are already disjoint and both occur in
THM-1253's full chronological invoice.  Numerically, each quantum in (9) is
at most `1/(14d_j)`, consistently with

```text
1/[14 lcm(d_j,h_-)]+1/[14 lcm(d_j,h_+)]
 <=1/(7d_j)=|J|.                                     (10)
```

Thus every coherent blocker edge exports two oppositely oriented,
phase-located lcm needles.  On an ascent they are additionally protected
third-support seams, strengthening the one first-boundary seam recorded in
THM-1248.

There is an independent fallback that does not require a common subcover.
At each spoke choose the least-speed dangerous label.  If a wall-provider
tooth covered `J`, its width would force its speed below `d_j`, while it would
also contain `t_i`, contradicting that least-speed choice.  The two proper
seams still satisfy (7)--(10); if their raw intervals overlap, trim inward
quantum subneedles as in (10).  The coherent gauge is stronger because no
trimming is needed and both seams are actual entries of the same tooth word.

## 3. The exact fork/toothpick dichotomy

There are two branches.

### 3.1 Distinct wall providers

If `h_-!=h_+`, the two pair edges

```text
{d_j,h_-},                    {d_j,h_+}               (11)
```

form a rank-two star.  Both credits in (9) may therefore be entered
simultaneously in one forest-Hunter inequality.  Chronologically they are the
nonbacktracking turn

```text
h_- -> d_j -> h_+,
```

and both occurrences enter THM-1253's full lcm and weighted functional
invoices.  On the full slow gap, the two-edge forest also extends through
THM-1250's connected chronological support to a five-edge spanning tree.
When (3) is the outgoing edge of `d1`, it is automatically an ascent; on the
protected slowest-spoke component `K=K_1`, the star is already a rank-two
positioned forest admissible in THM-1244's Hunter inequality.

### 3.2 One wall label on both sides

Suppose `h_-=h_+=h`.  The two provider teeth cannot be the same tooth: one
open interval containing both `x_-` and `x_+` would contain all of `J`, which
was excluded above.  Hence they are distinct `h`-teeth and the two raw seam
intervals themselves are disjoint.  In the coherent word they form the
immediate backtrack

```text
h -> d_j -> h,
```

so THM-1253 charges both occurrences rather than collapsing them to one
unordered pair edge.

Let their addresses be `M_-<M_+` and put

```text
r=M_+-M_- >=1.                                        (12)
```

At least the first complete `h`-safe gap between these teeth lies strictly
inside `J`, immediately giving

```text
6/(7h)<1/(7d_j),                  hence h>6d_j.        (13)
```

Keeping the whole address difference gives the sharper exact functional
law.  The interval between the right wall of the left provider tooth and the
left wall of the right provider tooth has length `(7r-1)/(7h)`.  Therefore

```text
omega_-+omega_+
 =1/(7d_j)-(7r-1)/(7h)
 =[h-(7r-1)d_j]/(7d_jh).                             (14)
```

Define the positive detuning

```text
E=h-(7r-1)d_j.                                        (15)
```

Then

```text
E>0,                   gcd(d_j,h) divides E,
E>=gcd(d_j,h),
omega_-+omega_+=E/(7d_jh)
                 >=1/[7 lcm(d_j,h)].                 (16)
```

This is the `H`-drift's exact toothpick functional form.  The same-provider
failure does not erase the second seam: it converts the two-wall fork into a
multiplicity-two credit on one lcm edge, with the entire seam mass equal to
the detuning from the resonant multiplier `7r-1`.

THM-1233 makes this toothpick alphabet finite.  Since `h<=d6<2345c` and
`d_j>c`, equations (14)--(15) imply

```text
1<=r<=335.                                             (17)
```

The coarse first rung `h>6d_j` also triggers THM-1238's nonzero pair-sum beat
for `(d_j,h)`: at a separate phase in `G`, both defining labels are safe and
one of the other four labels is a literal third blocker.  Thus the toothpick
branch emits an alternate mixed-clock obligation rather than closing into a
two-label local trap.

The branch is locally sharp and cannot be deleted by a purely centered-spoke
argument.  Take

```text
(c;k;d_i,d_j,h)=(5;2;7,8,49),       t_i=1/2.          (17a)
```

The minimum dangerous speed at `t_i` among the compact packet
`(7,8,49,50,51,52)` is `d_j=8`.  Its target tooth is

```text
J=(55/112,57/112).
```

The `49`-teeth of addresses `24` and `25` contain the left and right walls,
respectively, and each seam has the exact quantum

```text
omega_-=omega_+=1/5488=1/[14 lcm(8,49)].             (17b)
```

Here `r=1`, `E=49-6*8=1`, and (14) reads
`omega_-+omega_+=1/2744`.  The row is not asserted to be a six-cover; it is
an exact guardrail showing that minimum-blocker geometry, the compact ratio
box, and the centered affine spoke are all compatible with the same-label
branch at the smallest possible detuning.

Nor does a binary THM-1248 digit improve the rung bound.  For every
`1<=r<=335`, take an even target `j`, put

```text
e=1 if r is odd,                 e=2 if r is even,
h=(7r-1)j+e,
N=j/2,          M_-=(h-r)/2,    M_+=(h+r)/2.          (17c)
```

Then `gcd(j,h)=e`, both wall seams have the sharp length
`e/(14jh)`, and the address return is exactly `r`.  Choose odd `c<j<i` and
the half-circle carrier gap.  The centered numerator of `j` has the two
nearest choices

```text
P_j=k+N                  or                  P_j=k+N+1,
```

so the relative digit is respectively `0` or `1`.  Taking
`(c,i,j)=(10001,10003,10002)` makes this a speed descent and keeps
`h<2345c` through `r=335`.  Thus even on the exact binary-descent edge used
by the coherent drift theorem, binary address data leaves the entire finite
toothpick alphabet (17) alive.

## 4. Full-word and positioned Hunter invoices

Let `(s_a)` be the owner word of the coherent irredundant cover.  THM-1253
charges every handoff, so the two wall seams appear literally in

```text
H_fast>=1/c+(7/12)sum_(a=1)^(A-1)1/lcm(s_a,s_(a+1)), (18)
```

and in the weighted functional invoice

```text
F_6>= (c/16)sum_(a=1)^(A-1)1/lcm(s_a,s_(a+1)).       (19)
```

In the distinct-provider branch, the two displayed terms are

```text
1/lcm(d_j,h_-)+1/lcm(d_j,h_+).                       (20)
```

In the same-provider branch they are two separate occurrences of
`1/lcm(d_j,h)`.  Thus no Cayley one-third averaging and no spanning-tree
selection is needed to consume the second wall.  Equations (14) and (18)
also identify the exact raw mass behind those two occurrence quanta.

Let `g_B=gcd(d2,d3,d4,d5,d6)`.  Apply the result to the minimum blocker of
the slowest spoke `d1`—or, in the coherent gauge, to its selected blocker.
All seams lie in THM-1244's protected component `K`.

If the two wall labels differ, (11) is a rank-two forest and

```text
sum_(e in (11)) mu(K intersect D_e)
 >=1/[14lcm(d_j,h_-)]+1/[14lcm(d_j,h_+)]
 >=g_B/(7d6^2).                                       (21)
```

If the labels agree, the single pair has two disjoint occurrences and

```text
mu(K intersect D_(d_j) intersect D_h)
 >=2/[14lcm(d_j,h)]>=g_B/(7d6^2).                    (22)
```

THM-1244 proves that the overlap support on `K` has graphic rank at least
two.  In the same-label branch choose any second chronological edge distinct
from `{d_j,h}`.  The resulting two-edge forest has the stronger invoice

```text
sum_(e in F) mu(K intersect D_e)>=3g_B/(14d6^2),      (23)
```

one and a half times THM-1244's generic two-quantum floor.  This local forest
statement remains useful inside the protected component, while (18)--(20)
are the stronger global consumers on the whole slow gap.

## 5. Mixed-circuit and address interpretation

The two wall determinants in (7)--(8) are invariant under integral
translation of the time lift.  They attach opposite oriented germs to the
same target address `N`, exactly where THM-1248's finite relative blocker
word previously retained only the first boundary met on the path to `t_j`.
The proof-bearing local object is therefore

```text
(centered edge i->j; target tooth (d_j,N);
 left wall (h_-,M_-,omega_-); right wall (h_+,M_+,omega_+)). (24)
```

If the wall labels differ, (24) is a typed rank-two circuit fragment.  If
they agree, the address return `M_- -> M_+` has rung `r`, positive detuning
`E`, and two disjoint seam occurrences.  This is precisely the alternative
needed by the mixed chronological/centered identity of THM-1250: a failed
single-edge germ continuation is no longer anonymous.  It either branches
to a second label or records a finite positive address return whose failure
mass is `E/(7d_jh)`.

The coherent choice gives a direct same-edge composition with THM-1254.
Choose its speed-descent cycle edge.  The marked target tooth `J` has the two
adjacent selected wall handoffs proved here.  The same edge's THM-1248 digit
is binary.  If its next marked tooth is earlier in the word, the original
coherent factor is `k+delta>=0`; if it is later, reflection gives the factor
`c-k-delta>=0`.  Consequently every putative six-cover has one marked edge
carrying simultaneously

```text
two adjacent full-word wall occurrences at its target tooth,
and R>=n_0 n_r Delta>0 on its original/reflected mixed circuit.   (25)
```

This removes the former choice incompatibility between an abstract blocker
cycle and an independently selected overlap tree.  What remains is a
quantitative comparison between the positive reciprocal-centre invoice in
(25) and the complete overlap invoice (18).

The wall providers on this descent edge may include its source; that does not
affect either occurrence credit.  A speed-ascent edge gives the stronger
source-safe statement (5b), so its providers are genuine third labels and
its fork lies in `K_i`.  The family (17c) shows that even the exact binary
descent used in (25) does not further shorten the toothpick rung alphabet.

The theorem does not yet transport a private owner germ all the way around a
blocker cycle, and therefore does not prove six-comb noncoverage or LRC(14).
It does close the local “second seam” gap: every coherent marked edge—and,
independently, every minimum-speed blocker edge—carries two oriented and
disjoint lcm credits.  Ascents additionally carry source-safe third-support
placement.

The remaining operation-level question is whether the distinct-provider
fork closes with the centered affine circuit, or whether iterating the
finite detuned rungs in (14)--(17) forces an alternate-gap descent.

## 6. Tournament and alternate-carrier audit

On runner vertices, speed order remains the transitive tournament with score
histogram `(0,1,2,3,4,5)`, no directed triangles, singleton SCCs, and one
Hamiltonian path.  It loses both walls and the detuning `E`.

On centered-spoke obligations, the observable is strict blocker incidence.
The primary switch chooses a tooth from the common irredundant word and uses
its chronological position, with speed order as the tie Hamiltonian path.
The minimum-speed switch is the independent fallback.  Either output is a
loopless functional digraph, not a tournament.  Every marked edge has a
two-wall fork; every directed cycle also has a binary speed descent carrying
the unconditional original/reflected mixed invoice.

On wall events as vertices, chronological side is the switch.  The two
oriented vertices over one target tooth either project to the rank-two star
(11), or to two parallel occurrences of one pair edge with address return
`r`.  Quotienting them to runners preserves neither disjointness nor (14).

We also challenged gaps, fixed sections, individual teeth, section
boundaries, residues, cover arcs, Fourier modes, lcm clocks, Fano lines, and
matroid circuits.  The smallest faithful carrier is (24), the entire
irredundant tooth word, its marked cycle positions, and the current slow-gap
address.  In the same-provider branch the skeleton `7r-1` is the negative
mod-seven toothpick shape, while the positive detuning `E` is its literal
metric seam mass.  A bare `chi_7` or Fano colour forgets exactly that
decisive coordinate.

## 7. Formalization scope

The dependency-free exact referee checks all `15,625` loopless blocker maps
and `409` directed cycles, `695,520` all-edge centered target-tooth margin
rows, `36,872` same-label two-wall address rows (`11,648` spanning rows
excluded by irredundancy and `25,224` proper returns), `1,205,660` compact
gcd-compatible rung rows, `216,000` two-quantum packing rows, and all `335`
binary-digit-compatible sharp rungs.  It replays the exact `(5,2,7,8,49)`
guardrail and the full-word coefficient.  Normal and optimized outputs are
byte-identical.

The Lean arithmetic core checks the all-edge target-tooth/carrier-margin
inequality (5), that two half-tooth quanta fit disjointly in one target
tooth, the same-label `h>6d_j` implication, the exact rung/drift identity
(14), the compact rung ceiling `r<=335`, the detuning gcd quantum, the
three-quantum positioned invoice (23), binary original/reflected signs, and
the terminal compact rung.  The nearest-integer spoke, strict open-cover wall
continuation, irredundant adjacency, minimum-blocker fallback, and forest
extension are explicit paper topology providers.  There are no proof
placeholders or `native_decide` calls.

Frozen artifact hashes are

```text
source         2227fbdc7fc7f101313ccf4fea6c78d0643b99e99fc97359376240d9f4793680
output         fb29e9e6d9f315303a0d333f08a825ca397c669b5fb5298e56530bc46e714d47
formalization  3f4ca6f473f5345a2f1ff3012addf53625f1d36b6dd97f4577996b67afa3d88d
```
