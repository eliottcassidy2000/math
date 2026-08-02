# LRC14 reflected levels: the three-halves cone closes

**Proof candidate + frozen exact referee, 2026-08-01.**  This note is not a
canon promotion.  It gives a complete candidate proof that the sufficient
reflected `k=1` family of
[THM-2941, the critical-seven-slot scalar wall and balanced boundary](../01-canon/theorems/THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary.md)
closes throughout

```text
D >= 6,                       2m >= 3D,
m=min q_e,                    D=max q_e-min q_e.
```

Together with THM-2941's audited robust-edge-eight bank, this confines the
remaining reflected *certificate-failure* locus to the same `561` bodies and

```text
D >= 6,                       1 <= m < 3D/2.
```

This is neither a physical-survivor census nor a proof of LRC(14).

## 1. Inheritance and the new object

The closest proved mechanism is the accepted full-`C=2` cone proof in
[`lrc14-reflected-c2-full-cone-selector-proof-codex-20260801.md`](lrc14-reflected-c2-full-cone-selector-proof-codex-20260801.md).
It combines the exact primitive phase fibre

```text
F(P,Q)=1/49+c(P mod 14,Q mod 14)/(PQ),       c >= -12/49,
```

with oriented cross-determinant transport

```text
eta_g = g(Qa-Pb)/(PgL-a).
```

The hostile inherited example is the saturated zero-channel path
`32--24--36--27`: zero-channel graphs are not always matchings.  The
least-used relevant sidecar is the projective component scale, which was
present in the `C=2` proof but becomes the natural primary coordinate here.

The new cone allows

```text
max q_e/min q_e <= 5/3.
```

The right object is therefore not one uniform quality graph.  It is a
**body-specific finite rational constraint graph**: every label edge carries
the complete set of reduced ratios whose primitive phase is too small for
that edge's exact worst-corner loss.

## 2. The infinite cone has one exact corner

For labels `a<b`, write `delta=b-a`.  If `p,q` lie between `m` and `m+D`,
then

```text
|qa-pb|/(pL-a) <= B(D,m):=[m delta+D b]/[mL-b].          (1)
```

For fixed `D`, the exact finite difference is

```text
B(D,m)-B(D,m+1)
 = b(DL+delta)/[(mL-b)((m+1)L-b)] > 0.                  (2)
```

On the real boundary `m=3D/2`, put `A=(3/2)delta+b`.  Then

```text
B(D,3D/2)=DA/[(3D/2)L-b],
B(D,3D/2)-B(D+1,3(D+1)/2)
 = Ab/[((3D/2)L-b)((3(D+1)/2)L-b)] > 0.                (3)
```

Using the real boundary is important: it proves the integer condition
`2m>=3D` uniformly, including odd `D`, without pretending that `3D/2` is an
integer.  Singleton debt decreases with every level.  Hence the whole cone
is bounded by the exact corner `D=6,m=9`, with all six singleton levels set
to nine.

Across all `3,003` bodies and all label pairs, the resulting homotopy guard is

```text
|eta| <= 57/500
```

at `H=(1,2,3,4,6,12)`, labels `(1,12)`.  Every integer slope is at least
nine, so every moving slope stays at least `4443/500>1`.

## 3. Why every low edge has a finite exact alphabet

For a body `E`, ruler `L=14*lcm(E)`, and label pair `(a,b)`, define

```text
t_E(a,b)
 = 2[9(b-a)+6b]/[9L-b]
   + sum_(e in E) e/[7(9L-e)].                           (4)
```

If `t_E(a,b)>=1/49`, the coarse asymptotic phase gives no restriction.  If
`t_E(a,b)<1/49`, failure of the pair certificate forces

```text
F(P,Q) <= t_E(a,b).
```

The correction floor makes this finite:

```text
PQ <= 12/[1-49t_E(a,b)].                                (5)
```

There are `7,612` constrained body edges.  The largest bound in `(5)` is

```text
1434779927914204558290 / 149447587740929953 < 9601,
```

on body `(1,2,5,7,10,14)`, slots `(3,4)`.  The complete primitive bank under
that global bound has `1,492` unordered channels.  Thus the edge alphabets
are exact consequences of `(5)`, not a height-cutoff experiment.

Repeated levels have already been discharged by THM-2941's universal
same-level graph, and neither of its two exceptional bodies lies among the
`561`.  An uncertified word therefore has six distinct levels.

## 4. Exact projective CSP and its component theorem

For each constrained edge `ij`, require

```text
max(x_i,x_j)/min(x_i,x_j) in A_ij,
1 <= x_i <= 5/3,
x_i distinct.                                             (6)
```

Here `A_ij` is the complete phase alphabet from `(4)--(5)`.  In every
connected component, normalize one vertex to one and propagate along an
exposed edge.  A new value is its old neighbour times or divided by one of
finitely many ratios in `A_ij`; every already exposed edge is checked.  This
exhausts the component.

Only two existence bits matter: does the component have a realization of
span strictly below `5/3`, and does it have one of span exactly `5/3`?  A
short component has a nondegenerate interval of relative scales.  Its scale
can avoid every collision with previously placed components, since a
component of size `r` meets at most

```text
r(6-r) <= 9
```

old/new collision equations.  The referee uses `102` exact rational scale
candidates and then clears denominators.  Two forced full-span components
cannot coexist: both must use the global endpoints `1,5/3` and would violate
distinctness.  Conversely, at most one forced full component plus any number
of short components embeds by the finite avoidance argument.

The exact result is

```text
561 residual bodies;
540 have no projective realization of (6);
21 remain as coarse CSP traps.                            (7)
```

The search is replayed with the opposite vertex and candidate order; all
`561` verdicts agree.  The `21` are over-approximating coarse traps, not
physical packets and not asserted failures of every ignored pair edge.

## 5. Nineteen oriented located selectors

Every trap except

```text
H =(1,2,3,4,6,12),
H2=(1,3,4,6,8,12)
```

has one selected constrained pair.  Slots `(0,1)` work on all nineteen.
For a permitted unordered ratio, both orientations are retained.  This gives
`138` exact oriented controls.

For channel `(P,Q)`, let

```text
g_0=ceil(9/min(P,Q)).
```

At a selected body-safe cell, the referee computes the primitive skeleton,
the exact `eta_(g_0)`, the singleton debt at nine, and the direct reflected
overlap.  Every row first verifies

```text
skeleton-2|eta_(g_0)| > debt_9 > 0,                     (8)
```

and only then drops the favourable factor

```text
c^-1,                     c=1-a/(Pg_0L) in (0,1).
```

For all later common scales,

```text
g/(PgL-a)-(g+1)/(P(g+1)L-a)
 = a/[(PgL-a)(P(g+1)L-a)] > 0,                          (9)
```

so `|eta_g|` decreases.  Debt also decreases.  The weakest of the `138`
strict margins is

```text
112500079753279838083 / 13371962875471688608120
```

on body `(2,3,4,6,8,12)`, slots `(0,1)`, oriented channel `(13,8)`,
`g_0=2`, cell `323`.

If the actual channel is absent from the selected edge alphabet, the coarse
phase inequality already closes it.  Thus these located controls close all
nineteen bodies uniformly.

## 6. The two analytic tails

For `H`, select labels `(1,2)` and put `s=min(P,Q)`, `r=Q/P`.  On
`r in [3/5,5/3]`, convexity gives

```text
2|r-2| <= 14/5.                                          (10)
```

Using `PQ>=s^2`, `Pg>=s`, and the phase-correction floor gives

```text
1/49 - 12/(49s^2)
 - (14/5)/(168-1/s) - debt_H(9).                        (11)
```

At `s=16`, this is

```text
1372088014812922259 / 11438740346075160664000 > 0.
```

It is nonpositive at `s=15`.  Both loss terms decrease strictly thereafter;
the first increment is `3810023831/34763034555200>0`.  There are exactly
`96` oriented primitive channels with `s<16`.  Enumeration in the square of
side `2*16` agrees with an oversized square of side `4*16`, and every channel
has a located strict control.  The weakest is channel `(23,14)`, cell `155`,
with margin

```text
303724266295587869 / 256953933039218608375.
```

For `H2`, select labels `(3,4)`.  The endpoint bound is

```text
2|3r-4| <= 22/5,                                         (12)
```

and the corresponding envelope is

```text
1/49 - 12/(49s^2)
 - (22/5)/(336-3/s) - debt_H2(9).                       (13)
```

It first becomes positive at `s=7`, with value

```text
14634270155616164201 / 21275256513920852115855.
```

The complete head has `16` oriented channels.  Its weakest located row is
channel `(8,5)`, cell `311`, with margin

```text
117875406536401545877 / 16457827196385809800020.
```

This closes the final two bodies and proves the candidate theorem.

## 7. Fixed-cell 155 lemma and its exact boundary

The hostile audit suggested that the final safe singleton cell `155` of `H`
is intrinsically special.  It is: for labels `(1,2)`, every coprime oriented
primitive channel in the whole ratio interval `[3/5,5/3]` has skeleton

```text
S_155(P,Q) >= 1/126.                                     (14)
```

For `PQ>=20`, `(14)` follows immediately from

```text
F(P,Q) >= 1/49-12/(49PQ) >= 1/126.
```

The complete smaller bank consists of six orientations:

```text
(2,3):47/1008,  (3,2):1/126,
(3,4):29/1008,  (4,3):19/2016,
(3,5):1/35,     (5,3):1/105.
```

Equality is unique at `(3,2)`.  This is a genuine all-channel fixed-cell
lemma, not merely the finite scan that inspired it.

It does not by itself replace the located-head proof.  At debt nine the
coarse transport bracket at cell `155` fails on exactly

```text
(3,2), (4,3), (5,3), (5,4).
```

Their direct reflected overlaps at the least scales are positive, but no
all-scale direct monotonicity is claimed here.  The variable-cell located
proof above is the load-bearing argument.  This sharp failure anatomy says
what a future fixed-cell theorem would still have to control.

## 8. Exact referee and audit boundary

Run

```bash
python3 04-computation/lrc14_j7_reflected_three_halves_cone_closure_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_three_halves_cone_closure_thm2941.py
```

The ordinary and optimized transcripts are byte-identical.  The referee
pins the accepted `C=2` source/output/semantic image, rebuilds the `561`
body universe, freezes the complete `1,492`-channel primitive bank, checks
both CSP search orders, constructs exact component embeddings for every
coarse trap, verifies all `138+96+16` located controls against direct
reflected arcs, audits both finite heads in oversized squares, and checks the
fixed-cell lemma's complete exceptional bank.

```text
source:   15bcd2d67b4867d1ac6380a10a3b9618d677b031ad4354729371ba77501b50e6
output:   7d7ad662745d170c2492e93bf296bdd0b15372ef2ca83aed82bda41770a1089d
semantic: 6d2e760f3a929ad2eaa9af9eba11db1389ead59a383e11d8e8737c722b67aaaa
located:  b2bbe5ca1b44031030755c209f8f4c41ba40b3b7c5f0db4bd8cd75dcccb5b8a2
heads:    c59b68a9eeb4b8cefbee771827f055abc01b10ff2afb8f2226c184fa42fff6f7
```

The remaining mathematical boundary is below `3D/2`.  The decisive feature
of this proof is not a new zero channel—the zero phase set remains
`{2:3,3:4}`—but the replacement of a uniform quality tier by an exact
edge-labelled ratio alphabet plus an oriented location sidecar.
