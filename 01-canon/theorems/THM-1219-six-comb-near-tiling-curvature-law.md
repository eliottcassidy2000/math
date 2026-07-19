---
id: THM-1219
title: THE SIX-COMB NEAR-TILING CURVATURE SPECTRUM — consecutive q=1 stalks have only O(a^-2) survivor, whose mod-seven charges are exact positive active-pair slacks
status: PROVED (elementary exact interval and active-pair theorem; dependency-free Fraction referee; Lean rational arithmetic core)
source: codex-2026-07-19-S82
depends_on: [THM-1176, THM-1198, THM-1205, THM-1216]
related: [THM-1094, THM-1156, THM-1178, THM-1192, THM-1217, THM-1218, HYP-7855, HYP-7865]
script: 04-computation/lrc14_six_comb_near_tiling_curvature_referee_codex_S82.py
output: 05-knowledge/results/lrc14_six_comb_near_tiling_curvature_referee_codex_S82.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCSixCombNearTiling.lean
script_sha256: 8558229ab73b8952c879c48069a074491fdd8b4868f86800e011308dc79268ee
output_sha256: 33f81fdf369e24e02dd712a8bd2e20d3bec3ccba5e66906310f58eae9a1b6c70
formalization_sha256: 5979c911de614cc612f227c8359d9a00ee6ac41cdcc5ebcf28c345fd39ecd070
---

# THM-1219 — the six-comb near-tiling curvature spectrum

## Statement

Let `m>=6`, put

```text
a=7m,                 k=m,
b_r=a+r,              r=1,...,6,                       (1)
```

and consider the complete closed safe gap of the `a`-comb

```text
G=G_m(a)=[(14m+1)/(14a),(14m+13)/(14a)].              (2)
```

At radius `1/14`, write

```text
D_b={t:||bt||<1/14}.                                  (3)
```

Then each `D_(a+r)` meets `G` in exactly one tooth, centred at the common
integer `m+1`.  In chronological order the six teeth are disjoint and leave
exactly six safe intervals.  Five are the adjacent handoff gaps

```text
delta_r=(13-2r)/[14(a+r)(a+r+1)],    r=1,...,5,        (4)
```

and the right boundary gap is

```text
delta_0=13/[14a(a+1)].                                 (5)
```

Consequently the exact survivor length is

```text
S(a)=|G minus union_(r=1)^6 D_(a+r)|
    =13/[14a(a+1)]
      +sum_(r=1)^5 (13-2r)/[14(a+r)(a+r+1)].           (6)
```

It obeys the sharp asymptotic sandwich

```text
24/[7(a+6)^2] <= S(a) < 24/(7a^2),                    (7)
a^2 S(a) -> 24/7.                                     (8)
```

Since `|G|=6/(7a)`, the surviving fraction satisfies

```text
4a/(a+6)^2 <= S(a)/|G| < 4/a,                         (9)
S(a)/|G| ~ 4/a.                                       (10)
```

Thus these six open combs asymptotically cover the whole slow gap, while
never covering it.  In particular there is no constant `C>0` such that every
six-comb survivor in every `a`-slow gap has length at least `C/a`.  For every
multiple `a=7m>=175`, this family already has

```text
S(a)<1/(49a),                                         (11)
```

so THM-1198's universal five-comb floor cannot extend unchanged after adding
a sixth comb.

This is not an LRC counterexample.  Equations (4)--(6) exhibit six literal
safe intervals.  It is a sharp obstruction to a class of proof methods and
an exact model for the remaining six-comb geometry.

## One tooth per row

For `t` running across `G`, multiplication by `b_r=a+r` gives the endpoint
range

```text
b_r min(G)=m+(1+2r)/14+r/(14a),
b_r max(G)=m+(13+2r)/14+13r/(14a).                    (12)
```

For `1<=r<=6` and `a>=42`, the first quantity lies strictly after the tooth
centred at `m`, while the second lies strictly before the tooth centred at
`m+2`:

```text
b_r min(G)>m+1/14,
b_r max(G)<m+27/14.                                   (13)
```

The range crosses the danger interval

```text
(m+13/14,m+15/14),                                    (14)
```

so (14) is the unique dangerous tooth it meets.  Its preimage is

```text
I_r=((14(m+1)-1)/[14(a+r)],
     (14(m+1)+1)/[14(a+r)]).                           (15)
```

The threshold `m>=6` is a convenient exact uniform range.  It guarantees the
last inequality in (13) even for `r=6`; no asymptotic approximation is used.

## The curvature defects

As `r` increases, the reciprocal centre `(m+1)/(a+r)` moves left, so the
chronological order is

```text
I_6,I_5,...,I_1.                                      (16)
```

The gap between the right endpoint of `I_(r+1)` and the left endpoint of
`I_r` is exactly

```text
 [14(m+1)-1]/[14(a+r)]
-[14(m+1)+1]/[14(a+r+1)]

=((14(m+1)-1)(a+r+1)
  -(14(m+1)+1)(a+r))
  /[14(a+r)(a+r+1)]

=(13-2r)/[14(a+r)(a+r+1)],                            (17)
```

because `a=7m`.  Its numerator runs through

```text
11,9,7,5,3                                             (18)
```

for `r=1,...,5`, so every adjacent gap is positive.  Read from left to right,
these are `3,5,7,9,11`.

The left endpoint of `I_6` lies before the left endpoint of `G` by

```text
6/[14a(a+6)],                                         (19)
```

and `G`'s left endpoint lies strictly inside `I_6`.  There is therefore no
left boundary survivor.  At the other end, direct subtraction gives

```text
max(G)-max(I_1)=13/[14a(a+1)],                        (20)
```

which is the sixth safe interval.  Equations (17)--(20) prove (4)--(6), with
the open-tooth convention retaining every listed endpoint as safe.

The numerators in (18) plus the terminal numerator `13` sum to

```text
3+5+7+9+11+13=48.                                    (21)
```

Every denominator in (6) is strictly larger than `14a^2` and at most
`14(a+6)^2`.  Summing (21) proves (7).  Equations (8)--(10) follow by
squeezing.  Finally, `a>168` gives

```text
24/(7a^2)<1/(49a),                                    (22)
```

which proves (11) on the stated multiples of seven.

## The full mod-seven curvature spectrum

The phenomenon is not confined to multiples of seven.  Let

```text
a=7m+s,        0<=s<=6,        k=m,        n=m+1.    (23)
```

For every `m>=6`, the six speeds `a+r`, `1<=r<=6`, still contribute exactly
the tooth centred at `n`.  The chronological tooth order is unchanged.  The
signed left-boundary, internal-handoff, and right-boundary numerators are

```text
12s-6,
13-2s-2r,                     r=5,4,3,2,1,
13-2s.                                                   (24)
```

Write `[x]_+=max(x,0)`.  The exact survivor is therefore

```text
S_s(a)=[12s-6]_+/[14a(a+6)]
       +sum_(r=1)^5 [13-2s-2r]_+/[14(a+r)(a+r+1)]
       +(13-2s)/[14a(a+1)].                              (25)
```

Negative values in (24) are overlaps rather than holes.  The resulting
positive curvature charges and component counts are

```text
s                         0   1   2   3   4   5   6
N_s=sum positive slacks  48  42  43  46  51  58  67
# safe components          6   7   6   5   4   3   2.    (26)
```

Every positive denominator lies strictly above `14a^2` and at most
`14(a+6)^2`, so

```text
N_s/[14(a+6)^2] <= S_s(a) < N_s/(14a^2),
a^2 S_s(a) -> N_s/14.                                   (27)
```

The smallest charge occurs at `s=1`, not `s=0`.  Thus the family
`a=7m+1`, `k=m` has seven literal safe components but

```text
S_1(a)~3/a^2.                                           (28)
```

This strengthens the scale obstruction and shows that the boundary phase
participates in the second-order law; even the asymptotic coefficient is not
visible from `q=1` or the speed-order tournament.

## Active-pair slack factorization

The new active-pair formulation of THM-1205 does not merely resemble (24):
it is exactly the same arithmetic.  If two tooth owners `x<y` straddle their
pair-sum critical time with determinant `D`, their common distance is

```text
g=D/(x+y).
```

The signed handoff length between their radius-`1/14` danger intervals is

```text
(x+y)/(xy) * (g-1/14)
 = [14D-(x+y)]/[14xy].                                  (29)
```

For the five internal edges, take

```text
x=a+r,      y=a+r+1,      D=n=m+1.
```

For the right boundary use `(x,y,D)=(a,a+1,m+1)`.  For the left boundary
use `(x,y,D)=(a,a+6,a-6m)=(a,a+6,m+s)`.  Their active-pair slacks are
respectively

```text
14D-(x+y)=13-2s-2r,
14D-(x+y)=13-2s,
14D-(x+y)=12s-6,                                      (30)
```

which proves (24) and factors every term of (25).  Hence the faithful local
object is the seven-edge handoff path

```text
a_left -- b6 -- b5 -- b4 -- b3 -- b2 -- b1 -- a_right, (31)
```

with each edge labelled by `(D,x+y,xy)` and retained exactly when
`14D-(x+y)>0`.  The two appearances of `a` are boundary roles, not duplicate
runners.

This connection sharpens the global obstruction.  A positive edge is a
literal pair-sum witness for its two owners; a third runner can destroy it
only by supplying a tooth there, replacing an edge by a three-owner blocking
obligation of the type isolated in THM-1156.  What remains open is an
all-packet theorem preventing every positive active edge from being so
blocked.  A tournament on runners forgets the slack and blocker; the natural
vertices are pair-sum critical points or their blocking obligations.

## Why this is the hard stalk

This family survives the coarse geometry that precedes the exact interval
calculation:

```text
b1/a=(a+1)/a<13/6,
a sum_(r=1)^6 1/(a+r)>1.                              (32)
```

Moreover

```text
q=b6-b5=1.                                            (33)
```

Thus it realizes THM-1216's unique reduced clock with no beat residue hole:
the `q=1` beat sees only residue zero.  THM-1217 cannot manufacture a cyclic
run escape on a one-point full clock.  Nevertheless the *intervals between*
beat points retain the reciprocal-curvature defects (17).  This cleanly
separates two levels of information:

```text
beat mask at q=1                 full / no information,
continuous endpoint chronology  positive O(a^-2) active-pair slacks.       (34)
```

The missing global invariant must therefore be sensitive to second-order
endpoint drift, exact phase address, or a different beat denominator.  A
scale-free density margin alone cannot distinguish (1) from a perfect
six-interval tiling.

## Tournament and assumption audit

Taking runners as vertices and chronological tooth order as the pairwise
observable gives the transitive path

```text
b6 -> b5 -> ... -> b1,                                 (35)
```

with score histogram `(0,1,2,3,4,5)`, no directed cycles, six singleton
SCCs, and one Hamiltonian path.  This order is real but insufficient: it is
the same at every scale and forgets the vanishing gaps.

The faithful vertices are instead the seven boundary/handoff active-pair
obligations in (31), carrying their signed slacks and endpoint owners.  The
positive-part path preserves exact noncoverage and the curvature scale; the
unweighted tournament destroys them.  Its switch/gauge is the choice of
which integer each pair straddles, while ties `14D=x+y` are exact tooth
abutments.  The challenged assumption is therefore not merely that
tournament vertices should be runners; here the right vertices are the
*pair-sum defects of an almost Hamiltonian tiling*.

## Reproducibility and formal boundary

The dependency-free referee reconstructs the danger intervals directly and
independently compares their merged complement with (6) for every
`6<=m<=2000`.  It additionally checks (24)--(30) for all seven residues and
every `6<=m<=1000`, totalling 6,965 exact spectrum rows.  It checks the six
common tooth centres, component counts, active-pair factorizations, both
sandwiches, the coarse gates (32)--(33), and (11) for every tested `m>=25`.
Normal and optimized Python outputs are byte-identical.

`LRCSixCombNearTiling.lean` kernel-checks the symbolic adjacent-gap,
left-boundary, terminal-gap, and active-pair rational identities, the charge
table, and the exact arithmetic implications used in (7), (9), and (11).  The mapping
from real danger combs to the six intervals (15), their chronological order,
and the finite-union measure identity remain the explicit paper geometry.
No global six-comb or LRC(14) noncoverage conclusion is claimed.
