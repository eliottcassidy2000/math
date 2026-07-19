---
id: THM-1244
title: THE SLOWEST CARRIER SPOKE FORCES A POSITIONED HANDOFF DEBT — its safe component requires three blocker labels, two rational seam credits, and a macroscopic private stalk
status: PROVED (exact slowest-spoke component floor; uniform two-label exclusion; rank-two irredundant handoff forest; hybrid positioned-Hunter and spanning-tree debts; private-owner mass; dependency-free exact referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S78 continuation with debt-spoke-bridge agent
depends_on: [THM-1198, THM-1233, THM-1237, THM-1240]
related: [THM-1156, THM-1239, HYP-7870]
script: 04-computation/lrc14_slowest_spoke_component_handoff_debt_thm1244.py
output: 05-knowledge/results/lrc14_slowest_spoke_component_handoff_debt_thm1244.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCSlowestSpokeHandoffDebt.lean
script_sha256: 0044b9693b8af861a81f517ef19a83a235e15522ca635c7013c949051d68cbbe
output_sha256: c08cb1ff7e3fd8b80441eb4f1c3e43e80e2299399cdd5c8541eaeff1ddb6efe5
formalization_sha256: 36f1a179dfb24d4d10bcb0f924b0d284e3aeea8a6a9f362d60a381eb9012902e
---

# THM-1244 — slowest-spoke component handoff debt

## 1. The protected slowest-spoke component

Let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)]                      (1)
```

be a complete closed safe gap of speed `c`, and suppose the six strict
danger combs of

```text
c<d1<d2<d3<d4<d5<d6                                  (2)
```

cover `G`.  Let `t0=(2k+1)/(2c)`, put `q=c+d1`, choose an integer `p`
nearest `qt0`, and set `t1=p/q`.  This is THM-1240's centered carrier spoke.

Let `h=floor(d1t1)`.  Since

```text
delta=||d1t1||=||ct1||>1/4,                           (3)
```

the closed `d1`-safe component containing `t1` is

```text
S1=[(14h+1)/(14d1),(14h+13)/(14d1)].                  (4)
```

Define the protected component

```text
K=G intersect S1.                                     (5)
```

Then

```text
|K| >=(6d1-c)/[7d1(c+d1)]>5/(14d1).                  (6)
```

## 2. Proof of the component floor

Write `x=|t1-t0|`.  Nearest-integer selection and the centered spoke identity
give

```text
x<=1/[2(c+d1)],                 delta=1/2-cx.          (7)
```

The distance from `t1` to each endpoint of `S1` is at least

```text
r=(delta-1/14)/d1=(3/7-cx)/d1.                        (8)
```

Its distance to each endpoint of `G` is at least `3/(7c)-x`, and

```text
3/(7c)-x-r
 =(d1-c)[3/(7c)-x]/d1>0.                              (9)
```

Thus `[t1-r,t1+r]` lies in `K`.  Using (7),

```text
2r >=(6d1-c)/[7d1(c+d1)],                             (10)
```

and the last comparison in (6) clears exactly to `d1>c`.

Because `d1` is safe throughout `K`, the five remaining combs
`d2,...,d6` cover this connected closed interval.

## 3. Three blocker labels are unavoidable

Resolve that five-comb cover into its individual open teeth.  Compactness
gives a finite subcover; choose one irredundant under deletion and order its
teeth by their left endpoints:

```text
I_a=(alpha_a,beta_a),          a=1,...,N.              (11)
```

Irredundancy excludes containment, so the right endpoints increase too.
Connected coverage gives strict handoffs

```text
alpha_(a+1)<beta_a,                                   (12)
```

and every consecutive overlap in (12) lies in `int(K)`.  Since each tooth
has width less than `1/(7d1)`, (6) already implies `N>=3`.

In fact at least three distinct speed labels occur.  Suppose only two labels
`d1<u<v` covered `K`.  Every component of `D_u union D_v` meets at most one
`u`-tooth because a `v`-tooth is shorter than the safe gap between consecutive
`u`-teeth.

If `v<=6u`, the gap between consecutive `v`-teeth is at least the width of a
`u`-tooth, so at most one `v`-tooth meets it and the component span is

```text
<1/(7u)+1/(7v)<2/(7d1)<5/(14d1).                     (13)
```

If `v>6u`, all `v`-teeth attached to the unique `u`-tooth extend it by less
than one `v`-width on either side, so the span is

```text
<1/(7u)+2/(7v)<4/(21u)<5/(14d1).                     (14)
```

Both contradict the length of the connected interval `K`.  Hence the
projected handoff graph on blocker labels is connected on at least three
vertices and has graphic rank at least two.

## 4. Two located gcd quanta

Choose two chronological handoffs whose two distinct projected edges form a
forest `F0`.  If one handoff crosses from the tooth of speed `u`, address
`n`, to the tooth of speed `v`, address `m`, its overlap length is

```text
omega=[v(14n+1)-u(14m-1)]/(14uv)>0.                   (15)
```

The positive integer numerator is divisible by `gcd(u,v)`, hence

```text
omega>=gcd(u,v)/(14uv).                               (16)
```

For the two selected handoffs put `Q0=omega1+omega2` and

```text
gB=gcd(d2,d3,d4,d5,d6).                               (17)
```

Their overlap intervals lie inside `K`, so

```text
sum_(e in F0) mu(K intersect D_e)>=Q0
                                      >=gB/(7d6^2).   (18)
```

This is phase-bearing credit, not a global pair average.  A farthest-right
choice with speed/address tie-breaking makes `F0` canonical, although only
existence is needed.

## 5. Hybrid positioned-Hunter transfer

For a blocker edge `e={u,v}`, write

```text
rho_e=mu_T(D_u intersect D_v),
e_e=rho_e(1-rho_e)/gcd(u,v),
w_e(L)=L rho_e-e_e,
HB=sum_(j=2)^6 1/dj.                                  (19)
```

For every forest `F` on the five blocker labels extending `F0`, the
five-active-wall restriction of THM-1237 gives

```text
0 >=2|K|/7-6HB/49+Q0
      +sum_(e in F minus F0) w_e(|K|).                 (20)
```

Indeed forest Hunter starts from `|K|-5|K|/7=2|K|/7`;
the singleton discrepancy costs `6HB/49`; the two located edges use their
actual credit (18); and only the remaining edges pay periodic positioning
errors.  Thus the exact max-extension debt is

```text
Q0+max_(F superset F0) sum_(e in F minus F0)w_e(|K|)
 <=6HB/49-2|K|/7.                                    (21)
```

Taking `F=F0` and inserting (6), (18) yields the clean scalar consequence

```text
HB >=(6d1-c)/[3d1(c+d1)]+7gB/(6d6^2).                (22)
```

## 6. Only two ordinary tree debts remain

Extend the rank-two forest `F0` to a spanning tree `T` on the five blocker
labels.  It adds exactly two edges.  THM-1237 gives on each such edge

```text
rho_e>=1/91,                    e_e<=13/[196gcd(e)].   (23)
```

Define

```text
Gamma(F0)=min_(T superset F0 spanning)
             sum_(e in T minus F0)1/gcd(e).           (24)
```

Equations (20), (23) give the sharper protected debt

```text
HB+(13/24)Gamma(F0)
 >=(98/39)|K|+(49/6)Q0
 >=14(6d1-c)/[39d1(c+d1)]+7gB/(6d6^2).               (25)
```

The slowest spoke has therefore paid two of the four spanning-tree edges by
exact local geometry.  Only two unlocated gcd-period errors remain.

## 7. Macroscopic private stalk, microscopic seams

THM-1198 gives `d1/c<13/6`.  The normalized component coefficient

```text
phi(r)=(6r-1)/[7r(1+r)]                               (26)
```

is strictly decreasing for `r>=1`, so (6) implies

```text
|K|>phi(13/6)/c=432/(1729c).                          (27)
```

The same theorem bounds the multiply-covered part of all `G` by `4/(21c)`.
Since `d1` is inactive on `K`, more than

```text
432/(1729c)-4/(21c)=44/(741c)                         (28)
```

of `K` is uniquely covered by one blocker.  Therefore some label among
`d2,...,d6` uniquely owns more than

```text
44/(3705c).                                           (29)
```

The forced proof object is a macroscopic private-owner stalk joined by two
microscopic gcd-quantized seams.  It is not a macroscopic pair overlap.

## 8. Tournament and carrier audit

Speed order on the five blockers is the transitive tournament, with score
histogram `(0,1,2,3,4)`, no directed triangles, singleton SCCs, and one
Hamiltonian path.  It loses `K`, the chronological tooth word, both overlap
numerators, and private ownership.  The chronological handoff graph need not
be a tournament and can revisit labels; its invariant is connected graphic
rank, from which `F0` extracts two proof edges.

We challenged runners, carrier spokes, safe components, individual teeth,
wall crossings, overlap components, gcd channels, private owners, and proof
obligations as vertices.  The faithful carrier is

```text
(K; irredundant labelled tooth chain; rank-two F0;
 exact omega1,omega2; private-owner regions).          (30)
```

This also corrects the interpretation of THM-1240.  An external label may
block all six sampled spoke phases, but it cannot cover the protected
slowest-spoke component alone or with only one other label.  Continuous
handoff rank, not sampled outdegree, supplies the stable obligation.

## 9. Verification and scope

The dependency-free exact referee checks `22,400` centered spoke/component
rows, `943,200` two-label span rows, `123,962` positive gcd-quantized tooth
overlaps, all `291` forests on five labels, and all `9,312` forest/active-set
Hunter rows.  It also replays every constant in (22), (25), and (27)--(29).
Normal and optimized outputs are byte-identical.

The Lean module kernel-checks the centered radius floor, strict component
constant, both two-label span branches, positive-divisor overlap quantum,
max-extension/scalar/spanning-tree Hunter rearrangements, and private-mass
constants.  Compact irredundant subcover extraction and its chronological
handoffs remain the explicit topological paper provider; there are no proof
placeholders or `native_decide` calls.

Frozen hashes are

```text
source         0044b9693b8af861a81f517ef19a83a235e15522ca635c7013c949051d68cbbe
output         c08cb1ff7e3fd8b80441eb4f1c3e43e80e2299399cdd5c8541eaeff1ddb6efe5
formalization  36f1a179dfb24d4d10bcb0f924b0d284e3aeea8a6a9f362d60a381eb9012902e
```

THM-1244 does not yet make (21) contradictory throughout the compact ratio
box or prove LRC(14).  It materially narrows HYP-7870's remaining transfer:
two located edges and one macroscopic owner are now forced on a canonical
phase component.  The next closing target is to correlate that owner stalk
with the positive-holonomy blocker cycle or a transverse address escape.

