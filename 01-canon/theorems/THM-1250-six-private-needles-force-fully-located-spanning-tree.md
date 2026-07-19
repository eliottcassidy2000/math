---
id: THM-1250
title: SIX PRIVATE NEEDLES FORCE A FULLY LOCATED SPANNING TREE — every Hunter edge is paid by an actual tooth handoff inside the slow gap
status: PROVED (six-label essentiality; irredundant chronological chain; five gcd/lcm-quantized interior handoffs; multiplicity-averaged scale-covariant Hunter debt; uniform per-owner private interval stalk and oriented-germ dichotomy; dependency-free exact referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S78 continuation with owner-reuse and located-tree audits
depends_on: [THM-1166, THM-1178, THM-1198, THM-1233]
related: [THM-1237, THM-1240, THM-1244, THM-1248, HYP-7870]
script: 04-computation/lrc14_six_private_needles_fully_located_tree_thm1250.py
output: 05-knowledge/results/lrc14_six_private_needles_fully_located_tree_thm1250.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCSixPrivateLocatedTree.lean
script_sha256: 1973d9af1b28b053e6cd0e7b21fb059763983693da8bbe23bc40b66c380ec552
output_sha256: 454e04bfb558bad9cc935875daca9c7535d7c212cef013905faec487eaaf2a9f
formalization_sha256: e3c6de04437c51e74d5737e066828cc194f8af6515286f826eb826c3812d5181
---

# THM-1250 — six private needles force a fully located tree

## 1. Six labels are genuinely present

Let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)]                     (1)
```

be a complete closed safe gap of the positive integer speed `c`, and suppose
the six strict danger combs of

```text
c<d1<d2<d3<d4<d5<d6                                 (2)
```

cover `G`.  Put `D_i={t:||d_i t||<1/14}` and

```text
H=sum_(i=1)^6 1/d_i.                                  (3)
```

For every label `i`, THM-1198 applied after deleting `D_i` gives a survivor
of the other five combs in `G` of measure at least

```text
1/(49c).                                              (4)
```

Since all six strict combs cover, after discarding the finite set of tooth
boundaries this survivor lies strictly in `D_i` and in no other danger comb.
Thus every label has a positive private-provider region `Q_i`, with

```text
|Q_i|>=1/(49c).                                       (5)
```

In particular no label may disappear from any tooth subcover of `G`.

## 2. The chronological handoff graph has rank five

Only finitely many open teeth of the six combs meet `G`.  Choose a
deletion-minimal subcover and order its teeth by their left endpoints:

```text
I_a=(alpha_a,beta_a),                a=1,...,N.        (6)
```

Minimality excludes containment, so both endpoint sequences are strictly
increasing.  Connected coverage forces

```text
alpha_(a+1)<beta_a.                                  (7)
```

Moreover every overlap in (7) lies in `int(G)`.  If `alpha_(a+1)` were at or
to the left of the left endpoint of `G`, the preceding tooth would be
redundant because `beta_a<beta_(a+1)`; the right-end argument is symmetric.

Consecutive teeth cannot have the same speed, because distinct teeth of one
danger comb are disjoint.  Record the speed labels along (6), and join two
labels whenever they occur consecutively.  This projected transition graph
is connected: its chronological word visits all six labels by (5), and every
newly visited label is joined to its predecessor.  Hence it contains a
spanning tree

```text
T on {d1,...,d6},                   |T|=5.             (8)
```

For each edge of `T`, retain one chronological tooth handoff realizing it.
All five edges are actual overlaps inside the selected carrier gap; none is
a global pair average or an independently positioned period.

## 3. Every tree edge has a gcd/lcm quantum

Suppose one retained handoff runs from the tooth of speed `u` and address
`n` to the tooth of speed `v` and address `m`.  In chronological order its
overlap length is

```text
omega_uv
 =[v(14n+1)-u(14m-1)]/(14uv)>0.                      (9)
```

The numerator is a positive integer divisible by `gcd(u,v)`.  Therefore

```text
mu(G intersect D_u intersect D_v)
 >=omega_uv
 >=gcd(u,v)/(14uv)
 =1/[14 lcm(u,v)].                                   (10)
```

This is the improvement over THM-1178's phase-free seam envelope.  An
arbitrary component of `(G intersect D_u) intersect (G intersect D_v)` can
have a carrier endpoint owned by `c`, so THM-1178 safely kept only
`1/(14uv)`.  The irredundant *individual-tooth* chain makes both selected
endpoints fast-tooth endpoints and promotes the numerator by `gcd(u,v)`.

## 4. Fully located Hunter debt

Let `C(t)` be the number of active fast danger combs.  At every phase, the
forest induced by the active labels has at most `max(C(t)-1,0)` edges.  Thus
the ordinary forest-Hunter inequality on `G` gives

```text
0 >= |G|-sum_i mu(G intersect D_i)
       +sum_(uv in T)mu(G intersect D_u intersect D_v). (11)
```

The sharp singleton discrepancy from THM-1166 is

```text
mu(G intersect D_i)<=|G|/7+6/(49d_i),                (12)
```

and `|G|=6/(7c)`.  Substituting (10)--(12) yields

```text
0 >=6/(49c)-6H/49+sum_(uv in T)omega_uv.             (13)
```

Consequently

```text
H >=1/c+(49/6)sum_(uv in T)omega_uv
  >=1/c+(7/12)sum_(uv in T)1/lcm(u,v).               (14)
```

Equivalently, for the harmonic slack `delta=cH-1`,

```text
delta >=(7c/12)sum_(uv in T)1/lcm(u,v)>0.            (15)
```

All five tree credits are located.  There are no pair-period positioning
errors left in (13)--(15).  The correction is also invariant under a common
dilation of `c,d1,...,d6`, because `c/lcm(u,v)` is unchanged.

Two useful consumers are immediate.  Every selected edge satisfies

```text
lcm(u,v)>=7c/(12delta),                               (16)
```

so a low-slack cover carries a connected tree of uniformly large exact
clocks.  If `g0=gcd(d1,...,d6)`, then

```text
delta>=35c g0/(12d6^2).                              (17)
```

The tree in (14) is the one realized by chronology, not an arbitrarily
optimized numerical tree.  This distinction preserves phase truth.

### Repeated owners strengthen the debt

Let `n_i` be the number of selected teeth carrying label `d_i`.  Those teeth
cover `Q_i`, while one complete `d_i`-tooth has length `1/(7d_i)`.  Hence

```text
n_i>=ceil(d_i/(7c)).                                  (17a)
```

For an unordered pair `{u,v}`, let `m_uv` count its consecutive occurrences
in the chronological word.  The corresponding handoff intervals are
pairwise disjoint.  Two occurrences either use different teeth of both
speeds, or share one middle tooth in a word `u-v-u`; in the latter case their
two `u`-teeth are disjoint.  Therefore

```text
mu(G intersect D_u intersect D_v)
 >=m_uv/[14 lcm(u,v)].                                (17b)
```

Give edge `uv` of `K6` weight `m_uv/lcm(u,v)`.  Every edge occurs in exactly
`432` of the `6^4=1296` labelled spanning trees, so the average tree weight
is one third of the total.  A maximum-weight tree can be chosen inside the
connected positive support.  Applying Hunter once with that fixed tree gives
the stronger owner-reuse inequality

```text
H >=1/c+(7/36)sum_(u<v)m_uv/lcm(u,v).                 (17c)
```

Since

```text
sum_(u<v)m_uv=N-1,
N>=sum_i ceil(d_i/(7c)),                              (17d)
```

and `1/lcm(u,v)>=g0/d6^2`, (17c) has the scalar consumer

```text
H >=1/c+[7g0/(36d6^2)]
          [sum_i ceil(d_i/(7c))-1].                  (17e)
```

Unlike a sampled blocker outdegree, (17c) becomes stronger when an external
label repeatedly interrupts the owner word.  Reuse is a charge, not a way to
erase obligations for free.

## 5. A literal private interval stalk for every owner

THM-1233 gives `d_i/c<2345`.  A tooth of speed `d_i` meeting `G` has its
integer centre address in an interval of length

```text
d_i|G|+1/7 <(6*2345+1)/7=2010+1/7.                  (18)
```

Hence at most `2011` teeth of each label meet `G`, and at most `12066` fast
teeth occur altogether.  Pigeonholing (5), every label has one tooth carrying
private mass at least

```text
1/(49*2011*c)=1/(98539c).                            (19)
```

The other five labels contribute at most `5*2011=10055` teeth.  Their
endpoints split the selected owner tooth into at most `10056` private
components.  Therefore every one of the six labels owns a literal interval
stalk of length at least

```text
1/(98539*10056*c)=1/(990908184c).                    (20)
```

At least one endpoint of each such stalk is internal unless the stalk covers
all of `G`, which is impossible because the other five private regions are
nonempty.  Thus (20) supplies an oriented owner germ at an actual wall event,
not merely positive measure spread among unboundedly many cells.

There is also an exact two-branch law for any private component `J` whose two
endpoints lie in `int(G)`.  At its left endpoint some label `u` exits; at its
right endpoint some label `v` enters.  The private owner `b` remains strictly
inside one tooth across both endpoints—otherwise openness of the covering
label would intrude into `J`.  If `u!=v`, then

```text
|J|>=gcd(u,v)/(14uv)=1/[14 lcm(u,v)],                 (20a)
```

and the adjacent `u-b` and `b-v` overlaps form a rank-two V rooted on the
same private stalk.  If `u=v`, the two walls are consecutive teeth of `u`,
so `J` is its complete safe gap and lies strictly inside the `b`-tooth:

```text
|J|=6/(7u)<1/(7b),                 hence u>6b.         (20b)
```

Thus every fully interior private germ either exports two differently
labelled seams plus an lcm-quantized intervening cell, or enters the precise
toothpick/self-similar ratio branch.  At most two owner labels can have all
their private components touching the two endpoints of `G`; at least four
labels admit this dichotomy.  In addition, (17a) guarantees three owner
teeth—and hence an interior owner address—as soon as `d_i>14c`.

## 6. Relation to the blocker cycle and tournaments

The speed-order tournament remains transitive, with score histogram
`(0,1,2,3,4,5)`, no directed triangles, singleton SCCs, and one Hamiltonian
path.  It loses all five handoff positions.  A more relevant pairwise
observable is

```text
W_uv=mu(G intersect D_u intersect D_v),               (21)
```

switched at `W_uv>0`; its positive support is connected and has graphic rank
five.  Completing missing pairs by speed is only a tournament gauge and adds
no proof content.  The chronological transition graph can revisit labels and
is properly treated as a weighted handoff graph rather than forced into a
tournament.

We challenged runners, gaps, fixed sections, individual teeth, tooth
boundaries, private components, overlap components, lcm clocks, centered
spokes, and proof obligations as vertices.  The faithful carrier produced
here is

```text
(G; six private interval stalks; irredundant tooth word;
 rank-five T; five addresses and omega_uv; lcm clock on each edge). (22)
```

THM-1240/1248's centered blocker holonomy instead uses the affine clocks
`Q_i=c+d_i`.  Formula (15) lives on the speed clocks `d_i`.  There is no
general identity transporting `gcd(d_i,d_j)` to
`gcd(c+d_i,c+d_j)`.  This **affine c-drift** is now the exact remaining
coupling obstruction.  The tree is acyclic and cannot manufacture holonomy
by itself; it must be glued to the centered directed cycle through the
oriented wall germs in (20).

## 7. Verification and scope

The dependency-free exact referee checks `987,496` positive tooth-to-tooth
handoff numerators, all `1,296` labelled six-vertex trees against all `64`
active sets (`82,944` Hunter checks), and `112,320` surjective chronological
owner words with repeats.  It verifies the `432/1296` Cayley incidence,
all `112,320` multiplicity-weighted maximum-tree inequalities, `196,950`
private-owner recurrence rows, the tooth ceilings, and every scalar constant.
Normal and optimized outputs are byte-identical.

The Lean module kernel-checks the positive-divisor quantum, natural-number
gcd/lcm identity, fully located Hunter and multiplicity-averaged lcm debts,
scale-covariant form, low-slack clock consumer, tooth-address ceiling, both
private pigeonholes, owner recurrence, the same-label toothpick branch, and
common-dilate invariance.  Compact irredundant subcover extraction and
the interval chronology remain explicit paper topology providers; there are
no proof placeholders or `native_decide` calls.

Frozen hashes are

```text
source         1973d9af1b28b053e6cd0e7b21fb059763983693da8bbe23bc40b66c380ec552
output         454e04bfb558bad9cc935875daca9c7535d7c212cef013905faec487eaaf2a9f
formalization  e3c6de04437c51e74d5737e066828cc194f8af6515286f826eb826c3812d5181
```

THM-1250 does not exclude six-comb covers or prove LRC(14).  It replaces the
two located edges available on the protected slowest-spoke component by a
fully located five-edge tree on the complete slow gap, charges every repeated
owner handoff, and extracts a uniform private interval germ for every label.
The highest-leverage next step is to transport one such germ from the `d_i`
lcm tree through the `c+d_i` relative-address/positive-holonomy cycle, or
charge the first failed transport as another occurrence in (17c).
