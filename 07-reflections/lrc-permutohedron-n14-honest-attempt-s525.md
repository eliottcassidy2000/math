---
source: oracle-2026-06-01-S525
status: reflection + framework + HONEST n=14 attempt (NOT a proof; obstruction characterized)
tags: [LRC, n14, permutohedron, braid-arrangement, alcove-walk, circle-covering, round-tournament, SCC, CRT, THM-381, HYP-2000]
---

# The permutohedron geometry of LRC, and an honest n=14 attempt

**Prompt (user):** understand the deep permutohedron geometry; attempt an n=14
proof.

**Honesty up front.** LRC for 14 runners is open (the conjecture is proved only
for <=7 runners). I did not prove it and I am not claiming to. What follows is the
permutohedron framework, a reformulation that I believe is the right one, a
*verified* structural characterization of the n=14 lonely configuration, and the
precise obstruction. The repo's n=14 history (S514, S519, opus-S524) is a chain of
honest negatives; this continues it with sharper geometry.

## 1. The permutohedron picture (synthesizing the trinity + S522-S524)

Per `the-polygon-simplex-staircase-trinity.md`: a tournament on `n` vertices is a
point in the permutohedron `Pi_n` = the zonotope of the `A_{n-1}` root system;
each arc = a root `alpha_{ij}`, the score sequence = the projection to the base
direction, vertices of `Pi_n` = the `n!` linear orders = transitive tournaments.

LRC adds the *dynamics*. For speeds `s = (s_1,...,s_m)` (here `m=13`), the runner
positions `(s_i t mod 1)` trace the **speed line** on the torus `T^m`. Two
arrangements live on `T^m`:
- the **half-turn walls** `s_i t - s_j t in Z + {0,1/2}` (where the runner-runner
  tournament changes) — the affine braid arrangement `A_{m-1}` (= the permutohedron's
  normal fan, periodised);
- the **slab walls** `||s_i t|| = 1/n` (where runner `i` enters/exits the danger
  zone around the observer).

So the LRC movie `t -> T_S(t)` (THM-381) is an **alcove walk**: the speed line
crosses walls one at a time, each crossing = an adjacent transposition = an edge
of `Pi_m`; it is a closed walk (period 1). **Transitive tournaments = vertices of
`Pi_m` = the runners packed in a semicircle (THM-374).** From S522-S524:
- the runner sub-tournament is always **round** (half-turn of circle points), a
  vanishing `A000016` slice of `A000568` (S523);
- round `<=>` `#SCC in {1, m}` (HYP-2000): a lonely config is either **transitive**
  (semicircle, a permutohedron vertex) or **one pure strong block** (wrap-around);
- the most symmetric strong block = the **regular m-gon = roots of unity** (S522).

## 2. The reformulation: LRC@14 is a circle-covering problem

Each runner `i` has measure `|B_i| = 2/n = 1/7` of blocking time
`B_i = { t : ||s_i t|| < 1/14 }` (a union of `s_i` slab-intervals). Then:

> **LRC@14  <=>  the 13 blocking arcs `B_1,...,B_13` do not cover `[0,1)`.**

Total measure `sum |B_i| = 13/7 ~ 1.857 > 1`, so covering is a priori possible;
LRC asserts the overlaps always leave a gap. The arc endpoints are exactly the
slab walls = permutohedron alcove boundaries; a **gap is a lonely alcove**, and the
**tight case is a gap of measure zero (a single wall)**.

opus-S524's CRT view is the mod-7 grouping of these arcs (6 pair-arcs + 1
singleton-arc); the "7-way correlation" is the inclusion-exclusion of the union
measure. Same problem, two coordinates.

## 3. What is verified (`lrc_n14_permutohedron_covering_s525.py`, 127 sets)

- **0 LRC failures** over the initial segment + 120 random + 6 adversarial sets.
- **The initial segment `{1,...,13}` is the UNIQUE wall-only set**: `d_open = 12/13`
  (never 13 safe in an open alcove), lonely exactly at `t = 1/14` on a wall. Every
  other set is lonely in an open alcove. (Confirms opus-S524.)
- **At every lonely time, the runner sub-tournament has `#SCC in {1, 13}` only**
  (121 single-block, 6 transitive; never `2..12`). This is forced — it is a
  corollary of THM-381 + HYP-2000 (lonely => runner sub-tournament round =>
  `#SCC in {1,m}`) — and seeing it hold across 127 n=14 sets confirms the chain.

So the n=14 lonely configuration is pinned: a **round** tournament at a
block-structure **extreme**, the tight one being the **regular-polygon wall**
(initial segment at `t=1/14`). The permutohedron alcove walk reaches loneliness
either at a vertex (semicircle/transitive) or by threading to the regular-polygon
wall (single strong block).

## 4. The attempt, and exactly where it fails

**Strategy.** Two ingredients should combine into a proof:
(i) *Extremality*: the AP `{1,...,13}` is the hardest speed set (tightest gap), and
it is lonely (at the wall `t=1/14`). (ii) *Monotonicity*: any non-AP set is
"looser" and lonely in an open alcove. If both held, LRC@14 would follow.

**Why it does not close.**
- (i) "the AP is extremal" is the **Lonely Runner extremal-configuration
  conjecture** — itself open. The data (AP uniquely wall-only) is strong evidence,
  not a proof.
- (ii) There is no monotone scalar ("gap width") that provably decreases from the
  AP to every other set; the alcove walk's gap is a non-monotone function of the
  speeds. The covering union `|∪ B_i|` is `= 1` (gap measure 0) exactly at the AP
  and `< 1` elsewhere in the data, but I cannot bound `|∪ B_i| < 1` for all
  non-AP sets — that bound *is* LRC@14.
- opus-S524's "musical chairs" handoff has the same gap: proving that consecutive
  last-blockers do not **overlap** (so a free instant exists) is exactly proving
  the union has a gap. The block-extreme structure says the free instant, if it
  exists, lands on a round/regular config — but not *that* it exists.

**Net.** The permutohedron lens converts LRC@14 into a clean statement — *13
measure-`1/7` arcs, walls on the braid arrangement, do not cover the circle* — and
localizes the entire difficulty to a measure-zero neighborhood of the
**regular-polygon wall** (the AP). That is progress in *understanding* (the
obstruction is now a single covering gap at a single geometric point), not a
proof. The honest crux: **bound the covering union below 1 for the AP-adjacent
regime**, equivalently prove AP-extremality.

## 5. A concrete sub-target that might be provable (→ HYP-2001)

Not all of LRC@14, but a real bite: **prove LRC@14 for every speed set whose
lonely alcove is OPEN** (i.e. `#SCC = 13` transitive *or* a strong block reached
in an open cell), leaving only the wall-only sets. The data says only the AP is
wall-only at n=14; if "wall-only `=>` AP (up to scaling/symmetry)" can be proved
(a finite-type statement about which speed sets force `d_open < m`), then LRC@14
reduces to the single explicit check of the AP — which is lonely at `t=1/14`.
That conditional ("wall-only => AP") is the cleanest next computation/proof.

## Anchor
`04-computation/lrc_n14_permutohedron_covering_s525.py` (+ `.out`): 0 failures/127;
AP uniquely wall-only; lonely `#SCC in {1,13}`. Builds on opus-S524 (CRT),
THM-381, HYP-1998 (round=A000016), HYP-2000 (round=block extremes), HYP-1995
(regular polygon), and `the-polygon-simplex-staircase-trinity`.
