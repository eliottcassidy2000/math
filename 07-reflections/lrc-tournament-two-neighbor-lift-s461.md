---
source: codex-2026-05-31-S461
status: exploratory synthesis with exact endpoint audit
tags:
  - lonely-runner
  - tournaments
  - nearest-neighbor
  - endpoint-handoff
  - zeckendorf
---

# LRC Tournament Structure and the Two-Neighbor Lift

The cleanest answer to the prompt is:

```text
Tournament structure probably enters LRC at handoffs, not at the scalar
nearest-runner minimum.
```

The usual LRC objective is

```text
max_t min_v ||v t||.
```

That collapses every time slice to one number and one blocking runner.  The
endpoint program says this is too little data.  At a dangerous endpoint, the
important event is not "who is nearest?" but:

```text
one owner speed is exactly at threshold,
some protector speed is strictly inside,
and the left/right nearest runners decide whether the boundary point survives.
```

That is a handoff.  Handoffs are the same kind of object as tournament
good-cut protection: a cut or endpoint is safe only because a directed
protector crosses it.

## Three Tournament Lifts

### 1. Round tournament from runner positions

Include the stationary runner as speed `0`.  At a generic time `t`, the runner
positions on the circle define a cyclic order.  Orient `i -> j` if `j` lies in
the forward open semicircle from `i`.

This gives a round/circular tournament, not an arbitrary tournament.  It
changes only at collision and antipodal events

```text
(v_i-v_j)t = 0 or 1/2 mod 1.
```

This lift is useful for decomposing the time circle into order cells.  Inside
one cell the adjacent-gap vector is linear.  Loneliness of the stationary
runner is then a two-sided adjacent-gap condition, not just a nearest-neighbor
condition.

### 2. Blocker-rank tournament

At a fixed time, compare nonstationary speeds by distance to `0`:

```text
i beats j if ||v_i t|| < ||v_j t||.
```

This is transitive at a single time, so it is not a deep tournament by itself.
The useful object is the time-varying schedule of rank changes.  A
counterexample would be a full cover by first-place blockers.  The endpoints
where first place changes are precisely where the second nearest runner must
be remembered.

### 3. Labelled endpoint-handoff tournament

For every endpoint owned by speed `a`, draw a directed labelled handoff

```text
a -> b
```

when speed `b` strictly protects that endpoint.  This is the closest LRC
analogue of a backward tournament arc protecting Hamiltonian-path cuts.

The speed-level collapse is not enough.  The labels matter: endpoint value,
side, owner sign, protector slack, denominator depth, and which other rows the
same protector is already paying.

The right tournament question is therefore:

```text
Can every all-protected labelled handoff core be peeled, triangularized, or
assigned positive slack divergence?
```

That is the LRC version of good-cut/SCC and endpoint-transfer rank technology.

## S461 Exact Audit

The script

```text
04-computation/lrc_tournament_two_neighbor_s461.py
```

keeps, for each forbidden endpoint, the two-sided nearest neighbors of the
stationary runner and the second distinct nearest distance.  It also builds a
speed-level handoff graph and a pairwise protection-dominance tournament
shadow.

A companion artifact already present in the workspace,

```text
04-computation/lrc_pairwise_tournament_s470.py
```

pushes the same idea away from the stationary vertex: at selected times it
builds semicircle orientations from all runner positions and blocker-pressure
graphs from deletion relief in first/second-nearest distances.  Its synthesis
agrees with the endpoint audit: second-nearest distance is the missing datum
that distinguishes a rigid two-sided moat from a fragile one-blocker moat.

Main rows:

```text
initial n=14:
  endpoints = 176, witnesses = 6
  endpoint states = 42 redundant handoffs, 128 single handoffs, 6 witnesses
  handoff SCC sizes = (12,1)
  pairwise protection majority = complete transitive tournament

n14 seven-ladder:
  endpoints = 1150, witnesses = 84
  endpoint states = 356 redundant handoffs, 710 single handoffs, 84 witnesses
  handoff SCC sizes = (13)
  pairwise protection majority = complete transitive tournament

n14 S380 gate ladder:
  endpoints = 2298, witnesses = 168
  endpoint states = 710 redundant handoffs, 1420 single handoffs, 168 witnesses
  handoff SCC sizes = (13)
  pairwise protection majority = complete transitive tournament

initial n=16:
  endpoints = 232, witnesses = 8
  endpoint states = 60 redundant handoffs, 164 single handoffs, 8 witnesses
  handoff SCC sizes = (14,1)
  pairwise protection majority = complete transitive tournament
```

The important negative lesson is that the compressed pairwise speed tournament
is too orderly: it becomes transitive in these samples.  The real obstruction
does not live in majority comparisons between speeds.  It lives in the
labelled endpoint handoff graph, where a dense SCC can still have exposed
boundary witnesses.

That also explains why "nearest only" is misleading.  In the S380 gate ladder,
the overwhelming majority of endpoints are single handoffs:

```text
owner at threshold + exactly one strict protector inside.
```

Those are fragile tournament cut-protection events.  They should be treated as
private-pivot rows or labelled backward arcs, not as scalar distance samples.

## Reframed Methodologies

### Replace scalar height by two-sided moat data

Track

```text
left nearest distance,
right nearest distance,
nearest speed,
second distinct nearest distance,
threshold owner set,
strict protector set.
```

The LRC witness condition is the two-sided statement

```text
left >= 1/n and right >= 1/n.
```

The scalar minimum forgets which side failed and whether a boundary handoff is
single or redundant.

### Use endpoint handoff cores instead of interval volume alone

Volume and maximum gap measure how much cover exists.  The endpoint handoff
core measures whether cover can close without boundary residue.  This matches
the tournament warning: count and score shadows are weak unless incidence and
cut protection survive.

### Build an LRC Omega graph

The tournament OCF graph `Omega(T)` has vertices equal to odd-cycle packets and
edges equal to incompatibility.  The LRC analogue should have vertices equal to
labelled handoff packets:

```text
(owner endpoint, protector, side, slack, denominator depth).
```

Edges record incompatibility: same speed budget, conflicting side, impossible
arithmetic slack, or exported debt collision.  Independent sets would be
compatible repair packets.  This is the natural place to import hard-core
polynomial and tournament packet technology.

### Use all pairwise runner distances as a moving round tournament

Distances between all runners are

```text
||(v_i-v_j)t||.
```

They define pair labels on the round tournament of the current circular order.
This gives a richer object than the stationary-runner star.  A proof might
look for a frame runner whose two adjacent gaps become large, then translate
back to the stationary frame by a difference-speed operation.  The danger is
that arbitrary pairwise data can forget the special stationary vertex, so this
should be used as a cell-decomposition and deletion tool, not as a replacement
for endpoint incidence.

## Zeckendorf and Non-Consecutivity

The deep Zeckendorf point is not Fibonacci numerology.  It is a local carry
normal form.

Zeckendorf digits live on a path.  If two adjacent Fibonacci digits are active,
the recurrence supplies a carry move:

```text
F_k + F_{k+1} = F_{k+2}.
```

So adjacent `1`s are not stable normal-form data.  The canonical form is the
independent-set form on the Fibonacci path.

The general version is Ostrowski numeration.  Continued-fraction denominators
satisfy

```text
q_{k+1} = a_{k+1} q_k + q_{k-1}.
```

The digit rule is local: digits are bounded by the partial quotients, and when
a digit reaches its maximum, the previous digit must vanish.  When every
`a_k=1`, this becomes Zeckendorf's no-consecutive-ones rule.

So the fundamental structure creating non-consecutivity is:

```text
rank-two recurrence + local carry confluence + canonical greedy normal form.
```

The "two strands braided like DNA" image is useful if it is made precise.  In
continued fractions, the two strands are:

```text
denominator scales q_k,
signed closest-return errors.
```

They alternate left/right around the circle.  In LRC endpoint language, the
finite analogue is:

```text
left-neighbor handoff strand,
right-neighbor handoff strand.
```

A repair on one side changes which adjacent layer can be repaired next.  If
the endpoint core quotients to a path, then adjacent repairs should carry into
the next layer instead of coexisting.  That is exactly the non-consecutivity
condition.

## Working Picture

The two-neighbor lift suggests this proof route:

```text
1. Build labelled endpoint handoff graph.
2. Peel private single-handoff rows as tournament good cuts.
3. Collapse remaining compatible handoff packets into an LRC Omega graph.
4. Look for a path/Fibonacci-cube quotient of the residual debt.
5. Use Zeckendorf/Ostrowski normal form to show the residual debt is nonzero.
```

This reframes the usual LRC attack.  The target is not just to find a large
gap or optimize `min_v ||vt||`; it is to prove that every attempted full cover
has unavoidable nonzero boundary debt after all labelled handoffs are
accounted for.
