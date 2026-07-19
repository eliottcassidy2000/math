---
id: THM-1181
title: The six-box standoff proposal is refuted; the exact object is a cyclic-gap polytope
status: REFUTED PROPOSAL plus PROVED counterexample.  The flow/sojourn viewpoint is valid, but B is not proved to be six boxes and nonproportional directions can hit a balanced centre exactly.  Direction (3,6,13) at u=1/4 is an explicit counterexample and is bad on an interval of width 5/637.  The uniform 2/21 ceiling and r=5 remain OPEN
source: kind-pasteur-2026-07-18-S128 (cont.78; owner: prove d ∝ (1,2,3) is the maximiser — second request)
depends_on:
  - THM-1177    # exact cyclic-gap carrier and correction of the balance sketch
  - THM-1174, THM-1173
related: [MISTAKE-171]
script: 04-computation/tube_geometry_kps_S128c78.py, six_boxes_kps_S128c78.py (+ .out)
---

# THM-1181 — the six-box standoff proposal is refuted

## 1. What survives from the proposal

Badness depends only on

```text
g(u)=(frac(-d_2u),frac(-d_3u),frac(-d_4u)) in T^3,
```

so the bad measure is indeed the sojourn time of the closed integer geodesic
`u -> g(u)` in a fixed region `B`.  The sampled tube picture is refuted by its
own distance data.  The stronger replacement “`B` is six boxes” was never
derived: sampled bad points reached sup-distance `0.0917` from the proposed
centres, whereas the volume-matched box half-width was only `0.0412`.

THM-1177 gives the exact replacement.  After sorting the three coordinates,
`g` lies in `B` if and only if the four cyclic gaps of

```text
{0,g_1,g_2,g_3}
```

are all at most `2/7`.  Thus `B` is a union of six permuted rational
polytopes, one for each cyclic order.  It is centered near the permutations
of `(1/4,1/2,3/4)`, but it is not a union of the sampled cubes.

## 2. Exact counterexample to centre-hitting rigidity

The congruences in the original argument contain independent integer lifts;
substitution does not force equality of rational speed ratios.  Take

```text
(d_2,d_3,d_4)=(3,6,13),             u=1/4.             (1)
```

Then

```text
frac(-3u)=1/4,   frac(-6u)=1/2,   frac(-13u)=3/4.      (2)
```

The direction is not proportional to `(1,2,3)`, yet its geodesic hits the
balanced centre exactly.  This refutes both the claimed iff criterion and
the proposed positive standoff for all nonproportional directions.

The counterexample is robust, not a measure-zero contact.  For

```text
12/49 <= u <= 23/91,                                  (3)
```

the cyclic order of `{0,{3u},{6u},{13u}}` is

```text
0 < 13u-3 < 6u-1 < 3u < 1,
```

with gap word

```text
13u-3,       2-7u,       1-3u,       1-3u.             (4)
```

Every entry of (4) is at most `2/7` precisely throughout (3).  Hence this
single bad interval has exact width

```text
23/91-12/49=5/637.                                    (5)
```

Reflection supplies another interval of the same width.

## 3. Consequence for the continuum frontier

The balanced centre is not arithmetically rigid; many residue patterns can
hit it through different integer lifts.  Therefore neither centre avoidance
nor distance from six sample centres can prove the `2/21` extremal bound.
The remaining statement is the full polytope-sojourn inequality from
THM-1177:

```text
measure{u: max cyclic gap of {0,{d_2u},{d_3u},{d_4u}} <=2/7}
  <=2/21.                                               (6)
```

Any proof of (6) must sum the lengths of all cyclic-order chambers; it cannot
classify only which directions pass through their centroids.  The scripts
cited above remain useful negative telemetry for the failed tube model, but
their printed centre-hitting conclusion is superseded by (1)--(5).

## Carrier/tournament audit

The six cyclic orders are directed Hamilton cycles on the four labels
`{0,d_2,d_3,d_4}`.  Their edge observables are the four directed fractional
increments and the switch is membership in `(0,2/7]`.  Completing one order
to a pair tournament loses the gap lengths; retaining the Hamilton cycle and
its four weights preserves the bad predicate exactly.  The challenged
vertex choice is again decisive: gap events, not runners alone, are the
faithful vertices.

## Honest status

The flow identity, cyclic-gap polytope, counterexample (1), and interval
(3)--(5) are exact.  The uniform bound (6), finite-scale transfer, and uniform
`r=5` theorem remain open.
