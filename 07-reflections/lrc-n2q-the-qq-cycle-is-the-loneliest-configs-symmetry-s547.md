---
source: oracle-2026-06-01-S547
status: result (verified) + honest negative — the (q,q) cycle type IS the loneliest config's maximal odd symmetry; the apex is the reflection co-observer
tags: [lonely-runner, doubled-prime, qq-cycle, apex, symmetry, regular-polygon, burnside, n2q]
---

# At n=2q the (q,q) cycle type IS the loneliest config's symmetry — and the apex is the co-observer

**Question (HYP-2044, vigorously):** at `n = 2q` (a doubled prime), does the `(q,q)`
cycle type control the loneliest LRC configuration, with apex speed `=` the repeated
cycle length `q`?

**Answer: yes, and precisely.** The `(q,q)` cycle type is a genuine **automorphism**
of the loneliest configuration; the apex (speed `q = n/2`) is the second fixed point
of the *marked* reflection; and the value `q` is shared everywhere because `q = n/2`
is the order-2 pivot. Verified for `q = 3,5,7,11` (`n = 6,10,14,22`).

## The loneliest config is the regular 2q-gon

At `n = 2q` the AP is lonely at `t* = 1/(2q)`, where the `2q` points sit at `j/(2q)`,
`j = 0..2q-1` — the **`2q`-th roots of unity**, observer at `0`, collar `= 1/n`
exactly (verified). It is the maximally symmetric configuration; its symmetry group
is the dihedral `D_{2q}`. The question is *which* part of that symmetry is the
"doubled prime."

## (A) The apex is the co-observer (the marked reflection's 2nd fixed point)

The *marked* symmetries (those fixing the observer at `0`) are the reflections
through `0`. The reflection `v <-> 2q - v` (i.e. position `j -> -j`) fixes exactly

```
{ 0, q }  =  the OBSERVER and the APEX
```

(verified: fixed points `[0,3],[0,5],[0,7],[0,11]`). So **the apex is the unique
runner fixed by the same `Z_2` reflection that fixes the observer** — the antipodal
co-observer (`q = n/2`, the order-2 element, at position `1/2`). This is *why* the
apex is special (S530), in one line.

## (B) The (q,q) cycle type IS the loneliest config's maximal-odd automorphism

The `(q,q)` cycle type is realised by **rotation-by-2** of the `2q`-gon:
`j -> j+2 (mod 2q)`. Verified for every `q`:

- its cycle type is exactly `(q,q)` — two `q`-cycles = the **even** vs **odd**
  positions (the rotation-by-2 orbits);
- **it preserves the half-turn relation** — it is a genuine *automorphism* of the
  loneliest (unmarked) tournament;
- the **observer (0) lies in the even `q`-cycle**, the **apex (`q`, odd) in the odd
  `q`-cycle**, and the **bracketing runners `±1` (odd) are in the apex's cycle**.

So **the loneliest configuration's symmetry contains a `(q,q)` automorphism**, and
`(q,q)` is the *maximal odd-order rotation* of the `2q`-gon (rotation-by-1 is a
single `2q`-cycle — even — which is a Burnside `0`). The configuration is extremal
*because* it carries this maximal symmetry; the `(q,q)` cycle type is precisely the
symmetry that makes the loneliest config tight. In that sense the `(q,q)` cycle type
**does control the loneliest configuration.**

## The Burnside tie (this is the same fact as S546)

S546: among all-odd cycle types, the equal pair `(q,q)` gives the maximal between-
cycle Burnside term `gcd(q,q) = q`. Here that very `(q,q)` is the automorphism of the
loneliest config. **The doubled-prime cycle's Burnside dominance and the loneliest
config's symmetry are one fact**: the `(q,q)` rotation is where the symmetry of the
`2q`-vertex object concentrates, both in the count (`gcd = q`) and in the LRC
extremiser (the regular `2q`-gon). Everything carries the value `q` because

```
apex speed = n/2 = the order-2 element = the (q,q) cycle length = gcd(q,q) = q.
```

"Doubling is pairing" (S546): `n = 2q` makes `q = n/2` the self-paired pivot, and it
shows up simultaneously as the apex, the reflection axis, and the `(q,q)` cycle.

## Honest negative (recorded)

One tempting form is **false**: the apex is **not** the runner the cascade traps.
Clearing runners `1..2q-1` in order, the zero conditional clearance (S545) falls on
the **last** runner `2q-1`, never the apex `q` (verified `q=3,5,7,11`). So the apex's
role is **geometric** (the reflection's co-fixed point / the antipodal diameter), not
the cascade obstruction. The `(q,q)` cycle "controls" the loneliest config as its
*symmetry*, not by making the apex a bottleneck of the clearance cascade.

## Net

> At `n = 2q`, the loneliest configuration is the regular `2q`-gon, whose maximal
> Burnside-live automorphism is exactly the doubled-prime `(q,q)` cycle (rotation-by-2,
> even/odd `q`-cycles); the apex (speed `q = n/2`) is the marked reflection's
> co-observer. The shared value `q = n/2` unifies the apex, the reflection, the cycle
> length, and the Burnside `gcd` — the self-paired pivot of the doubling. The cascade
> bottleneck is a *different* runner (the last), an honest caveat.

## Open (→ HYP)
- Does requiring the `(q,q)` automorphism + the marked reflection single out the
  regular `2q`-gon as the *unique* tight observer-source config at `n = 2q`? (A
  symmetry-forces-extremiser statement — would pin AP-tightness at the doubled primes.)
- The cascade trap being the last runner `2q-1 = -1 mod 2q`: is `2q-1` distinguished
  (the maximal speed, the reflection-image of `1`)? Relate to the apex via the
  reflection.

## Anchor
`04-computation/lrc_n2q_qq_cycle_control_s547.py` (+ `.out`): loneliest = regular
`2q`-gon; reflection fixes `{0,q}` (observer+apex); rotation-by-2 = `(q,q)`
automorphism (even/odd `q`-cycles); cascade trap = last runner (not apex). Builds on
S546 (doubled primes/Burnside), S530/S541 (apex), S522 (regular polygon), S545
(cascade).
