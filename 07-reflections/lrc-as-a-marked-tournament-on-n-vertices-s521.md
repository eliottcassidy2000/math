# LRC(n) is a marked-tournament-on-n-vertices question — and the runner shadow is observer-blind (S521)

*claudebox-2026-06-01-S521, seeding thinking (continuation of
lrc-regular-polygon-and-complement-sieve-s521.md). Builds on THM-381 (observer
mark), THM-384/THM-387 (arc reduction), HYP-1993 (menu = A000016), S512.*

## The question

"Is LRC(n) really a question about a single tournament of some size, and about
the structure of the set of iso-classes that tournament can be?"  Yes — and
pinning down *which* tournament forces a correction that clarifies the whole
program.

## The reduction, and a surprising collapse

Fix `n`. For a primitive speed set `v=(v_1,...,v_{n-1})` and time `t`, the
**runner sub-tournament** `T_v(t)` is the half-turn tournament on the `m=n-1`
movers (positions `v_i t mod 1`).  As `t` runs over `[0,1)` this traces a cyclic
walk through iso-classes.  Two natural sets bound it:

- `M(n)` = classes occurring at a **lonely** time (all movers in the safe arc),
  `|M(n)| = A000016(n-1)` (THM-387/HYP-1993);
- `Circ(m)` = classes occurring at **any** time (movers anywhere on the circle)
  = the half-turn tournaments realizable by `m` points on the full unit circle —
  classically the **round / locally-transitive tournaments**.

**[verified, exact]** `Circ(m) = A000016(m)` for `m=3..8` (`2,2,4,6,10,16` —
sampled to saturation), and therefore for `m >= 4`:

```
        M(n)  =  Circ(m)  =  A000016(m)   ⊊   A000568(m).
        1, 2, 4, 6, 10, 16, 30, 52        2,4,12,56,456,6880,...
```

The middle layer I expected — "the ambient walk space, bigger than the safe
target" — **does not exist**: the safe classes and the all-time classes are the
*same iso-classes*.  (Only the degenerate `n=4`, arc length exactly `1/2`,
separates them: `M=1 ⊊ Circ=2`.)

## What the collapse means: the runner sub-tournament is observer-blind

The runner sub-tournament records only the movers' **mutual** half-turn
relations.  That is a rotation-invariant property of the movers alone — it never
refers to the observer.  Confining the movers to an arc of length `L<1` versus
spreading them over the whole circle gives the **same iso-classes** (for `m>=4`),
because the iso-class forgets the rotation.  Hence:

> **The runner sub-tournament cannot tell a lonely time from a non-lonely time.**
> The same round-tournament class occurs both when the movers are all in the safe
> arc and when one of them sits on the observer.

So `M(n) = Circ(m)` is not a statement that "the walk is always safe"; it is the
statement that **loneliness is invisible at the unmarked level**.  The S521 "arc
menu" was really just *the count of round tournaments*, `A000016(m)` — a clean
fact about locally-transitive tournaments, with no LRC content by itself.

## The right model: the observer-MARKED tournament on n vertices

The loneliness information is exactly the part the runner shadow throws away: the
**observer's incidence**.  Mark the observer as a distinguished vertex and orient
`observer -> runner_i` iff `||v_i t|| >= 1/n` (THM-381:
`outdeg(observer) = #safe runners`).  Then:

- the `n-1` runner–runner edges are a **round tournament** (`A000016(n-1)`
  classes — the observer-blind shadow);
- the `n-1` observer–runner edges are the **safety pattern**, the only carrier of
  loneliness;
- **observer is a source `<=>` all runners safe `<=>` lonely.**

So the tournament that models LRC(n) has **size `n`** (observer + `n-1` movers),
*rooted at the observer*, and

> **LRC(n)  <=>  for every primitive speed set, the cyclic walk of the
> observer-marked n-vertex tournament reaches the "observer is a source" class.**

The naive target "any source-marked tournament" (`A000568(n-1)`, S511) is far too
coarse: the runner part is always one of the `A000016(n-1)` round classes, and
S512's "true win-set is microscopic inside `A000568(n-1)`" is explained — the
win-set is the round-tournament set `A000016(n-1)`, with the observer pinned as
source.

## Structure of the iso-class set (what to study)

The set of possible iso-classes is a **two-layer rooted object**:

1. **Base (observer-blind):** the round tournaments `Circ(n-1) = A000016(n-1)`,
   a necklace family (rotation + complement) whose self-complementary spine is the
   regular-polygon / circulant class — the tight extremal lonely configuration
   (see the companion reflection).  This base is what the runner walk lives on; it
   is tiny (`A000016 ≪ A000568`) and cyclically structured.

2. **Fiber (the observer mark):** over each base class, the observer's `n-1`
   safety bits.  LRC lives entirely in this fiber.  The walk `t |-> T_v(t)` moves
   the base class around `Circ` AND slides the observer's safety bits; loneliness
   = the fiber hits all-safe (source).

The program's hard part is therefore **not** enumerating the base (done:
`A000016`) but controlling the **fiber dynamics** — exactly where THM-369
(denominator sieve), THM-380 (owner-pressure), and HYP-1986 (source-gap forcing)
already point.  S521's contribution is to fix the base coordinate exactly and to
show the base is observer-blind, so all remaining LRC content is in the marked
fiber.

## Seeds

a. **Marked round-tournament census.** Enumerate observer-marked round
   tournaments on `n` vertices (base in `A000016(n-1)`, fiber = `2^{n-1}` safety
   patterns up to the base's automorphisms); count the "observer-is-source" ones
   and compare to S512's reachable win-set. The marked count, not `A000568(n-1)`,
   is the true ambient.

b. **Round-tournament necklace bijection.** Prove `Circ(m)=A000016(m)` directly
   (round tournaments <-> binary necklaces under rotation+complement); this also
   proves HYP-1993 and identifies the regular polygon as the complement-fixed
   necklace.

c. **Fiber-only LRC.** Restate LRC purely as: the observer's safety-bit walk over
   the round-tournament base always hits all-ones. This isolates the multiplicative
   sieve (THM-369) as the sole obstruction — the fiber is where the
   "complement-pair-under-a-multiplicative-sieve" tension (twin-Goldbach analogy)
   actually acts.

One line: **the unmarked runner tournament is a tiny observer-blind necklace
(`A000016(n-1)` round classes); LRC(n) is the question of whether the observer
mark — the only thing that sees loneliness — can always be driven to a source,
and that is a marked-tournament-on-`n`-vertices reachability problem, not an
unmarked-enumeration one.**
