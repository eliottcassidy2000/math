---
source: oracle-2026-06-01-S544
status: reframe + verified computation (the hole is guaranteed; LRC is defect transport / coverage, not existence)
tags: [lonely-runner, reframe, pigeonhole, three-gap, defect-transport, apex, coverage, creative]
---

# Global spread guarantees local emptiness: LRC is defect TRANSPORT, not existence

**Prompt (user):** get creative about reframing the problem and the choices I made,
so that global spread can guarantee local emptiness.

Here is the reframe, and it inverts the difficulty. **Local emptiness is free.** The
hard part was never *creating* a gap — it was *steering* a guaranteed one.

## The free lunch: pigeonhole makes emptiness automatic

The `n` points (observer at `0` + runners at `v_i t`) have gaps summing to `1`, so
the **largest gap is `>= 1/n` at every instant** — a region of local emptiness
(a "hole") of width `>= 1/n` *always exists*, generated for free by the global
spread of the points. Verified (`lrc_wandering_hole_transport_s544.py`): `min_t G >=
1/n` for every system tested (AP and generic, n=5,7,14).

So the LRC quantity splits cleanly into a **guarantee** and a **transport**:

```
G(t) = largest gap among all n points   = the HOLE   >= 1/n  ALWAYS   (guarantee)
O(t) = observer collar = min_i ||v_i t|| = the COLLAR  >= 1/n = LRC    (transport)
```

`O(t) >= 1/n` means a gap of width `>= 2/n` sits *around the observer*. The hole `G`
exists for free; LRC asks the hole to **arrive at `x=0`**.

Verified split:

| system | `min_t G` (hole) | `max_t O` (collar) |
|---|---|---|
| **tight AP** (n=5,7,14) | `= 1/n` exactly | `= 1/n` exactly |
| generic | `> 1/n` | `> 1/n` (often `1/2`) |

The AP is the unique critical case where the hole is *minimum size* and *only just*
reaches the observer. Generic systems have holes to spare, transported with room.

## The choices, reconsidered (this is where the reframe lives)

Each modelling choice I made was hiding the transport view; flipping it exposes it.

- **Choice: fix the observer, move the runners.** → Flip to the **hole's frame**:
  the observer wanders relative to a (largely) standing hole; LRC = the observer
  *enters* the hole. The motion is relative; put it on the defect.
- **Choice: "find a lonely time."** → Reframe as **"transport the guaranteed defect
  to `x=0`."** Existence (`∃` empty region) is automatic; the verb is *transport*,
  not *exist*.
- **Choice: the danger arcs `B_i` (cover the time-line).** → Their complement is the
  hole's visits to the observer; the primary object is the **hole**, not the danger.
- **Choice: the half-turn tournament / the winding braid.** → The hole is the
  **apex / source-sink arc** (S530) and the **fat collar** (S541); the
  tournament/braid is just the bookkeeping of where the hole is.
- **Choice: a specific observer.** → By the center/frame-shift (S541), *which* point
  is the observer is a free choice. So LRC `=` **the wandering hole sweeps EVERY
  point** (every potential observer gets its moment) — a *coverage* statement for
  the defect's trajectory.

## Three creative reframes of "global spread => local emptiness"

1. **Defect transport (quasiparticle).** The hole is a particle of width `>= 1/n`,
   guaranteed to exist, with a worldline `θ_hole(t)`. LRC = its worldline (widened
   to `2/n`) covers `x=0`. The runners are the *medium* whose spread *generates and
   drives* the defect; the conjecture is a transport/reachability statement. The AP
   is critical confinement: the defect is minimal and barely reaches the observer.

2. **Non-adaptive linear guards.** Read the runners as guards trying to keep the
   observer *covered* (someone always within `1/n`). Their total covering budget is
   `Σ 2/n = 2 - 2/n > 1`, so by measure they *could* cover all time — but they move
   at **constant speeds (non-adaptive, linear)**. LRC says linear guards **cannot**
   cover every instant: their global (equidistributing) spread forces an uncovered
   moment. The obstruction to covering is exactly the *linearity* of the motion —
   "global spread" is the flow filling the torus, and it leaves the observer line
   un-coverable as a set.

3. **Hole-sweep coverage.** As `t` runs, the hole's location sweeps an arc of the
   circle. LRC = the sweep covers `x=0`. By frame-symmetry, LRC `=` *the hole visits
   every point* = the defect trajectory is "space-filling enough." This recasts LRC
   as an **equidistribution/coverage** of the wandering-hole process — and explains
   why the AP (most rigid, most symmetric driving) gives the *least* sweeping (the
   hole barely moves past `x=0`), hence the tight case.

## Why this is progress

The reframe **removes half the problem**: the existence of local emptiness is a
one-line pigeonhole, not a conjecture. Everything hard is concentrated in the
**transport** of the guaranteed defect to the observer — which is exactly where the
AP sits critically (`G = O = 1/n`). It also unifies the thread: the hole `=` apex
(S530) `=` fat collar (S541); transport `=` the alcove walk / affine-symmetric-group
motion (S525/S541); "the hole barely sweeps" `=` the entropy-0 critical line (S543);
"sweeps every point" `=` the covering (S526). A proof should now target the
*transport law of a width-`>= 1/n` defect under a linear flow*, with the only thing
to rule out being the AP-style critical confinement where the defect is minimal and
the sweep is degenerate.

## Open (→ HYP)
- Make the hole worldline `θ_hole(t)` and its width `w(t) >= 1/n` explicit (a
  piecewise-constant/continuous process with jumps at the gap-reorderings), and
  state LRC as: the set `{ x : some t has w(t) >= 2/n and x in the hole }` `= ` the
  whole circle. Is its complement provably empty unless the system is AP-like?
- The "non-adaptive linear guards can't cover" form: is there a clean
  flow/equidistribution lemma that a measure-`(2-2/n)` union of `n-1` linear tubes
  can never cover a single rational fiber as a SET (only up to the AP boundary)?

## Anchor
`04-computation/lrc_wandering_hole_transport_s544.py` (+ `.out`): `min_t G >= 1/n`
always (the hole is guaranteed); AP critical at `G = O = 1/n`. Builds on S530 (apex
= largest gap), S541 (collar/frame-shift), S543 (critical line), S526 (covering).
