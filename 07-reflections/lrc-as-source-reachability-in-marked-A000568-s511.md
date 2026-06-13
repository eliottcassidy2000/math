---
source: oracle-2026-06-01-S511
status: exploratory synthesis (LRC <-> A000568, marked source-reachability)
tags:
  - lonely-runner
  - a000568
  - tournament-analysis
  - marked-tournaments
  - paley
  - addition-multiplication
---

# LRC as Source-Reachability in the Marked A000568 Quotient

The user's hypothesis: the Lonely Runner Conjecture is *completely analogous to a
problem on tournament isomorphism classes and A000568*. codex (S507/S509/S510)
established the bridge but hit a wall — **LRC safety is not a function of the
half-turn tournament class** — and retreated to a sheaf/fiber picture. This
session removes the wall by choosing the right arcs, and turns LRC into an exact
**marked-class reachability** statement whose target is literally `A000568(n−1)`.

## Why codex's lift lost the information

The half-turn comparator `i→j iff frac((s_i−s_j)t) ∈ (0,½)` meters the **½-gap**
(S26), not the LRC **1/n-gap**. Its walls are at `m/(2(s_i−s_j))`, the wrong
places. So the half-turn class genuinely cannot see loneliness — e.g. for speeds
`(2,3,5,7)` codex found a safe `t=3/28` and an unsafe `t=67/560` sharing one
regular class `H=15`. The fix is to use the **LRC walls** and **anchor on the
observer**.

## The observer-source construction (LRC walls + a mark)

Vertices: observer (speed 0) + runners `v_1,…,v_{n−1}`; threshold `1/n`.

```
observer — runner i :  observer → i   iff   ‖v_i t‖ ≥ 1/n   (runner i is FAR;
                       in the safe arc [1/n, 1−1/n]); else  i → observer
runner i — runner j :  half-turn      (i → j iff frac((v_i−v_j)t) ∈ (0,½))
```

The observer's edges flip exactly at the LRC endpoint walls
`t = (mn ± 1)/(n v_i)`. And then, **identically**,

> observer is LONELY at `t`  ⇔  observer beats every runner  ⇔  observer is a
> **SOURCE** (out-degree `n−1`).

Computed (`lrc_observer_source_tournament_s511.py`): `observer-source = LRC-safe`
with **0 mismatches** across `n=4,5,6` families. So:

> **LRC(observer)  ⇔  the observer-marked tournament clock reaches a cell whose
> marked class has the observer as a source.**

That is a pure statement about **marked tournament isomorphism classes**
(A000568 with one distinguished vertex). Loneliness *is* a marked-class property
— the thing codex found false for the half-turn class becomes true once the arcs
are the LRC arcs and the observer is marked.

What the data shows about the regimes:

```
LRC-easy (positive gap):  >=1 open SOURCE cell   (n4 (0,1,2,5): 4; n5 (0,2,3,5,7): 6)
tight:                    0 open source cells, source only on the boundary
                          (n4/n5/n6 initial; the sporadic {1,3,4,5,9}: all 0)
counterexample:           NO source cell, open OR boundary
                          = the marked walk AVOIDS every source class.
```

So a Lonely-Runner counterexample is exactly **a marked clock-walk in `G_n` that
never visits a source class** — a concrete, finite, combinatorial target.

## The target is A000568(n−1)

A source-marked class has the observer beating all runners; the only remaining
data is the runner–runner sub-tournament, an arbitrary tournament on the `n−1`
runners. Hence

> **# observer-source marked classes = A000568(n−1).**

Verified: `n=4→2, 5→4, 6→12, 7→56, 8→456, 9→6880, 10→191536` = `A000568(3..9)`.
The loneliness *target set* is, on the nose, an A000568 number. LRC says the
runner walk must hit this `A000568(n−1)`-element target (counting the boundary).
This is the cleanest sense in which "LRC is a problem on A000568": its **win
condition is sized by A000568(n−1)**, and its dynamics are a walk in the marked
A000568 quotient.

Not all `A000568(n−1)` source classes are reachable: at a source cell the runners
all lie in the safe arc `[1/n, 1−1/n]` (length `1−2/n < 1`), so their half-turn
sub-tournament is a **circular tournament of `n−1` points in a bounded arc** — the
small "circular menu" of S24. So the reachable source target is a tiny, arithmetic
subset of `A000568(n−1)`, and LRC is the claim that this subset is always
non-empty for primitive speeds.

## The two faces: addition and multiplication meet in A000568

This frames the addition/multiplication duality the user pointed at.

- **Addition** drives the *dynamics*: runner positions add (`v_i t`), the LRC
  walls are additive endpoint crossings, the staircase of differences (S25) sets
  the clock, and "observer is a source" is an additive-geometry event (everyone
  in the safe arc). This is the *walk* in the marked quotient.
- **Multiplication** sets the *menu and the base count*: A000568(n) itself comes
  from the **odd-cycle Burnside sum** (`Fix(σ)≠0` only for all-odd-cycle `σ`) —
  pure multiplicative/group structure — and the canonical regular classes are the
  **Paley (quadratic-residue) tournaments**: orienting speed differences by
  `(i−j) ∈ QR mod p` (`p ≡ 3 mod 4`) gives the regular tournament
  (`p=3,7,11 → H=3,189,95095`). The multiplicative residue structure of the
  speeds picks out *which* A000568 classes are even reachable.

So the marked A000568 quotient is exactly the place the two operations meet:
**addition moves the walk, multiplication shapes the space and the target.** The
`x+2 / x·2` grid is the same split one level down — the additive `+1/+2` chain of
differences (the staircase) versus the `·2` doubling tower (the dyadic n=16 debt,
the Cayley–Dickson ladder), with A000568's odd-Burnside core sitting on the odd
(multiplicative) axis.

## The convoluted analogy, assembled

```
runner speed set
   │  (additive dynamics: LRC walls)
   ▼
observer-marked tournament clock  ──►  closed walk in marked A000568 quotient
   │                                        │
   │  observer SOURCE ⇔ LONELY (exact)      │ menu/reachability set by
   ▼                                        ▼ multiplication (QR/Paley, odd-Burnside)
LRC  ⇔  walk hits the source target  (target size = A000568(n−1))
counterexample ⇔ walk avoids all A000568(n−1) source classes.
```

This is a faithful, *exact* (mism=0) realization of the user's hypothesis: LRC is
a reachability problem for a marked walk in the A000568 quotient, with an
A000568(n−1)-sized win set, additive dynamics, and a multiplicatively-determined
reachable menu.

## Next

1. Characterise the *reachable* source classes = circular `(n−1)`-tournaments in
   an arc of length `1−2/n`; count them vs `A000568(n−1)` (the real LRC target).
2. Translate the tight/counterexample dichotomy into a property of the marked
   walk's homology in `G_n` (extremal set = boundary-only source-touching).
3. Tie the Paley/odd-Burnside menu to the sieve (THM-369): which residue classes
   make a source class reachable.

## Artifacts

```
04-computation/lrc_observer_source_tournament_s511.py
05-knowledge/results/lrc_observer_source_tournament_s511.out
05-knowledge/results/lrc_qr_paley_speed_tournament_s511.out
```
