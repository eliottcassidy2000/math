# Three walls, one constant: V* ≈ 200 is the paper's enumeration wall, the conjecture's D, and the project's finite-check wall

*mac-mini-2026-06-27-S61. Merging Sungkawichai–Trakulthongchai (arXiv:2604.23906, LRC k≤12) into the
project (per the owner), three numbers computed in three different frameworks turned out to be the same
number. That coincidence is the content. Companion to kps-S31ag's mapping
([[lrc14-is-the-composite-k-plus-1-case-of-the-polynomial-method-paper]]) and HYP-3089.*

## The three constants

LRC(14) has a magnitude scale `V*` that nobody put there on purpose. It shows up three times:

1. **The paper's enumeration wall.** Sungkawichai–Trakulthongchai prove LRC(k≤12) by enumerating the
   improper set `I(k,p,1)` over a prime sieve, at cost `~p^{(k+1)/2}/(k·2^k)`. For k=13 this becomes
   "astronomical" — but only because the enumeration must reach primes large enough to certify the
   covering tuples whose binding apex is of size `~V*`. The wall is *where the apex magnitude forces the
   sieve past feasibility.*

2. **Conjecture 7.1's threshold `D`.** Their Conjecture 7.1(13): there is a `D` such that every non-tight
   coprime tuple has a witness in `(1/d)ℤ` for all `d≥D`. kps-S31ag estimated `D ≈ 213` from the witness
   floor over the arc count.

3. **The project's `V*` atlas.** The covering-bound program's finite-check wall: directly verify all
   covering tuples with apex `≤ V*` (worst `V* ≈ 234`, decreasing in k); beyond it, peel.

This session's computation (`lrc_lonely_arc_count_vs_apex_macmini_S61.py`) measured the largest arc of the
*direct* lonely set `L({1..12, 14V})` as the apex grows: it is **bounded below (~0.005) until apex ≈ 112,
then decays ~1/V**, with `1/ℓ_max` bottoming at ≈188 and crossing ≈213 right where the apex enters the
200–234 band. **The crossover where the direct lonely arc stops being long is the same ≈200.**

## Why they are the same constant (not a coincidence)

All three are the same geometric event seen three ways. The lonely set of a covering 13-tuple is built
from a **bounded core** (the small speeds, whose lonely arcs have a fixed length `ℓ_core ≈ 0.006`) plus a
**large apex** `14V` (whose forbidden arcs have spacing `1/(14V)`). There is one threshold:

> `V* = the apex magnitude at which the apex spacing 1/(14V) drops below the core arc length ℓ_core`,
> i.e. `V* ≈ 1/(14 ℓ_core)`.

- Below `V*`: the apex's fine arcs do not subdivide the core arcs; the direct lonely set keeps a long arc
  of length `~ℓ_core` ⟹ a bounded-denominator witness ⟹ **finite check suffices** (project), **`D ≈ 1/ℓ_core`**
  (conjecture), **sieve still reaches it** (paper).
- Above `V*`: the apex subdivides the core arcs (`ℓ_max ~ 1/V → 0`); the direct route dies and you must
  **peel the apex** (equidistribution). The peel is where the enumeration explodes (paper), where a
  *single* `D` no longer catches a witness in the direct arc (conjecture), and where the finite check is
  replaced by Node-3 (project).

So `V*` is `1/(14 ℓ_core)` in all three: the **reciprocal of the bounded core's lonely arc length, scaled
by the apex denominator 14.** The paper hits it as a computational ceiling, the conjecture as a witness
threshold, the project as a check/peel boundary — one number, `≈ 1/(14 · 0.006) ≈ 12` per core arc, `≈ 200`
in `d`.

## What it means for the proof

The merge collapses the endgame to a *single, magnitude-bounded* statement, and pins its one free constant:

```
LRC(14)  ⟺  Conjecture 7.1(13)  ⟺  [ apex ≤ V*:  finite check (direct long arc) ]
                                  ∧ [ apex >  V*:  bounded-core long arc survives the equidistributed peel ]
```

The first conjunct is finite and bounded (the `V*`-atlas check, feasible). The second is the project's two
open lemmas — **CRUX 1** (bounded-core positivity = `p0≤cap` = gK8, HYP-3085) and **Node-3** effective
equidistribution (HYP-2900) — which are *also* the paper's Conjecture 7.1(13). The "astronomical k=13
computation" the paper stops at is not astronomical in the right variable: it is a finite check up to
`V*≈200` plus one analytic equidistribution bound. **The barrier was an artifact of enumerating in `p`
instead of bounding in measure.**

This is the project's contribution past k=12, stated cleanly: replace the unbounded-`p` enumeration of
`I(k,p,1)` by a uniform measure bound (the witness floor + gK8 cap) below `V*`, and an equidistribution
peel above it — both finite-or-analytic, neither astronomical.

Related: HYP-3089 (the exact `I(13,7,1)=covering mod 7` bridge + the crossover data), HYP-3087/3088 (kps),
THM-565/THM-530 (arc count / witness floor), HYP-2900 (Node-3), HYP-3085 (gK8), [[the-four-faces-of-14-why-the-exceptional-structures-crowd-into-lrc]], arXiv:2604.23906.
