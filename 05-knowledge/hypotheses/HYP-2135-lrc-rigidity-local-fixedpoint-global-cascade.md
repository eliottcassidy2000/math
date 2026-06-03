---
id: HYP-2135
status: EXPLORATION + EVIDENCE — rigidity as one lens: the LRC witness is pinned by a reflection
  ±-orbit of binding runners symmetric about the apex (the fixed point). Local rigidity = the apex
  fixed point; global rigidity = the ±-orbit cascade. The pair-sum multi-sieve clears the ±-orbit
  it cannot clear with one corrector. Conceptual synthesis; the orbit-structure claim is verified
  on tight families, the proof bridge is a conjecture.
source: claudebox-2026-06-03-S593
related:
  - HYP-2130
  - HYP-2063
  - HYP-2075
  - HYP-2120
  - THM-398
---

# HYP-2135: rigidity is local (apex fixed point) + global (the ±-orbit cascade)

A unifying reading of "rigidity," which the repo uses in two opposite senses, and what it says
about where the Lonely Runner argument works and fails.

## The word points two ways; symmetry reconciles them

- **Pinning rigidity** (THM-398: a small owner pins the endpoint to the arc centre, no slack; the
  tight witness is an isolated lonely time): a degree of freedom is *forced*.
- **Automorphism rigidity** (HYP-2130: trivial Aut, a *generic* tournament): no symmetry.

These are opposite — one means "forced," the other "generic." The reconciler is that **symmetry
converts one into the other**: a symmetric config is automorphism-*flexible* but witness-*rigid*
(its lonely time is pinned/tight); a generic config is automorphism-*rigid* but witness-*flexible*
(an open set of lonely times). Verified: the tight families ({1..6}, {1..7}, {1..13}) have positive
3-term/symmetry and `G=δ` (rigid witness); the asymmetric/translated/Sidon families have `3term≈0`
and `G≫δ` (flexible witness). The two rigidities are anti-correlated because symmetry is their
common hinge.

## The geometry: a reflection, its fixed point, and its orbits

Run the maximin witness `t* = argmax_t min_i ‖v_i t‖`. The **binding runners** (those at distance
`=G`) pin `t*`. Across *every* tight config tested — including the **sporadic** (non-AP) ones — the
binding pins are a **±-reflection pair `{v, n−v}`** (summing to `n`):

```
{1..6}    (n=7)  binders {3,4}   3+4=7      {1,3,4,7}   (n=5) binders {2,3}  2+3=5
{1..13}   (n=14) binders {5,9}   5+9=14     {1,3,4,5,9} (n=6) binders {1,5}  1+5=6
```

The relevant symmetry is the reflection `ι : x ↦ −x (mod n)`. It has two fixed points, `0` and
**`n/2` — the apex** (HYP-2063: at the tight witness the apex sits at `frac = 1/2`, distance `1/2`,
the safest spot; it is its own mirror image). So:

- **Local rigidity around a fixed point** = the apex. It is the `ι`-fixed runner, parked at the
  non-smooth maximum `1/2`. The single corrector's algebra trips exactly here (`0 ≠ 0`, HYP-2063) —
  the fixed point is where the *local* obstruction concentrates.
- **Global rigidity that cascades** = the binding pair `{v, n−v}`, an `ι`-orbit of size 2. The pin
  at runner `v` is the mirror image — the *isomorphic copy* — of the pin at `n−v`; they bind `t*`
  from opposite sides (one `+v` slope, one `−v` slope), and **cannot be cleared independently**
  because `ι` ties them. This is the rigidity that "permeates from isomorphic to another."

Safe configs differ: their binding set is not closed under the *central* reflection (at most under
an off-centre affine reflection `x↦−x+b`), so the pins are generic and pickable one at a time.

## Why this is exactly the single-corrector / multi-sieve line

- A **single corrector** (one modulus, Sungkawichai–Trakulthongchai Prop 4.1) acts on one runner at
  a time. It can pin/handle the **apex fixed point** (set it aside, as the tight-tuple repair does),
  but it **cannot clear an `ι`-orbit** `{v, n−v}`: clearing `v` does not clear `n−v`, and the
  symmetry regenerates the obstruction. This is the LRC-*hard* (witness-rigid, symmetric) regime.
- The **pair-sum multi-sieve** (HYP-2075: "multi-sieving has no apex") uses **pair** moduli — it
  acts on `{v, n−v}` *as a unit*, matching the orbit. That is *why* it dissolves the apex
  obstruction: the obstruction was never a single runner, it was an `ι`-orbit, and a pair-based
  sieve is exactly the tool sized to a 2-element orbit.

> **Claim.** "Where the LRC works" = where the loneliness obstruction is a single local pin (the
> apex fixed point), clearable by one corrector. "Where it doesn't" = where reflection symmetry
> makes the binding an `ι`-orbit `{v, n−v}` that cascades — clearable only by a method whose
> granularity matches the orbit (the pair-sum sieve). The transition is the local→global rigidity
> transition: from a fixed point to its orbit.

This refines HYP-2130 (rigidity = observer-coupling; single vs multi corrector) with the precise
symmetry (`ι : x↦−x mod n`), its fixed point (apex), and its orbit (the binding ±-pair), and it
*explains the form of the cure*: the multi-sieve is pair-based because the cascade is a 2-orbit.

## Honest nuances / open

- The binding pins are ±-paired even for sporadic tight sets whose *full* speed set is not
  `ι`-symmetric — so the relevant symmetry is of the *binding subset / the witness*, not always the
  whole config. Characterise exactly which symmetry of the witness is forced at tightness.
- Larger orbits: at `n` with more prime structure (n=18=2·3²), the relevant symmetry group is
  bigger (units mod n), so the cascade is a larger orbit needing a higher-arity sieve. Predict the
  sieve arity from the orbit size of the witness symmetry — test against the n=18 gate ladder
  (HYP-1992) and the `k+1=2·prime` frontier.
- Make "the multi-sieve granularity = the witness-symmetry orbit size" a theorem.

**Artifacts:** `04-computation/lrc_rigidity_local_global_s593.py` (+`.out`),
`07-reflections/lrc-rigidity-the-fixed-point-and-its-orbit-s593.md`. Builds on HYP-2130
(perspectives = rigidity), HYP-2063 (apex), HYP-2075 (pair-sum multi-sieve), THM-398 (pinning),
HYP-2120 (circuit-free/3-term).
