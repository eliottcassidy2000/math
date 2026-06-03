---
source: claudebox-2026-06-03-S593
status: REFLECTION — rigidity as a unifying lens: the two opposite senses, the fixed-point/orbit
  duality, and the LRC single-corrector/multi-sieve line as a local→global rigidity transition.
tags: [rigidity, fixed-point, reflection-symmetry, apex, orbit, cascade, deformation, LRC,
  single-corrector, multi-sieve, observer, tournaments]
---

# Rigidity: the fixed point and its orbit

**Prompt (human):** see where rigidity appears; be abstract, rigorous, free; look for local
rigidity around a fixed point, and global rigidity that cascades, permeating from isomorphic to
another in symmetric objects.

## Rigidity is one word for two opposite things

The repo says "rigid" constantly, and it means opposite things in different rooms. THM-398: a small
owner makes "a rigid pin" — a degree of freedom *forced* to one value. HYP-2130: an
"automorphism-rigid" tournament is *generic*, no symmetry. The directed cycle is "too rigid" to
change `H` (HYP-1239); the regular tournaments are "rigid classes" (HYP-665); the eigenvalue
spacing is "rigid," meaning deterministic (HYP-704); the witness matrix has a "24-entry rigidity"
(HYP-850). One word, two poles: *forced* and *generic*.

The reconciliation is the whole point. **Symmetry is the hinge that turns one rigidity into the
other.** A symmetric object is automorphism-*flexible* (it has many automorphisms) but its
distinguished problems are *pinned*: the symmetry leaves the answer no room. A generic object is
automorphism-*rigid* (trivial Aut) but its problems are *loose*: nothing forces the answer, so it
floats in an open set. Rigidity is conserved across this trade — it moves from the structure to
the dynamics. In the Lonely Runner this is exact: the symmetric configs (the arithmetic
progressions) are the tight ones, their lonely time an isolated point with no slack; the generic
configs are safe, their lonely times an open interval you can slide around in.

## Local rigidity is a fixed point; global rigidity is its orbit

The cleanest place to watch both at once is the maximin witness `t* = argmax_t min_i ‖v_i t‖`. The
runners that bind there — at distance exactly `G` — are what pin `t*`. And the reflection
`ι : x ↦ −x (mod n)` organizes them. `ι` has two fixed points, `0` and `n/2`. The second is the
**apex**, and at the tight witness the apex parks itself precisely there, at `frac = 1/2`, the
non-smooth top of `‖·‖` — the safest possible spot and, not coincidentally, its own mirror image.

That is *local rigidity around a fixed point* in the most literal sense the human asked for: the
apex is the `ι`-fixed runner, pinned at the reflection's fixed point, and it is exactly where the
single-corrector algebra dies (`0 ≠ 0`, HYP-2063). The fixed point is where the local obstruction
lives, and a local tool can sit on it.

But the *binding* danger — the runners actually at distance `G` — is never the apex (it's safe at
`1/2`). It is an `ι`-**orbit**: a pair `{v, n−v}` summing to `n`. I checked every tight family,
including the sporadic non-AP ones: `{1..6}→{3,4}`, `{1,3,4,7}→{2,3}`, `{1,3,4,5,9}→{1,5}`,
`{1..13}→{5,9}`. The witness is pinned from both sides by a runner and its mirror — one with `+v`
slope, one with `−v` — and *that* is the rigidity. It is *global* and it *cascades* exactly as the
human described: the pin at `v` is the isomorphic image, under `ι`, of the pin at `n−v`; you cannot
clear one without clearing the other, because the symmetry copies it over. Rigidity permeates from
one runner to its mirror.

## The cure has to be shaped like the orbit

Here is the part I find genuinely beautiful, because it explains not just *that* the single
corrector fails but *what replaces it*. A single corrector is a one-runner tool. It can handle the
apex — the fixed point, the orbit of size one — which is why the tight-tuple repair works by
setting the apex aside. It cannot handle the `±`-pair, the orbit of size two: clearing `v`
regenerates at `n−v`. The fix the repo found empirically (HYP-2075, "multi-sieving has no apex") is
the **pair-sum** multi-sieve — and now one sees it was forced: the obstruction is a 2-element
orbit, so the tool that resolves it must act on pairs. *The granularity of the sieve must match the
arity of the cascade.* The apex obstruction was never about the apex; it was about insisting on a
one-runner tool against a two-runner orbit.

This reframes "where and why the LRC works and doesn't" as a single statement about rigidity:

> The loneliness obstruction is an orbit of the witness symmetry. Where that orbit is a fixed point
> (the apex), a local one-runner corrector clears it — LRC easy. Where reflection symmetry makes it
> a `±`-pair, the obstruction cascades and only a pair-arity sieve, sized to the orbit, clears it —
> LRC hard. The proof must change phase exactly when the rigidity goes from local (fixed point) to
> global (orbit).

## Seeing other problems through the lens

The lens generalizes the moment one asks "what is the symmetry, and what are its fixed points and
orbits?" At `n = 18 = 2·3²` the witness symmetry is larger (units mod `n`), so the cascade is a
bigger orbit and one predicts a higher-arity sieve — matching the gate-ladder proliferation
(HYP-1992). The whole `k+1 = 2·prime` frontier is the family where the central reflection has the
apex as a zero-divisor fixed point. And on the tournament side it closes the loop with HYP-2130:
the "perspectives" count is the rigidity of the structure (vertex-orbits), and a marked vertex —
the observer — is rigid (clearable) exactly when it is alone in its orbit, blind (uncatchable by
one move) when symmetry gives it company. Local rigidity is a vertex you can name; global rigidity
is an orbit you must name all at once.

## The transcending pattern

Hardness is symmetry wearing an arithmetic mask. A problem is *locally* rigid — pinned, isolated,
hard to perturb — precisely when a *global* symmetry has cascaded its constraints into a single
fixed configuration with no slack. The apex is the fixed point you can stand on; the orbit is the
ring of mirrors around it that no single step escapes. Every tool that has worked on the Lonely
Runner — the corrector at a fixed point, the pair-sum sieve on a pair-orbit — has secretly been a
tool whose symmetry matches the obstruction's. The conjecture, in this light, is the search for a
tool whose arity dominates the largest orbit the witness symmetry can produce.

**Artifacts:** `04-computation/lrc_rigidity_local_global_s593.py` (+`.out`); new **HYP-2135**.
Builds on HYP-2130, HYP-2063 (apex), HYP-2075 (pair-sum sieve), THM-398 (pinning), HYP-2120.
