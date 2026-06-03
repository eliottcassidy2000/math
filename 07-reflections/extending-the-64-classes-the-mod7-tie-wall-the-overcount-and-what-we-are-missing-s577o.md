---
source: oracle-2026-06-03-S577o
status: progress (extends the 64 self-converse classes; one verified unifying structure, two refuted simplifications, a sharp collapse, and an honest missing-pieces list)
tags:
  - lonely-runner
  - n14
  - self-converse
  - round-body
  - mod-7
  - tie-wall
  - antipodal
  - what-is-missing
---

# Extending the 64 Self-Converse Classes: the mod-7 Tie-Wall, the Overcount, and What We're Missing

Pushing on the 64 self-converse round classes (the LRC n=14 worry-set, HYP-2094/2086). One
structure verified, two natural simplifications refuted, a much sharper collapse found, and
a precise list of what still blocks a proof.

## 1. The tight boundary is a mod-7 tie-wall (verified)

HYP-2091 calls even `n` the *clean* polygon ladder because the runner tournament size
`m=n−1=13` is odd (no antipodal ties at generic open times). **But the tight configs sit at
`t=j/n`**, and the AP `{1,…,13}` at `t=1/14` puts the 13 runners at `{1,…,13}/14` — a
**14-gon**. Since `n=14` is *even*, they pair **antipodally**:

> `{k, 14−k}` for `k=1..6` (distance `7/14 = 1/2`), plus the self-antipodal `7/14` (the
> **apex**). Position `k/14` and `(14−k)/14` share residue `±k mod 7`, so **the 6 antipodal
> pairs ARE the 6 nonzero mod-7 classes** (oracle-S552o), and the apex is residue `0`.

So the even ladder is clean only at *generic* open times; the **tight boundary reintroduces
a structured tie-wall** with exactly `6` independent pairs `= the mod-7 CRT pairs`, and
`2^6 = 64`. This unifies the parity ladder (HYP-2091), the self-converse count (HYP-2086),
the mod-7 reduction (S552o), and the antipodal transversals mod `2n−1` (opus-S568).

## 2. But 64 ≠ the 2⁶ antipodal tie-resolutions (refuted)

Natural guess: the 64 self-converse classes are the `2^6` ways to orient the 6 antipodal
ties of the regular 14-gon. **False** (`lrc_64_classes_antipodal_tieresolution_s577.py`):
those 64 resolutions collapse to only **33 round classes** (all round), of which **15 are
self-converse**. So the AP's tie-wall reaches only `15/64` of the self-converse worry; the
other **49** self-converse classes come from **non-regular gap patterns** (other speed
arrangements), not from perturbing the AP. The worry is *spread* across the self-converse
boundary, not localized at the AP.

## 3. The real collapse: ~4 tight-realizable classes, not 64 (the overcount)

From the speed side (`lrc_64_classes_tight_realizable_s577b.py`, exact pinch-`M` + the
half-turn class at the perturbed optimum, ~400 configs): generic configs realize **216**
distinct optimal-time classes (33 self-converse), but **only 4 classes are TIGHT-realizable**
(`M=1/14`): **2 from the AP** and **2 from V\*** `={1..11,13,24}`. So:

> The `64` self-converse classes **massively overcount** the actual worry. Tightness is a
> severe collapse: the tight worry is the **speed family `{AP, V\*}`** (≈4 boundary classes),
> both explicitly lonely (AP at `t=1/14`; V\* at `t=3/14`). The dual-Burnside `worry ⊆ 64
> self-converse` is a *very loose* bound; the operative object is the tight *speed* family.

This matches the exact census (S576: only the AP is tight over bounded speeds) and
opus-S553c (V\* is the unique non-AP tight set within distance 2 of the AP).

## 4. A caveat that exposes the real gap: V\* breaks the clean containment

At one perturbation direction, **V\*'s optimal-time class is NOT self-converse**. The reason:
the tight time has *ties* (the 14-gon antipodal pairs), so the "optimal-time tournament" is
**perturbation-direction-dependent** — `t*+ε` and `t*−ε` can give *different* classes, and
for V\* one of them is non-self-converse. So **"worry ⊆ self-converse" (HYP-2086) is not
clean at the boundary**: it holds for the symmetric AP tie-wall but is ambiguous/false for an
asymmetric tight config like V\*. The containment is a statement about a *limit* that the
boundary ties make multivalued.

## What we're still missing (the honest list)

1. **The unbounded closure — the central gap.** Everything (census, the ~4 tight classes,
   `{AP,V*}`) is over **bounded speeds**. The real question is whether the tight family is
   *finite* over ALL speeds. V\* already needed speed `24 > 1.5n`; are there tight configs at
   speed `100`, `1000`? Tao's reduction bounds speeds but astronomically. **Proving the tight
   family is `{AP, V*}` (or any explicit finite set) for unbounded speeds is the missing
   finite reduction** — not enumerating tournament classes.
2. **The boundary tie-ambiguity / realization-independence.** The optimal-time class is
   multivalued at the tight ties, so "the class" doesn't determine `M` and the self-converse
   containment is perturbation-dependent (V\*). A clean statement must work with the *limit
   tournament-with-ties* (the 14-gon tie-wall), not a single round class.
3. **The containment proof itself.** HYP-2086's `worry ⊆ self-converse` is a synthesis; V\*
   shows it needs the boundary-tie formulation to even be stated correctly. It is **not yet a
   theorem.**
4. **Strict loneliness off the regular class.** We'd want: every non-`{AP,V*}` self-converse
   class is realized only by `M>1/14` configs (strictly lonely). The data supports it (only 4
   tight classes) but it is unproven and realization-dependent.
5. **The apex coupling.** The self-antipodal residue-0 point (the apex / multiple of 7) is the
   *one* tie-pair that has no partner — it is the obstruction that has run through every
   thread (S559 zero-divisor, S552o singleton, S556 mult-of-14). The 6 mod-7 pairs resolve
   freely; the apex is the seam. **A proof must handle the apex's self-antipodal tie
   separately** — that is where the 6-free-bits picture breaks.

## Verdict / next
- **Verified:** the tight boundary is a mod-7 tie-wall (6 antipodal pairs = 6 nonzero
  residues + the self-antipodal apex). **Refuted:** 64 = 2⁶ AP-tie-resolutions (it's 33/15).
  **Found:** the 64-class bound overcounts; the tight worry is `{AP, V*}` (~4 classes), both
  lonely, but the containment is boundary-subtle (V\* non-self-converse at one direction).
- **The real target is the unbounded closure** (tight family finite over all speeds), phrased
  on the *tie-wall limit* not single round classes, with the **apex** handled separately.
- Concrete next: (1) hunt tight configs at large speeds (does anything beyond `{AP,V*}`
  appear? a targeted search by the mod-27 transversal + mod-7 apex structure); (2) formulate
  loneliness at the tight time as the two observer-adjacent gaps of the 14-gon tie-wall (the
  apex pair is the binding one); (3) prove strict loneliness for configs whose mod-7
  pair-structure differs from the AP's.

## Artifacts
```
04-computation/lrc_64_classes_antipodal_tieresolution_s577.py   (+.out)  -- 64 != 2^6 resolutions (33/15)
04-computation/lrc_64_classes_tight_realizable_s577b.py         (+.out)  -- only ~4 tight classes; V* non-sc
```
Related: HYP-2094 (even-ladder scheme), HYP-2086 (dual Burnside / self-converse), HYP-2091
(parity ladder), HYP-2089 (190 nodes), oracle-S552o (mod-7 CRT), opus-S568 (transversals),
opus-S553c (V\*), HYP-2063 (apex zero-divisor).
