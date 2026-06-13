---
source: opus-2026-06-02-S565 (remote-control)
status: SYNTHESIS — the LRC's two regimes ARE Burnside's two sides (orbit=round/A000016/open; fix=self-converse/boundary); verified
tags: [LRC, Burnside, A000016, round-tournament, self-converse, orbit-fix-duality, worry-ignore, n14, S547, HYP-1998, S563, S564]
---

# Dual Burnside for the LRC: open=orbit side (round), boundary=fix side (self-converse)

**Prompt (user):** try dual Burnside regarding the LRC.

Burnside is the repo's tournament-counting engine (A000568, A000016). Applied to
the LRC, its **two sides — orbit-counting and fixed-points — turn out to be the
LRC's two regimes** (the worry/ignore split, S564). That is the "dual Burnside."

## The two sides

Burnside: `#orbits = (1/|G|) Σ_g |Fix(g)|`. The lemma has an *orbit* face (count
classes) and a *fixed-point* face (the symmetric configurations). The LRC realizes
both:

### ORBIT side — open times → ROUND tournaments = A000016 (the IGNORE regime)

At a generic (open) time the runner sub-tournament is **round** (out-set = a
clockwise arc), and round classes are the **necklace Burnside count** A000016 (only
odd divisors, same odd-cycle structure as A000568) — HYP-1998. Verified:

| m (runners) | A000016 (round) | A000568 (all) | round / total |
|---|---|---|---|
| 3 | 2 | 2 | 1 |
| 5 | 4 | 12 | 0.33 |
| 7 | 10 | 456 | 0.022 |
| 9 | 30 | 191 536 | 1.6e-4 |
| 11 | 94 | 903 753 248 | 1.0e-7 |
| **13** | **316** | **48 542 114 686 912** | **6.5e-12** |

So at n=14 the lonely sub-tournament is **one of only 316 round classes out of 48
trillion** tournaments. This is the IGNORE regime: open times have positive measure,
the orbit is generic/round, loneliness is automatic (S563/S564).

### FIX side — boundary times → SELF-CONVERSE, maximally-symmetric (the WORRY regime)

The hard cases live on the **boundary** (antipodal-tie times), where the tournament
need not be round (HYP-1998). The loneliest configuration — the regular round
tournament (roots of unity) — is exactly a **Burnside fixed point**:

- **Self-converse:** `T = T^op` via the reflection `i ↦ -i` (verified for all
  m=3..13). It is fixed by the complement/reversal involution.
- **Rotation-symmetric:** full cyclic automorphism `C_m` (order m, verified).
- It carries the **(q,q) cycle-type automorphism** at `n=2q` with Burnside
  `Fix = 2^{n-1}` (S547).

These are precisely the configs the *fixed-point* side of Burnside picks out: the
ones with large automorphism group. **The WORRY set = the Burnside-fixed
(maximally symmetric, self-converse) boundary configs** = the resonance-maximal AP
family (S563) = the tight family (S553).

## The duality, stated

> **LRC's two regimes ARE Burnside's two faces.** The easy, positive-measure
> majority lives on the **orbit side** — its sub-tournament is one of the few round
> classes (A000016), a closed-form Burnside *count*. The hard, measure-zero worry
> set lives on the **fixed-point side** — the self-converse, rotation-symmetric,
> `(q,q)`-automorphic configurations that Burnside's `Fix` terms isolate.

This unifies, under one lemma:
- round = A000016 (HYP-1998) — the orbit count;
- (q,q), `Fix=2^{n-1}` (S547) — the fixed-point count of the extremiser;
- worry/ignore (S564), resonance-maximal (S563) — the two sides named dynamically.

## Why it is useful (and the "dual" payoff)

- **It bounds the worry-set by symmetry.** A counterexample's sub-tournament cannot
  be a generic round class (those are open/lonely); it must be a **boundary,
  self-converse** class — a Burnside-fixed object, of which there are *few*
  (self-converse tournaments are a thin, Burnside-counted family; cf. the repo's
  A002785/A005639 self-converse work). So the LRC worry-set is contained in the
  **self-converse boundary classes** — a small, enumerable target.
- **The lonely menu sits between the two sides:** `A000016(m) ⊆ lonely menu ⊆
  A000016 ∪ {self-converse boundary classes}` (HYP-1998 open-Q-B). The "dual
  Burnside" says the menu's surplus over round is exactly the *fixed-point* layer.
- **Proof shape:** show every **self-converse boundary** runner tournament (the fix
  side) is still source-realizable (lonely). The orbit side (round/open) is free.
  This is the Burnside translation of "prove LRC only for the resonance-maximal
  tight family."

## Honest status

The orbit side (round = A000016) and the fixed-point anchors (regular config
self-converse + C_m + (q,q) Fix=2^{n-1}) are computed/verified and known
(HYP-1998, S547). The **new** content is the *duality framing*: that LRC's
worry/ignore split is literally Burnside's fix/orbit split, which localises the
worry-set to the self-converse boundary classes and gives the proof its
Burnside-shaped target. Not a proof; an organising principle plus a concrete
containment (worry ⊆ self-converse boundary).

**Artifacts:** `04-computation/lrc_dual_burnside_s565.py` (+`.out`). Builds on
HYP-1998 (round=A000016), S547 ((q,q)), S563 (resonance/orbit), S564
(worry/ignore), and the repo's self-converse counts (A002785). New: **HYP-2086**.
