---
source: opus-2026-07-11-S250
status: A REFRAMING (owner-proposed), developed + verified. "The problem is a change of bases" maps exactly
  onto the project's clean-ruler/clearing framework: base = modulus q, digit = residue v_i·a mod q, clean base
  = clearing modulus. The base-landscape re-derives WHY AP and V* are extremal: they are RESONANT (best=0) in
  every base 2..13 and hit EXACTLY 1/14 only at base 14 — "base-rigid". The reframing is correct as far as it
  goes; the obstruction it makes vivid is that a base change cleaning the runners' MOTION distorts the danger
  ZONE — LRC couples the archimedean place (the 1/14 band) with the finite places (the moduli), asking for a
  finite certificate of an archimedean fact. Actionable direction: prime-local (p-adic) decomposition, 14 = 2·7.
tags:
  - lrc14
  - change-of-base
  - modulus-as-base
  - archimedean-vs-finite-places
  - base-rigidity
  - observer-lens
  - reframing
---

# LRC(14) as a change of bases: archimedean vs finite places

**opus-2026-07-11-S250.** The owner proposed: a runner's position is an infinite base-b expansion; each runner
has a base in which its constant change is "the unit," so calculations are clean from that perspective; **the
problem is a change of bases.** Developed and verified — it is a genuine and accurate reframing, and it
re-derives the crux.

## The dictionary (owner's language → the machinery)

| owner | project |
|---|---|
| runner position `{v_i t}` (digits after the point) | point on the circle — the **archimedean** place |
| "base" | **modulus `q`** (the clean-ruler / clearing modulus) |
| digit of runner `i` in base `q` | residue `v_i·a mod q` |
| a base where runner `i`'s change is "the unit" | a `q` where runner `i` is **non-resonant** (spread, not stuck at 0) |
| "the problem is a change of bases" | hunt for a modulus `q` (+ multiplier `a`) that cleans **all** runners |

Formally `best_q(v) = max_a min_i ‖v_i a/q‖` is the floor achievable in base `q`, and `M(v) = max_q best_q(v)`.
LRC(14) ⟺ `M(v) ≥ 1/14` for every family — i.e. **some base cleans the configuration to margin ≥ 1/14.** The
owner's reframing is exactly right: the project has been doing "change of base" all along, under the name
*clearing at a modulus*.

## What the reframing correctly captures — and re-derives

**Per-runner a clean base always exists** (trivial: any `q` coprime to `v_i` spreads runner `i`). The whole
difficulty is that you need **one base clean for all runners simultaneously** — and the CRT-combination of the
per-runner clean bases is exactly the **covering-system** structure the project studies.

Verified base-landscapes re-derive why the extremizers are extremal:

- **AP `{1..13}` and V\* `{1..11,13,24}` are BASE-RIGID.** `max_q best_q = 1/14`, peaked at base `q* = 14`. Both
  are **resonant** (`best_q = 0`, a runner stuck at the origin) in **every base `2..13`** — because each
  contains a multiple of every `d ∈ {2..13}` (divisor-complete). Their *only* clean base is `14`, where they
  hit **exactly** `1/14` with **zero margin**. That is precisely why they are the hard extremizers: they
  defeat every base below 14 and only marginally clear at 14.
- **AP and V\* are base-TWINS** — identical coarse landscapes (same top bases `14, 85, 71, 57`, same resonant
  set). The mod-14 fine structure (S249: full residues vs vacate-8/double-2) is the *only* distinguisher,
  invisible at coarse resolution. V\* is a twin of AP because base `14 = 2·7` is composite (the doubling
  `12→24`). The change-of-base picture thus reproduces the composite-14 signature.
- A **clearable** family has a *small* clean base beating `1/14` (e.g. `5/38` at `q=38`) — easy. A **divisor-
  complete adversarial** family is resonant/marginal in every small base; its clean base is pushed **large**
  (MISTAKE-110's unbounded-modulus phenomenon). "No small base works" is the residual, in base language.

## The obstruction the reframing makes vivid

Why "just a change of base" does not *close* it: a base change that makes a runner's **motion** the unit (a
unimodular relabel of the speed lattice `ℤ^k`) simultaneously **distorts the danger zone**, which is defined
per-runner in the **original** coordinates (the band `‖·‖ < 1/14`). You can clean the **dynamics** (finite
places / the choice of `q`) *or* keep the **danger metric** simple (the archimedean place), **not both**.

LRC **couples the archimedean place** (the `1/14` band — a real-distance fact) **with the finite places** (the
bases/moduli). The rational-time reduction says: to certify the archimedean fact "a lonely point exists,"
choose a **finite** place (base `q`) and check residues. The conjecture asserts a **finite certificate for an
archimedean fact.** This is exactly why pure finite-place (p-adic) reasoning cannot alone settle it — the band
lives at infinity — and why the extremizers are the divisor-complete near-AP families that are *bad at every
finite place at once* yet only *marginally* good at the one base (14) tied to the archimedean scale `1/(k+1)`.

## Net — and the actionable direction

The reframing is correct and clarifying: LRC(14) = "**some base cleans the configuration to margin ≥ 1/14**,"
and the hard families (AP, V\*) are the **base-rigid** ones — resonant in every base `2..13`, exactly-`1/14` at
base 14. It doesn't hand us a proof (the archimedean/finite-place coupling is the genuine content), but it
points cleanly at the **prime-local (p-adic) decomposition**: work base = prime power, splitting `14 = 2·7`
into the 2-adic and 7-adic pieces, where the clean ruler is strongest (THM-712) and the apex prime 7 is the
seam. That is the natural place to try to turn "no *single* small base is clean for the extremizers" into a
prime-by-prime certificate — the concrete next experiment this frame suggests.

→ opus-S249 (mod-14 fine structure = the base-twin distinguisher), opus-S248/S247 (empty window, `{AP,V*}`),
THM-712 (prime clean ruler — the p-adic direction), THM-366 (divisor-complete = resonant in every base 2..14),
MISTAKE-110 (unbounded/adaptive clean base), the observer-lens thread (per-runner frames). Files:
`lrc14_base_change_landscape_opus_S250.py` (+`.out`).
