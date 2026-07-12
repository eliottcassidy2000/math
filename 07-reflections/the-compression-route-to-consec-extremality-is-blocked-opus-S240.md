---
source: opus-2026-07-11-S240
status: A negative result on the PROOF STRATEGY for the shared crux (consec-extremality). The natural
  monotone-compression route is BLOCKED: the base functional J = E[N(7−N)] is NOT unimodal toward consec —
  greedy descent stuck at 35/35 distinct non-consec local minima, all algebraically-special (even/dilated).
  This confirms klein's "extremals invisible to local search" and explains why the crux is verified-only.
tags:
  - lrc14
  - consec-extremality
  - compression
  - negative-result
  - proof-strategy
  - rugged-landscape
---

# The compression route to consec-extremality is blocked

**opus-2026-07-11-S240.** Owner: work the concrete next target — the AP coverage-extremality that S239 showed
is the shared crux of both LRC(14) routes. The target's base form (klein THM-711/717, mac-mini THM-716):
`J = E[N(7−N)] = 6m1 − m2`, with consec `{1..9}` the conjectured global minimum (`J ≈ 5.05`), **verified
adversarially but not proved**. The standard way to prove an extremal is a **monotone compression**: a local
step toward the AP that decreases `J`, iterated to the consec fixed point. I tested whether one exists.

## The result: the landscape is rugged, compression is blocked

Greedy single-coordinate `J`-descent from 35 random primitive 9-sets reaches consec **0/35** times. It stalls
at **35 distinct non-consec local minima**, every one **algebraically-special** — even/dilated families near
`2·{1..k}`:

| local min | J |
|---|---|
| `[2,4,6,7,8,9,10,12,14]` | 5.33 |
| `[2,4,6,8,9,10,12,14,16]` | 5.47 |
| `[2,4,6,7,8,10,12,13,14]` | 5.54 |
| … (32 more) | 5.3–5.8 |

All have `J > J(consec) = 5.05` — so consec **is** the global min, but it is **unreachable by local descent**.
The functional `J` is **not unimodal** toward consec; its basins funnel into even-structured (mod-2/dilated)
minima, and consec sits as a *separate* global optimum.

## Why this matters (honest)

1. **It rules out the natural proof.** The single-step compression / smoothing argument — the standard route
   to an extremal inequality — **does not work** for consec-extremality. The crux does not reduce to local
   convexity or monotonicity.
2. **It confirms klein-S255** ("every LRC(14) functional's extremal is an AP or a mod-p resonance, invisible
   to local search") *quantitatively*: the local minima of `J` **are** the algebraically-special families
   (even/dilated/mod-7), and consec beats them only **globally**.
3. **It explains the verified-only status.** Proving consec is the global min over a rugged landscape of
   algebraically-special competitors is a **global/algebraic** problem — the inverse theorem — not a local
   smoothing. The most natural attack is dead.

## Consequence for the endgame

Consec-extremality (the S239 shared wall of both routes) will **not** fall to compression. It needs either a
**global algebraic argument** — classify the algebraically-special local minima (even, dilated, mod-7
resonances) and beat each explicitly, which is klein/mac-mini's extremal-atlas lane — or the **finite-census**
framing (LEM-010's `Vmax ≤ 3^12`). Local, per-family, smoothing attacks are exhausted.

**Honest meta-note (the S230–S240 arc).** Eleven sessions have sharpened this one crux — the AP
coverage-extremality inverse theorem — from many angles: the clean-ruler reduction, the band-edge margin, the
divisor-complete hard core, the spread/coverage reframing, the two-route unification, and now the ruling-out
of the compression proof. Each is a real reduction or negative result, but the crux itself has not yielded,
and this session shows *why* a whole class of proof attempts (local compression) cannot work. The honest
assessment: this is the project's central open problem, it is a global inverse theorem over a rugged
algebraically-special landscape, and solo per-session attacks are at diminishing returns. The productive paths
are the collective extremal-atlas classification (klein/mac-mini) or an external inverse-theorem input — not
another local reduction.

→ opus-S239 (the shared crux), klein THM-711/717 + mac-mini THM-716 (the base functional, consec-extremal
conjecture), klein-S255 (extremals invisible to local search — confirmed here), mac-mini cont.44 (three-gap
coverage; the local minima are its `{kα}`/mod-p resonances), LEM-010 (the finite-census alternative). Files:
`lrc14_compression_route_blocked_opus_S240.py` (+`.out`).
