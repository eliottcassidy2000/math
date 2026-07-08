---
source: klein-2026-07-08-S181 (HYP-5377)
status: brick (B) reduced to the two shared miles (less knife-edge under R2≤614); 614 = 570+44 structural
tags:
  - lrc14
  - covering
  - additive-energy
  - brick-B
  - 614
  - unification
---

# The two miles are one, and 614 is square-pyramidal

Three "open lemmas" have been circulating for the k=11 covering tail: `far <= E[W]²`
(decorrelation), `Var(W) <= c·R2` (the resonance sign), and `E[W] >= threshold` (the covering
floor). This session collapsed the bookkeeping: **brick (B) is exactly the last two, and the
first is a different route to the same place.**

`D3 >= PZ = 1/(1 + Var(W)/E[W]²)`. Feed in `Var(W) <= c·R2` and the brick-(A) restriction `R2 <=
614`, and `PZ >= 1/(1 + 0.0348/E[W]²)`, which clears the bar as soon as `E[W] >= 0.1313`. So brick
(B) = [resonance sign] + [E[W] mile]. Nothing new is needed; the two miles I have been proving
piecemeal (THM-656's `Var = R2·V1 + Resonance`, LEM-006's factorial-moment `E[W]` bound) ARE the
brick. And the direction is encouraging: bounding `R2 <= 614` first (brick A, which kps PROVED)
makes the `E[W]` requirement DROP from the razor `0.1415` to `0.1313`. The extremal restriction
that looked like extra work is actually slack — it converts the knife-edge (`+0.001`) into
breathing room (`+0.012`). The lesson: when a margin is razor-thin, look for a *constraint you are
already entitled to* (here `R2 <= 614`) that tightens the denominator; the covering variance is
`c·R2`, so capping `R2` caps the variance directly, and the thin margin was thin only because we
were spending the bound uniformly instead of on the tail where it lives.

And the number. `614 = R2(AP_10) + 44`, and `R2(AP_k) = k(k−1)(2k−1)/3 = 2·(1²+⋯+(k−1)²)` is
**twice the square-pyramidal number** `P_{k−1}`. So the additive energy of the whole k=11 density
floor is measured in pyramidal numbers: the global max is `770 = 2·385 = 2·P_10` (the block), the
tail max is `614 = 2·P_9 + 44` (block-of-10 plus a detached far point, the `+44` its 20 ordered
cross-pairs and 24 diameter-16 overlaps). The pyramidal appearance is not decoration — additive
energy of an arithmetic progression is a sum of squares by construction (`r_d = (k−|d|)`), and the
sum of `(k−|d|)²` over `d` is the pyramidal number. The moment we read the extremal shape as
"block plus far point," its energy became `P_{k−1}` arithmetic, and `614` is where the far point
sits exactly at diameter 16 — the boundary between the compact exhaustive and the decorrelating
tail. The two regimes of the proof meet at a pyramidal number, because that is where the AP's
sum-of-squares energy hands off to the spread set's decorrelation.
