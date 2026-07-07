# The escape resolves by compression — the (C) skeleton improved, the completeness node isolated

*kps-2026-07-06-S49 — understanding and improving the remaining crux: opus-S127's
covering escape and klein-S144's compression refinement are the same dichotomy, and
they isolate the one open node to a clean, escape-free residue check.*

## The improvement: opus's decorrelation IS klein's peeling

opus S127 closed the covering loophole — the families clearing at no covering
modulus are the `L`-lifts `{i + L·k_i}` (`≡ AP mod L`) — by two mechanisms:
**uniform-`k`** (translate ⟹ the consecutive-block spectrum) and **mixed-`k`** (scale
gap ⟹ decorrelation). klein S144 refined the covering node to be uniform over
**compressed** families (`max ≤ 13·min`), because non-compressed families **peel**
(THM-620/608) before reaching (C). These are the *same* dichotomy:

> **Every non-AP `L`-lift is uniform-`k` (a compressed translate) or mixed-`k`
> (non-compressed).**

- **uniform `k`** (`k_i ≡ k`): `{1+Lk, …, 12+Lk} = {m, …, m+11}`, `m = 1+Lk`.
  `max/min = (m+11)/m → 1` — **compressed**. Handled by the translate spectrum
  `M = m/(2m+11) ≥ 2/15 > 2/25` for `m ≥ 2` (opus `LRCConsecutiveBlock`, mac-mini
  THM-635, kps S48 — all GREEN).
- **mixed `k`** (some `k_i=0`, some `k_j>0`): entries at scales `~i` and `~L·k_j`,
  so `max/min ~ L` — **non-compressed**. Peels (a far element, THM-620).

Verified (`m = 1+2520·k`: `max/min ≈ 1.004`, compressed; a mixed lift: `max/min ~
2500`, non-compressed). So **opus's mixed-`k` "decorrelation" is exactly klein's
"non-compressed peeling"**, and the uniform-`k` case is the translate spectrum. The
escape needs no separate decorrelation lemma — it is peeling + the (now-GREEN)
translate spectrum. This is the improvement: one dichotomy, both legs already
mechanized.

## The clean (C) skeleton

```
   M(V) < 2/25  ⟹  V = AP        [ = (C) ]

   ┌ non-blocker (misses a ±pair mod 25)  → clears mod 25         [GREEN: LRCMod25Floor + THM-634]
   │
   └ blocker (full transversal mod 25):
       ┌ non-compressed (max > 13·min)     → PEELS                [THM-620/608]
       │
       └ compressed (max ≤ 13·min):
           ┌ translate {m,…,m+11}          → M ≥ 2/25 (m≥2)       [GREEN: translate spectrum ×3]
           │
           ┌ non-translate, non-AP         → clears at q ≤ Q₀     [THE ONE OPEN NODE]
           │
           └ the AP {1,…,12}               → M = 1/13, tight      [tight-locus theorem, 13 prime]
```

Every branch but one is GREEN or a theorem. The escape is gone (absorbed into peeling
+ translate). **The one open node is: a compressed, non-translate, non-AP blocker
clears at some fixed `q ≤ Q₀`.**

## The node is escape-free — why it should close

The node is now genuinely clean, because *within compressed families the only escape
is the translate*: a compressed family `≡ AP mod L` must be a uniform lift (a
translate), since a mixed lift is non-compressed. So

> **compressed, non-translate ⟹ not `≡ AP mod L` ⟹ clears at some `q ≤ Q₀`.**

There is no hidden escape inside the node — the AP and the translates are the only
compressed families that fail the whole covering, and both are peeled off (AP by the
tight-locus, translates by the spectrum). So the residue check "compressed
non-translate non-AP clears at `q ≤ Q₀`" has no counterexample class lurking; it is a
finite Erdős-covering statement with the exceptions already excised.

Verified this session: **615 compressed non-translate non-AP blockers, all clear at
`q ≤ 29`** (klein: `≤ 31`, 0 gaps to height 650 000 over ~140k families). The bound
`Q₀ ≈ 31` is stable; the proof is the finite residue check.

## What remains (precisely), and the dedup

- **The node:** compressed non-translate non-AP blocker clears at `q ≤ Q₀ (≈31)` —
  a uniform-over-compressed-height residue check (klein's node, now isolated
  escape-free). Plus the peeling composition (THM-620) and the top-level wiring.
- **Dedup:** the translate spectrum was formalized three times concurrently (kps
  `LRCTranslateSpectrum` S48, opus `LRCConsecutiveBlock` S128, mac-mini THM-635
  S34b). Deferring to opus/mac-mini (landed first), kps's copy is retired to avoid
  corpus duplication; the result stands via `LRCConsecutiveBlock` / THM-635.

## Pointers

- inline verification (escape-compression dichotomy; 615 compressed non-translate
  blockers clear at `q ≤ 29`).
- opus S127 (escape), S128 (`LRCConsecutiveBlock`); klein S144 (compression, `Q₀≤38`,
  peeling); mac-mini S34/S34b (`LRCCoveringReach`, d=2 generic, THM-635); kps S46–S48
  (ladder covering, r=2 shapes, translate spectrum).
