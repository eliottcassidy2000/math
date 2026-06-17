---
id: THM-525
title: "The easy-dominates-hard covering reduction of LRC(14): hard covering 13-sets are an easy LRC(12)-core C plus a runner parked in section 0; centering the parked runner ('perfect middle of section 0') turns OPEN-Q-108's uniform-fattening crux into a transversality statement about C's gap-1/14 lonely set"
status: STUB / RESERVED (kind-pasteur-2026-06-17-S2). Being built by a workflow (structure investigators + adversarial verify). Fill with PROVED/VERIFIED/CONJECTURE marks per part once results land.
source: kind-pasteur-2026-06-17-S2 (user reframe: hardest configs park a runner in section 0 forever; center it in the perfect middle; easy & hard cases come hand in hand; use the easy cases' structure to kill the hard cases)
depends_on:
  - THM-524   # binding-pair reduction + regions/sections reframe (covering = off-grid, parked runner)
  - THM-523   # covering needs a multiple of every q in 2..14 (the hard core IS the covering family)
  - THM-522   # measure-side scale-invariance + quantization
related:
  - OPEN-Q-108  # the uniform fattening lemma = THE crux; this attacks it constructively via easy/hard
  - HYP-2573
  - HYP-2574
  - HYP-2575
external: proven LRC(12) and LRC(13) (Sungkawichai–Trakulthongchai 2026) — the EASY cores; Goddyn–Wong tight family.
---

# THM-525 — Easy dominates hard: the centered-parked-runner covering reduction (STUB)

**Claim under construction.** Every "hard" 13-speed config for LRC(14) (covering set: contains a
speed `w ≡ 0 mod 14`, so on-grid `gridM=0`, the runner sits in section 0 forever) decomposes as
`S = C ∪ {w}` with `C` a 12-speed **easy core** that LRC(12) already proves lonely at gap `1/13`
(hence positive gap-`1/14` measure `meas(G_C)`). The reduction: bound `M(S) ≥ 1/14` by **centering
the parked runner** — choosing the lonely time `τ` inside `G_C` so that `‖wτ‖` is pushed toward the
**middle of its safe band** (the "perfect middle of section 0"), and showing the easy core's
gap-`1/14` set is fat enough to contain such a `τ`. This converts OPEN-Q-108 (uniform
`meas(G_C) ≥ c`) from existence into the **transversality** of `w`'s danger comb against `G_C`.

Baseline (VERIFIED): covering core `{1..11,13,84}` has `M=7/89` at `τ=37/89`, where 84 is a
BINDING runner (dist `7/89=M`, at frac `82/89 ≈ −7/89`) — i.e. the parked runner sits at the
EDGE of section 0, not the middle. So the naive optimizer does NOT center `w`; the centering is a
constructive LOWER-BOUND device, not the optimizer's behavior. Quantifying the gap between
"centered-`w` safe time" and "optimal time" is the content.

(Parts A–E + status to be filled from the S2 workflow.)
