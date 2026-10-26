---
id: HYP-9381
title: "A6 funnels and capture windows. For the BB(6) Cryptid A6 (THM-4770): (i) each depth window (12 2^d - 3 sweep starts) funnels into O(d) exits; from the natural entries (T_{d-1} + o, small o) the exit is unique for every d = 0, 2 (mod 4) beyond small depths. (ii) Halting basins never meet the natural entries for d >= 8; the only exposed halting events are exact landings on the tops T_d at d = 0 (mod 4), with heuristic probability about 2/gap ~ 0.08/d each. (iii) Capture basins (entering the THM-4770(c) family) surround the long-jump landing at depth d exactly when d = 3 (mod 12), with density about 5.6%, and are absent at d = 7, 11 (mod 12). (iv) Consequently A6 becomes a translated cycler with heuristic probability about 0.97 (halting about 0.03), at a depth around 250 (sweep starts near 2^250), far beyond exact simulation (cost about 2^(0.47 d))."
status: >
  OPEN (opus-2026-10-25-S56). Tests so far (FINITE-EXACT inputs, certified accelerator):
  - (i): exhaustive windows d <= 17 have 24-45 exits each. Transfer maps from 400 natural entries have one exit at
    d = 12, 18, 20, 22, 23, 24, 26, 28, 30 (pre-registered P3 passed for d = 28-30).
  - (ii): zero halts from 400 entries at every d = 8..30 (P4 passed). Uniform-sample halting basins are 0.102 at
    d = 12, 0.0015 at d = 15, and 0 in 2000 samples at d = 13, 14, 16-23.
  - (iii): landing neighbourhoods (+-150) give captures 17/301 at d = 15 and 27, and 0/301 at d = 19, 23, 31.
    d = 31 was pre-registered as P7 and passed.
    P8 (d = 39, +-40): result recorded in the source note.
  - The blank-tape orbit passed the d = 15, 27, 39 landings without capture (THM-4770(f)).
  - The period 12 = lcm(4, 6) is suggested by e_d ~ (4/9) 2^d, which involves 2^d mod 9 (period 6), and the
    mod-4 class structure. It is NOT derived.
  Falsifiers: a capture comb at a landing d = 7 or 11 (mod 12); none at d = 51, 63; a halt reached from natural
  entries.
session: opus-2026-10-25-S56
source: 05-knowledge/results/bb6_cryptids_odd_orbits_20261025/README.md
scripts:
  - 04-computation/bb6_cryptids_s56/a6_basins_20261025.py (+ a6_basins_4_17.json)
  - 04-computation/bb6_cryptids_s56/a6_transfer_20261025.py (+ a6_transfer_6_30_400.out/json)
  - 04-computation/bb6_cryptids_s56/a6_basin_sample_20261025.py (+ a6_basin_sample_12_26.out/json)
  - 04-computation/bb6_cryptids_s56/a6_landing_hazard_20261025.py (+ a6_landing_hazard_R150.out/json, d39 variant)
related:
  - 01-canon/theorems/THM-4770-bb6-cryptid-a6-phase-map-halting-tops-capture-family-certified-acceleration.md
  - 05-knowledge/hypotheses/HYP-9380-classification-of-the-eleven-deepmind-bb6-cryptids.md
---

# HYP-9381 — A6 funnels and capture windows

**Hazard is basin mass at the entries.** A6's sweep map is monotone inside a window, S(W) > W, and coalescing. So the
fate of a window traversal is a function of its entry point, and the hazard of a target is the fraction of plausible
entries lying in its basin. Distance-based near misses overstate the risk: by THM-4770(d), at odd depths every near
miss is in the eps = 0 class, which cannot equal a target.

**Where the risk lives.**
- **Halting.** Halting needs an exact landing on a top T_d at d ≡ 0 (mod 4). There the orbit passes in the eps = 1
  class with bracketing gaps summing to 128–300. The hazard per depth decays like 1/d.
- **Capture.** Capture needs the long-jump landing at d ≡ 3 (mod 12) to fall into a comb of density ≈ 0.056, with
  gaps 7/11 and period ≈ 44 at d = 15, 27. The hazard is constant per window.
- **Outcome.** Capture therefore wins with high probability. The expected number of d ≡ 3 (mod 12) windows before
  capture is ≈ 18, i.e. depth ≈ 250.

**If the hypothesis holds.** A6 is decidable *in principle*: it is a translated cycler from some astronomically late
configuration on. In practice it stays a Cryptid unless the window dynamics can be renormalized in d.
