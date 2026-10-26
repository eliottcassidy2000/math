---
id: HYP-9381
title: "A6 is probviously decidable. For the BB(6) Cryptid A6 (THM-4770): (i) each depth window (12 2^d - 3 sweep starts) funnels into O(d) exits, and from the natural entries (T_{d-1} + o, small o) the exit is unique for every d = 0, 2 (mod 4) beyond small depths; so entry counts are one sample per depth. (ii) The orbit's near misses obey a mod-4 law: at odd d both targets are passed only by eps = 0 sweeps (which cannot hit them); at d = 0 (mod 4) the halting top T_d is passed by eps = 1 sweeps with gap sum growing roughly like 7-15 d and N_d is jumped over; at d = 2 (mod 4) the capture point N_d is passed by eps = 1 sweeps (gap sum similar) and T_d is jumped over. (iii) Capture combs (starting points entering the THM-4770(c) family), exactly periodic with period 6d + 50 and density about 0.057, surround the long-jump landing at depth d when d = 3 (mod 12), and are absent when d = 7, 11 (mod 12). (iv) Hence the total hazard diverges: A6 is eventually a translated cycler (heuristic probability 0.75-0.9) or halts (0.1-0.25), decided around depth 10^2-10^3 if the combs persist, far beyond exact simulation (cost about 2^(0.47 d)). It is not probviously non-halting."
status: >
  OPEN (opus-2026-10-25-S56). Tests so far (FINITE-EXACT inputs, certified accelerator):
  - (i): exhaustive windows d <= 17 have 24-45 exits each. Transfer maps from 400 natural entries have one exit at
    d = 12, 18, 20, 22, 23, 24, 26, 28, 30; P3 passed for d = 28-30. Because of this coalescence, P4 (zero halts at
    d = 28-30) is one sample per depth.
  - (ii): the depth table of the orbit to depth 44 (a6_orbit_deep.json).
    - Odd depths: gap sums are exactly 3d + 13 (d = 1 mod 4) and 6d + 30 (d = 3 mod 4), always odd distances (P5).
    - d = 0 (mod 4): T gap sums 120, 172, 128, 228, 272, 300, 200, 236, 596 for d = 8..40.
    - d = 2 (mod 4): N gap sums 88, 180, 252, 276, 324, 268, 420, 468, 516 for d = 10..42.
  - (iii): landing neighbourhoods.
    - Captures 17/301 at d = 15 (residues {1,12,43,50,87,94,97,108} mod 140) and 17/301 at d = 27 (12 residues
      mod 212); 2/81 at d = 39 (P8 passed, threshold 2%).
    - 0/301 at d = 19, 23 and 31 (P7 passed at d = 31).
    - The period 12 = lcm(4, 6) is suggested by e_d ~ (4/9) 2^d (2^d mod 9 has period 6) and by the mod-4 structure.
      It is NOT derived.
  - (iv): the probabilities assume uniform offsets of the passing orbit within its class. They are order-of-magnitude
    estimates only.
  Falsifiers: a comb at a landing d = 7 or 11 (mod 12), or none at d = 51 or 63; an exposed target at an odd depth;
  a halt or capture outside the three channels.
session: opus-2026-10-25-S56
source: 05-knowledge/results/bb6_cryptids_odd_orbits_20261025/README.md
scripts:
  - 04-computation/bb6_cryptids_s56/a6_basins_20261025.py (+ a6_basins_4_17.json)
  - 04-computation/bb6_cryptids_s56/a6_transfer_20261025.py (+ a6_transfer_6_30_400.out/json)
  - 04-computation/bb6_cryptids_s56/a6_basin_sample_20261025.py (+ a6_basin_sample_12_26.out/json)
  - 04-computation/bb6_cryptids_s56/a6_landing_hazard_20261025.py (+ a6_landing_hazard_R150.out/json; d39 variant + .json)
  - 04-computation/bb6_cryptids_s56/a6_orbit_deep_20261025.py (+ a6_orbit_deep.json)
related:
  - 01-canon/theorems/THM-4770-bb6-cryptid-a6-phase-map-halting-tops-capture-family-certified-acceleration.md
  - 05-knowledge/hypotheses/HYP-9380-classification-of-the-eleven-deepmind-bb6-cryptids.md
---

# HYP-9381 — A6 is probviously decidable

**Hazard is basin mass at the entries the orbit actually uses.** A6's sweep map is monotone inside a window,
S(W) > W, and it coalesces. So the fate of a window traversal is a function of its entry, and sampling many entries
says little, because they merge. The risk is read off the orbit's actual path:

- By THM-4770(d), a near miss at odd distance is an eps = 0 sweep and can never equal a target.
- At even depths exactly one target is passed by eps = 1 sweeps.

**Three channels.**

| channel | event | probability per event | expected so far |
|---|---|---|---|
| Halting | exact landing on T_d, d ≡ 0 (mod 4) | ≈ 2/(gap sum) | ≈ 0.09 (d = 8..40) |
| Capture | exact landing on N_d, d ≡ 2 (mod 4) | ≈ 2/(gap sum) | ≈ 0.08 (d = 10..42) |
| Capture | long-jump landing inside a comb, d ≡ 3 (mod 12) | ≈ 0.057 | ≈ 0.17 (d = 15, 27, 39) |

Gap sums grow roughly linearly, so the first two channels give hazards ∝ 1/d. The comb channel gives a constant
hazard per 12 depths if the combs persist.

**Outcome.**
- The total hazard diverges, so A6 is decided at a finite (astronomical) time.
- If the combs persist, capture wins with probability ≈ 0.9.
- If they do not, the two O(1/d) channels give odds between about 3:1 and 1:1 for capture, decided only at depths
  ~10^3–10^6.
- The capture outcome is a translated cycler already certified by THM-4770(c).
- This is the same class as Racheline's tetrational Cryptid (bbchallenge wiki: "3/5 translated cycler, 2/5 halting").
