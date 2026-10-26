---
id: HYP-9381
title: "A6 is probviously decidable. For the BB(6) Cryptid A6 (THM-4770): (i) each depth window (12 2^d - 3 sweep starts) funnels into O(d) exits (19-45 for d <= 17), and from the natural entries (T_{d-1} + o, small o) the exit is unique for d = 0, 2 (mod 4) beyond small depths, so entry counts are one sample per depth. (ii) The orbit's near misses obey a mod-4 law: at odd d the targets are passed only by sweeps that cannot hit them (eps = 0 sweeps of depth d, and the first start of depth d+1); at d = 0 (mod 4) the halting top T_d is passed by eps = 1 sweeps with gap sum growing roughly like 7-15 d and N_d is jumped over; at d = 2 (mod 4) the capture point N_d is passed by eps = 1 sweeps and T_d is jumped over. (iii) Capture combs (starting points entering the THM-4770(c) family) surround the long-jump landing at depth d when d = 3 (mod 12), with density about 0.057 (8 residues mod 140 at d = 15; 12 residues consistent with period 212 at d = 27; a two-point fit gives period 6d + 50), and are absent when d = 7, 11 (mod 12). (iv) Hence the total hazard diverges and A6 is eventually captured (a certified infinite macro chain; on the tape its period grows quadratically, so not a bbchallenge translated cycler) or halts: capture about 0.9-0.95 and halting 0.05-0.1 if the combs persist; if they do not, halting is slightly favoured (about 0.55 of decided outcomes) with median decision depth about 10^5. The probability that the combs persist is not estimated."
status: >
  OPEN (opus-2026-10-25-S56; independently audited, AUDIT.md, MISTAKE-669; (iv) restated per scenario after the audit found the
  first odds internally inconsistent). Tests so far (FINITE-EXACT inputs, certified accelerator):
  - (i): exhaustive windows d <= 17 have 19-45 exits each. Transfer maps from 400 natural entries have one exit at
    d = 12, 18, 20, 22, 23, 24, 26, 28, 30 (P3 passed for d = 28-30). Uniform samples at d = 24: 1988/2000 share the
    main-stem exit.
  - (ii): depth table of the orbit into depth 44 (a6_orbit_deep.json).
    - Odd depths: gap sums are exactly 3d + 13 (d = 1 mod 4) and 6d + 30 (d = 3 mod 4).
    - d = 0 (mod 4): T gap sums 120, 172, 128, 228, 272, 300, 200, 236, 596 for d = 8..40.
    - d = 2 (mod 4): N gap sums 88, 180, 252, 276, 324, 268, 420, 468, 516 for d = 10..42.
    - Pre-registered P5 holds for the three approaches inside depth d. The fourth, the first start past T_d, lies in
      depth d+1 and is eps = 1 there (audit). P5 as literally stated therefore fails for it, harmlessly.
  - (iii): landing neighbourhoods.
    - Captures 17/301 at d = 15 and at d = 27.
    - 0/301 at d = 19, 23 and 31 (P7 passed at d = 31).
    - 2/81 at d = 39 (P8 passed at threshold 2%).
    - The comb gaps at d = 15 are 11, 31, 7, 37, 7, 3, 11, 33. The windows cover at most 2.15 periods.
    - Period 12 = lcm(4, 6) is suggested by e_d ~ (4/9) 2^d (2^d mod 9 has period 6) and by the mod-4 structure.
      It is NOT derived.
  - (iv): uniform-offset model for the passing orbit. The audit simulated it with gap slopes 11 and 12.
    - With the combs: halting 0.05-0.10, capture 0.90-0.95.
    - Without them: halting 0.27, capture 0.25, and 0.48 still undecided at depth 2e5. Survival goes like
      (D/44)^-0.09, so 10% is undecided beyond about 1e13.
  Falsifiers: a comb at a landing d = 7 or 11 (mod 12), or none at d = 51 or 63; an exposed target at an odd depth;
  a halt or capture outside the three channels.
session: opus-2026-10-25-S56
source: 05-knowledge/results/bb6_cryptids_odd_orbits_20261025/README.md
scripts:
  - 04-computation/bb6_cryptids_s56/a6_basins_20261025.py (+ a6_basins_4_17.json)
  - 04-computation/bb6_cryptids_s56/a6_transfer_20261025.py (+ a6_transfer_6_30_400.out/json)
  - 04-computation/bb6_cryptids_s56/a6_basin_sample_20261025.py (+ a6_basin_sample_12_26.out, a6_basin_sample_12_26_2000.json)
  - 04-computation/bb6_cryptids_s56/a6_landing_hazard_20261025.py (+ a6_landing_hazard_R150.out/json; d39 variant + .json)
  - 04-computation/bb6_cryptids_s56/a6_orbit_deep_20261025.py (+ a6_orbit_deep.json)
related:
  - 01-canon/theorems/THM-4770-bb6-cryptid-a6-phase-map-halting-tops-capture-family-certified-acceleration.md
  - 05-knowledge/hypotheses/HYP-9380-classification-of-the-eleven-deepmind-bb6-cryptids.md
---

# HYP-9381 — A6 is probviously decidable

**Hazard is basin mass at the entries the orbit actually uses.** A6's sweep map increases (S(W) > W, proved by the
Q-monotonicity of THM-4770(e)) and coalesces. So the fate of a window traversal is a function of its entry, and
sampling many entries says little, because they merge. The risk is read off the orbit's actual path:

- By THM-4770(d), an eps = 0 start of depth d can never equal a target.
- At even depths exactly one target is passed by eps = 1 sweeps.

**Three channels.**

| channel | event | probability per event | expected so far |
|---|---|---|---|
| Halting | exact landing on T_d, d ≡ 0 (mod 4) | ≈ 2/(gap sum) | ≈ 0.09 (d = 8..40) |
| Capture | exact landing on N_d, d ≡ 2 (mod 4) | ≈ 2/(gap sum) | ≈ 0.08 (d = 10..42) |
| Capture | long-jump landing inside a comb, d ≡ 3 (mod 12) | ≈ 0.057 | ≈ 0.17 (d = 15, 27, 39) |

Gap sums grow roughly linearly. So the first two channels give hazards ≈ 1/(2κd) per depth, with κ ≈ 11–12, and halting
is marginally ahead of N_d capture. The comb channel gives a constant ≈ 0.0047 per depth if the combs persist.

**Outcome.**
- The total hazard diverges, so A6 is decided at a finite (astronomical) time.
- If the combs persist, capture wins with probability ≈ 0.9–0.95.
- If they do not, the outcome is roughly even (halting slightly favoured), and the decision is typically beyond depth
  10^5.
- A capture means entering THM-4770(c)'s family. That is a certified infinite chain on the macro level, hence a proof
  of non-halting from the paper's Lean-verified rules. On the tape its period grows quadratically, so it is not a
  bbchallenge translated cycler.
- This is the same "probviously decidable" class as Racheline's tetrational Cryptid (3/5 translated cycler, 2/5
  halting), with a different non-halting mechanism.
