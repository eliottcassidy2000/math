---
id: HYP-9380
title: "Heuristic classification of the eleven new BB(6) Cryptids of Dos Santos et al. (arXiv:2610.22603v2, A.1-A.11), by halting mechanism. (i) A11 is probviously HALTING: each countdown end (THM-4771(c)) halts with probability 1/4 under uniform end parities, and countdown lengths grow tetrationally; the blank-tape orbit has survived 4 ends and the 5th lies about 2^336 rounds ahead. (ii) A1, A3, A4, A7, A8, A9 are probviously non-halting with residual halting probability below 2^-(current bit size): their halting conditions are exact coincidences between exponentially growing parameters. (iii) A2, A5, A10 are probviously non-halting with residual about c/L (L = current scale, about 2^26 after 6e7 macro steps): halting needs a parity-breaking episode (about one per doubling of L, THM-4771(d)) followed by an exact hit inside it (probability about c/L). (iv) A6 is probviously decidable, not probviously non-halting: three channels with divergent total hazard (halting at the tops T_d, d = 0 mod 4; capture at N_d, d = 2 mod 4; capture combs at landings d = 3 mod 12), so it most likely becomes a translated cycler (0.75-0.9) but may halt (0.1-0.25) (HYP-9381)."
status: >
  OPEN / HEURISTIC (opus-2026-10-25-S56). Classification claims about unsolved halting problems; the structural
  inputs are PROVED (THM-4770, THM-4771). Measured inputs:
  - A11 end-parity uniformity: halting fraction 0.257/0.243 (seed 5) and 0.262/0.233 (seed 77, pre-registered P9)
    over random countdowns.
  - A2/A5/A10: 11, 5 and 10 dangerous episodes in 6e7 macro steps each, at log2-scales spread over 3..24, with
    smallest halting margins of order L/8 after the first.
  - A6: exact orbit to depth 44 with a mod-4 law of exposed targets; capture combs of density about 0.057 at the
    landings d = 15, 27, 39; none at d = 19, 23, 31. The 400-entry transfer samples coalesce, so they are one sample
    per depth.
  - Growth of the fast machines: long_run_60s.out; A1 by the certified engine, 6,711 bits at 17,241 segments (2400 s).
  - A11 directly on the rule system: from random countdown configurations the halting fraction at the first end is
    0.253 / 0.272 / 0.268 / 0.256 for (c mod 2, d) = (0,3), (0,6), (1,3), (1,6) (a11_rule_check.out).
  Falsifiers:
  - A halt of A1, A2, A3, A4, A5, A7, A8, A9 or A10 in any longer run.
  - A non-uniform end-parity law for A11 (a structural invariant forcing a0 = 1 or the wrong b parity at every end
    would make A11 non-halting).
  - For A6, see HYP-9381.
  Precedent for "probviously decidable/halting" tetrational Cryptids: Racheline's
  1RB1LE_0LC0LB_1RD1LC_1RD1RA_1RF0LA_---1RE (bbchallenge wiki).
session: opus-2026-10-25-S56
source: 05-knowledge/results/bb6_cryptids_odd_orbits_20261025/README.md
scripts:
  - 04-computation/bb6_cryptids_s56/a11_countdown_20261025.py (+ a11_countdown_P9.out), a11_rule_check_20261025.py (+ .out)
  - 04-computation/bb6_cryptids_s56/parity_episodes_20261025.py (+ parity_episodes_*.out)
  - 04-computation/bb6_cryptids_s56/long_run_20261025.py (+ long_run_60s.out), uncapped_residues_20261025.py
  - 04-computation/bb6_cryptids_s56/rs_deep_20261025.py (+ rs_deep_A1.out)
related:
  - 01-canon/theorems/THM-4770-bb6-cryptid-a6-phase-map-halting-tops-capture-family-certified-acceleration.md
  - 01-canon/theorems/THM-4771-odd-orbit-and-countdown-forms-of-bb6-cryptids-antihydra-a4-a11-a2-a5-a10.md
  - 05-knowledge/hypotheses/HYP-9381-a6-funnels-capture-windows-and-probviously-decidable.md
---

# HYP-9380 — classification of the eleven

| machine | mechanism | verdict (heuristic) |
|---|---|---|
| A1 | d = b + 2 at a cascade end; exponential growth per segment | non-halting, residual < 2^-6000 |
| A2, A5, A10 | double parity protection, slow growth | non-halting, residual about 2c/L, about 1e-7 now |
| A3 | Euclid-like race, K = M or small (K, M) | non-halting, residual < 2^-199000 |
| A4 | a = 2^(m+2) - m - 4 at visits (Antihydra core) | non-halting, residual < 2^-2000 |
| A6 | funnels; exposed tops T_d and points N_d; landing combs | **probviously decidable**: translated cycler 0.75-0.9, halt 0.1-0.25 |
| A7 | coincidence, then a parity | non-halting, residual < 2^-186000 |
| A8, A9 | coincidence of two exponentially large counters | non-halting, residual < 2^-270000 |
| A11 | parity check with probability 1/4 at tetrationally spaced countdown ends | **halting**, at a tetrational time |

"Non-halting" means probviously non-halting throughout.

**Why A11 differs.** In every other machine the halting test compares two large quantities exactly. In A11 it compares
two *bits* (the parity a0 and the parity of b) that a Hydra-type map scrambles over the whole countdown. So each
countdown end is a fair 1/4 chance. The countdowns themselves get tetrationally longer, which makes the halting time
astronomically large but, heuristically, finite.
