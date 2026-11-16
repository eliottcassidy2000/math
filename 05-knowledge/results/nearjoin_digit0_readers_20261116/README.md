# Why the near-join exceptions are the counts without a strong digit-0 vertex; readers, divisors and overshoot; the field with one element

Session klein-2026-11-16-nearjoin (2026-11-16), on the owner's three questions: can the chosen reader break the
unwanted divisor without overshooting the absolute target (the design rule of `metagraph_flow_observers_20261116` §5);
why the counts with no digit-0 vertex are exactly the near-join exceptions; the field with one element. It begins with
an erratum to the previous packet (MISTAKE-776: a padded sink always has `B = H`, so the digit-0 statement is about
strong realisers). Canon: THM-5263 (new); HYP-9453 updated; THM-5260 (C) corrected.

## Status

- **PROVED** ([NEARJOIN.md](NEARJOIN.md), THM-5263). Three identities: the digit `σ_u = B_u − H` counts the two-path
  covers with a backward junction at `u`; for a sole exit `z` with port `t_0`, `σ_z = e_{t_0}(S − z)`; the near-join
  `T → M_k → z` (`z → t_0`) is the rooted append at a dead-end port, so its `z` is a digit-0 vertex. Hence every prime
  near-join value has a strong realiser with a digit-0 vertex, and **11, 13, 47, 67 are never near-join values from any
  base**: HYP-9453's seven exceptions are permanent. FINITE-EXACT: every digit-0 vertex of a strong tournament of order
  ≤ 8 has out-degree 1 (344 of 344), so strong digit-0 counts are the dead-end appends `h(B′) + B_{t_0}(B′)` with
  `e_{t_0}(B′) = 0`, which from bases ≤ 7 miss exactly `1, 7, 11, 13, 21, 47, 67` below 125.
- **Answered** ([READERS.md](READERS.md)). Homogeneous readers are module insertions (THM-5260 (A)); a homogeneous reader can break a divisor only if it does not divide the core's readout gcd, and in the +2 window it never does without overshooting: on every strong core of orders 3–7 the cheapest breaking reader costs at least +6 over the midpoint 3H and no single-cycle reader lands on 3H + 2 or 3H + 4; the +2 needs the inhomogeneous reader (the digit vertex).
- **Typed** (NEARJOIN.md §5). The `F_1`-tournament is the linear order; its sink is the universal digit-0 vertex; a
  strong digit-0 vertex is an internal sink; the near-join orbit of the `F_1`-point is the line `L_1`.
- Erratum: MISTAKE-776 recorded; FORCING.md §2–3 and THM-5260 (C) of the previous packet corrected in place.
- Pre-registration: [PREREGISTERED.md](PREREGISTERED.md) (P-1–P-4 confirmed).

## Files

- `scripts/`: `nearjoin_check.py`, `digit0_structure.py`, `probes.py`.
- `receipts/`: `nearjoin_check.out` / `.optimized.out` (identical), `digit0_structure.out`, `probes.out`.
- `prediction_nearjoin.txt`.
