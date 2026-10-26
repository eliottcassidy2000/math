# Pre-registration (opus-2026-10-25-S56), committed before reading the runs it refers to

Context: certified accelerators `a6_accel_20261025.py` (A6) and `rs_accel_20261025.py` (A1, A2, A5, A10);
three-point lemma for translated cycles. Runs already launched but not yet read: `rs_deep_{A1,A2,A5,A10}`
(2400 s each), `a6_orbit_deep_3000s`, `a6_transfer_6_30_400` beyond d = 27.

Predictions (each is falsifiable by the named output file):

- **P1 (no halts).** None of A1, A2, A5, A10 halts within its 2400 s accelerated run (`rs_deep_*.out`).
- **P2 (A6 depth).** The A6 orbit reaches depth >= 43 within 3000 s without HALT or CAPTURE
  (`a6_orbit_deep_3000s.out`).
- **P3 (A6 coalescence).** For d = 28, 29, 30 the transfer maps from entries 1..400 have at most 3 distinct exits
  for d = 0 or 2 (mod 4) and at least 10 distinct exits for d = 1 (mod 4) (`a6_transfer_6_30_400.out`).
- **P4 (A6 halting basins unreachable).** No entry 1..400 halts at any d in 28..30.
- **P5 (A6 parity protection).** In `a6_orbit_deep_3000s.out`, at odd depths >= 29 the orbit's
  nearest approaches to T_d and N_d come from sweeps whose first phase has eps = 0 (w odd); targets have
  eps = 1, so these near misses carry zero hazard. To be checked by recomputing the near-miss sweeps.
- **P6 (capture is the only non-returning mechanism).** Every certified capture found anywhere in this session
  (A6 or the generic engine) is a period-2 depth-1 chain of the form P(v0, 4 v0 + 29) / P(v0, 4 v0 + 33) for A6;
  for the other machines no capture is found within the runs.
