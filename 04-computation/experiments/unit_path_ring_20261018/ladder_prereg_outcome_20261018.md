# Outcome of the d = 8, 9, 10 clock-graded ladder pre-registrations (mac-mini-2026-10-18-unitpath)

Registered in `ladder_d8_prereg_20261018.json` (commit ab45d7d9a1) and `ladder_d9_d10_prereg_20261018.json` (commit
a23e1be37c), before any d >= 8 line was printed. Computed by the C port `clock_graded_ladder.c`. The port reproduces the
Python run on d = 1..7 digit for digit, and S41's exact E_8, E_9, E_10 to 10 digits. Outputs:
`clock_graded_ladder_c_d9.out` (d = 1..9, 56 s, 0.52 GB) and `clock_graded_ladder_c_d10.out` (d = 10, 123 s, 1.5 GB).
The Python run was stopped before printing d = 8.

| d | cross (x2 units) | lag 0 (twins) | lag 1 | lag 2 | lag 3 |
|---:|---:|---:|---:|---:|---:|
| 7 | -0.059204 | -0.101833 | +0.040117 | +0.003866 | -0.000596 |
| 8 | -0.043218 | -0.087732 | +0.038274 | +0.006118 | +0.000011 |
| 9 | -0.034613 | -0.078642 | +0.037308 | +0.008185 | -0.001904 |
| 10 | -0.028655 | -0.071510 | +0.036366 | +0.007725 | -0.001911 |

- **P1 FAILED**: T_8 = -0.0877 is outside [-0.106, -0.098].
- **P2 FAILED**: lag 1 at d = 8 is +0.0383, below d = 7's +0.0401.
- **P3 FAILED**: the other lags sum to +0.0445, outside [+0.055, +0.063]. P3 was implied by P1, since cross_8 was known in
  advance (audit G1, N7), so it is not an independent test.
- **Q1 FAILED**: T_9 = -0.0786 and T_10 = -0.0715.
- **Q2 FAILED**: lag 1 decreases, +0.0373 then +0.0364.
- **Q3 PASSED**: E_9 = 1.556295496813 and E_10 = 1.527640161230 agree with S41 a4.

**Reading.** The brief-2 heuristic, repeated in this session's note before the run, was that "the twin part does not
decay; positive odd lags grow against a persistent same-clock dissipation". It is refuted from d = 8 on:
- the same-clock part shrinks by a factor of about 0.9 per depth after d = 7;
- the lag-1 part peaks at d = 7 and then declines slowly.

The same-clock sign law (twin part < 0) still holds at every computed depth, d = 2..10.

The cross term remains a difference of two parts of comparable size, but both decay. Its d^-1.9 law (S41 values) is
not a cancellation of non-decaying parts. Notes updated accordingly; HYP-9325(B1) is untouched, since it concerns
finite-N first arrivals.
