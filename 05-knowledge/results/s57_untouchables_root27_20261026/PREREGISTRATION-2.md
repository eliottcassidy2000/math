# S57 pre-registration 2: the 10^10 extension (committed before the 10^10 array is analysed)

opus-2026-10-26-S57. The exact census to `N = 10^10` (`untouchable_sieve_20261026.py 1e10`) was launched before this
file was written; none of the statistics below has been computed from it. The only number it prints on its own is
the count at `10^10`, already public in OEIS A057978. Windows are `W_e = (2^(e-1), 2^e]`, even `k` only unless stated.
Values from `W_14` to `W_29`, measured to `10^9`, are recorded in `pp_residue_drift_1000000000.json`,
`pp_xbin_drift.json` and `pp_truncated.json`.

| id | statistic in `W_33 = (2^32, 2^33]` | 10^9 trend (last value) | predicted range |
|---|---|---|---|
| U-C | untouchables `<= 10^10` (control) | — | `1,681,871,718` (A057978) |
| U-P6 | untouchables in `W_33` divided by `2^32` (local density, all integers) | `+0.0009` per octave (`0.16601` at `W_29`) | `[0.1685, 0.1700]` |
| U-P7 | untouchable share, `x = k/s(k)` in `[1.8, 2.0)` (deficient end) | flat (`0.2256` at `W_30`) | `[0.215, 0.235]`; Pollack–Pomerance per-`x` value `0.141` |
| U-P8 | untouchable share, `x` in `[0.2, 0.4)` (abundant end) | `+0.005` per octave (`0.4909` at `W_30`) | `[0.495, 0.520]`; PP value `0.690` |
| U-P9 | max/min of the shares over the residues `k = 1,2,3,4 (mod 5)` | `-0.015` per octave (`1.500` at `W_29`) | `[1.42, 1.48]`, argmax `1`, argmin `4` |
| U-P10 | max/min over the residues `1..12 (mod 13)` | `-0.0075` per octave (`1.222`) | `[1.17, 1.215]` |
| U-P11 | share(`k = 2 mod 3`) / share(`k = 1 mod 3`) | `+0.003` per octave (`1.084`) | `[1.085, 1.11]` |

Mechanism being tested (to be written up in the S57 note): Pollack–Pomerance's class model assigns the same
untouchability probability to every residue `r != 0 (mod l)`, because it depends on `k` only through `gcd(k, A_y)`.
The prime-extension identity `s(ap) = s(a)p + sigma(a)` instead makes each core `a` with `l` not dividing `s(a)`
remove the residue `r = sigma(a) (mod l)` (there `l` divides the would-be prime), and each core with `l | s(a)` serve
only `r = sigma(a)`. A fixed core's weight decays like `1/log k`, and the cores with `l` not dividing `sigma(a)` have
counting function `x/(log x)^(1/(l-1))` (Wirsing). The predicted spreads therefore decrease toward 1, first like
`1/log k` and eventually like `(log k)^(-1/(l-1))`. A spread that stops decreasing (U-P9 or U-P10 above its range)
counts against the mechanism. A deficient-end share that falls into PP's range (U-P7 below `0.2`) counts against the
claim that the per-`x` discrepancy is stable.
