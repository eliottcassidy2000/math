# Pre-registration (opus-2026-11-11-S92): the class law of F1 modulo odd prime powers on an unseen window

Registered and pushed before `odd_modulus_census.py eval 28 27 prereg_window_28_27.json` and
`odd_modulus_census.py cum 28` are run. Nothing of the window `[2^27, 2^28)` has been enumerated by this session for
any odd prime-power or composite modulus (the audit of S91 enumerated it for prime moduli only, and this session has
not seen class-level numbers from it).

## Object

F1: `N = x1 x2 + 2^s y1 y2 - [s = 0]` over ordered triples `(h1, h2, s)`, one-apex states `(x, y)`, `g = gcd(x, y)`.
Sheet 0 is `s = 0` (`c_0 = -1`), sheet 1 is `s >= 1` (`c_1 = 0`). For an odd modulus `m` and a divisor `d` of `m`, a
triple with `h1 != h2` is of type `d` if its state of smaller word has `gcd(g, m) = d`; `F_d` is the share of the
triples of type `d`. Statement under test (THM-5134 (B), draft): the class `c` modulo `m` has share
`sum over d | gcd(m, c - c_sigma) of (d/m) F_d`, up to a power-saving error. So all classes with the same
`g0 = gcd(m, c - c_sigma)` have the same class factor `m * share = Phi(g0) = sum_{d | g0} d F_d`.

## A-priori prediction

`prereg_window_28_27.json` (committed with this text) holds, for the moduli `9, 27, 81, 25, 125, 49, 45, 75` and both
sheets, the predicted `F_d` and `Phi(g0)` on the window of odd `N` in `[2^27, 2^28)`. They are computed without
enumerating the window: for a smaller state `(h1, s)` the number of second states `h2 > h1` with value in the window
is modelled by `W = int_{1/2}^{1} 2 d rho (1/2) max(0, B/f - max(A/f, y1))`, `f = x1 rho + 2^s y1` (the second state
has `x2/y2` uniform on `[1/2, 1]`, THM-5128 (B3)), and `F_d` is the `W`-share of the smaller states of type `d`.
Predicted class factors, sheet 0 / sheet 1:

    m = 9:   g0 = 1: 0.8577 / 0.8552    g0 = 3: 1.2143 / 1.2156    g0 = 9: 1.4249 / 1.4377
    m = 27:  g0 = 1: 0.8577 / 0.8552    g0 = 3: 1.2143 / 1.2156    g0 = 9: 1.4198 / 1.4335    g0 = 27: 1.4350 / 1.4461
    m = 25:  g0 = 1: 0.9418 / 0.9393    g0 = 5: 1.2308 / 1.2412    g0 = 25: 1.2409 / 1.2496
    m = 49:  g0 = 1: 0.9815 / 0.9837    g0 = 7: 1.1098 / 1.0967    g0 = 49: 1.1188 / 1.1051

(limits of THM-5134 (D): `0.8889, 1.1852, 1.2963` at 9; `0.9600, 1.1520, 1.1920` at 25; `0.9796, 1.1195, 1.1399` at 49).

## In-sample calibration

Windows `[2^23, 2^24)`, `[2^24, 2^25)`, `[2^25, 2^26)`, `[2^26, 2^27)` (`window_insample_*.out`). Largest
`|group mean - prediction|` over the four windows: `0.031` (9), `0.035` (27), `0.035` (81), `0.008` (25), `0.0095` (125),
`0.151, 0.129, 0.118, 0.107` (49, decreasing), `0.054` (45), `0.057` (75). Largest deviation of a class from its group
mean: `0.023` (9), `0.026` (27), `0.031` (81), `0.0053` (25), `0.0098` (125), `0.071` (49), `0.047` (45), `0.046` (75).
In all four windows the measured factor of every group with `g0 > 1` at `m = 9, 27, 25` is closer to the prediction
than to the limit (14 of 14). The ordering of the factors along divisor chains FAILS in sample (for instance at 25,
sheet 1, the classes with `v = 1` have `1.271` and the class `v >= 2` has `1.261`): it is not registered.

## Registered tests (window `[2^27, 2^28)`)

- **P1.** For each of the 8 moduli and each sheet, `|group mean - prediction| <= tol_m` for every group, with
  `tol = 0.04 (9), 0.04 (27), 0.045 (81), 0.008 (25), 0.009 (125), 0.115 (49), 0.045 (45), 0.035 (75)`.
- **P2.** For each modulus and sheet, every class is within `flat_m` of its group mean, with
  `flat = 0.025 (9), 0.027 (27), 0.03 (81), 0.005 (25), 0.008 (125), 0.075 (49), 0.045 (45), 0.04 (75)`.
- **P3.** For every group with `g0 > 1` at `m = 9, 27, 25` (14 comparisons) the measured group mean is closer to the
  a-priori prediction than to the limit of THM-5134 (D).
- **P4 (cumulative census at `X = 2^28`).** The exact remainder of THM-5134 (B),
  `m max_c |#{N = c} - sum_{d | gcd(m, c - c_sigma)} (d/m) #type d| / #T`, is smaller at `X = 2^28` than at `X = 2^26`
  for at least 8 of the 12 pairs (modulus, sheet) with `m` in `9, 27, 81, 25, 125, 49`. Values at `2^26` (sheet 0 /
  sheet 1): `0.02356 / 0.02029`, `0.02934 / 0.02216`, `0.03395 / 0.02309`, `0.00367 / 0.00495`, `0.00470 / 0.00657`,
  `0.13044 / 0.12524`.

The tolerance at 49 (`0.115`) is below two of the four in-sample values: it tests the decrease observed in sample.
A failure of any test is reported as a failure. What the tests do not establish: any statement beyond the window, or
the limits themselves (at these sizes the factors are far from their limits, as THM-5128 (D2) explains for primes).
