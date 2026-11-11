# S93 pre-registration (opus-2026-11-11-S93), written and pushed before any lane result of the wave was read

The wave of S93 (five lanes: counting exponent, two-apex family, digit systems, seed slopes, thin-orbit bridge) was
launched at the commit that reserves THM-5141, THM-5142, HYP-9671, HYP-9672. This file is pushed while the lanes are
running and before the coordinator has seen any of their output. The Chan-Pak conjecture is OPEN and nothing below
is a claim about it.

## What is already known to the coordinator (in sample; scratch probes in `04-computation/experiments/s93_chanpak_wave_20261111/seed/`)

- Two-apex family `G(2,0)` (spine with apexes `u -> v`, 36-state transfer count): the words of length at most 14
  attain every odd `N < 2^21` except 7 and 21; minimum multiplicity over the odd `N` of `[2^18, 2^19)`: 262, of
  `[2^19, 2^20)`: 86, of `[2^20, 2^21)`: 6 (lengths at most 14; maximal value at length 14 about `2^24`).
- One-apex words: colliding pairs per word with common prefix of length at least `t`, for `n = 18, ..., 26`:
  `0.96` (`t = 1`), `0.29` (`t = 3`), `0.148` (`t = 4`), `0.073` (`t = 5`), `0.037` (`t = 6`), `0.019` (`t = 7`),
  `0.009` (`t = 8`), `0.0046` (`t = 9`); share of distinct (prefix, value) classes at `t = ceil(n/4)`: `0.93` to `0.98`.
- `G(tau) = sum over ordered pairs (T, T')` of `tau`-bit strings with different last bits of
  `gcd(alpha, beta)/max(|alpha|, |beta|)`, `alpha = rev(T) - rev(T')`, `beta = T - T'`: `G(tau)/2^tau = 3.06, 4.42,
  5.39, 7.24, 8.48, 11.08, 12.73, 16.48, 18.63` for `tau = 4, ..., 12`, so `G(12)/6^6 = 1.636`.

## Predictions (out of sample)

- **W1 (two-apex completeness).** The two-apex words of length at most 16 attain every odd `N` in `[2^21, 2^23)`.
  Refuted by a single odd `N` of that range with no word of length at most 16.
- **W2 (two-apex multiplicity).** Over the odd `N` of `[2^21, 2^22)`, the minimum multiplicity among the words of
  length at most 16 is at least 100. Refuted if it is below 100.
- **W3 (prefix families).** For `n = 28` and `n = 30`, the colliding pairs per word with common prefix of length at
  least `t` are within the intervals `[0.118, 0.185]` (`t = 4`), `[0.015, 0.024]` (`t = 7`), `[0.0037, 0.0058]`
  (`t = 9`) (that is, within a factor `1.25` of the values at `n = 26`), and the share of distinct (prefix, value)
  classes at `t = ceil(n/4)` is at least `0.95`. Refuted if any of the eight numbers falls outside.
- **W4 (the gcd sum).** `G(13)/6^(13/2)` and `G(14)/6^7` lie in `[1.2, 2.2]`. Refuted otherwise. (If `G(tau)` grew
  like `2^tau` times a polynomial the ratios would fall below `1.2`; if it grew like `4^tau` they would exceed `2.2`.)
- **W5 (the counting exponent; a prediction about a proof, not a number).** The prefix-family argument (the block
  identity, the linear collision equation `beta c - alpha P = gamma` in the lattice point `(P, c)`, and a count of
  the window pairs `(T, T')`) yields a rigorous bound `A(X) >= X^theta/(log X)^C` with some `theta > 2/3`. This
  fails if the lane and its verifier find a gap that cannot be closed with an explicit exponent above `2/3`; it is
  then recorded as FAILED, with the gap.

No tolerance will be changed after the fact. W1-W4 are evaluated with the coordinator's own scripts or with a
lane's output, whichever comes first, and the receipt is named in the session note.
