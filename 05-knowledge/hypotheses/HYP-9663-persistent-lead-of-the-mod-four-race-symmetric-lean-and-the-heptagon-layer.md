---
id: HYP-9663
title: "Open statements about the mod-4 race of Hamiltonian path counts. (G) Persistent lead: D_n = P(H = 1 mod 4) - P(H = 3 mod 4) is positive at every order n >= 3 (proved for n <= 6 by the closed forms of THM-5132, exact for n <= 11, ten standard errors at n = 12; the stronger form D_n ~ c 4^(-n) was pre-registered and FAILED at n = 12: 989 +- 99 against about 2 * 10^3). (A) Symmetric lean: at every order n >= 4 the non-self-converse tournaments with a non-trivial automorphism have more labelled members with H = 3 than with H = 1 mod 4. (H) Heptagon layer (questions): do the exact identities E[{7}] = E[{3,5}] E[{3,5,7}] and E[{3,7}] = E[{3,5}] E[{5,7}] of order 10 hold at order 14; is there a description of the parity of c_7 that extends THM-5132; the fraction of bases of order 1 mod 4 with zero triangle-pentagon cube sum has a limit near 0.2206. The triangle-pentagon law first recorded here as (K) is PROVED (THM-5132 (5))"
status: >
  OPEN.
  (G): PROVED for n <= 6 (THM-5131 (T), THM-5132 (5)); the sign is FINITE-EXACT for n <= 11 and NUMERICAL at n = 12 (ten
  standard errors). The rate is OPEN: the pre-registered plateau D_12 4^12 in [1500, 2500] FAILED (989 +- 99 against
  1706, 2722, 1972, 2541, 1988 at n = 7..11).
  (A): FINITE-EXACT for 4 <= n <= 10 (signed sums in units of 4^-n: -64, -80, -216, -7.0, -104.2, -73.0, -54.1).
  (H): the two product identities are FINITE-EXACT at n = 10 only (one order; their analogues with c_7 = 0 are false
  at n = 6), so they are recorded as a question, not as a conjecture for all n = 2 mod 4; NUMERICAL for the zero
  fraction (0.2206 +- 0.0002 at m = 17, 21, ..., 37; 4,000,000 bases each).
  RESOLVED: the triangle-pentagon law (K) of the first draft of this file (E (-1)^(c_3 + c_5) = 2^(-(n-2)/2) and
  E (-1)^(c_5) = 2^(-(n-2)) for n = 2 mod 4), pre-registered as R7 and confirmed by sampling, is PROVED: THM-5132 (5).
  Two earlier forms were pre-registered and FAILED at n = 10 (R5, R6 of the packet).
  AUDIT: AUDITSTATUS
session: mac-mini-2026-11-11-race
source: 05-knowledge/results/mod4_race_square_term_20261111/README.md
related:
  - THM-5131
  - THM-5132
  - THM-5127
  - HYP-9659
  - HYP-9355
---

# HYP-9663 - the persistent lead of the mod-4 race, the symmetric lean, the heptagon layer

Definitions, proofs of what is proved and all tables:
[source note](../results/mod4_race_square_term_20261111/README.md), THM-5131, THM-5132.

## (G) Persistent lead

**Statement.** `D_n > 0` for every `n >= 3`: among labelled tournaments of every order, `H = 1 (mod 4)` is more
frequent than `H = 3 (mod 4)`.

**What is proved.** `|D_n| <= 1/2` (THM-5127 (C) (6)); for `n <= 4` the race is the triangle race and
`D_n = 2^(-floor((n-1)/2))` (THM-5131 (T)); for `n <= 6` it is the triangle-pentagon race and `D_5 = 13/32`,
`D_6 = 1/4` (THM-5132 (5)); for even `n` the self-converse tournaments contribute exactly
`SC_n/2^C(n,2) <= (n-1)!! 2^(-n(n-2)/4)`, all on the side of `H = 1` (THM-5131 (S)). Whether `D_n` tends to 0 is OPEN
(HYP-9659 (B'')).

**Evidence.** Exact: `D_n = 1/2, 1/2, 13/32, 1/4, 853/8192, 1361/32768, 63103/2^23, 650471/2^28, 260546929/2^39` at
`n = 3, ..., 11`. Monte Carlo with exact averaging over the last vertex (14.4 million bases of order 11):
`D_12 = (5.89 +- 0.59) * 10^(-5)`. The strong tournaments alone have a positive lead for `4 <= n <= 11`; the rigid,
non-self-converse tournaments for `5 <= n <= 10`.

**The rate (no conjecture).** In units of `4^(-n)` the lead is `1706, 2722, 1972, 2541, 1988` at `n = 7, ..., 11` and
`989 +- 99` at `n = 12`. The plateau was registered as a prediction for `n = 12` (R1, interval `[1500, 2500]`) and
FAILED. The decay is faster than `4^(-n)`; nothing more is known. The second-order quantities registered with it
held (rank law at order 11 within `3 * 10^(-4)` of the random alternating law; `rho_11 = 1.017`).

**Why the lead is not explained.** Forced coincidences - a source or sink (`chi(T) = chi(T - v)`), a gate
(`chi(Gamma R) = chi(R)`), twins - account for the second-order correlations of `chi` (THM-5127 (C) (5)), but they
contribute to `D_n` only through `D_(n-1)`, `D_(n-2)`, which are small: sources and sinks carry about `16 n 2^(-n)` of
`D_n`. The self-converse tournaments are the only family found on which `chi` is forced, and they die like
`2^(-n^2/4)`.

**Falsifier.** `D_n <= 0` at some order (order 13 is reachable by the Monte Carlo of the packet with about
`6 * 10^7` samples of order 12 if `D_13` is near `10^(-5)`).

## (A) Symmetric lean

**Statement.** For every `n >= 4`, among labelled tournaments on `n` vertices that are not self-converse and have a
non-trivial automorphism, those with `H = 3 (mod 4)` outnumber those with `H = 1 (mod 4)`.

**Evidence.** Exact for `4 <= n <= 10` (`out/race_tables.out`). No mechanism is known: an automorphism of odd prime
order `p` controls `c_odd` only mod `p`.

**Falsifier.** A positive signed sum at `n = 11` or `12` (needs the classes with a non-trivial automorphism only).

## (H) Heptagon layer

**Background.** THM-5132 solves the race for the cycle lengths 3 and 5: the two parities are weights mod 4 of the
vectors `s = A1` and `t = As`. The same construction fails for length 7 (`tr(A^3 + A^5 + A^7)` is not
`C(|A^3 1|, 2) + floor(n/4)`; `c_7` is not a trace, THM-500).

**Statements.**

1. (Question, one order of evidence.) At `n = 10`, exactly, `E[{7}] = E[{3,5}] E[{3,5,7}]` and
   `E[{3,7}] = E[{3,5}] E[{5,7}]`: the sign `(-1)^(c_3 + c_5)` is uncorrelated with `(-1)^(c_3 + c_5 + c_7)` and with
   `(-1)^(c_5 + c_7)`, and with no other sign involving `c_7` or `c_9`. Does this hold at `n = 14`? What is proved
   (THM-5132, remark to (5)): for `n = 2 (mod 4)` the conditional mean of `(-1)^(c_3 + c_5)` given the score-parity
   vector is the constant `2^(-(n-2)/2)`, so this sign is uncorrelated with every function of the score parities.
   The order-10 identities would follow from the same constancy given a finer invariant that determines
   `c_5 + c_7 (mod 2)`; none is known. At `n = 6`, where `c_7 = 0`, the analogous products are false.
2. (Conjecture.) For bases `R` of order `m = 1 (mod 4)`, the probability that the cube sum `G(R)` of THM-5132 (3)
   vanishes (that is: `A^T u = 1` is solvable over `F_2` and every solution has weight `0 mod 4`) tends to a limit,
   near `0.2206`.
3. (Question.) Is `c_7 (mod 2)`, or `c_3 + c_5 + c_7 (mod 2)`, a weight mod 4 of a vector built from the tournament?
   The direct continuation of THM-5132 (1) is false (`tower_test.py`).

**Evidence.** (1) `35881/2^24 = 2^(-4) * 35881/2^20` and `9019/2^22 = 2^(-4) * 9019/2^18` (`out/cycle_tables.out`).
(2) `0.2177` at `m = 13` and `0.2209, 0.2207, 0.2207, 0.2206, 0.2206, 0.2205` at `m = 17, 21, ..., 37`
(`out/k35_summary.out`); exact class fractions at `m = 5, 9`: 2 of 12, 43993 of 191536.

**Falsifier.** (1): the exact table at `n = 14` (not reachable by census; a sampler with the exact vertex-cube sum
of THM-5127 restricted to lengths at most 7 reaches it). (2): a drift of the fraction at `m = 41, ..., 61`.
