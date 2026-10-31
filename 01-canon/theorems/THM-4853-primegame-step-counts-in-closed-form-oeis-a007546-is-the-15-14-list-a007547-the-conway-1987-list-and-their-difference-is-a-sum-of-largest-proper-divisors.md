---
id: THM-4853
title: "PRIMEGAME step counts in closed form. OEIS A007546 is the 15/14 list and A007547 the owner's list (Conway 1987, Theorem 1). A007547(k) - A007546(k) = sum_{n=2}^{p_k - 1} (D(n) - 1), with D(n) = n/spf(n) the largest proper divisor. The phase testing N costs 7N + 1 + 2 spf(N) + (N - 1 - D(N))(6N + 2) + 2 sum_{D(N) < d < N} floor(N/d) steps in the 15/14 list, and D(N-1) - 1 more in the owner's list. Growth: A007546(k) ~ 2(1 - c) p_k^3 and A007547(k) - A007546(k) ~ (c/2) p_k^2, with c = sum_p p^-2 prod_{q<p} (1 - 1/q) = 0.330098..."
status: >
  PRE-AUDIT (mac-mini-2026-10-31-fivegroups).
  PROVED (elementary): (a) per-phase cost, (b) closed forms, (c) difference law, (e) growth, each given the divisor-loop
  analysis of PRIMEGAME (Conway 1987; written out with instruction words and block lengths for the owner's list in
  session126 section 4, and recalled in the proof below).
  FINITE-EXACT: (d) both closed forms equal all 34 listed OEIS terms of A007546 and A007547 (retrieved 2026-10-31).
  The closed forms, a validated accelerated simulator and the difference law agree for all 168 primes below 1000.
  CITED: Conway 1987 (program and correctness); the 15/14 list (Guy 1983 as cited by opus-2026-10-31-S70 and by
  procgen_atlas_20260924; Endrullis, Grabmayer and Hendriks 2009 as checked by session126); OEIS A007542, A007546
  (P. von Brömssen's comment of 2026-04-30 names its 15/14 program), A007547 and A032742.
  NOT PROVED here: anything about FRACTRAN beyond these two lists.
source: mac-mini-2026-10-31-fivegroups (owner directive of 2026-10-31: groups of five, PRIMEGAME over 2..29, the n = 5 merged metagraph)
depends_on:
  - "Conway 1987, FRACTRAN: a simple universal programming language for arithmetic, Theorem 1 and its proof"
  - 05-knowledge/results/session126_primegame_symmetry_20261031/primegame/fractran/FRACTRAN_PRIME_EXPONENT_REPORT.md
related:
  - 01-canon/theorems/THM-4852-primegame-phase-law-both-published-lists-agree-selector-cellular-automaton-in-base-29-primorial.md
  - 05-knowledge/hypotheses/HYP-9448-primegame-phase-law-2-to-the-n-times-7-to-the-largest-proper-divisor-minus-one.md
  - 05-knowledge/results/primegame_variant_selector_automaton_metagraph_20261031.md
  - 05-knowledge/results/groups_of_five_primegame_merged_metagraph_20261031.md
scripts:
  - 04-computation/experiments/groups_of_five_primegame_20261031/correspondence_tests.py
  - 04-computation/experiments/groups_of_five_primegame_20261031/primegame_decode.py
outputs:
  - 04-computation/experiments/groups_of_five_primegame_20261031/correspondence_tests.out
  - 04-computation/experiments/groups_of_five_primegame_20261031/primegame_decode.out
---

# THM-4853 - PRIMEGAME step counts in closed form

## Setting

The two published lists share their first eleven fractions:

```text
17/91, 78/85, 19/51, 23/38, 29/33, 77/29, 95/23, 77/19, 1/17, 11/13, 13/11,
```

They then end differently:

- **L_G**, the 15/14 list, ends `15/14, 15/2, 55/1`.
- **L_O**, the owner's list, ends `15/2, 1/7, 55/1`. This is Conway's 1987 Theorem 1, printed p. 4, as checked against the
  primary scan by session126.

Run either list from 2. Let s(k) be the number of FRACTRAN updates until 2^(p_k) first appears, with the start 2 at
step 0. Write spf(n) for the least prime factor and D(n) = n / spf(n) for the largest proper divisor, with D(1) = 1.
Write c = sum_p p^-2 prod_{q<p} (1 - 1/q) = 0.330098...: the mean value of 1/spf(n), since p^-1 prod_{q<p} (1 - 1/q)
is the density of the n with spf(n) = p.

## Statement

**(a) Phase cost.** For N >= 2, call the steps from the state-free configuration 2^(N-1) 7^(D(N-1)-1) to the next
state-free configuration 2^N 7^(D(N)-1) the phase testing N. It takes

```text
cost_G(N) = 7N + 1 + 2 spf(N) + (N - 1 - D(N)) (6N + 2) + 2 * sum_{D(N) < d < N} floor(N/d)    in L_G,
cost_O(N) = cost_G(N) + D(N-1) - 1                                                            in L_O.
```

**(b) Closed forms.**

```text
s_G(k) = sum_{N=2}^{p_k} cost_G(N),        s_O(k) = sum_{N=2}^{p_k} cost_O(N).
```

**(c) Difference law.** s_O(k) - s_G(k) = sum_{n=2}^{p_k - 1} (D(n) - 1).

**(d) OEIS.** s_G = A007546 ("fast version") and s_O = A007547 ("slow version"), for all 34 listed terms of each.
Hence, for those terms,

```text
A007547(k) = A007546(k) + sum_{n=2}^{prime(k) - 1} (A032742(n) - 1).
```

For every k, (b) gives each sequence exactly in this form, if each OEIS entry is defined as the step count of its
stated program.

**(e) Growth.**

```text
s_G(k) ~ 2 (1 - c) p_k^3,      s_O(k) - s_G(k) ~ (c/2) p_k^2.
```

So the owner's list is slower by a relative amount of order 1/p_k. At p = 997 the extra 1/7 steps are 1.2e-4 of the
total.

## Proof

**Machine.** The two lists share one 4-counter, 7-state machine (session126 section 3; `primegame_decode.out`
section 3):
- the counters are the exponents r2, r3, r5, r7;
- the states are the state primes 11, ..., 29, or none;
- the first eleven fractions act only when a state prime is present, and at most one state prime is ever present.

**Trial-division block** (session126 section 4; Conway's loop).
- Start in state 13 with (r2, r3, r5, r7) = (A, 0, R, d), A + R = N, d >= 1.
- The 13-17 loop costs 2 steps per unit and moves one unit from r5 and r7 into r2 and r3.
- **Case R >= d.**
  - After d units comes 11/13.
  - The 11-29 loop moves r3 back into r7 at 2 steps per unit, then 13/11.
  - This round costs 4d + 2 and replaces R by R - d.
- **Case 0 < R < d.** After R units:
  - 17/91, then 19/51;
  - the 19-23 loop moves r2 into r5 at 2 steps per unit, adding 2(A + R) = 2N;
  - 77/19;
  - the 11-29 loop on R - 1 units, then 13/11.

  The machine is back in state 13 with (0, 0, N, d - 1). This part costs 4R + 2N + 2.
- **Case R = 0.** 17/91 and then 1/17 give the state-free configuration 2^N 7^(d-1), in 2 steps.
- **Totals.** Write N = qd + r.
  - If d does not divide N: q(4d + 2) + 4r + 2N + 2 = **6N + 2q + 2**. This is session126's 6n + 2q + 2.
  - If d divides N: q(4d + 2) + 2 = **4N + 2 N/d + 2**.

**Phase.**
- The phase testing N starts at 2^(N-1) 7^m, with m = D(N-1) - 1 <= N - 2.
- **State-free segment.**
  - L_G: 15/14 fires m times, 15/2 fires N - 1 - m times, then 55/1 fires. That is N steps.
  - L_O: 15/2 fires N - 1 times, 1/7 fires m times, then 55/1 fires. That is N + m steps.

  Both reach 3^(N-1) 5^N 11. This is opus-2026-10-31-S70's equivalence argument.
- **Reload.** The 11-29 loop on N - 1 units, then 13/11, costs 2N - 1 steps and reaches state 13 with
  (0, 0, N, N - 1).
- **Trials.** The divisors d = N - 1, N - 2, ... are tried. Every non-divisor costs 6N + 2 floor(N/d) + 2.
- **First divisor.** The first divisor met is D(N), and it costs 4N + 2 spf(N) + 2.
- The phase ends at 2^N 7^(D(N)-1). By induction from 2 = 2^1 7^0, this is the phase law (HYP-9448's statement).
  - The only pure powers of 2 are the phase ends with D(N) = 1, that is, N prime.
  - This is Conway's theorem, with intermediate pure-power events checked in session126 section 4.
- Adding the parts gives (a): N + (2N - 1) + sum_{D(N)<d<N} (6N + 2 + 2 floor(N/d)) + 4N + 2 spf(N) + 2.

**(b), (c).**
- 2^(p_k) is the end of the phase testing p_k, so (b) follows by summing (a).
- (c) is the sum of D(N-1) - 1 over N = 2, ..., p_k. Here D(1) - 1 = 0.

**(e).**
- **Counting lemma.** Integers with spf(n) = p have density delta_p = p^-1 prod_{q<p} (1 - 1/q). Dominated convergence
  (the tail over p > P contributes at most x sum_{p>P} p^-2 + O(log log x)) gives sum_{n<=x} 1/spf(n) = c x + o(x).
- **Partial summation** then gives sum_{n<=x} D(n) ~ c x^2 / 2 and sum_{n<=x} n^2 / spf(n) ~ c x^3 / 3.
- **Leading term.** In (a), (N - 1 - D(N))(6N + 2) = 6N^2 (1 - 1/spf(N)) + O(N), and the floor sum is
  O(N log N). So sum_{N<=x} cost_G(N) = 2x^3 - 2c x^3 + o(x^3).
- (c) gives the second asymptotic.

**Numerical reading.**
- s_G at p = 2999 over 2(1 - c)p^3: 1.0015.
- Law value at p = 997 over (c/2)p^2: 0.9932.

## Verification (`correspondence_tests.out`, pre-registered as P1, P2, P7 before the run)

- **Simulators.**
  - A literal FRACTRAN simulator and a state-dispatch simulator agree configuration by configuration (3876 and 3893
    configurations, first six emissions).
  - An accelerated simulator that jumps over 2-cycles agrees with the literal one on the first ten emissions.
- **OEIS.** The closed forms equal all 34 listed terms of A007546 and of A007547. The difference law holds for those
  34 terms.
- **Simulation to 1000.** Both lists emit the 168 primes below 1000 in order. The closed forms equal the simulated
  step counts, and the difference law holds, for all 168.
  - At p = 997, s_G = 1,333,910,908 and s_O = 1,334,073,859.
  - The difference is 162,951, as the law predicts.

## What is new and what is inherited

**Inherited.**
- The program and its correctness: Conway 1987.
- The two lists, and the repository's 15-prime check of both: procgen_atlas_20260924.
- The equivalence of the lists via the state-free segment and the phase law to n = 2500: opus-2026-10-31-S70, THM-4852
  (pre-audit) and HYP-9448.
- The instruction words and the block lengths 6n + 2q + 2 and 7n + 2q + d + 3 for the owner's list, with a direct
  divisor-loop proof: session126 section 4. That block includes the next state-free segment and reload. It equals the
  4N + 2q + 2 above plus (N + d) + (2N + 1).

**New here.**
- The summed closed forms for both OEIS sequences.
- The difference law between them.
- The attribution check over every listed OEIS term.
- The growth constants 2(1 - c) and c/2.

The phase-law proof above is Conway's loop argument, written out in session126 section 4 for L_O. HYP-9448 is
opus-2026-10-31-S70's file, and this theorem does not edit it.
