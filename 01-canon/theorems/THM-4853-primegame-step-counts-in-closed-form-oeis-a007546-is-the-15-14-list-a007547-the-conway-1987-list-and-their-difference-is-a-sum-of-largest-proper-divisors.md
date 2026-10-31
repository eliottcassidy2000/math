---
id: THM-4853
title: "PRIMEGAME step counts. (a) A proof of A. P. Heinz's 2009 per-phase formula for OEIS A007546 (Guy's 15/14 list). (b) The matching closed form for A007547 (the owner's list, Conway 1987 Theorem 1): each phase costs D(N-1) - 1 more steps. (c) Hence the difference law A007547(k) - A007546(k) = sum_{n=2}^{p_k - 1} (D(n) - 1), with D(n) = n/spf(n) the largest proper divisor. (d) Both closed forms equal the full OEIS b-files (5000 and 1000 terms). (e) Growth: A007546(k) ~ 2(1 - c) p_k^3 with c = sum_p p^-2 prod_{q<p} (1 - 1/q) = 0.330098...; the difference is ~ (c/2) p_k^2, an OEIS corollary (Eldar 2022, A032742)."
status: >
  PROVED (elementary, given Conway's divisor loop) + FINITE-EXACT. Independently audited 2026-10-31: SOUND WITH
  CORRECTIONS, all applied (04-computation/experiments/groups_of_five_primegame_20261031/audit/AUDIT.md; MISTAKE-698).
  After the audit the session found that the per-phase formula for A007546 is A. P. Heinz's Maple program on that OEIS
  entry (2009). It is now CITED; this file supplies its proof. The credit is recorded in MISTAKE-698.
  - CITED:
    - Conway 1987 (program, Theorem 1, the loop proof in section 5 Figure 2);
    - Guy 1983 (the 15/14 list, "fast version"; per opus-2026-10-31-S70's audited THM-4852 and OEIS; not read here);
    - OEIS A007546, including Heinz's 2009 per-phase formula and P. von Brömssen's 2026-04-30 comment naming its
      program;
    - OEIS A007547 (Hasler 2026: A007542(a(n)) = 2^prime(n)) and A007542;
    - OEIS A032742 (Eldar 2022: sum_{k<=n} a(k) ~ 0.165049... n^2, which is c/2);
    - session126 section 4: instruction words and block lengths for the owner's list.
  - New here:
    - the written proof of Heinz's formula;
    - the closed form for A007547;
    - the difference law, as an identity between the two OEIS sequences;
    - the constant 2(1 - c).
    The A007547 closed form and the difference law are one-line corollaries of the two reset words (Conway 1987
    section 5; opus-2026-10-31-S70 THM-4852 (1)) once the trial costs are known.
  - FINITE-EXACT:
    - the session checked the 34 OEIS data terms of each sequence and all 168 primes below 1000 by validated
      simulation;
    - the independent audit checked all 5000 terms of the A007546 b-file (p <= 48611) and all 1000 terms of the A007547
      b-file (p <= 7919), with the difference law on those 1000.
  - NOT PROVED here: anything about FRACTRAN beyond these two lists.
source: mac-mini-2026-10-31-fivegroups (owner directive of 2026-10-31: groups of five, PRIMEGAME over 2..29, the n = 5 merged metagraph)
depends_on:
  - "Conway 1987, FRACTRAN: a simple universal programming language for arithmetic, Theorem 1 and section 5"
  - 05-knowledge/results/session126_primegame_symmetry_20261031/primegame/fractran/FRACTRAN_PRIME_EXPONENT_REPORT.md
related:
  - 01-canon/theorems/THM-4852-primegame-phase-law-both-published-lists-agree-selector-cellular-automaton-in-base-29-primorial.md
  - 05-knowledge/hypotheses/HYP-9448-primegame-phase-law-2-to-the-n-times-7-to-the-largest-proper-divisor-minus-one.md
  - 05-knowledge/results/primegame_variant_selector_automaton_metagraph_20261031.md
  - 05-knowledge/results/groups_of_five_primegame_merged_metagraph_20261031.md
scripts:
  - 04-computation/experiments/groups_of_five_primegame_20261031/correspondence_tests.py
  - 04-computation/experiments/groups_of_five_primegame_20261031/primegame_decode.py
  - 04-computation/experiments/groups_of_five_primegame_20261031/audit/audit_primegame.py
outputs:
  - 04-computation/experiments/groups_of_five_primegame_20261031/correspondence_tests.out
  - 04-computation/experiments/groups_of_five_primegame_20261031/primegame_decode.out
  - 04-computation/experiments/groups_of_five_primegame_20261031/audit/audit_primegame.out
---

# THM-4853 - PRIMEGAME step counts in closed form

## Setting

The two published lists share their first eleven fractions:

```text
17/91, 78/85, 19/51, 23/38, 29/33, 77/29, 95/23, 77/19, 1/17, 11/13, 13/11,
```

They then end differently:

- **L_G**, Guy's 1983 list ("fast version"), ends `15/14, 15/2, 55/1`.
- **L_O**, the owner's list, ends `15/2, 1/7, 55/1`. This is Conway's 1987 PRIMEGAME (Theorem 1, printed p. 4), the
  "slow version".

**Notation.**
- Run either list from 2. s(k) is the number of FRACTRAN updates until 2^(p_k) first appears; the start 2 is step 0.
- A **phase entry** is the start 2, or a configuration produced by 1/17. The state-free configurations 15 and 30 on
  the orbit are not phase entries.
- spf(n) is the least prime factor, and D(n) = n / spf(n) is the largest proper divisor, with D(1) = 1.
- c = sum_p p^-2 prod_{q<p} (1 - 1/q) = 0.330098... is the mean value of 1/spf(n), since p^-1 prod_{q<p} (1 - 1/q) is
  the density of the n with spf(n) = p.

## Statement

**(a) Phase cost.** For N >= 2, the steps from the phase entry 2^(N-1) 7^(D(N-1)-1) to the next phase entry,
2^N 7^(D(N)-1), number

```text
cost_G(N) = 7N + 1 + 2 spf(N) + (N - 1 - D(N)) (6N + 2) + 2 * sum_{D(N) < d < N} floor(N/d)    in L_G,
cost_O(N) = cost_G(N) + D(N-1) - 1                                                            in L_O.
```

cost_G(N) equals Heinz's `f(N) = N - 1 + (6N + 2)(N - D(N)) + 2 * sum_{d=D(N)}^{N-1} floor(N/d)` (OEIS A007546, Maple,
2009): split off the term d = D(N).

**(b) Closed forms.**

```text
s_G(k) = sum_{N=2}^{p_k} cost_G(N),        s_O(k) = sum_{N=2}^{p_k} cost_O(N).
```

**(c) Difference law.** s_O(k) - s_G(k) = sum_{n=2}^{p_k - 1} (D(n) - 1).

**(d) OEIS.** s_G = A007546 and s_O = A007547.
- The OEIS entries attribute the programs (P. von Brömssen's comment on A007546).
- The match was verified here on the 34 data terms of each. The audit verified it on all 5000 terms of the A007546
  b-file and all 1000 of the A007547 b-file.
- So, on those ranges,

```text
A007547(k) = A007546(k) + sum_{n=2}^{prime(k) - 1} (A032742(n) - 1).
```

**(e) Growth.**

```text
s_G(k) ~ 2 (1 - c) p_k^3,      s_O(k) - s_G(k) ~ (c/2) p_k^2.
```

- The second asymptotic is Eldar's A032742 asymptotic, sum_{k<=n} a(k) ~ 0.165049... n^2, applied through (c).
- The owner's list is slower by a relative amount of order 1/p_k. At p = 997 the extra 1/7 steps are 1.2e-4 of the
  total.

## Proof

**Machine.** Both lists drive the same 4 counters r2, r3, r5, r7 and the same 7 control states (the state primes
11, ..., 29, or none). Their transition tables are identical except for the state-free row (session126 section 3;
`primegame_decode.out` section 3).
- The first eleven fractions act only when a state prime is present.
- At most one state prime is ever present.
- When one is present, only the first eleven fractions fire (the converse step, opus-2026-10-31-S70 THM-4852 (1)).

**Trial-division block** (Conway 1987 section 5; session126 section 4).
- Start in state 13 with (r2, r3, r5, r7) = (A, 0, R, d), where A + R = N and d >= 1.
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
- **Case R = 0.** 17/91, then 1/17, give the phase entry 2^N 7^(d-1) in 2 steps.
- **Totals.** Write N = qd + r.
  - If d does not divide N: q(4d + 2) + 4r + 2N + 2 = **6N + 2q + 2**. This is session126's 6n + 2q + 2.
  - If d divides N: q(4d + 2) + 2 = **4N + 2 N/d + 2**.

**Phase.**
- The phase starts at the phase entry 2^(N-1) 7^m, with m = D(N-1) - 1 <= N - 2.
- **State-free segment** (Conway's reset word; THM-4852 (1)).
  - L_G: 15/14 fires m times, 15/2 fires N - 1 - m times, then 55/1 fires. That is N steps.
  - L_O: 15/2 fires N - 1 times, 1/7 fires m times, then 55/1 fires. That is N + m steps.

  Both reach 3^(N-1) 5^N 11.
- **Reload.** The 11-29 loop on N - 1 units, then 13/11, costs 2N - 1 steps and reaches state 13 with
  (0, 0, N, N - 1).
- **Trials.** The candidates d = N - 1, N - 2, ... are tried in turn. Each non-divisor costs 6N + 2 floor(N/d) + 2.
- **First divisor.** The first divisor met is D(N), and it costs 4N + 2 spf(N) + 2.
- The phase ends at the phase entry 2^N 7^(D(N)-1).
  - By induction from 2 = 2^1 7^0, this is the phase-entry form of the phase law (Conway 1987 section 5 Figure 2;
    HYP-9448, CLOSED as KNOWN by MISTAKE-697).
  - The only pure powers of 2 are the phase entries with D(N) = 1, that is, N prime (Conway's theorem).
- Adding the parts gives (a): N + (2N - 1) + sum_{D(N)<d<N} (6N + 2 + 2 floor(N/d)) + 4N + 2 spf(N) + 2.

**(b), (c).**
- 2^(p_k) is the phase entry that ends the phase testing p_k, so (b) follows by summing (a).
- (c) is the sum of D(N-1) - 1 over N = 2, ..., p_k. Here D(1) - 1 = 0.

**(e).**
- **Counting lemma.** Integers with spf(n) = p have density delta_p = p^-1 prod_{q<p} (1 - 1/q). Dominated
  convergence (the tail over p > P contributes at most x sum_{p>P} p^-2) gives sum_{n<=x} 1/spf(n) = c x + o(x).
- **Partial summation** then gives sum_{n<=x} D(n) ~ c x^2 / 2 (this is Eldar's A032742 formula) and
  sum_{n<=x} n^2 / spf(n) ~ c x^3 / 3.
- **Leading term.** In (a), (N - 1 - D(N))(6N + 2) = 6N^2 (1 - 1/spf(N)) + O(N), and the floor sum is
  O(N log N). So sum_{N<=x} cost_G(N) = 2x^3 - 2c x^3 + o(x^3).
- (c) gives the second asymptotic.

**Numerical reading.**
- s_G over 2(1 - c)p^3: 1.0015 at p = 2999, and 1.00006 at p = 48611 (audit).
- Law over (c/2)p^2: 0.9932 at p = 997.

## Verification

- **Session** (`correspondence_tests.out`, pre-registered as P1, P2, P7 before the run).
  - **Simulators.**
    - A literal FRACTRAN simulator and a state-dispatch simulator agree configuration by configuration.
    - An accelerated simulator that jumps over 2-cycles agrees with the literal one on the first ten emissions.
  - **OEIS data.** The closed forms equal the 34 OEIS data terms of each sequence.
  - **Simulation to 1000.** For the 168 primes below 1000, the closed forms equal the simulated counts, and the
    difference law holds.
    - At p = 997, s_G = 1,333,910,908 and s_O = 1,334,073,859.
    - The difference is 162,951, as the law predicts.
- **Independent audit** (`audit/audit_primegame.out`, own code).
  - Every block cost, in 3540 (list, N, d) cases with N <= 60.
  - The phase entries, by literal integer FRACTRAN through 25 primes.
  - The closed forms against all 5000 terms of the A007546 b-file and all 1000 of the A007547 b-file.
  - The difference law on those 1000 terms.
  - The owner's orbit against the A007542 b-file.

## What is inherited and what is new

**Inherited.**
- The program and its correctness, including the phase law and the reset word: Conway 1987.
- The 15/14 list: Guy 1983.
- The per-phase formula for A007546: Heinz 2009 (OEIS Maple program, given there without proof).
- The attribution of the two programs to A007546 and A007547: OEIS.
- The asymptotic of sum D(n): Eldar 2022 (A032742).
- The block lengths for the owner's list (6n + 2q + 2 and 7n + 2q + d + 3, the latter including the next reset and
  reload): session126 section 4.
- The two-list reduction: opus-2026-10-31-S70 THM-4852 (1).

**New here.**
- A written proof of Heinz's formula.
- Its companion closed form for A007547.
- The difference law between the two OEIS sequences.
- The constant 2(1 - c) in the growth of A007546.
