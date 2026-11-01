---
id: THM-4856
title: "Three halts-iff-counterexample FRACTRAN programs built verbatim on PRIMEGAME's divisor loop. GOLDBACHGAME (37 fractions, 30 primes, start 43) halts iff some even N >= 4 is not a sum of two primes. FERMATGAME (36 fractions, 32 primes, started at m = 5) halts iff a sixth Fermat prime exists. LEGENDREGAME (33 fractions, 29 primes) halts iff some interval (n^2, (n+1)^2) contains no prime. Each consists of PRIMEGAME's first eleven fractions (exit 1/17 redirected) and a 22-26-fraction controller."
status: >
  Independently audited 2026-11-01: SOUND WITH CORRECTIONS, all applied
  (04-computation/experiments/circulant_fractran_goldbach_20261101/audit/AUDIT.md; MISTAKE-700). The audit
  re-ran all three programs literally and with its own accelerator, reproducing every step count.
  PROVED: the three equivalences, given the PRIMEGAME kernel contract (Conway 1987; THM-4853's proof; session126
  section 4).
  FINITE-EXACT (pre-registered).
  - GOLDBACHGAME:
    - no halt for every even N <= 1000, accepting the least-prime split each time (2.99e9 steps);
    - literal FRACTRAN = state dispatch on 4,119 configurations.
  - FERMATGAME:
    - started at m = 0..4 it halts after testing exactly F_m;
    - its inverted control passes F_0..F_4 and reaches the F_5 test;
    - literal = dispatch on 1,911 configurations.
  - LEGENDREGAME:
    - no halt for windows n = 1..120, accepting the least prime above n^2 (1.39e11 steps);
    - literal = dispatch on 1,396 configurations.
  FAILED AS STATED (reported, post-hoc anatomy in posthoc_controls.out):
  - G4: the control halted during N = 12, as predicted, but the check misread the target as 11, because a FRACTRAN
    test consumes the unit it tests.
  - L3: the control was mis-specified; truncating W also shortened the window jump. A corrected control halts at
    n = 5, as predicted.
  CREDIT:
  - That such encodings EXIST is classical and generic: every Pi_1 statement is a statement about some FRACTRAN
    program or generalized Collatz map (Conway 1972, 1987; Kurtz-Simon 2007). The repository's procgen atlas
    (2026-09-24) already lists Goldbach and Legendre this way.
  - New here are explicit, small, verified programs on PRIMEGAME's own kernel.
  - The audit's priority search (bbchallenge Fractran and Cryptid pages, Esolang Fractran) found no published
    FRACTRAN program for Goldbach, Legendre or Fermat primes.
  - Concurrent: session135 (`session135_fractran_tm_structure_20261101/`, checkpoint 1 at 08:23, two minutes after
    this file's pre-audit checkpoint) treats PRIMEGAME basins and priority-aware acceleration, not these
    encodings.
  NOT CLAIMED:
  - minimality of any program;
  - novelty outside the literature searched.
  - Description sizes are not invariants of the predicates (session65's warning about the 25-state machine applies).
source: mac-mini-2026-11-01-fractranTM (owner directive of 2026-11-01: FRACTRAN and Turing-machine connections with more repository questions; PRIMEGAME and the 25-state Goldbach machine)
depends_on:
  - 01-canon/theorems/THM-4853-primegame-step-counts-in-closed-form-oeis-a007546-is-the-15-14-list-a007547-the-conway-1987-list-and-their-difference-is-a-sum-of-largest-proper-divisors.md
related:
  - 01-canon/theorems/THM-4852-primegame-phase-law-both-published-lists-agree-selector-cellular-automaton-in-base-29-primorial.md
  - 01-canon/theorems/THM-4847-fermat-primes-barriers-to-proving-exactly-five-the-function-field-analogue-and-why-tournament-counts-are-blind-while-three-state-counts-are-pepin.md
  - 05-knowledge/results/session65_cryptid_tournament_interfaces_20261026/session65-twentyfive/REPORT.md
  - 05-knowledge/results/fractran_turing_circulant_goldbach_20261101.md
scripts:
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/goldbachgame.py
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/primekernel_programs.py
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/posthoc_controls.py
outputs:
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/goldbachgame.out
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/primekernel_programs.out
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/posthoc_controls.out
---

# THM-4856 - halts-iff-counterexample FRACTRAN programs on PRIMEGAME's divisor loop

## The kernel

PRIMEGAME's first eleven fractions are

```text
17/91, 78/85, 19/51, 23/38, 29/33, 77/29, 95/23, 77/19, X/17, 11/13, 13/11,
```

with Conway's exit 1/17 replaced by X/17 for a controller state prime X. They are a primality subroutine with this
contract (Conway 1987 section 5; THM-4853 proof; session126 section 4):

> From state 11 with (r2, r3, r5, r7) = (0, M-1, M, 0), M >= 2, the machine reaches state X with
> (r2, r3, r5, r7) = (M, 0, 0, D(M) - 1), where D(M) = M/spf(M). So M is prime iff r7 = 0.

## The programs

All three use fixed primes for their registers and states, and every fraction carries exactly one state prime in
its denominator. FRACTRAN cannot express a self-loop, since s·a/(s·b) reduces to a/b and loses the guard. So each
loop s -> s goes through a helper state s' with the single fraction s/s'. Every program is listed in its `.out` file.

**GOLDBACHGAME.** 37 fractions; registers 2, 3, 5, 7, P' = 31, Q' = 37, F = 41; 23 states; start integer 43 (the
state TP).

```text
17/91, 78/85, 19/51, 23/38, 29/33, 77/29, 95/23, 77/19, 53/17, 11/13, 13/11, 705/1333, 825/43, 43/47, 79/2173,
59/53, 1891/118, 67/56699, 59/61, 97/469, 71/67, 1095/2627, 33825/71, 71/73, 3071/158, 89/108151, 79/83, 97/623,
107/89, 101/679, 103/97, 97/101, 1333/3811, 113/103, 4033/3317, 58867/107, 107/109
```

The controller:
- holds the current pair as P' = P - 2 and Q' = Q - 2, with target N = P + Q;
- tests P, then Q, with the kernel, using the flag F to tell the two returns apart;
- on a composite, steps to the next split (P+1, Q-1);
- on success, moves to N + 2 with P = 2 (Q' := P' + Q' + 2);
- after the split (N-2, 2) fails, enters the halt state H, which has no fractions.

**FERMATGAME.** 36 fractions; registers K = 2^m, CK, K2, E, T; started at m = 5 from 53 · 31^32 · 43 (state CPK = 53, K = 2^5, E = 1).
- It doubles E K times to get 2^(2^m), tests M = E + 1, and halts at the first prime.
- Otherwise K := 2K and it continues.

**LEGENDREGAME.** 33 fractions; registers C1 = C - 1, W (candidates left after C), L = 2n, T.
- It tests C. A prime moves to the next window (C1 += W + 2, L += 2, W := L - 1).
- A composite moves to the next candidate, or halts if W = 0.

## Statement

1. GOLDBACHGAME, started at 43, halts iff some even N >= 4 is not the sum of two primes. Before halting, it accepts
   N = 4, 6, 8, ... in order, each with the least prime P such that N - P is prime.
2. FERMATGAME, started at m = m0, halts iff F_m is prime for some m >= m0. With m0 = 5, it halts iff there is a
   sixth Fermat prime.
3. LEGENDREGAME, started at 47 · 31 · 37 · 41^2 (state TEST = 47; C1 = 1, W = 1, L = 2, i.e. n = 1), halts iff some n >= 1 has no prime in (n^2, (n+1)^2).

## Proof

**One-hot control.** Every fraction has exactly one state prime in its denominator and one in its numerator, and the
start integer has exactly one state prime. By induction, every reachable integer has exactly one state prime, to
exponent 1. So the only fractions that can apply are those of the current state, in their listed priority order,
and FRACTRAN executes the compiled Minsky machine step for step. The helper states only return control.

**GOLDBACHGAME.** Invariant at each entry to TP:
- P + Q = N with P, Q >= 2;
- every split (P_0, N - P_0) with 2 <= P_0 < P has a composite member;
- every N_0 < N with 4 <= N_0 is a sum of two primes.

Each step preserves the invariant:
- **TP** sets (r5, r3) = (P, P-1) with r2 = r7 = 0, and enters the kernel.
- **At X with F = 0**, the controller moves r2 = P back into P' and removes 2.
- **If P is composite** (r7 > 0), it clears r7 and goes to NEXT.
- **If P is prime**, TQ tests Q the same way with F = 1.
- **NEXT** advances to (P+1, Q-1) if Q' > 0. Otherwise it is at P = N-2, so every split of N has failed, and it
  halts.
- **SUCC** (both prime) sets Q' := P' + Q' + 2 and P' = 0. That is N := N + 2 and P := 2.

So the program halts exactly at the first even N that is not a sum of two primes, and runs forever if there is
none. Each test terminates by the kernel contract, so there are no other non-halting traps.

**FERMATGAME.** At CPK with K = 2^m and E = 1:
- CPK sets CK := K and K2 := 2K;
- the OUT/DA/DB loops double E once per unit of CK, so E = 2^K;
- TM sets (r5, r3) = (E+1, E);
- the kernel decides whether F_m is prime;
- on a composite, the controller clears r7 and r2, sets K := K2 = 2^(m+1) and E := 1, and returns to CPK;
- on a prime, it halts.

**LEGENDREGAME.** Invariant at TEST:
- C lies in the window (n^2, (n+1)^2) with L = 2n;
- W = n^2 + 2n - C;
- every candidate of the window below C is composite;
- every earlier window has a prime.

On a prime, the next window starts at n^2 + 2n + 2 = C + W + 2 and has W := 2(n+1) - 1 candidates after its first.
On a composite, the controller moves to the next candidate, or halts when W = 0, i.e. when the window has no prime.

## Verification (pre-registered; see the status for G4 and L3)

| program | literal = dispatch | positive run | negative control |
|---|---|---|---|
| GOLDBACHGAME | 4,119 configurations | N = 4..1000, least-prime splits, 2.99e9 steps | P restricted to {2, 3}: halts during N = 12 (G4's check misread N; post-hoc confirmed) |
| FERMATGAME | 1,911 configurations (m = 2) | from m = 0..4: halts after F_m (F_4: 2.58e10 steps) | inverted: passes F_0..F_4, reaches F_5 = 4294967297 |
| LEGENDREGAME | 1,396 configurations (n <= 3) | n = 1..120, least primes, 1.39e11 steps | L3 mis-specified; corrected budget control halts at n = 5 (post-hoc) |

## Structure and size

**Kernel and controllers.** All three programs share the kernel, PRIMEGAME's 11 fractions and 6 states. Their
controllers have 26, 25 and 22 fractions. A repository question that asks "for every x, some y in a finite list is
prime" compiles to the kernel plus a controller that enumerates the list.

**Comparison with Leng's 25-state Goldbach Turing machine** (session65; `goldbach_tm25_structure.out`).
- **Size.** The machine has 49 instructions on 25 states and 2 symbols. GOLDBACHGAME has 37 fractions on 23 control
  states and 7 counters.
- **Shared design, no map (DICTIONARY).** Both search by trial division, but in different orders:
  - the machine tests divisors in ascending order and tries only splits x <= y;
  - PRIMEGAME's kernel tests downward from M - 1.

  No predicate-preserving map between them is constructed. The structural link is GOLDBACHGAME's literal embedding
  of PRIMEGAME's kernel, not anything about the Turing machine.
- **Cost per target.** The Turing machine spends about n^3.26 steps per stage (n <= 74). GOLDBACHGAME spends
  8.70-8.97 N^2 steps per target (8.97e6 at N = 1000), so about 9N^2. Its main term is the test of a prime Q near N,
  sum_d (6Q + 2 floor(Q/d) + 2) ≈ 6Q^2 (THM-4853 (a)); the remaining tests in the target account for the rest.
- **Sizes are not invariants.** These are description sizes only (session65): reachable padding changes them, and a
  21-state machine exists (DeepMind cryptids, 2026).
