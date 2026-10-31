---
id: HYP-9448
title: "PRIMEGAME phase law: Conway's PRIMEGAME (and the owner's variant with 15/14 split into 15/2, 1/7) passes, at the end of testing each n >= 1, through exactly 2^n * 7^(D(n)-1), D(n) = n / spf(n) the largest proper divisor; hence the two programs emit the same primes forever; pre-registered verification to n = 2500"
status: OPEN for all n; pre-registered P1-P3 CONFIRMED to n = 2500 (2026-10-31)
session: opus-2026-10-31-S70
date: 2026-10-31
related: [THM-4852]
note: 05-knowledge/results/primegame_variant_selector_automaton_metagraph_20261031.md
script: 04-computation/experiments/primegame_variant_phase_law_20261031.c
---

# HYP-9448 - the PRIMEGAME phase law

**Statement.** Run PRIMEGAME (17/91, 78/85, 19/51, 23/38, 29/33, 77/29, 95/23, 77/19, 1/17, 11/13, 13/11, 15/14,
15/2, 55/1) from 2. Every time the configuration contains none of the state primes 11, 13, 17, 19, 23, 29, it equals
2^n * 7^(D(n)-1) for the next n = 1, 2, 3, ... in turn. Here D(n) = n / spf(n), with D(1) = D(p) = 1.

The same holds for the owner's variant, in which 15/14 is replaced by the two fractions 15/2 and 1/7 in that order.

**Why it matters.**
- From 2^n 7^(D(n)-1) the standard program applies 15/14 D(n)-1 times, 15/2 the remaining n - D(n) + 1 times, and then
  55/1.
- The variant applies 15/2 n times, 1/7 D(n)-1 times, and then 55/1.
- Both reach 3^n 5^(n+1) 11. The rules before 15/14 are identical and need a state prime.
- So the phase law makes the two orbits coincide at every phase boundary. In particular they emit the same primes.
  (The invariant actually needed is only 7-count <= 2-count at phase entry.)
- The law expresses Conway's trial division downward from n - 1: D(n) is the first divisor met (Conway 1987; Guy 1983).
  This file treats its all-n form as a hypothesis to be checked.

**Exploratory evidence (before this pre-registration).** It holds for both programs for every n <= 1000, and the 168
primes below 1000 appear in order with identical outputs (`primegame_variant_phase_law_20261031.c`, NMAX = 1000).

## Pre-registered predictions

Confirmatory run: `primegame_variant_phase_law 2500`.

- **(P1)** The phase law holds for both programs for every n <= 2500.
- **(P2)** The standard program emits exactly 2^p for the 367 primes p < 2500, in order.
- **(P3)** The variant emits the identical sequence.

Any failure is a refutation at that n, to be reported with the configuration.

## Results (confirmatory run, 2026-10-31, after the pre-registration commit adfb264520)

`primegame_variant_phase_law 2500` (9 min 40 s) printed `NMAX=2500 standard_ok=1 variant_ok=1 phases=2500/2500
powers_of_2=367/367 primes_exact=1 same_powers=1`.
- **(P1) CONFIRMED:** the phase law holds for both lists at all 2500 state-free configurations.
- **(P2) CONFIRMED:** Guy's list emits exactly the 367 primes below 2500 in order.
- **(P3) CONFIRMED:** the owner's list emits the identical sequence.

**Naming correction (found while writing the note).** The repository's PRIMEGAME check
(`procgen_atlas_20260924_primegame.py`) calls the owner's list (ending 15/2, 1/7, 55/1) the Wikipedia/OEIS form, and
the list ending 15/14, 15/2, 55/1 Guy's 1983 form. This file's 'standard' is Guy's form and its 'variant' is the
owner's. The simulator labels are unchanged.
