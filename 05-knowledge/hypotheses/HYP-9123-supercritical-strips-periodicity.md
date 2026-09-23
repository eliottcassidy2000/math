---
id: HYP-9123
title: "Supercritical strips: no rational has a non-eventually-periodic parity vector of bounded discrepancy around a supercritical slope"
status: >
  OPEN HYPOTHESIS, the smallest open instance of the HARD class. For every
  slope alpha in (log_3 2, 1) and every width C, no rational with odd
  denominator has a non-eventually-periodic 3x+1 parity vector with
  |a_s - alpha s| <= C for all s. Width < 1 is PROVED (Theorem S,
  Sturmian words, audited). The critical slope alpha = log_3 2 with any
  width is PROVED for positive integers (in-house discrepancy theorem)
  and for rationals (Proposition B, which relies on an in-house sketch).
  No existing mechanism reaches width >= 1 at supercritical slopes: the
  orbit grows exponentially, so capacity fails, and the words have
  positive entropy, so repetition fails.
source: collatz-procgen-20260922 (mac-mini), transversality lane candidate C2
depends_on:
  - 05-knowledge/results/collatz_procgen_20260922_transversality_foundry.md
  - 05-knowledge/results/collatz_guards_20260921_discrepancy.md
  - 05-knowledge/results/collatz_procgen_20260922_foundry.md
---

# HYP-9123 -- supercritical strips (candidate C2)

**Refutation form:** a rational `x` (odd denominator), a slope
`alpha > log_3 2` and a width `C` such that the parity vector of `x` is
not eventually periodic and stays in the strip.

**Why this is the smallest instance.** Foundry v4 gives every all-orbits
target of the family the same missing mechanism: 2-adic non-integrality
or irrationality of Bernstein numbers on positive-entropy, supercritical,
non-repetitive words. The targets are Collatz divergence, the Periodicity
Conjecture, both E-SCC halves, Mahler and Erdős. Among the HARD words, the
bounded-discrepancy ones around a fixed supercritical slope are the most
regular that still have positive entropy (strip width `>= 1`).

**Wave-4 update ([HARD-class lane](../results/collatz_procgen_20260922_hard_class.md), audited).**
* For integers, C2 is **equivalent** to a *coupled* Z-number statement:
  no `L != 0` and strip word `w` with
  `frac(L 3^(a_s)/2^s) = frac(E_s(w))` for all `s`, where
  `E_s = -Phi_R(tail) > 0` (Proposition M). The decoupled Mahler relaxation
  implies C2, but it is vacuous for wide strips.
* C2 is PROVED on every strip word with `Dio > mu` (Theorem D) and on the
  square-swap words (Theorem Y).
* Random strip words have `Dio = 1` almost surely.
* "Smallest instance" is superseded: a single explicit zero-entropy word
  (HYP-9127) is already open.
