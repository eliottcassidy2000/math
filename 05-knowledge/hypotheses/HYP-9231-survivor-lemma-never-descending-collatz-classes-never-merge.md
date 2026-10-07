---
id: HYP-9231
title: "Survivor lemma: two distinct never-descending Collatz classes (parity words of equal length s and equal weight a with 3^(a_j) > 2^j for all prefixes) never collide mod 3^a, i.e. never merge with each other; consequently translation joins are always subsumed by backward-branch certificates in the maximal sieve (THM-4594)"
status: >
  OPEN (CONJECTURE). FINITE-EXACT: 0 collisions among all survivors for every s <= 30 (12.8M words at s = 30); joins subsumed
  by branches for every K <= 30. HEURISTIC: the elementary collision move shifts the j-th odd step by 2*3^j, which needs an
  even run too long for a survivor.
source: mac-mini-2026-10-07-golden (core reader), 05-knowledge/results/golden_collatz_resonance_20261007.md
related:
  - 01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md
  - 01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md (collisions of words)
scripts:
  - 04-computation/experiments/golden_20261007_readers/core/survcollide.c
---

# HYP-9231 — never-descending classes never merge

**Statement.**
* Call a parity word `w` of length `s` a **survivor** if every prefix has `3^(a_j) > 2^j`, i.e. the class never descends.
* For two distinct survivors `w ≠ w′` of equal length and weight `a`: `c(w) ≢ c(w′) (mod 3^a)`. Equivalently, their classes never merge at equal time.

**Why it matters.**
* If `n` merges with `n − d` at time `s`, and `n − d` descends at some `j < s`, then `T^j(n − d)` is a branch certificate for `n`.
* So a join that is not subsumed would need two survivors to merge.
* The hypothesis is exactly the statement that the classical translation joins (Roosendaal, Barina) are redundant once all backward branches are used.

**Evidence.** Exhaustive for `s ≤ 30`: 0 collisions among 12.8M survivors at `s = 30`.
