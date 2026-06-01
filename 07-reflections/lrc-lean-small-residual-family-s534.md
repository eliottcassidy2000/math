---
source: codex-2026-06-01-S534
status: formalization milestone
tags:
  - lonely-runner
  - lean
  - n3
  - residual-family
  - formalization
---

# Lean Small Residual Family for LRC at `n = 3`

This session pushed the current Lean LRC frontier one step past the denominator
sieve.  The existing `three_lonely_sieve_cover` proves `n=3` whenever one of the
two denominator traps is absent:

- no speed divisible by `3`, witnessed at `t = 1/3`;
- no speed divisible by `2`, witnessed at `t = 1/2`.

The new theorem `three_one_three_mul_lonely` handles an infinite family inside
the residual zone: `{1, 3r}` for every `r > 0`.

## The Witness

The proof uses the explicit time

```text
t = 1/3 + 1/(9r).
```

Runner `1` is safely inside `[1/3, 2/3]`, while runner `3r` lands at
fractional part `1/3`:

```text
3r * t = r + 1/3.
```

This is small, but it is the right kind of small: it attacks the part not seen
by the unit-denominator sieve.

## What Remains

To close full LRC at `n=3` in Lean, the cleanest path still appears to be S522o:

1. prove a scaling/time-change lemma for `Lonely`;
2. reduce two positive speeds to a primitive ordered pair `a < b`;
3. formalize the center-grid pigeonhole: the points
   `(a(2k+1))/(2b)` for `k = 0,...,b-1` form a full `1/b` grid;
4. show any closed interval of length at least `1/b` meets that grid.

THM-388 supplies a useful test case for that route: the witness is exactly a
controlled perturbation off the `1/3` wall, and the proof stays entirely in the
central-box language of THM-386.

## Artifacts

- Canon theorem: `01-canon/theorems/THM-388-lrc-n3-one-three-multiple-lean-family.md`
- Lean module: `04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean`
- Lean output: `05-knowledge/results/lean_lrc_small_residual_s533.out`
