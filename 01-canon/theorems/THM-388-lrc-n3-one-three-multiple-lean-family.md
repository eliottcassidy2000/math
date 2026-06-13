---
id: THM-388
name: lrc-n3-one-three-multiple-lean-family
status: PROVED (Lean-checked)
date: 2026-06-01
session: codex-2026-06-01-S534
depends_on:
  - THM-369
  - THM-386
lean: 04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean
---

# THM-388: A checked residual `n = 3` LRC family

For every positive natural number `r`, the two-speed family `{1, 3r}` has a
`3`-lonely time.  In the Lean module this is the theorem
`three_one_three_mul_lonely`.

## Lean statement

In words:

```text
0 < r -> exists t, Lonely 3 (1, 3r) t.
```

The explicit witness is

```text
t = 1/3 + 1/(9r).
```

## Proof

The witness lies in the central band for the first runner:

```text
1/3 <= 1/3 + 1/(9r) <= 2/3.
```

For the second runner,

```text
3r * (1/3 + 1/(9r)) = r + 1/3,
```

so its fractional part is exactly `1/3`.  THM-386 identifies this central-band
condition with the standard nearest-integer Lonely Runner predicate.

The Lean proof is axiom-clean: the audit output for
`three_one_three_mul_lonely` reports only `propext`, `Classical.choice`, and
`Quot.sound`, matching the existing THM-369/THM-386 foundation audit.

## Significance

The earlier machine-checked `n=3` result, `three_lonely_sieve_cover`, proves the
case whenever no speed is divisible by `3` or no speed is divisible by `2`.
The family `{1, 3r}` lies directly beyond that sieve cover because it contains a
multiple of `3`.  THM-388 is therefore a first checked step into the residual
`6`-entangled kernel that the paper proofs S522o/S526 already close
informally.

This is not the full `n=3` formalization.  The next Lean target is a scaling
lemma plus the one-dimensional center-grid pigeonhole from S522o, which should
turn this explicit-family proof into the complete two-runner theorem.

## Verification

- Lean file: `04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean`
- Audit output: `05-knowledge/results/lean_lrc_small_residual_s533.out`
