---
id: THM-389
name: lrc-trienerment-tie-conventions
status: PROVED (observer equivalence Lean-checked; circle packing proof elementary)
date: 2026-06-01
session: codex-2026-06-01-S542
depends_on:
  - THM-381
  - THM-383
  - THM-386
  - HYP-2029
lean: 04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean
---

# THM-389: The LRC trienerment needs strict observer ties and weak global ties

Let `x_0=0,x_1,...,x_{n-1}` be the observer and runner positions on
`R/Z`, and let `d(i,j)` be circular distance.  There are two natural LRC
trienerment tie conventions:

```text
strict tie: d(i,j) <  1/n
weak tie:   d(i,j) <= 1/n
```

Then:

1. The observer is lonely in the usual closed LRC sense,
   `d(0,i) >= 1/n` for every runner `i`, iff the observer has strict
   tie-degree `0`.
2. Every `n`-point configuration has at least one weak tie.
3. A strict globally tie-free configuration exists iff the `n` points are a
   regular `n`-gon, i.e. the cyclic gaps are all exactly `1/n`.

Thus the S539 slogan must be read as a two-layer convention.  Strict ties give
the local observer equivalence.  Weak ties give the global pigeonhole fact.
With strict ties, the pure-tournament slice is not empty: it is exactly the
regular wall slice, including the AP/unit-wall witnesses.  This is why the
compactified symbolic target in HYP-2029 needs both open targets `G` and wall
targets `W`.

## Proof

For (1), a strict observer tie to runner `i` means

```text
exists m in Z such that |v_i t - m| < 1/n.
```

Negating this for every runner is exactly

```text
forall i, forall m in Z, 1/n <= |v_i t - m|,
```

which is the repo's Lean predicate `Lonely n v t`.  The Lean theorem
`lonely_iff_observerTieFree` proves this equivalence directly from
`not_lt`.

For (2), sort the `n` points cyclically and write the consecutive gaps as
`g_0,...,g_{n-1}`.  They are nonnegative and sum to `1`, so some gap is at
most `1/n`.  The two endpoints of that gap form a weak tie.

For (3), first suppose there is no strict tie.  Then no two points coincide,
and every consecutive cyclic gap is at least `1/n`; otherwise the adjacent
points would have strict distance below `1/n`.  Since the `n` gaps sum to `1`,
all gaps are exactly `1/n`, so the configuration is a regular `n`-gon.
Conversely, in a regular `n`-gon every nonzero circular distance is at least
`1/n`, so no strict tie exists.

## Consequences

- The LRC trienerment should keep a boundary flag, not only a ternary edge
  state.  The open edge state detects `dist < 1/n`; the compactified boundary
  state detects `dist = 1/n`.
- The observer-source tournament of THM-381 is the strict-observer-tie-free
  quotient of the same data.
- HYP-2029's `G/W` split is forced: `G` is an open strict-tie-free observer
  chamber, while `W` is a regular-wall or boundary witness that would be
  invisible if one only counted strict open chambers.

## Verification

Lean:

```text
04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean
```

Audit/computation:

```text
04-computation/lrc_trienerment_tie_conventions_s542.py
05-knowledge/results/lrc_trienerment_tie_conventions_s542.out
05-knowledge/results/lean_lrc_trienerment_s542.out
```
