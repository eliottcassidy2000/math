---
source: codex-2026-06-01-S513
status: formalization synthesis
tags:
  - lonely-runner
  - a000568
  - tournament-analysis
  - formalization
---

# Formalizing the LRC/A000568 results from S506-S512

This pass separated theorem-grade statements from the larger speculative proof
program that accumulated across the recent LRC Tournament Analysis sessions.

## Canonized now

**THM-381: observer-source reachability.**  The exact part of HYP-1981 is now a
theorem.  If the observer-to-runner arcs use the LRC threshold
`||v_i t|| >= 1/n`, then

```text
observer lonely at t  <=>  observer is a source.
```

The moving-runner arcs can still carry half-turn/A000568 structure, but they do
not affect the source test.  Deleting the marked source gives a bijection from
source-marked classes on `n` vertices to ordinary tournament classes on `n-1`
vertices, so the target set has size `A000568(n-1)`.

This is the clean formal repair to the S509 projection-defect result.  Raw
half-turn class loses the LRC bit.  The observer-source marked class recovers it
exactly.

**THM-382: bounded threshold-fiber separation.**  S512's small audit is now a
finite theorem, not just a hypothesis paragraph.  In the exact bounded state
spaces

```text
total n=3, moving speeds <= 16
total n=4, moving speeds <= 10
```

raw phase, observer-marked phase, gap-rank-only, and pair-deficit-rank-only
fibers are mixed and certify no sampled speed set.  Threshold-decorated
gap/pair fibers have no mixed classes and certify every sampled speed set.

The scope is intentionally finite.  The theorem does not claim an all-speed
classification.  It gives a reliable proof template: keep the A000568 base, but
put the `1/n` threshold colors in the fiber.

**THM-383: boundary compactification audit.**  S512 also becomes a finite
theorem about the half-turn class menu.  In the bounded probe, open phase
classes are small:

```text
n=3: 2/2
n=4: 2/4
n=5: 4/12
```

Adding equality-wall samples and the fixed tie Hamiltonian path expands the
menu:

```text
n=3: 2/2
n=4: 4/4
n=5: 11/12
```

This makes the boundary issue formal enough to stop hand-waving.  Tight LRC
states can live on equality walls, so an A000568 proof language needs either a
boundary compactification or a separate argument eliminating boundary-only
witness dependence.

## Left as hypotheses

The following claims remain open or empirical and should not be cited as canon.

The S506-S507 gauge work says LRC needs a bundle of metrics rather than one
tournament.  That is well supported computationally, but the useful gauge vector
is not yet a theorem.  It should remain HYP-1972/HYP-1973/HYP-1974/HYP-1975
until a finite certificate or general invariant is isolated.

The S509-S511 operation-grid work says static `x+2`/`x*2` labels become moving
metrics when weighted by current danger.  This remains HYP-1976/HYP-1980.  The
data are strong enough to drive experiments, not strong enough to become a
canon theorem.

The reachable-source-class claim remains the main unproved A000568 statement.
THM-381 proves the target size `A000568(n-1)`, but not which source classes are
reachable by arithmetic runner clocks, nor that every primitive clock reaches
one.

## Next formal targets

1. Define the reachable source menu: source-marked classes whose runner
   sub-tournament is realizable by `n-1` points inside an arc of length
   `1 - 2/n`.
2. Turn the S512 threshold-fiber audit into an exhaustive theorem for all
   primitive total `n=3` and `n=4` speed systems, or prove a reduction that
   makes the bounded windows complete.
3. Build the n=14 decorated-fiber corridor object:

```text
phase class
gap-threshold class
pair-deficit-threshold class
endpoint-pressure class
```

and prove a forced exit into a good-only class, or identify the exact mixed
fiber that blocks that route.
