# HYP-1816: Fourteen-runner CRT-gate descent

**Status:** EXPLORATORY proof-search hypothesis.

## Statement

In the reduced Lonely Runner Conjecture for `k=13` speeds, every
counterexample must pass the `14`-gate: it must contain at least one speed
divisible by `14`.

The hypothesis is that passing this gate forces a descent in the
endpoint-protection graph.  More precisely, any nonempty all-protected
endpoint core for a primitive `k=13` speed set projects along the
divisible-by-`14` channel

```text
W = {v/14 : 14 divides v}
```

to a smaller all-protected Bohr-boundary core.  Since `|W| <= 12` in a
primitive `k=13` set, this projected core should contradict the currently
proved finite frontier.

## Evidence

`04-computation/lonely_runner_fourteen_runner_s363.py` checks the first exact
signatures:

- The initial segment `(1,...,13)` is boundary-only with the unit skeleton
  modulo `14`.
- Replacing the final speed by the mandatory gate speed `14` creates a
  positive gap rather than a tighter obstruction:

```text
(1,2,...,12,14): gap/thresh = 0.045455
```

- In the exhaustive primitive `k=13`, `max_speed=16` box, all `455` sets that
  pass the `14`-gate are positive-gap:

```text
full_measure=0
open_cover=0
nonempty_cores=0
```

- The tightest gated positive example in that box has a long but empty peel:

```text
(1,2,3,4,5,7,8,9,10,11,12,13,14)
gap/thresh = 0.037879
peel_depth = 27
core_E = 0
```

Random `14`-gated candidates up to speed `80` also stayed positive-gap and
core-empty.

## Why It Matters

The latest published finite frontier proves the reduced conjecture through
`k=12`.  The next target is `k=13`, where `k+1=14` is composite.  This makes
the odd-prime polynomial sieve from the recent thirteen-runner work less
directly applicable, but it gives a CRT structure:

```text
14 = 2 * 7.
```

The quotient colorings `t=1/2`, `t=1/7`, and `t=1/14` are exact filters.  A
counterexample must defeat all of them.  HYP-1816 says defeating the strongest
one, `t=1/14`, opens a descent channel instead of creating a stable
obstruction.

## Next Tests

- Optimize endpoint peeling so larger `k=13` gated boxes can be scanned.
- Print residue data for the late peel layers of the tightest gated example.
- Define the projected core map from a `14`-protected endpoint core to the
  divided speed set `W`.
- Compare projected cores with Sungkawichai-Trakulthongchai improper residue
  tuples.

## See

- `04-computation/lonely_runner_fourteen_runner_s363.py`
- `05-knowledge/results/lonely_runner_fourteen_runner_s363.out`
- `07-reflections/lonely-runner-fourteen-runner-crt-gate-s363.md`
- HYP-1810
- HYP-1811
- HYP-1813
