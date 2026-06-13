---
source: codex-2026-05-31 S363
status: exploratory
tags:
  - lonely-runner
  - fourteen-runners
  - distance-graphs
  - circular-coloring
  - crt
  - endpoint-protection
---

# Fourteen Runners: The CRT Gate

## Frontier

The current finite-checking frontier is the reduced Lonely Runner Conjecture
through `k=12`, i.e. thirteen runners including the stationary runner.  The
next natural case is therefore:

```text
k = 13 speeds, threshold = 1/14.
```

The recent Sungkawichai-Trakulthongchai paper proves `k in {10,11,12}` and
also gives a polynomial-method result for tuples congruent to
`(1,2,...,k) mod p` when both `k+1` and a large `p` are odd primes.  The next
case is different:

```text
k+1 = 14 = 2 * 7.
```

So the prime-denominator mechanism does not apply directly.  That is a hint,
not a bug: the fourteen-runner case should have a composite CRT proof shape.

## Distance-Graph Coloring View

For a speed set `D` with `|D|=13`, the problem asks for a multiplier circular
`14`-coloring of the integer distance graph `G(Z,D)`:

```text
c_t(x) = tx mod 1,
||d t|| >= 1/14 for every d in D.
```

The three obvious quotient colorings are:

```text
t = 1/2
t = 1/7
t = 1/14
```

To be a counterexample, `D` must defeat all three:

```text
t=1/2  forces at least one even speed,
t=1/7  forces at least one speed divisible by 7,
t=1/14 forces at least one speed divisible by 14.
```

The last condition is exactly THM-360's unit-boundary filter.  In endpoint
language: the unit witnesses `a/14` can only be protected by speeds divisible
by `14`.

## S363 Probe

I added:

```text
04-computation/lonely_runner_fourteen_runner_s363.py
05-knowledge/results/lonely_runner_fourteen_runner_s363.out
```

It reuses:

```text
S356 exact forbidden interval union
S362 endpoint/interval protection-core peeling
```

The script audits:

1. the tight initial segment `(1,...,13)`;
2. hand-built `14`-gated variants;
3. the exhaustive primitive `k=13`, `max_speed=16` box after the mandatory
   `14`-gate;
4. random `14`-gated candidates up to speed `80`.

## Findings

The initial segment is tight, exactly as expected:

```text
speeds=(1,2,...,13)
classification=boundary_only
unprotected=6
unit_skeleton=True
peel_depth=15
core_E=0
```

The six unprotected endpoints are the unit residues modulo `14`.

But inserting the mandatory `14`-gate breaks the tight skeleton instead of
making it more dangerous:

```text
speeds=(1,2,...,12,14)
classification=positive_gap
max_gap=1/308
gap/thresh=0.045455
unprotected=24
core_E=0
```

The exhaustive gated box was clean:

```text
k=13, max_speed=16
total_primitive=560
pass_14_gate=455
positive_gap=455
full_measure=0
open_cover=0
nonempty_cores=0
```

The tightest positive gated example in that box was:

```text
speeds=(1,2,3,4,5,7,8,9,10,11,12,13,14)
gap/thresh=0.037879
unprotected=8
peel_depth=27
core_E=0
```

This is the most interesting object from the run: not a counterexample, but a
long-leaking endpoint system.  It may be a good toy model for proving that the
`14`-gate creates a leak rather than an all-protected core.

Random `14`-gated candidates up to speed `80` also stayed positive-gap and
core-empty.  The best random gap ratio was about `0.047619`.

## Creative Proof Angle

The possible proof shape is:

```text
1. Unit skeleton:
   The initial tight family has witnesses at units a/14.

2. Gate:
   A counterexample must include at least one speed 14w.

3. Descent split:
   Speeds divisible by 14 protect the unit layer, but they also expose a
   lower-scale speed set W={w : 14w in D}.

4. Leak:
   Non-divisible speeds can disrupt some lifted residues, but the endpoint
   protection graph seems to peel anyway.  The computational signature is
   long finite peeling, not a stable core.

5. Target theorem:
   In the n=14 endpoint system, any self-supporting all-protected core must
   project to a smaller all-protected core for W.  Since |W|<=12 and the
   reduced cases through k=12 are known, the core cannot exist.
```

This is not yet a proof.  But it is sharper than a blind finite search: it
uses the exact place where fourteen differs from the previous prime-style
methods.

## New Hypothesis

HYP-1816:

```text
Fourteen-runner CRT-gate descent.

Every primitive k=13 counterexample must pass the 14-gate, but passing that
gate forces an endpoint-protection descent along the divisible-by-14 channel.
The terminal core would project to a forbidden core for a <=12-speed instance,
contradicting the proved frontier.
```

## Next Computation

1. Optimize S362 peeling for larger `k=13` boxes by avoiding `Fraction` in the
   inner containment loop.
2. For the tightest gated example
   `(1,2,3,4,5,7,8,9,10,11,12,13,14)`, print all 27 peel layers and identify
   the late surviving residues.
3. Implement the explicit projection from a hypothetical `14`-core to the
   divided channel `W={v/14 : 14|v}`.
4. Compare this with the Sungkawichai-Trakulthongchai projection/lifting
   language: their "improper residue tuple" should correspond to a nonempty
   endpoint core here.

## Sources

- Sungkawichai and Trakulthongchai, `arXiv:2604.23906`.
- Rosenfeld, `arXiv:2509.14111`.
- Jensen, `arXiv:2605.27941`, for mixed-threshold / safe-product inspiration.
- Distance-graph coloring surveys and regular chromatic number literature.
