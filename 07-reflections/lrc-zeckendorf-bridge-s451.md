---
source: codex-2026-05-31-S451
status: integration note
tags:
  - lonely-runner
  - zeckendorf
  - anti-bohr
  - ostrowski
  - endpoint-debt
---

# LRC and Zeckendorf: The No-Adjacent-Carry Lens

This note merges the Zeckendorf decomposition thread into the current LRC
analogy thread.  The shortest honest version:

```text
Zeckendorf is not LRC.
Zeckendorf is a canonical boundary-normal-form model for recursive debt.
```

The repo's earlier Zeckendorf work already found the important bridge:

```text
Zeckendorf: independent sets in P_infty, evaluated at x=1.
Tournament OCF: independent sets in Omega(T), evaluated at x=2.
```

HYP-1901 says exact LRC analogues must preserve the protected anti-Bohr
boundary.  Zeckendorf becomes relevant when that boundary recursion is
path-like: a repair at one layer forbids or changes what can happen at the
adjacent layer.  That is exactly the no-consecutive-digits rule.

## What Is Actually Shared?

The shared object is an independence/carry system.

```text
Zeckendorf:
  chosen Fibonacci digits form an independent set in a path.

Tournament OCF:
  chosen odd cycles form an independent set in a conflict graph.

LRC endpoint debt:
  chosen protectors try to cover endpoint rows without creating adjacent
  exported debt layers.
```

So the right LRC import is not "the numbers are Fibonacci."  It is:

```text
Can the endpoint-protection hypergraph be quotiented to a path or Fibonacci
cube in which uncovered rows have a unique no-adjacent normal form?
```

If yes, Zeckendorf uniqueness would turn a messy repair branch into a canonical
certificate.

## Golden Rotations

There is also a sharper Diophantine interpretation.  For a circle rotation with
irrational slope alpha, Ostrowski numeration gives the canonical expansion of
return times using the continued-fraction denominators of alpha.  When
alpha is the golden slope, every continued-fraction digit is `1`, and
Ostrowski numeration becomes Zeckendorf decomposition.

That makes the golden rotation the clean model of anti-Bohr boundary recursion:

```text
continued-fraction denominators = Fibonacci layers
no adjacent Ostrowski carries = Zeckendorf
closest-return gaps = anti-Bohr boundary events
```

LRC uses rational speed vectors, not one irrational slope, but the recursive
structure is the same kind of question: what happens to boundary debt when a
near-miss is repaired and pushed to the next denominator layer?

## S451 Computation

`04-computation/lrc_zeckendorf_bridge_s451.py` records two tables.

First, the path fugacity bridge:

```text
I(P_m,1) = Fibonacci/Zeckendorf regime
I(P_m,2) = Jacobsthal/tournament path-conflict regime
```

The notable old obstruction is still visible:

```text
I(P_4,2)=21
```

Second, the S450 LRC boundary quantities get Zeckendorf normal-form labels:

```text
initial n=14 debt  = 6   = F1 + F4
n=14 lower speeds  = 13  = F6
seven-ladder debt  = 84  = F5 + F7 + F9
S380 exported debt = 168 = F3 + F7 + F11
```

The seven-ladder debt is a tight gap-2 Zeckendorf chain.  The S380 gate ladder
doubles exposed endpoints but spreads them to wider Fibonacci-index gaps.  I
would not call that evidence of a Fibonacci theorem yet.  I would call it a
good coordinate label for the exact S450 moral: gates move debt; they do not
annihilate it.

## Proof Direction

S452 adds a sharper local version of the same idea.  If the actual runner
positions are ordered on the circle, then a runner is lonely exactly when the
two adjacent gap edges are both safe.  Thus no-lonely safe-gap masks are
independent sets in the cycle `C_n`; after cutting at an unsafe gap they become
path independent sets, the native Zeckendorf graph.  This does not replace the
endpoint-debt automaton below, but it gives a concrete place to test whether
the endpoint layers really collapse to a path/cycle normal form.

The bridge suggests a new kind of LRC certificate:

```text
endpoint rows -> recurrence graph layers
protectors -> carry moves
unprotected endpoints -> nonzero Zeckendorf normal form
counterexample -> impossible zero normal form after all carries
```

For `n=14`, combine this with HYP-1910:

```text
unit rows force a 14-gate,
the 14-gate pays a fan tax,
the repair exports product-depth debt,
the exported rows should have a no-adjacent normal form.
```

For `n=16`, the same idea should be dyadic first and Zeckendorf second:
dyadic depth gives the tree, while Zeckendorf labels the no-adjacent repair
schedule along a chosen path quotient.

## Practical Next Step

Build a graph from an endpoint incidence matrix:

```text
vertices = denominator/product-depth layers
edge i--j = a single protector can transfer debt between i and j
```

Then test whether the exposed debt subgraph in known near-counterexamples is
path-like, a forest, or a Fibonacci cube quotient.  If it is path-like, compute
the Zeckendorf normal form of the exposed row multiset and see whether every
repair is a legal carry or creates adjacent debt.
