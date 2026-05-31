---
id: HYP-1856
status: EXPLORATORY
source: codex-2026-05-31-S388
related:
  - THM-365
  - THM-366
  - HYP-1841
  - HYP-1853
---

# HYP-1856: Labelled endpoint cycles carry a positive slack potential

## Statement

Every labelled endpoint-protection cycle required by THM-365 carries a
strict-slack vector.  For an endpoint

```text
e = (n*m + eps)/(n*u)
```

protected by speed `p`, the integer label condition is

```text
|p*(n*m+eps) - a*n*u| < u.
```

Define the arrow slack

```text
s(e -> p) = u - |p*(n*m+eps) - a*n*u|  > 0.
```

The hypothesis is that for integer-realizable LRC endpoint systems, this
positive slack cannot close around a terminal protection cycle without leaking
one of the THM-366 small-denominator witnesses or exporting endpoint debt to a
higher denominator layer.

Equivalently, a true counterexample would need a labelled directed cycle whose
total slack is compatible with all denominator transitions.  The conjectural
proof object is a potential

```text
Phi(endpoint, owner speed, protector speed)
```

that strictly descends along every realized protection arrow unless the arrow
pays a small-denominator gate.  A closed cycle would then force an unpaid gate
or a positive gap.

## Motivation

S384 and THM-365 show that bare endpoint cycles are necessary for
counterexamples, but abstract circular-arc cycles are easy mirages.  S386
reframes labelled LRC cycles as circular versions of tournament good-cut
protection.  THM-366 adds a hard arithmetic boundary condition:

```text
small denominators m <= n require speeds divisible by m.
```

The next obstruction should combine these: a labelled cycle must not merely
exist; its slack must circulate through the denominator sieve without creating
a descending leak.

## Evidence

No direct theorem yet.  The circumstantial evidence is consistent:

```text
initial examples:        missing top denominator -> unit boundary skeleton
lcm repair examples:     sieve complete -> positive gap + endpoint debt
n14 seven-ladder:        sieve complete -> tiny gap, 84 unprotected, coreE=0
n14 14-ladder debt:      smaller gap -> 168 unprotected, coreE=0
```

Every sampled route toward closing the gap either misses a denominator gate or
exports endpoint debt before a nonempty endpoint core appears.

## Predictions

1. A cycle-first search should record arrow slack, owner denominator, protector
   denominator, and paid small-denominator gates.
2. Near-disproofs should contain almost-cycles whose slack sum fails by a small
   positive residue.  This residue should correlate with the critical-radius
   surplus from HYP-1843 better than with raw gap width.
3. The `n=14` seven-ladder's projected cycle-like behavior should fail because
   at least one slack edge points into a private endpoint leaf.
4. A formal proof may look like tournament SCC condensation: labelled
   protection arrows create a would-be core, but the slack potential induces a
   quotient-layer order that a directed cycle cannot respect.

## Test Plan

1. Extend `lonely_runner_endpoint_cycle_formal_s384.py` to enumerate short
   labelled protection paths and compute their slack vectors.
2. For `n=14`, search the seven-ladder and lcm-repair families for the longest
   slack-compatible almost-cycles before the first peel layer removes them.
3. Compare slack residue with endpoint private pivots from HYP-1842.
4. Try to prove a denominator-layer descent inequality for arrows whose
   protector speed is not paying a THM-366 gate.

## Sources

- THM-365
- THM-366
- HYP-1841
- HYP-1853
- `04-computation/lrc_small_denominator_sieve_s388.py`
- `05-knowledge/results/lrc_small_denominator_sieve_s388.out`
