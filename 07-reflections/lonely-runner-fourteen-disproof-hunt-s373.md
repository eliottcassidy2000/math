# Fourteen-Runner Disproof Hunt S373

The session tried to disprove the fourteen-runner reduced case directly:
`k=13` moving speeds, threshold `1/14`.  By THM-357, a true counterexample is
an exact full open cover, not merely a large forbidden-measure set.  By
THM-360, any primitive counterexample must contain a speed divisible by `14`
to protect the unit endpoint layer.

## Construction Routes Tried

The new script `lonely_runner_14_disproof_hunt_s373.py` combined a fast
`12012`-point grid cover heuristic with exact S356/S360 Fraction audits.  It
tested:

- single `14`-gate replacements of `(1,...,13)`;
- two- and three-gate overloads;
- CRT residue lifts, including S367 half-turn residue patterns;
- anchored greedy set cover;
- mandatory-gate mutation search.

No exact open-cover candidate appeared.  The final exact audit had

```text
open_cover_candidates=0
audited_boundary_only=1
audited_positive_gap=69
```

## Best Near-Miss

The best positive-gap object was not random.  It was the anchored seven-ladder

```text
(1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91).
```

Exact audit:

```text
forbidden_length = 142/143
max_gap = 5/12936
gap/thresh = 0.005411
boundary_witnesses = 84
unprotected = 84
first_unprotected = 9/98
```

This is disproof-shaped because the positive gap is tiny, but it is not
protection-shaped: it leaves many exposed endpoints.  It looks like a
quotient-ladder over-cover of the `7` layer.

The best simple single-gate replacement was

```text
(1,2,3,4,5,7,8,9,10,11,12,13,224),
```

with

```text
max_gap = 11/9408
gap/thresh = 0.016369
unprotected = 12
first_unprotected = 29/182
```

So the two useful pressures split:

- quotient ladders minimize max gap but expose many endpoints;
- near-initial gate replacements expose few endpoints but leak larger gaps.

## Lesson

Speed-first disproof search keeps paying a toll: protect the unit layer, and a
higher-denominator leak opens.  The next creative route should reverse the
direction.  First build a small directed endpoint-protection cycle involving
the unit layer and the first leak layers (`9/98`, `29/182`, `15/182`), then
solve the integer protection inequalities for speeds realizing that cycle.

In other words, search for the protection graph first and the speed set
second.  That is the new disproof hypothesis, HYP-1828.
