---
id: HYP-1839
status: EXPLORATORY
source: codex-2026-05-31-S382
related:
  - HYP-1810
  - HYP-1813
  - HYP-1828
  - HYP-1837
  - HYP-1838
---

# HYP-1839: Fourteen-runner tightness and counterexamples require opposite gate behavior

## Statement

At threshold denominator `n=14`, the tight stratum and the counterexample
stratum pull in opposite directions:

```text
tight / boundary-only behavior  wants no 14-multiple gate;
open-cover counterexample       requires at least one 14-multiple gate.
```

The second clause is structural.  If no speed is divisible by `14`, then all
six unit points

```text
1/14, 3/14, 5/14, 9/14, 11/14, 13/14
```

remain safe.  Therefore a full open cover must contain a `14`-divisible speed.

The first clause is still experimental: known and locally found tight examples
keep those unit points as boundary witnesses, while inserting a `14`-gate
protects the unit layer but creates positive gaps or descendant endpoint debt.

## Evidence

`lonely_runner_14_gate_tightness_s382.py` separates unit safety, endpoint
ownership, and open-cover classification.

Residue action on the unit skeleton is binary:

```text
speed residue 0 mod 14 covers all six unit points;
every other residue covers none of them.
```

Canonical audits:

```text
initial segment: boundary_only, has14=0, unit_safe=6, unit_unprotected=6
initial replace 6 by 14: positive_gap, has14=1, unit_safe=0
seven-ladder: positive_gap, has14=1, unit_safe=0
S380 14-multiple ladder: positive_gap, has14=1, unit_safe=0
```

One-swap scan around the initial tight set with new speed `<=112`:

```text
without a 14-multiple: 1 boundary_only, 1182 positive_gap
with a 14-multiple:    0 boundary_only, 104 positive_gap
```

The one additional local tight example was

```text
(1,2,3,4,5,6,7,8,9,10,11,13,24),
```

again with no `14`-multiple and with the six unit witnesses intact.

Pure gate replacements `remove r, add 14q` for `q<=12` produced

```text
156/156 positive_gap.
```

## Predictions

1. Any primitive `k=13` full-measure set with no `14`-multiple is
   boundary-only, witnessed at the unit skeleton or a quotient of it.
2. Boundary-only examples with a `14`-multiple, if they exist, are not in the
   initial-segment tight basin and should have nonunit boundary witnesses.
3. Every disproof search must pay the `14`-gate toll first, then prove it can
   close the descendant endpoint debt at denominators such as `98`, `182`,
   `196`, and higher lcm layers.
4. A useful branch-and-bound split is:
   `no 14-gate -> cannot be counterexample`;
   `has 14-gate -> unit layer closed, charge new endpoint debt`.

## Sources

- `04-computation/lonely_runner_14_gate_tightness_s382.py`
- `05-knowledge/results/lonely_runner_14_gate_tightness_s382.out`
- HYP-1810
- HYP-1813
- HYP-1828
- HYP-1837
- HYP-1838
