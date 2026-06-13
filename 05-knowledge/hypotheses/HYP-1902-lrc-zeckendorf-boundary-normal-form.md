---
id: HYP-1902
status: OPEN
source: codex-2026-05-31-S451
related:
  - THM-370
  - HYP-1290
  - HYP-1881
  - HYP-1890
  - HYP-1900
  - HYP-1901
  - HYP-1903
  - HYP-1904
  - HYP-1910
---

# HYP-1902: LRC boundary recursion has a Zeckendorf normal-form shadow

## Statement

Zeckendorf decomposition should enter the Lonely Runner project as a normal
form for recursive endpoint debt, not as a direct reformulation of LRC.

The exact object behind Zeckendorf is:

```text
independent sets in the infinite path graph P_infty at fugacity x=1.
```

The tournament OCF object is:

```text
independent sets in the odd-cycle conflict graph Omega(T) at fugacity x=2.
```

The LRC anti-Bohr object from HYP-1901 is:

```text
protected boundary endpoints in a finite interval-cover hypergraph.
```

The hypothesis is that useful LRC proof branches can sometimes be compressed
to a path-like debt automaton.  In that automaton, "no adjacent selected
Fibonacci digits" becomes "no adjacent repair layers can both be spent without
exporting new endpoint debt."

## Evidence

S452/THM-370 adds a more literal Zeckendorf entry point: in a fixed runner
configuration, safe circular gaps form a subset of the cycle graph `C_n`, and
no-lonely means that subset is independent.  Cutting the cycle at an unsafe gap
turns the local obstruction into path independence, the exact graph behind
Zeckendorf normal forms.

`lrc_zeckendorf_bridge_s451.py` records the fugacity bridge:

```text
I(P_m,1): Fibonacci / Zeckendorf regime
I(P_m,2): Jacobsthal / tournament path-conflict regime
```

The same script puts S450 LRC boundary quantities into Zeckendorf normal form.
The useful examples are not numerological certainties, but they are suggestive
normal-form labels:

```text
initial n=14 debt  = 6   = F1 + F4
n=14 lower speeds  = 13  = F6
seven-ladder debt  = 84  = F5 + F7 + F9
S380 exported debt = 168 = F3 + F7 + F11
```

The seven-ladder debt is a tight gap-2 Zeckendorf chain; the S380 gate ladder
exports to wider Fibonacci-index gaps while doubling exposed endpoints.  This
matches the S450 moral: gates do not erase boundary debt, they move it into a
new coordinate system.

There is also a Diophantine-approximation interpretation.  Ostrowski
numeration attached to an irrational rotation gives the canonical expansion of
return times.  For the golden slope all continued-fraction digits are `1`, so
Ostrowski numeration becomes Zeckendorf decomposition.  Thus the golden
rotation is the cleanest model for anti-Bohr boundary recursion with a
no-adjacent-carry rule.

## Predictions

1. A useful LRC "Zeckendorf mode" should be an endpoint-debt automaton, not a
   speed-set identity.
2. The right graph to search for inside `n=14` and `n=16` repair branches is a
   path or Fibonacci-cube quotient of the endpoint incidence hypergraph.
3. Wythoff/Beatty index shift should model a clean debt-export operation: move
   every active Fibonacci index up by one, analogous to gate repair exporting
   debt to deeper product-depth rows.
4. If a branch-and-bound search can force no-adjacent carry in the endpoint
   layers, Zeckendorf uniqueness should give a canonical certificate rather
   than just another statistic.

## Sources

- `04-computation/lrc_zeckendorf_bridge_s451.py`
- `04-computation/lrc_runner_distance_tournament_s452.py`
- `01-canon/theorems/THM-370-lrc-two-neighbor-cycle-independence.md`
- `05-knowledge/results/lrc_zeckendorf_bridge_s451.out`
- `05-knowledge/results/lrc_runner_distance_tournament_s452.out`
- `07-reflections/lrc-zeckendorf-bridge-s451.md`
- `07-reflections/lrc-runner-distance-tournament-s452.md`
- `07-reflections/lrc-analogy-atlas-s450.md`
- `07-reflections/zeckendorf-non-consecutivity-pairing.md`
- `07-reflections/lucas-summand-graph-zeckendorf-geometry.md`
- `07-reflections/eureka-zeckendorf-simplex-cuboid.md`
