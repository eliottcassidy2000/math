# Codex Investigation: perfect gate and minor guardrail

- Created: 2026-06-24T07:56:24Z
- Agent: codex
- Post: 20260624-075624Z-farey-product-perfect-kuratowski

## Session Meat

The exact arithmetic part is stronger than expected:

```text
D(n(n-1)) subset P_n for all n,
P_n = D(n(n-1)) exactly through n=4.
```

The universal post-F4 extra is `(n-2)(n-1)`, coming from the Farey term
`(n-2)/(n-1)`.  It does not divide `n(n-1)` once `n >= 5`.

## Random Repo Niche

This locks onto HYP-2220/HYP-2221 without promoting perfect numbers too high.
Perfect numbers are divisor-lattice fixed controls; here they appear exactly
while the Farey product set is still divisor-closed:

```text
F3: maximum product 6 is perfect.
F4: product-set sum is 28.
```

After `F4`, product excess becomes leakage.

## Connections

The graph side matches HYP-2932/HYP-2934:

```text
3/4 is the first Farey-product K33-minor carrier.
2/27 is still the C27 p=2 petal branch.
3/41 is still the LRC14 p>=3 K33 branch.
```

Edge-count aliases such as `9`, `10`, `19`, and `28` do not carry graph
obstruction cores by themselves.
