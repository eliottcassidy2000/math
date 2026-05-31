# Lonely Runner Natural Gates

Session: `codex-2026-05-31-S377-natural-gates`

The user asked to reconnect the old natural-number "incomplete tournament"
idea with the 14/15-runner feedback loop.  The old additive rule

```text
x -> z and y -> z iff x+y=z
```

collapses to the complete transitive order once the second input is forgotten.
The multiplicative analogue stays sparse as the divisor DAG.  So the live
object is the labeled two-input cospan, and the family

```text
x_1+...+x_k = x_1*...*x_k
```

is the critical-pair interface between the additive and multiplicative gates.

`04-computation/lonely_runner_natural_gate_feedback_s377.py` uses that
interface as a disciplined way to generate new LRC probes.

## Fourteen-Runner Route

For `n=14,k=13`, the exact normalized one-defect scan found:

```text
product-sum target coordinates <=13: 4,6,8,9,10,12
best collision target missed: 56
best non-target missed:      154
```

The best coordinate is `6` with residue `7`, the same coordinate-6 half-turn
from S367/S371:

```text
coord=6, residue=7, missed=56.
```

This is no longer just a residue coincidence.  The coordinate `6` is the first
visible distinct product-sum resonance:

```text
1+2+3 = 1*2*3 = 6.
```

The second-best coordinate is `12`, another product-sum target with several
small cores:

```text
coord=12, residue=7, missed=84.
```

Thus the 56-cell hard case sits exactly on the first natural-mode critical
target.

## Fifteen-Runner Route

For `n=15,k=14`, the target coordinates are:

```text
4,6,8,9,10,12,14.
```

Again target coordinates dominate the blocker search:

```text
best collision target missed: 120
best non-target missed:      260
```

The best one-defects split into a tied order-3 pair:

```text
coord=6,  residue=5 or 10, missed=120
coord=14, residue=5 or 10, missed=120
```

This sharpens the S377 torsion/CRT dichotomy.  The `15=3*5` route should first
prove the order-3 stencils at the product-sum coordinates `6` and `14`, then
only afterwards attempt the full `3x5` micro-staircase.

## Disproof Route

The script also generated speed sets from product-sum target modes and divisor
payloads, then ran exact forbidden-interval and endpoint-peel checks.

```text
k=13,n=14: 300 operation-critical probes, all positive-gap
k=14,n=15: 300 operation-critical probes, all positive-gap
```

The tightest `n=14` operation-critical probe replaced `13` by `26` and still
had a positive gap:

```text
gap/thresh = 0.038462.
```

The tightest `n=15` probe was the usual first gate replacement `14 -> 15`:

```text
gap/thresh = 0.030303.
```

Even the deepest product-sum overloads leaked.  The deepest `n=15` probe had
peel depth `31`, but still a positive complement gap and empty endpoint core.

## Synthesis

The natural-mode graph is now doing more than supplying analogies.  It gives a
coordinate priority for the LRC quotient problem:

```text
scalar line -> product-sum target coordinates -> all other coordinates.
```

That is a useful finite proof order.  The target coordinates are where the
cell system is most fragile, and they are exactly where the additive and
multiplicative natural-number gate systems collide.

This became HYP-1833: product-sum natural gates mark the fragile LRC quotient
coordinates.
