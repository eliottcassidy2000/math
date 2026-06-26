# LRC14 Active-Bottleneck Normal Fan

codex-2026-06-26-S182

## What Changed

The barcode carrier was useful, but it still hid the local constraints.  This
pass keeps the active bottleneck owners of the lower envelope `m_S(t)` at each
bar's left endpoint, peak, and right endpoint.

The key output is a local certificate schema:

```text
left threshold owners -> peak bottleneck owners -> right threshold owners
```

For example, K33 `12->36` has bars supported by `(5,36)`, while petal
`10->20` and `P10+GW` have `(7,20)`, and petal `13->26` has `(1,26)`.
Those support pairs are better proof handles than raw `M`, safe mass, or a row
family name.

## Why It Helps

HYP-3016 and HYP-3017 show automatic words and residue-language fibers mix
boundary and open rows.  The normal fan separates them at the local certificate
level:

```text
AP/GW: zero bars with boundary peak pairs summing to 0 mod 14
mixed open rows: explicit open bars with peak supports such as (5,7) or (5,96)
```

This makes the next proof target sharper: prove that a packet's normal-fan
class either has an open bar, is the AP/GW zero-bar boundary atom, or routes to
a named state-lift/Fejer residual.

## Assumption Challenge

The vertices are not runners.  They are proof carriers and active
bottleneck-support records.  Alternate vertices considered were runners, gaps,
fixed circle sections, section boundaries, wall crossings, residues, cover
arcs, Fourier modes, matroid circuits, persistence bars, active bottleneck
sets, and proof obligations.

The quotient preserves the LRC predicate by retaining strict bars and destroys
global runner identity only after endpoint owners, peak owners, exact peak
height/time, and support residue sums have been stored.

## Next Pull

Run this normal-fan sidecar over the full HYP-2963 bank and test whether the
HYP-3017 mixed automatic-word fibers become boundary/open pure after adding
peak support and endpoint residue sums.
