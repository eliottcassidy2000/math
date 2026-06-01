---
source: codex-2026-06-01-S540
status: computational probe plus proof-shape synthesis
tags:
  - lonely-runner
  - nowhere-zero-flow
  - support-flow
  - cut-flow
  - endpoint-protection
  - fourier
  - directed-graphs
---

# LRC Support-Flow and Cut-Flow Duality

The upstream HYP-2025 result gave the clean dictionary:

```text
full-support resonance = nowhere-zero flow on the speed-weighted dipole
runner positions       = tension/circular coloring on the observer star
```

That is a serious reframing, but it also names its own limit.  Full-support
inside debt is only one layer of the Fourier expansion.  LRC is about the
whole safe measure, so this session asks what happens when every support is
kept.

## The Support-Flow Expansion

For speeds `v_1,...,v_{n-1}`,

```text
|SAFE| = sum_{m : sum_i m_i v_i = 0} prod_i ghat_n(m_i).
```

The support `S={i:m_i!=0}` is a nowhere-zero integer flow on the speed-weighted
sub-dipole induced by `S`.  So the full measure is not one flow polynomial.  It
is a layered flow enumerator over all sub-dipoles:

```text
support 0: independence main term
support 2: pairwise flow corrections
support 3: first inside-debt layer
...
support n-1: full-support flow of HYP-2025
```

This matters because AP wall examples are not cancelled by a single magical
full-support object.  Their cancellation is distributed across many support
layers.

In the `n=14` AP wall, even the tiny truncation `M=3` already has the shape:

```text
support 0:  +0.13480057
support 2:  +0.07770106
support 3:  -0.36665706
support 4:  +0.59500528
support 5:  -0.82435649
support 6:  +0.86803411
support 7:  -0.68375725
support 8:  +0.40982896
...
```

That is not a parity bit.  It is an alternating support-current ledger.

## The Cut-Flow Side

There is a dual directed graph that is more geometric.  Each danger interval

```text
B_v = {t : ||v t|| < 1/n}
```

is a directed arc on the time circle.  The endpoints of all arcs cut the circle
into cells.  On each cell, the coverage count is a nonnegative flow value.

Then:

```text
positive lonely interval  = a zero-flow cut on the cell cycle
open full cover           = a nowhere-zero cover flow on every open cell
```

The AP wall has open full cover, so every open cell has positive cover flow.
But it is still not a counterexample, because it has a closed witness on the
wall.  The directed-graph sign of that fact is that the strict endpoint
protection core peels to empty.  There is no self-sustaining protected
circulation.

The S540 pattern was uniform in the selected audits:

```text
AP wall:
  open SAFE = 0, min coverage = 1, zero open cells = 0, strict core = 0

open witness:
  open SAFE > 0, min coverage = 0, zero open cells > 0, strict core = 0
```

So the cover-flow side says a genuine counterexample must be stronger than AP:
not just nowhere-zero coverage on open cells, but a nonpeeling endpoint core
while the best margin stays below `1/n`.

## The Combined Obstruction

The two sides now give a sharper counterexample shape.

On the Fourier side, a counterexample needs complete cancellation:

```text
main term + all support-flow layers = 0
```

On the cover side, it needs no zero-flow cut and no wall witness:

```text
cover flow positive on every open cell
best margin < 1/n
endpoint-protection core nonempty
```

This is a much narrower object than "many resonances exist."  Many resonances
are normal.  HYP-2025 already explains that bridgeless dipoles carry
nowhere-zero flows almost automatically.  The hard condition is arranging the
flows so they cancel the main term while simultaneously sustaining a
nonpeeling cover circulation.

## Proof Route

The flow proof route should not try to show that nowhere-zero flows are absent.
They are usually present.  It should try to show:

```text
every total support-flow cancellation forces the cover-flow core to peel
or lands on the regular-polygon wall class.
```

Equivalently:

```text
support-flow cancellation => zero cut or wall witness.
```

That is the graph-flow version of LRC.

## Relation to Existing Threads

HYP-2025 gives the full-support NZ-flow dictionary.  HYP-2026 extends it to all
support layers.

HYP-2022/HYP-2024 make the sector/section side explicit: the cover-cell cycle
is the real-space side of the same Fourier data.

THM-379/THM-380 and the endpoint-pressure work are the right next formal
target: the crude `strict_core_intervals=0` should become a labelled SCC or
Farkas certificate.  The flow language says the desired theorem should be a
cut/circulation alternative.

## Slogan

LRC is not "there are no nowhere-zero flows."  There are plenty.

LRC is:

```text
nowhere-zero support flows cannot assemble into a total cancellation
without exposing a zero cut or collapsing to a wall-only peelable cover.
```

That feels like the right directed-graph mindset: not flow existence, but
flow-cancellation versus cut exposure.
