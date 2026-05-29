# Blindspot Triads

**Date:** 2026-05-29
**Source:** `04-computation/meta_blindspot_probe_s95.py`

The useful question was not "what else connects?" but "what have we kept in
separate rooms?"

The prose scan found that some common topics almost never meet:

```text
cartan_attention + feedback_arc
dark_even_graph + vitali_gap
chromatic + feedback_arc
feedback_arc + magnitude_reachability
dark_even_graph + feedback_arc
```

That is a better flashlight than another broad synthesis. The repo has become
good at connecting Paley, OCF, Betti, Vitali, Krawtchouk, and tilings. It has
looked less at how order distance, graph-side sign, and continuous geometry
interact.

## Order and Entropy

The `H` versus min-FAS measurement says:

```text
n=6 exhaustive: corr(H,FAS)=0.8625
n=7 sample:     corr(H,FAS)=0.8482
```

The correlation is strong, but the functions are not equivalent. At `n=6`,
the same `H` can occur in different FAS strata:

```text
H=9  -> FAS 1 or 2
H=15 -> FAS 2 or 3
H=17 -> FAS 1 or 2
H=33 -> FAS 2 or 3
H=37 -> FAS 2 or 3
```

This suggests a cleaner language:

- `min-FAS` is the order parameter: distance from transitive structure.
- `H` is entropy inside an order stratum.
- The project may have been asking for a scalar formula where the right object
  is two-dimensional.

## Dark H Is Not Just Acyclic Orientations

The natural graph-side analog of `H` might have been the number of acyclic
orientations, `|chi_G(-1)|`. It fails as a separator:

```text
n=5: Royle-even/dark acyclic-orientation values overlap in 3 values
n=6: overlap in 22 values
```

At `n=6`, the median acyclic-orientation count is `94` on both sides.

So a dark analog of `H` probably needs the sign representation itself. The
candidate is not "count orientations"; it is "count orientations with edge-sign
character in the Burnside average."

## The Missed Room

The missed room appears to be:

```text
cycle-space quotient
  = tilings as GF(2) cycles
  + transitive distance / min-FAS filtration
  + Royle edge-sign representation
  + Krawtchouk low-pass constraints
```

The project has pieces of this:

- the cycle-space bijection from tilings to even graphs,
- the Krawtchouk/band-limited picture,
- the feedback-arc interpretation of metagraph diameter,
- the Royle-Praeger even graph equinumerosity.

But it has not yet made the feedback-arc filtration pass through the
tiling-to-even-graph bijection. That is the next natural measurement.
