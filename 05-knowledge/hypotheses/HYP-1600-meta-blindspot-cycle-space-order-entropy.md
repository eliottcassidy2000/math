---
id: HYP-1600
name: Cycle-Space Order/Entropy Blindspot
status: EXPLORATORY
source: codex-2026-05-29
related: [OPEN-Q-B10, dark-tournaments-and-the-sign-representation, what-we-can-and-cannot-know, meta_blindspot_probe_s95]
---

# HYP-1600: Cycle-Space Order/Entropy Blindspot

## Statement

Several under-bridged project threads may be projections of a single
cycle-space order/entropy structure:

1. `min-FAS(T)` is an **order parameter**: distance from the transitive class.
2. `H(T)` is an **entropy parameter**: number of Hamiltonian linearizations inside
   an order stratum.
3. Royle-even/dark graph parity is a **sign parameter**: the edge-space
   automorphism representation.
4. Krawtchouk band-limitedness is the **low-pass constraint** on all three.

The missing theory should live in the GF(2) cycle-space quotient rather than
in tournaments, graphs, or path complexes alone.

## Evidence

`04-computation/meta_blindspot_probe_s95.py` scanned project prose for common
but rarely co-mentioned topics. The highest under-bridged pair was:

```text
cartan_attention + feedback_arc: co=0, expected~1.2
```

Other under-bridged pairs included:

```text
dark_even_graph + vitali_gap
chromatic + feedback_arc
feedback_arc + magnitude_reachability
dark_even_graph + feedback_arc
```

This suggests that the project has studied "distance from order" and
"dark/even graph sign" mostly in separate channels.

## Small Measurements

### H versus min-FAS

All labeled tournaments at `n=6`:

```text
corr(H,FAS) = 0.8625
H range = 1..45
FAS range = 0..4
distinct H = 19
H values with multiple FAS = 5
score sequences with multiple (H,FAS) pairs = 11 / 22
```

Sample of 5000 labeled tournaments at `n=7`:

```text
corr(H,FAS) = 0.8482
H range = 1..189
FAS range = 0..7
distinct H = 77
H values with multiple FAS = 46
score sequences with multiple (H,FAS) pairs = 41 / 59
```

So `H` and `min-FAS` share a strong monotone component but neither determines
the other. Already at `n=6`:

```text
H=9  occurs with FAS in {1,2}
H=15 occurs with FAS in {2,3}
H=17 occurs with FAS in {1,2}
H=33 occurs with FAS in {2,3}
H=37 occurs with FAS in {2,3}
```

This argues against a direct formula `min-FAS = f(H)` and instead supports a
two-coordinate model: order stratum plus entropy within stratum.

### Dark-H Candidate

Royle-even graph counts were reproduced through `n=6`:

```text
n=3: Royle-even=2,  dark=2
n=4: Royle-even=4,  dark=7
n=5: Royle-even=12, dark=22
n=6: Royle-even=56, dark=100
```

The candidate graph-side analog `D(G) = number of acyclic orientations` is
not sufficient as a dark/light separator:

```text
n=5: acyclic-orientation value overlap count = 3
n=6: acyclic-orientation value overlap count = 22
```

At `n=6`, both Royle-even and dark classes have median acyclic-orientation
count `94`. Thus the dark analog of `H` should include the automorphism sign
representation, not just the unsigned chromatic-polynomial evaluation.

## New Hypotheses

### HYP-1600a

`H(T)` and `min-FAS(T)` are two low-frequency projections of the same
cycle-space structure. `min-FAS` is closer to the first Krawtchouk/order mode;
`H` contains the entropy residual inside each order mode.

### HYP-1600b

The graph-side "dark H" should be a signed acyclic-orientation count:

```text
D_sign(G) = signed Burnside average over acyclic orientations,
            weighted by the automorphism edge-sign character.
```

Unsigned acyclic orientations overlap too much between Royle-even and dark
classes to carry the distinction.

### HYP-1600c

The underexplored triad

```text
dark/even graphs + feedback arc distance + Krawtchouk band-limitedness
```

may be the right place to look for a natural Royle-Praeger bijection. The
cycle-space bijection already maps tilings to even graphs; the missing datum is
how the transitive-distance filtration transforms under that map.

## Next Tests

1. Compute `D_sign(G)` for graph classes at `n<=6`.
2. Push exact `H`/`min-FAS` class-level aggregation to unlabeled `n=7`.
3. Under the fundamental-cycle tiling-to-even-graph bijection, measure whether
   tournament `min-FAS` becomes a cut norm, signed chromatic statistic, or
   Krawtchouk mode on the graph side.
4. Test whether fixed `min-FAS` slices have narrower Krawtchouk bandwidth for
   `H` than the full tournament cube.
