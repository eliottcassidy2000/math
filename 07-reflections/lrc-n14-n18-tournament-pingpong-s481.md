# LRC n=14/n=18 Tournament Ping-Pong S481

This session treated the user's "back and forth many times" request as a
protocol: take an obstruction at `n=14`, translate it to `n=18`, inject one
forced-randomness idea when the scalar route stalls, then turn that into an
exact check.

The noise cards were not evidence.  They were prompts:

- Bruhat-Tits frontier: shrinking the real gap should push endpoint mass
  deeper, not erase it.
- Two-neighbor braid: nearest-runner distance is a projection of a left/right
  handoff pair.
- SC blowup: first-even rows create a twin sheet that hides inherited unit
  data.
- Basketball tie path: arbitrary labels are still a Hamiltonian path and
  should be declared.

The web search also nudged two directions.  The Quanta summary of recent LRC
progress emphasizes product-divisibility constraints and finite
computer-assisted sieves, while Jensen's May 2026 arXiv preprint on mixed
thresholds suggests studying unequal threshold vectors rather than only the
uniform `1/n` radius.  S481 used these only as inspiration: product rows became
bridge-fiber invoices, and mixed thresholds became a rank-like two-neighbor
control gauge.

## The Main Exact Comparison

The local gate invoice is unexpectedly parallel.

Cover endpoints owned by the `n`-gate using only lower columns `1..n-1`:

```text
n=14 forced=(1,3,5,7,9,11,13), exact=8, covers=6
n=18 forced=(1,5,7,9,11,13,17), exact=8, covers=2
```

Both rows require six unit residues, the half-gate, and one bridge.  But the
bridge fiber differs:

```text
n=14 free bridge: 2,4,6,8,10,12
n=18 free bridge: 6,12
```

So `n=18` is not a larger version of `n=14` locally.  It is locally more rigid.
The square core `9=3^2` collapses the bridge choice to `+-6 mod 18`.  That
suggests a proof route: solve a small bridge-fiber problem first, then push the
remaining difficulty into exported endpoint frontier mass.

## Gap-Debt Products

The row-parent, gate, and double-gate ladders preserve a product exactly:

```text
n=14:
  gap/th = 5/924, 5/1848, 5/3696
  debt   = 84, 168, 336
  product = 5/11 throughout

n=18:
  gap/th = 1/176, 1/352, 1/704
  debt   = 176, 352, 704
  product = 1 throughout
```

This is HYP-1866 in a cleaner comparative form.  The apparent scalar
improvement is exact endpoint export.  A counterexample would need both the
Archimedean gap and the exported endpoint frontier to vanish; these ladders do
the opposite.

The frontier mass is also stable inside each row:

```text
n=14 frontier mass = 66/7
n=18 frontier mass = 316/27
```

The `n=18` depth histogram is more informative than the raw debt:

```text
row-parent: {2:+0,3:+2}:96 {2:+0,3:+3}:16 {2:+4,3:+2}:64
gate:       {2:+1,3:+2}:192 {2:+1,3:+3}:32 {2:+5,3:+2}:128
double:     {2:+2,3:+2}:384 {2:+2,3:+3}:64 {2:+6,3:+2}:256
```

The `3`-depth is square-torsion payload; the `2`-depth moves down the row.  The
right statement is not just "debt doubles" but "the same square-core packet is
translated in the first-even row direction."

## Tournament Analysis Layer

S481 made four tournament gauges at the selected time for each row:

- semicircle orientation of positions;
- close-threshold switch over the base path;
- deletion-pressure tournament from two-neighbor relief;
- mixed-two-neighbor slack as a rank-like control.

The mixed two-neighbor gauge mostly collapses to a transitive order.  That is
useful negative evidence: simply giving each runner a two-neighbor slack score
throws away the cyclic data.  The edge-local gauges retain structure.

Examples:

```text
n=14 row-parent, close-threshold:
  close_flips=8, c3=38, scc=13,1, H=109371

n=14 row-parent, pressure:
  relief_ties=84, c3=14, scc=12,1,1, H=3601

n=18 gate, close-threshold:
  close_flips=10, c3=60, scc=17,1

n=18 gate, pressure:
  relief_ties=145, c3=22, scc=16,1,1
```

The proof object should therefore not be scalar nearest distance or scalar
two-neighbor slack.  It should be an edge-local tournament: either the
threshold-switch graph of close pairs or the deletion-pressure graph of who is
an irreplaceable blocker for whom.

## Working Lemma

First-even LRC rows `2m` should be attacked by a bridge-fiber lemma:

```text
unit/half fan rows force one bridge fiber,
and every bridge fiber either reopens a small-denominator row
or exports positive endpoint frontier mass under an edge-local tournament.
```

For `n=14`, the bridge fiber has six even choices and needs a fan-tax/product
depth certificate.  For `n=18`, the bridge fiber has only two choices, so the
local branch may be smaller.  The hard part is certifying that the square-core
endpoint packets cannot be made leafless under the pressure tournament.

## Next Checks

1. Enumerate all full speed sets containing an `n`-gate plus one local bridge
   fiber, with bounded filler columns, and rank them by endpoint frontier mass
   and pressure SCC size.
2. For `n=18`, branch separately on bridge `6` and bridge `12`; test whether
   one is just the antipodal/relabelled copy of the other after the base path.
3. Add a true mixed-threshold edge switch: each edge `{i,j}` uses a threshold
   derived from the two endpoint depths of `i` and `j`, not a vertex score.
4. Build the first bridge-fiber dual row matrix: rows are unit/half invoices,
   bridge invoices, small denominators, and exported endpoint-depth packets.

Sources used for outside-noise only:

- Quanta Magazine, "New Strides Made on Deceptively Simple 'Lonely Runner'
  Problem", March 6, 2026:
  https://www.quantamagazine.org/new-strides-made-on-deceptively-simple-lonely-runner-problem-20260306/
- Alathea Jensen, "Mixed thresholds in the Lonely Runner Conjecture",
  arXiv:2605.27941, submitted May 27, 2026:
  https://arxiv.org/abs/2605.27941
