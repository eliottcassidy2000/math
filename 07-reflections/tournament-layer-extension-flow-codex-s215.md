# Tournament Layer-Extension Flow

codex-2026-06-26-S215. Prompt: use the tiling model for the growth
`n -> n+1`, where a new diagonal layer is added and consecutive binary layers
of sizes `k` and `k+1` are connected by `k^2+k` complete-bipartite lines.

## Main Synthesis

The useful recursive object is not an unrooted tournament class. It is a
rooted extension state with a retained layer address.

In the fixed-Hamiltonian-path staircase, growing by one vertex adds a new
diagonal strip. If two consecutive strips have binary boundary words

```text
x in {0,1}^k,
y in {0,1}^{k+1},
```

then the complete bipartite line sheet between them should not be treated as
`k(k+1)` independent bits. The natural pair-potential model is

```text
e_ij = x_i XOR y_j.
```

That immediately forces the rectangle law

```text
e_ij + e_i'j + e_ij' + e_i'j' = 0 mod 2.
```

So the cross-sheet is a coboundary. Its `k^2+k` apparent edge values have only
`2k` degrees of freedom. The saved exponent is `k(k-1)`. This is the same
shape as the two-dimensional Haar product rule, rank-one sign matrices, and
Schur/observability sidecars from the matrix atlas: store the boundary
potentials, and treat rectangle defects as the meaningful residual.

## Extension States

Adding one tournament vertex gives the exact class-level version of the same
idea. For a parent class `[T]`, the new vertex is described by an incident
word on `V(T)`, but words related by `Aut(T)` are the same rooted extension
state. Hence

```text
E(n -> n+1)
  = sum_[T in A(n)] (1/|Aut(T)|) sum_{g in Aut(T)} 2^{cycles(g)}
  = R(n+1),
```

where `R(m)` is the number of rooted/perspective tournament classes on `m`
vertices.

This cleanly explains the old perspective observation:

```text
R(3)=4=A(4)
R(4)=12=A(5)
```

and also puts a warning sign exactly where the pattern stops:

```text
R(5)=48, A(6)=56
R(6)=296, A(7)=456
```

So "perspectives at one level count structures at the next" is true only for
the first few levels. The real recursive law is `E(n -> n+1)=R(n+1)`, followed
by a nontrivial unrooting quotient.

## Algorithmic Efficiency

The practical enumeration recipe is:

1. Keep parent class representatives and automorphism cycle indices.
2. Generate incident-word orbits under `Aut(parent)`.
3. Extend by one rooted vertex.
4. Cache the rooted child state before unrooting.
5. Merge rooted states to unrooted children and record the collision fiber.

This is the same lesson as the LRC sidecar work: do not flatten the carrier
before the proof-relevant coordinate has either been retained, reconstructed,
annihilated, descended, or routed to named residual debt.

For the tiling sheet, the analogous recipe is:

1. Store row and column boundary words.
2. Generate the rank-one/coboundary cross-sheet.
3. Track rectangle defects as residuals.
4. Only then quotient by reflection, complement, or isomorphism.

## Half-Tiling Reading

The proposed "unique binary half-tilings of a smaller size" should be sharpened
to:

> unrooted tournament classes are quotients of rooted/perspective half-tiling
> extension states.

Symmetry along diagonals is not just a visual feature. It is a group action:
parent automorphisms, root motion, and complement/reflection can all act on the
same layer word. The right count is therefore a Burnside count of fixed
boundary words and fixed rank-one sheets, not raw binary half-tilings.

## LRC Angle

For LRC packet tournaments, this suggests a new sidecar family:

```text
parent_class,
root_orbit,
incident_word_orbit,
layer_boundary_word,
rank_one_sheet_id,
rectangle_defect_rank,
unrooting_collision_fiber.
```

The preserved predicate is not immediately loneliness. It is controlled
extension identity: we know which local address was responsible for producing
a class or packet before a quotient erased it. The destroyed information is
the root/deletion address and the specific boundary potentials behind a
cross-layer interaction.

This is especially attractive for residual packets that currently look equal
after scalarization. A scalar may say two rows are the same; the rooted
extension carrier may say they came from different deletion addresses, just as
HYP-3041's AP-tail row needed the hidden `m mod 13` clock.

## Concrete Result

Script: `04-computation/tournament_layer_extension_flow_codex_s215.py`.

Stored output: `05-knowledge/results/tournament_layer_extension_flow_codex_s215.out`.

The scout verifies the layer rectangle law through `k=6`, computes unrooted
and rooted/perspective counts through `n=7`, and checks
`E(n -> n+1)=R(n+1)` through `n=6 -> 7`.

## Next Theorem Target

Prove a **rooted extension quotient theorem**:

> Any tournament class enumeration or LRC tournament quotient built by adding a
> layer may forget the new root only after the incident-word orbit and
> rectangle-defect sidecar are either retained, reconstructed from the child,
> annihilated by a symmetry certificate, or routed to a named residual fiber.

That would turn the tiling model from a picture into an enumeration and
proof-safety rule.
