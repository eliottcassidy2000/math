# Curvature and triangle-inequality capacity control

**PROVED / FINITE-EXACT / CITED standard metric inputs.**

The triangle inequality puts the images of r-step predecessors inside a
radius-Lr ball when every edge has length at most L. Separation or a bounded
fibre is a separate assumption. In curvature -kappa^2, separated-position
packing gives

    (N/K)(cosh(kappa delta/2)-1)
        <= cosh(kappa Lr+kappa delta/2)-1.

For a b-ary tree this forces kappa L>=log b. Negative curvature increases
ball capacity; it does not reverse or weaken the triangle inequality.

## Results

- **Exact Collatz replay:** 704 first-hit sources at radius 21 versus 694
  honeycomb positions, excess **10**; 2263 versus 2107 at radius 26 for the
  triangular lattice, excess 156. Full source/clock TSV is pinned and replayed.
- **Native unions:** the selected SF5 depth-8 tree has 511 odd nodes and 2001
  actual native states, not 2811 independent-subdivision states. Duplicate
  even tails account for 810. First-hit root stopping is stated precisely.
- **Dyadic tree:** exact outer radius
  acosh((5*2^k-2+5*2^(-k))/8)=k log2+log(5/4)+o(1). It is radially optimal
  but not quasi-isometric: some graph distances 2k have constant ambient
  distance acosh(3/2).
- **Gapped repair:** a two-parameter equal-geodesic-edge drawing has explicit
  coarse upper/lower metric bounds. It is proper, noncrossing, and uniformly
  separated for each parameter choice. It proves the binary edge threshold
  log2 is sharp as an infimum; the equality case is not settled. Small-width
  versions need not have minimum-distance pairs equal to the tree edges.
- **A distinct 20:** degree-5 spherical triangulations have 20 triangular
  faces by the Euler curvature sum. This is not the Collatz capacity gap.

See [PROOF.md](PROOF.md) for all quantifiers, proofs, failure boundaries, and
primary-source links. No Collatz convergence theorem, native hyperbolic
Collatz embedding, or regular hyperbolic hexagonal mesh is claimed.

## Replay

From any working directory:

```sh
python3 /path/to/curvature_triangle/verify.py
python3 -O /path/to/curvature_triangle/verify.py
```

The stdlib-only checker has no repository or network dependency. It verifies
five pinned source snapshots, exact predecessor enumeration and all 2263
source paths, lattice BFS through radius 105, native union counts through
odd depth 8, five-step supply, 56,388 rational vertex-pair checks, edge-circle
and noncrossing certificates, and exact dyadic radial controls. It raises
explicit exceptions, not optimization-sensitive assertions. Both saved
outputs are identical. Floating metric examples are illustrative only;
all finite algebraic geometry certificates use rational arithmetic.

`PUBLIC_FILES.json` selects the payload; `SHA256SUMS` seals the payload and
selector. `PORTABILITY.json` records an isolated-copy, different-cwd replay.
