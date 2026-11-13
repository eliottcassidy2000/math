# Arc-use shadows, two and three visits, and signed return gates

2026-11-13. **Chan–Pak and Collatz remain open.** The main construction
changes one actual arc (after an exterior relabelling), leaves the total
Hamiltonian count fixed, and changes both the visible support and the
higher path-cover response. The construction works for arbitrary common
internal tournaments, not just growing transitive chains.

## Main results

- [Reservoir polynomial](RESERVOIR_POLYNOMIAL.md): replace the common
  marked module by any nonempty tournament U. The two strong outputs
  have the same H = 9(H(U)+Q_U(2)) and the same Hamilton cycle count
  H(U), with identical internal and entry arc weights. For |U|>=2
  the switched side has strictly more ordered j-path covers for
  every 2<=j<=|U|+2, with no change at j=1. A common strong homogeneous
  continuation therefore separates the Hamilton counts; a later
  common homogeneous continuation cannot restore equality.
- [Dark arcs](RESERVOIR_DARK_ARCS.md): with |U|>=2 the two outputs have
  exactly one and three zero-use arcs, respectively, regardless of
  U's own support. Every internal arc of U is now used, even if it
  was dark in U alone. Thus the equal H values do not have equal
  weighted shadows. The underlying dark-arc classification and
  reconstruction from **directed** positive support are inherited;
  [arc gates](ARC_GATES.md) adds a direct quantitative cycle-cut test
  and an explicit repeatable gate on a marked arc's endpoint response.
- [Arc-use shadow](ARC_USAGE_SHADOW.md): the exact ordered-SCC law
  is internal scaling, rank-one adjacent seam, and zero distant seam.
  At order four, the three positive-weight multisets in the request
  are exactly the possibilities; the middle one belongs to the
  converse pair and the first/last to a single class each. The
  transitive class has three additional zeros if all six arcs are
  listed. At order 22, two tournaments that are neither isomorphic
  nor converses share the entire labelled **undirected** weighted
  matrix and the separate labelled start/end vectors. The anchored
  block construction gives 2^k converse-merged classes with one
  such shadow at order 21k+1. There is no claim this is the first
  order of a collision or that it already occurs in a strong class.
- [Two and three](TWO_AND_THREE.md): the 2^m and 3^m are surjection
  counts of transitive-chain vertices to two versus three visits.
  For the dyadic pair, the two-slot cover gap is exactly 4(3^m−3).
  The flipped arc itself carries 2/9 of the paths on each side.
  In a general bounded-visit transitive inflation, fixed-distance
  internal edge probabilities converge to the explicit geometric
  kernel for the maximal visit count; this local law does not
  determine the exterior response.
- [Signed clocks](SIGNED_CLOCK_RESERVOIR.md): the known three negative
  cycles and one positive cycle are separated by actual affine carry
  gates, not by a tournament-class bijection. Equal-clock blocks give
  an exact weighted hull and finite native inverse. Neither the two-
  phase nor seven-phase negative cycle allows a mixed return using
  its phase rotations as blocks. Positive shadows of the two-phase
  negative cycle keep gap two indefinitely, but their endpoints
  cannot be made equal by integer starts at a common depth.

The recurring principle is to distinguish **a preserved scalar, the
marked interface, and the physical moves that remain legal**. Counts
alone do not say which arcs stay available. Binary growth of one
observable can conceal a ternary contribution in another. Neither a
dense formal range nor homogeneous substitution supplies a native
exact return. The next Chan–Pak route needs a selective contact or
asymmetric continuation that uses the new response without imposing
this family's permanent factor of 9.

Replays: replay_shadow.py, cayley_certificate.py, replay_gates.py,
replay_layouts.py, replay_homogeneous.py, replay_occupancy.py,
replay_signed_clock.py. The layout and gate replays share arcscan.py
for direct subset path counts; the Cayley certificate checks only
the local fingerprints and uses the all-size SCC product proof.
