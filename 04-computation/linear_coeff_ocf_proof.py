"""
linear_coeff_ocf_proof.py — opus-2026-04-04-S2

PROVE: c_k = 2^(skip-1) for tile k = (x,y) with skip s = x-y.

From OCF: c_k = 2 * Σ_C [coeff of t_k in χ_C(t)]

The coefficient of t_k in χ_C is:
  +1 if C uses tile k in backward direction AND all other arcs of C are base-path arcs
  -1 if C uses tile k in forward direction AND all other arcs of C are base-path arcs
  (More complex if C has multiple tile arcs)

Actually, the linear coefficient in χ_C counts cycles that use tile k and have
all OTHER arcs evaluated at t=0 (base path), giving factor 1 from (1-0) for forward
tile arcs and 0 for backward tile arcs.

WAIT — at t=0, the indicator is χ_C(0) = 0 (no cycles in transitive tournament).
The linear coefficient of t_k in χ_C is: ∂χ_C/∂t_k evaluated at t=0.

For χ_C = Π factors:
  ∂χ_C/∂t_k |_{t=0} = (∂/∂t_k of the tile-k factor) × Π_{other factors at t=0}

The tile-k factor is either:
  (1-t_k): derivative = -1, evaluated at 0 gives -1 × Π{others at 0}
  t_k: derivative = +1, evaluated at 0 gives +1 × Π{others at 0}

The other factors at t=0:
  (1-0) = 1 for forward tile arcs (correct: arc present at t=0)
  0 for backward tile arcs (arc NOT present at t=0)
  1 for base-path arcs (always present)

So: ∂χ_C/∂t_k |_{t=0} is nonzero ONLY if all OTHER tile arcs are forward.
The cycle must use tile k (in either direction) and ALL other non-base-path arcs
in the forward direction.

For a 3-cycle (a,b,c) with a→b→c→a:
If it uses tile k forward, the derivative contribution is -1 × (product of other factors at 0).
If it uses tile k backward, the derivative contribution is +1 × (product of other factors at 0).
Other factors at 0 are nonzero only if all other tile arcs are forward (value 1).

So: the linear coefficient of t_k in χ_C is:
  +1 if tile k is backward AND all other arcs of C are either base-path or forward tile
  -1 if tile k is forward AND all other arcs of C are either base-path or forward tile

Now let's count how many such cycles exist for tile k = (x,y) with skip s = x-y.

CASE 1: Tile k is BACKWARD in C. The arc is (y-1) → (x-1) (lower → higher).
All other arcs must be base-path or forward tile (going from higher to lower).
The cycle visits some set of vertices including x-1 and y-1.

A cycle that goes from y-1 up to x-1 (backward) and then descends back to y-1
using only forward (high→low) arcs and base-path (k+1→k) arcs.

The descent from x-1 to y-1: at each step, we go from some vertex v to a lower
vertex w, either via base path (v=w+1) or forward tile (v>w+1). The descent must
pass through exactly the vertices between x-1 and y-1.

Wait, a 3-cycle has only 3 vertices. If the backward arc is y-1 → x-1, the third
vertex z must be between y-1 and x-1 (to allow descending arcs from x-1 to z to y-1).

For the 3-cycle (y-1, x-1, z) with z between y-1 and x-1:
  Arc y-1 → x-1: backward tile (CONTRIBUTES t_k)
  Arc x-1 → z: forward if x-1 > z+1 (tile arc) or base path if z = x-2
  Arc z → y-1: forward if z > y-1+1 (tile arc) or base path if z = y

For arcs x-1→z and z→y-1 to both be forward-or-base-path, we need z ∈ {y, y+1, ..., x-2}.
That's x-2-y+1 = s-2 choices (where s = x-y = skip).

Wait, but for longer cycles (5-cycles etc.), there are additional vertices between
x-1 and y-1, and the descent can take multiple steps.

CASE 2: Tile k is FORWARD in C. The arc is (x-1) → (y-1) (higher → lower).
All other arcs must be base-path or forward tile (also higher → lower).
But a cycle must go both up AND down — so if tile k provides the only non-descending
arc, the cycle is impossible (all arcs descend, can't form a cycle).

Wait, forward means going from higher to lower (x-1 > y-1). If ALL arcs in the cycle
go from higher to lower, the cycle can't close (would need to go up at some point).
So the linear coefficient from the forward case is 0!

NO — I think I confused the direction. Let me re-examine.

In the base-path model:
- Base-path arcs: k+1 → k (descending)
- Forward tile arcs: x-1 → y-1 where x > y+1 (also descending, since x-1 > y-1)
- Backward tile arcs: y-1 → x-1 (ASCENDING)

So at t=0 (all tiles forward), ALL arcs descend. The transitive tournament.
A cycle needs at least one ascending arc. So at t=0, no cycles exist (Redei). ✓

The linear coefficient of t_k comes from cycles that use tile k as their ONLY
ascending arc (backward direction), with all other arcs descending.

If tile k is used in FORWARD direction (descending), the derivative is -1 but the
other factors at t=0 require the other arcs to all be descending too. A cycle with
ALL arcs descending except one that's also descending (the forward tile) = still all
descending. No cycle possible. Contribution = 0 from forward case.

So: c_k = 2 × (number of directed odd cycles using tile k backward, with all other
arcs descending)
"""
import numpy as np
from collections import defaultdict

def count_descent_cycles(n, tile_x, tile_y):
    """
    Count directed odd cycles that:
    - Use the backward arc (tile_y-1) → (tile_x-1) (ascending)
    - All other arcs descend (base-path or forward tile)

    These are paths that go from tile_x-1 DOWN to tile_y-1 via descending arcs,
    then jump UP via the backward tile arc.

    A descending path from a to b (a > b) visits vertices a, z_1, z_2, ..., b
    where each step goes to a strictly lower vertex. The steps can be:
    - Base-path: v → v-1 (step size 1)
    - Forward tile: v → w where v > w+1 (step size ≥ 2)

    The number of descending paths from a to b equals the number of compositions
    of (a - b) into positive integers (each ≥ 1).

    Wait — the path visits intermediate vertices. From a to b, we need
    a > z_1 > z_2 > ... > b. Each step is a descending arc (base-path or tile).

    The number of descending paths from a to b (not passing through any specific
    set of vertices) is 2^(a-b-1).

    This is because at each intermediate vertex v (a > v > b), we either:
    - Pass through v (using arcs ... → v → ...)
    - Skip v (using a longer tile arc jumping over v)

    There are a-b-1 intermediate vertices, each independently included or excluded.
    So 2^(a-b-1) paths.

    But wait — we need an ODD-length cycle. The cycle has:
    - One backward arc: y-1 → x-1 (length 1)
    - A descending path from x-1 to y-1 of some length L
    - Total cycle length = 1 + L

    For the cycle to be odd-length, L must be even.

    Hmm actually, in our directed cycle definition, the cycle length equals the number
    of arcs, which equals the number of vertices (since it's a cycle).

    The cycle visits vertices x-1, z_1, ..., z_{L-1}, y-1 and then back to x-1.
    Total vertices = L + 1 (from x-1 to y-1, L steps + the backward arc gives L+1 vertices
    in the cycle... wait, the cycle visits x-1, then descends to y-1 through L-1 intermediate
    vertices, then jumps y-1 → x-1. So vertices are: x-1, z_1, ..., z_{L-2}, y-1.
    That's L vertices total. Arcs: L (descent) + 1 (backward) = L+1.

    No wait, a cycle on V vertices has V arcs. If V vertices are visited, V arcs are used.

    Let me recount: the cycle is (y-1, x-1, z_1, z_2, ..., z_{k-1}) for some vertices
    z_1 > z_2 > ... > z_{k-1} all between x-1 and y-1. Then:
    - Arc y-1 → x-1: backward tile (ascending)
    - Arc x-1 → z_1: descending
    - Arc z_1 → z_2: descending
    - ...
    - Arc z_{k-1} → y-1: descending

    Total arcs = 1 + k = cycle length. Total vertices = 2 + k.
    Cycle length = 2 + k must be odd, so k must be odd.

    The vertices z_1,...,z_k in the range (y-1, x-1) = {y, y+1, ..., x-2}.
    There are s-2 such vertices (s = x-y).

    For each SUBSET of these s-2 vertices, we get a valid descending path
    (by visiting them in decreasing order). The subset must have ODD size k
    for the cycle to have odd length.

    Number of odd-size subsets of {y, y+1, ..., x-2} (s-2 elements):
    = Σ_{k odd} C(s-2, k) = 2^(s-3) (if s ≥ 3)

    For s = 2: no intermediate vertices, k must be odd, k=0 is even. BUT k=0 gives
    cycle (y-1, x-1) which is a 2-cycle (even). No valid cycles? But we know c_1 = 2
    for skip-2 tiles...

    Wait, k=0 means 0 intermediate vertices: the cycle is (y-1, x-1) → y-1, a 2-cycle.
    That's even-length, not odd. So for skip-2 tiles, we need k=1 (one intermediate
    vertex), giving a 3-cycle. The intermediate vertex is... there are s-2 = 0 vertices
    between y-1 and x-1 when s=2 (since y-1 and x-1 are separated by 1).

    Hmm wait, x-y = 2 means x-1 = y+1. So vertices y-1 and x-1 = y are adjacent!
    The backward arc y-1 → y is... y-1 → (y-1+1) = y-1 → y = reverse base path.
    That's NOT a tile arc! Tile (x,y) = (y+2, y) has x-y=2, and the arc connects
    vertices y+1 and y-1 (0-indexed: x-1 = y+1, y-1 = y-1). So it's between vertices
    y+1 and y-1, skip = (y+1)-(y-1) = 2. These are NOT consecutive!

    OK I got confused with indexing. Let me be very explicit.

    Tile (x, y) with x-y >= 2 (1-indexed). In 0-indexed: vertices x-1 and y-1.
    x-1 - (y-1) = x-y = s >= 2. So x-1 > y-1 + 1 (non-consecutive).

    Forward: x-1 → y-1 (descending).
    Backward: y-1 → x-1 (ascending).

    Intermediate vertices between y-1 and x-1: {y, y+1, ..., x-2}.
    Count = (x-2) - y + 1 = x - y - 1 = s - 1.

    Wait, that's s-1, not s-2! Let me recount.
    y-1 < y < y+1 < ... < x-2 < x-1.
    Between y-1 and x-1 (exclusive): {y, y+1, ..., x-2}. Count = (x-2) - y + 1 = s - 2.

    Hmm: s = x - y. x - 2 - y + 1 = x - y - 1 = s - 1. So there are s-1 intermediate
    vertices?

    No: y-1 = y-1, the NEXT integer is y. The vertex BEFORE x-1 is x-2.
    So between y-1 and x-1 exclusive: {y-1+1, ..., x-1-1} = {y, ..., x-2}.
    Count = (x-2) - y + 1 = x - y - 1 = s - 1.

    So s - 1 intermediate vertices, not s - 2!

    For skip s=2: s-1 = 1 intermediate vertex. Subsets of size k (odd): k=1 gives C(1,1)=1.
    Total odd subsets = 1. Cycle count = 1. c_k = 2·1 = 2 = 2^(2-1). ✓

    For skip s=3: s-1 = 2 intermediate vertices. Odd subsets: k=1 gives C(2,1)=2.
    Total = 2. c_k = 2·2 = 4 = 2^(3-1). ✓

    For skip s=4: s-1 = 3 intermediate vertices. Odd subsets: k=1 gives C(3,1)=3,
    k=3 gives C(3,3)=1. Total = 4 = 2^2. c_k = 2·4 = 8 = 2^(4-1). ✓

    GENERAL: number of odd-size subsets of (s-1) elements = 2^(s-2).
    c_k = 2 · 2^(s-2) = 2^(s-1). ✓✓✓

    THIS IS THE PROOF! The number of cycles is the number of odd-size subsets of
    the s-1 intermediate vertices, which is 2^(s-2) by the standard identity.
    """
    s = tile_x - tile_y
    x0 = tile_x - 1  # 0-indexed upper vertex
    y0 = tile_y - 1  # 0-indexed lower vertex

    # Intermediate vertices: {y0+1, ..., x0-1}
    intermediates = list(range(y0+1, x0))
    num_inter = len(intermediates)

    # Count odd-size subsets
    total_odd = 0
    for k in range(1, num_inter + 1, 2):
        from math import comb
        total_odd += comb(num_inter, k)

    # Expected: 2^(s-2) = 2^(num_inter - 1) for num_inter >= 1
    expected = 2**(num_inter - 1) if num_inter >= 1 else 0

    return total_odd, expected, num_inter

def verify_formula():
    print("PROOF OF LINEAR COEFFICIENT FORMULA c_k = 2^(skip-1)")
    print("="*60)
    print()
    print("Tile (x,y) with skip s = x-y connects vertices x-1 and y-1 (0-indexed).")
    print("Between them: s-1 intermediate vertices.")
    print()
    print("A valid cycle: backward tile arc (y-1 → x-1) + descent through a SUBSET")
    print("of intermediates (all arcs descending). Cycle length = 2 + |subset|.")
    print("For ODD length: |subset| must be ODD.")
    print()
    print("Count of odd-size subsets of (s-1) elements = 2^(s-2).")
    print("OCF weight: 2^1 = 2 per cycle.")
    print("Total: c_k = 2 · 2^(s-2) = 2^(s-1). QED")
    print()

    # Verify for various s
    for s in range(2, 9):
        count, expected, num_inter = count_descent_cycles(10, s+1, 1)
        c_k = 2 * count
        target = 2**(s-1)
        status = "✓" if c_k == target else f"✗ (got {c_k})"
        print(f"  s={s}: intermediates={num_inter}, odd subsets={count}, "
              f"c_k={c_k}, 2^(s-1)={target} {status}")

    print()
    print("The key identity: #{odd-size subsets of n elements} = 2^(n-1) for n >= 1.")
    print("Proof: the map S → S Δ {1} (symmetric difference with {1}) bijects")
    print("  odd-size subsets ↔ even-size subsets. Total subsets = 2^n.")
    print("  So odd count = even count = 2^n / 2 = 2^(n-1).")

if __name__ == "__main__":
    verify_formula()
