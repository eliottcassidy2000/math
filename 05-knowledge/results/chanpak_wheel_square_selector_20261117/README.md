# First-25 prime wheel, square-sum paths through 36, and a free +2 seam

**New constructive result.** A literal marked core balances its unsigned
left/right reads and unbalances its signed reads by two. Reflecting the
two address roles therefore gives an exact +2 pair while both addresses
remain independent. Any prescribed odd residue can be the lower member,
even with frozen binary fringe ends. From an arbitrary large absolute
target N one can select both addresses and place such a pair within
O(N^(2/5)). This is not exact universal placement. See
[the free seam and selector](FREE_SEAM.md).

## What the two archives contributed

- **The first 25 primes:** the
  [first-25 exponent wheel](../first25_spread_target_returns_20261116/QUINTILES.md)
  for 99·2^k+23 is an exact periodic divisibility observer. Two of
  its filters are redundant, and a composite outside the wheel's
  prime support can pass all of them. The [99-slope rails](ray/RAY.md)
  now realize the tested sequence itself at *every* height, rather
  than inferring its physical existence from the wheel. Nine genuine
  cores share the slope; choosing the core places one of nine exact
  offsets (including a mod-7 selector and a +2 pair). These nine rails
  simultaneously fail the wheel at a common periodic phase. Any
  finite menu of such prime candidates with a common odd anchor
  has analogous repeatable simultaneous composite phases, possibly
  using larger primes.
- **The square-sum graphs:** the
  [arrival and chain-kernel archive](../session32_square_monster_20261023/square-sum/README.md)
  shows why merely knowing that a path exists at each order loses the
  endpoint and parity information used by
  [Gerbicz's all-order construction](https://oeis.org/A090461).
  The [protected-bridge analysis](square/PROOF.md) makes the loss
  measurable: orders 30–36 share exactly one parity state, forcing
  edge 19–30, while adding order 29 destroys the shared state. A
  different path at 35 admits the all-depth 25-fold lift carrying
  that edge to its scaled descendants. This is a marked transport,
  not the unchanged pair of numerical labels at arbitrary order.

The common design principle is to **put the exact move in a coordinate
that growth does not consume**. The square construction does this with
a parity-compatible marked bridge. In the tournament one obtains it by
cancelling the two unsigned reads and preserving the signed difference
under reflection. The [congruence theorem](CONGRUENCE.md) proves which
fringe constraints really survive a general modulus; it even permits
two independently interleavable blocks at every depth. The
[singleton control](REFLECTION.md) goes further locally: any fixed even
gap with any target congruence, yet only O(sqrt(X)) possible lower
values below X. This is an explicit counterexample to inferring
absolute coverage from complete local coverage for exact repairs.

For the next constructive step, seek a second physical change that
commutes with the free reflection and covers its remaining absolute
residual. Testing additional wheel phases would not address this
step: that observer is already saturated on the free seam.

## Reproduction

From this directory:

    python3 check.py
    python3 check_reflection.py
    python3 check_free.py
    python3 ray/replay.py
    c++ -O3 -std=c++20 square/replay.cpp -o /tmp/chanpak-square36
    /tmp/chanpak-square36
    python3 square/lift_edge.py

Python replays agree under normal and optimized execution. The literal
tournament checks count Hamiltonian paths by subset DP from the actual
arcs. Large approximate-selector controls use exact integer/rational
arithmetic and the proved interface law, not large path enumeration.
