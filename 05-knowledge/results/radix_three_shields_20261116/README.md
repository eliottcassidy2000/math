# Width-four windows, factor-three shields, and physical radix tests

- [Two independent dials](TWO_DIALS.md): 24 vertices,
  H=6,962,193+2i+6j for i,j=0,1,2,3. The sixteen choices
  cover thirteen consecutive odd counts. A second bank has
  an actual boundary weight 68, with an order-17 smaller
  test and a full order-28 rectangle. The factor three
  is an effective response ratio, not the standalone count
  of a second disjoint active module.
- [A literal count-three shield with the +2 window intact](COMPENSATION.md).
  All four counts 299,301,303,305 live next to the same
  rooted triangle, with identical shield incidences.
  Growing shields of rooted counts 9 and 17 work too.
  A full signature proof classifies every fixed shield
  in this architecture for these four seeds.
- [A cyclic triangle amplifier](TRIANGLE_WINDOW.md) scales
  a literal five-count window into a step-six window on
  nine vertices. [Obstructions](OBSTRUCTIONS.md) show that
  **no** width-four palette can scale by three twice with
  this same cyclic wrapper. A different uniform rooted
  amplifier does give exact 3^k, but erases the low digit
  modulo three from its continuation.
- [Strong reverse bridge and exact continuation](BRIDGE.md):
  a three-setting digit survives at any height with step
  2·3^k. The weighted endpoint factor is distinguished
  from an ordinary count-three block. [The factor-one
  bound](BRIDGE_BOUNDS.md) and [a repeatable independent
  deposit](DEPOSIT.md) identify a side interface: one may
  deposit a new bit while scaling all previous variations
  by 3, but the new deposit grows on a base-16 scale
  relative to the old digit, not the required base 3.
- [Symmetric and adjacent product barriers](NO_GO.md)
  and a [small-order independent-arc test](CUBE.md). These
  reject specific stacking architectures, not the two-dial
  construction or arbitrary coordinated edits.
- [17, 68, Fibonacci, and the primitive triple tree](TREE_AND_SEEDS.md).
  The redundant ternary representation count has two
  actual Stern shears, also used in a two-coordinate
  Berggren chart with a decreasing inverse. The arithmetic
  seed reduction is exactly odd primes, 7p, and 49,63,343.

What remains open: a physically reusable **bounded** deposit
with common continuation, so that an existing lower digit and
the next higher digit remain freely selectable. Exact
amplification, literal neighboring shields, and a finite
two-dial rectangle no longer need to be assumed. None supplies
the entire unbounded radix-three tower or proves Chan–Pak.
