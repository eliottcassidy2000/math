# Chan–Pak: choose an exact repair and place its baseline from the target

**Main result.** For every sufficiently large prescribed target N, a
deterministic construction produces one three-count window below N and
one above it. Each window consists of three consecutive odds, witnessed
by strong tournaments on O(log N) vertices. All six counts are within
**O(N^(3/8))** of N. This replaces the earlier coarse bound of about
N^0.8963 for this approach and, crucially, the +2 choices are retained
on the very baselines computed from N. It does not prove that a window
contains N.

## Choose the baseline and the repair together

- [Two-sided target bracketing](BRACKETS.md): hold the low bits of
  both addresses equal but charge their exact checksum contribution
  to the baseline. A common midpoint can then be chosen on the
  integer grid before opening a coarser search grid around it. A
  single integer square-root-like inversion picks each pair of
  addresses. About the target, the error from the high digits and the
  quadratic grid error balance at N^(3/8). Both windows
  contain all mod-3 classes.
- [Complete inverse for the three-count bank](INVERSE.md): a literal
  six-vertex center with three settings preserves the full interface
  to independently chosen transitive exteriors. For an exact request
  the residue mod 3 chooses its setting first; low-bit recursion then
  searches the two addresses in O(sqrt(N) log² N) arithmetic steps.
  For 35,099 it returns addresses 5 and 14 and a strong 14-vertex
  witness. It rejects 36,833 for this bank, even though the earlier
  marked-host construction reaches that target. A complete inverse
  to one construction is not a proof of existence.
- [Exact inverse for all one-apex sources](bounded/EXACT_INVERSE.md):
  the product estimate leaves only O(1) candidate opposite addresses
  per address tried. This reduces the exact test to O(sqrt(N) log N)
  arithmetic steps over the entire one-apex family. The literal
  contact word for 60,001 is `10000000001101000`. The same
  low-phase placement idea also gives an
  [O(N^(3/8)) bare one-apex approximation](bounded/TARGET_PLACEMENT.md).
- [Why not search only near squares?](IMBALANCE.md) For just the
  six-vertex center bank, bounding the aspect ratio by a fixed κ
  gives at most (log κ)/72 + o(1) of odd targets, even in each fixed
  fractional band and before any collisions are counted.

## Other exact repairs and obstructions

- [A four-count mirror window](windows/PORT_WINDOWS.md) starts at
  1279, 1281, 1283, 1285. It uses a nontransitive order-seven
  center, with matched exterior addresses. The same work shows
  that swapping two ordinal blocks cannot change just the apex
  count while preserving their entire two-sided interface.
- [Switch the root rather than just add a seam](seams/SINGLE_SEAM.md).
  A nontrivial isolated seam reversal costs at least half the
  baseline. But choosing between two backward arcs instead can
  give exactly `18Q−H` and `18Q+H` using the same module;
  a transitive module gives consecutive odds. There is a separate
  three-power scale and an exact factor-based inverse for a
  two-operand construction with singleton padding.
- [Two independent scales with a unit switch](seams/MIXED_SCALE.md).
  A strong six-vertex template with two inflated vertices produces
  exact adjacent odd counts with leading scale 18·3^m·4^n.
  Irrational log rotation moves this pair into any fixed relative
  interval around a sufficiently large target. The bank is sparse,
  and has a separate exact inverse; relative coverage is not exact
  coverage.
- [Small marked correction is not free](bounded/BOUNDED_REPAIRS.md).
  Every pair of marks contributes at least two cover choices.
  With bounded correction and bounded marked boundary components,
  a fixed source gives only a finite union of exact dyadic rays.
  The quotient test reconstructs any hit without a host-order scan.

## Next proof obligation

Use the same computed interfaces to absorb a general prescribed residual
of size O(N^(3/8)), or replace the two separate steps with an exact
selector that is guaranteed to find an address. Neither a fixed
three- or four-count window nor a separate approximation suffices.
The included graph-count replays, all-target inverse comparisons,
and bounded-state checks test the formulas; the uniform claims are
proved separately. Chan–Pak remains open.
