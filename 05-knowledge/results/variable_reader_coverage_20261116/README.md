# Variable tournament readers beyond sparse rays

- [Binary address, exact offsets, and absolute coverage](PROOF.md).
  Three strong rooted readers synchronize on the same continuation
  data. An arbitrary-length binary word then makes their common
  leading coefficient run through a full arithmetic progression
  while retaining every +2 step. A collision argument proves
  Ω(X^(2/3)) disjoint exact windows below X. A different, explicit
  decoder brackets every sufficiently large prescribed count
  within O(√X). Neither statement is just a residue census.
- [Five-count bank](FIVE.md). A second small seed bank increases
  each window to five consecutive odd counts, without restricting
  its binary address. All the same absolute lower bounds apply.
  The witnesses for a fixed address keep every arc incident to the
  growing appended portion.
- [Whole-block transfer](BLOCKS.md). Attaching a common block
  multiplies the low intercept by the block's Hamilton count
  and the high coefficient by its two-slot path-cover count.
  This exactly explains why common nontransitive blocks do not
  preserve the unit low step.

The width is still bounded, and neither positive density nor a
witness at an arbitrary exact target follows. The next concrete
problem is to match the selector's remaining signed error to
a compatible reader, or to prove stronger collision and local
gap bounds at the scale of single integers. Chan–Pak remains open.

Reproduce with Python 3 verify.py. For independent seed discovery,
compile seed_search.cpp and pass the order-eight CSV from the
[fixed-count census](../actual_hamilton_multiplicities_20261021/session24-multiplicities/).
