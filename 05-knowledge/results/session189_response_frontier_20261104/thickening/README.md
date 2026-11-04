# Every source admits a separated transitive thickening

For any nonempty tournament T, compute the Seidel exponent L_Q for every
prime quotient Q in its original strong-module tree. Let L be their maximum,
or1 if none exist. For any positive N divisible by2^L, replace every vertex
of T by TT_N.

The resulting source is separated from every distinct XOR target formed
using modules of that inflated source. New modules inside the cells are
included. The proof combines the geometric-profile fixed point, the all-order
ambient Seidel separator, preservation of prime quotient types under
transitive inflation, and the independently established highest-node
restoration argument.

This is source-relative and concerns the new inflated source. It does not
settle the uninflated source or all-pairs response injectivity.

`PROOF.md` specifies the finite constructive threshold and exact dependencies.
`verify.py` independently enumerates modules and strong nodes for155 literal
inflations through18 vertices, retaining every structural control. Those
small structural examples are not substituted for the larger dyadic profile
hypothesis. Normal and optimized receipts must match.
