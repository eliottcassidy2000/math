# LRC14 Triangular Affine-Operator Carrier

The equality from the prompt is best treated as a failed scalar identity and a
good structural warning.  It is equivalent to `x^2=1/2+log2(x)`, but the gap
`x^2-1/2-log2(x)` has positive minimum about `0.456964`, so there is no positive
real solution.  That matters: using this as an equality would smuggle a false
scalar into a proof route that is already sensitive to scalarization.

The useful part is the carrier.  `x*(2x-1)` is `T_{2x-1}` and becomes an even
perfect number when `x=2^(r-1)` and `2^r-1` is Mersenne prime.  This is exactly
the product lane: for LRC14, `p*q=p*(14p-1)` is the apex-14 version of the same
quadratic family.  It belongs beside `p*q`, `K_{p,q}`, and the product/coimage
ledger, not beside the binding denominator `q`.

The affine maps `a(x)=x/2` and `b(x)=x+1` give the cleaner proof metaphor.
Each `b` increment carries the number of future halvings as a depth label.
The staircase word `(ba)^n` has depths `n,n-1,...,1` and depth sum `T_n`;
the block word `b^n a^n` has depth sum `n^2`; the tail word `a^n b^n` has
depth sum `0`.  The same counts produce different carriers depending on order.
That is the same warning as HYP-2940's middle SCC: additive, product, and graph
packet labels cannot be collapsed to one scalar without losing the proof
predicate.

For LRC14, the proposed use is narrow.  First attach exact `M=p/q`, the C27
shell-transfer or two-swap splice label, and the Kpq/K33 owner label.  Only then
use the affine-depth packet as an order certificate: unit-visible entries should
route to the C27 petal/two-swap discharge; nonunit entries should route toward
the K33/octahedral/Clebsch state-lift packet and HYP-2908/THM-572.  This gives a
new way to phrase the remaining packet-construction obligation, not a completed
proof.
