# Order-eight attachment dimension and lifted bounds

This result builds on Session 145 `POOLED_PARITY_CARRIER.md` and the independently certified Session 146 order-seven result. All statements concern rational linear mean responses to all fixed source-oblivious attachments. No claim is made for a restricted S/B-only port language, source-dependent attachments, conditioning or shared-source products.

## Theorem

At eight named old vertices, the dimension of the span of fixed attachment mean responses is exactly **13028**. For normalized laws the constant coordinate is known, leaving exactly **13027** varying coordinates.

The pooled-parity carrier supplies the upper bound 13028. For the lower bound, 10795 coordinates have at least one isolated old vertex. Each such coordinate is individually in the response span: select seven old vertices containing its non-isolated support, and use the saturated-vertex lifting and rational interpolation proved in the order-seven packet. Distinct signatures have disjoint forest-monomial support and are independent.

There remain exactly 2233 signatures with no isolated old vertices. Project every response coefficient row onto just these 2233 columns. The accompanying 2233-context certificate gives an invertible 2233-by-2233 projected matrix with rank 2233 and determinant 70 modulo 101, and rank 2233 and determinant 94 modulo 103. A duplicated-row control has rank 2232 modulo 101. It follows that the original responses contain 2233 additional independent directions modulo the span of the 10795 individually observable coordinates. Thus the response dimension is at least 10795+2233=13028.

Projection is used only to prove independence modulo the already observable coordinates. The individual projected rows need not themselves be actual response rows. Exterior coefficient extraction is an exact rational combination of actual fixed attachment evaluations by Boolean multilinear interpolation. Nonzero modular determinant implies nonzero rational determinant because all coefficient denominators are powers of two.

## All-order bounds

Let d(n) be the rational response dimension and c(n) the previously proved pooled-parity carrier dimension. For n ≥ 8, the exact order-eight result lifts to every subset of eight named old vertices by the same two-fresh-leaves saturation and component-count interpolation used in the order-seven packet. Consequently every coordinate supported on at most eight old vertices is individually observable, so

    L(n) = 1 + binom(n,2) + 3 binom(n,3) + 9 binom(n,4)
      + 40 binom(n,5) + 120 binom(n,6) + 546 binom(n,7)
      + 2233 binom(n,8) ≤ d(n) ≤ c(n).

Here 1,0,1,3,9,40,120,546,2233 are the numbers of pooled-parity signatures with exactly k non-isolated old vertices for k=0,...,8, obtained by binomial inversion of c(0),...,c(8). In particular, 56296 ≤ d(9) ≤ 65908; 203581 ≤ d(10) ≤ 347626; 626066 ≤ d(11) ≤ 1902231; and 1685509 ≤ d(12) ≤ 10777900. These counts include the constant; subtract one for normalized laws. Equality at n ≥ 9 remains open.
