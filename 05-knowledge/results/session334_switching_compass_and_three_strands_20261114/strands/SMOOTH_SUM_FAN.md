# The internal additive fan of the smooth rail

**Proved.** Write M_i=2^(a_i)*3^(b_i). When is M_1+M_2 again smooth?
The only possibilities, after exchanging the input terms, are

| Exponent difference | Reduced addition | Output exponent change from smaller |
|---|---|---|
| (0,0) | 1+1=2 | (1,0) |
| (1,0) | 1+2=3 | (0,1) |
| (3,0) | 1+8=9 | (0,2) |
| (0,1) | 1+3=4 | (2,0) |

Factor the coordinatewise gcd. Strict incomparable exponent pairs
leave 2^d+3^e with d,e>0, which is >1 and coprime to6. Strict
comparable pairs leave 1+2^d*3^e, also >1 and coprime to6. On the
remaining axes, 1+2^d=3^j: d=1 works, d=2 fails. If d>=3 then
j must be even by mod8. Factor the difference of squares; the two
positive powers of two 3^(j/2)-1 and 3^(j/2)+1 differ by2, so they
are 2 and4, giving d=3. For 1+3^e=2^j, parity of e makes the result
2 or4 mod8, hence only e=1 works.

Geometrically the admissible input pairs are horizontal distances
0,1,3, or vertical distance1; genuine diagonal adjacency is forbidden.
In the translated H-sector H=1+M, sums are relative to one shared
origin: H_1+H_2-1. Beginning at M=8, repeated doubling and the rule
M+2M=3M generate the entire upper quadrant A>=3,B>=0 (all classified
pairs preserve that quadrant). The S333 strong rail gives literal
unmarked tournaments for each endpoint of this arithmetic fan. This
does not assert one physical tournament gluing operation for its
edges.

Stronger no-clean-regrouping rule: take a primitive smooth divisibility
antichain of at least two terms (a strictly increasing, b strictly
decreasing). The minimum a and minimum b of any subcollection of
size at least two are each unique, so they are exactly the 2- and
3-adic valuations of its sum. After dividing out both minima, the
sum is >1 and coprime6; it is not smooth. All higher-arity attempts
that regroup the given terms and return to the smooth rail at any
intermediate stage are thus excluded, not just binary regroupings.

The S333 examples 34=18+16,45=27+18,69=27+18+24 have normalized
non-smooth factors 17,5,23 respectively after stripping those
extremal valuations. The arithmetically known availability of
individual monomials, even with the same baseline and all exponent
pairs present, is not an absolute-target inversion for their
antichain sums. A genuinely new interface must support nonsmooth
intermediates, alter the given decomposition, or combine the terms
while paying the full absolute interaction correction.
