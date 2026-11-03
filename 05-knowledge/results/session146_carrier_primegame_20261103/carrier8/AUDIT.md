# Independent order-eight fixed-attachment quotient audit

The 13,028 pooled-parity signatures on eight old vertices consist of 10,795 with at least one isolated old vertex and 2,233 with no isolates. The order-seven result and lifting argument make the first 10,795 coordinate polynomials observable. Projecting fixed-attachment coefficient rows onto the remaining 2,233 coordinates, the supplied bank has:

| Field | Quotient rank | Determinant |
| --- | ---: | ---: |
| GF(101) | 2,233 | 70 |
| GF(103) | 2,233 | 94 |

Replacing the last row with the first gives rank 2,232 mod 101; deleting the last row gives rank 2,232 in each field. The context bank SHA-256 is `94df48ba7e434651698928b156d81ef253c39cf76f9df8c20a4c87b6adcad57b`. Over the rationals the projected rank is full, and combining the rows with the already-observable 10,795 coordinates proves dimension 13,028, matching the Session145 upper bound. On normalized laws this is 13,027 nonconstant coordinates. This concerns the full fixed source-oblivious attachment language, not S/B-only operations or arbitrary nonlinear protocols.

To lift an order-seven coefficient row to eight old vertices, use two fresh exterior leaves incident to the eighth old vertex. In any contributing old forest this outside vertex must be isolated, because its selected exterior degree is already two. The added length-two component changes the coefficient by a fixed sign times (c+1), where c is the number of components in the base forest. Adding t exterior isolates changes the ratio to the same sign times 2^(-t)(c+1)_(t+1). The possible positive integers c form a finite set, on which the polynomial factor c+1 is nonzero. A linear combination of these rising-factorial ratios therefore interpolates 1 over that set. This recovers each unweighted embedded order-seven row. Inverting the verified order-seven certificate and varying the seven-subset gives every signature with an isolate. The corresponding polynomials have disjoint Fourier monomial supports, so they are independent.

The independent verifier reconstructs the 2,233 no-isolate signatures from canonical simple paths and exact-cover recursion anchored at the least unused vertex. It checks representative masks and signs against `reps8.txt`, validates each external context, calculates coefficients by disjoint-set cycle detection and path traversal, and performs exact finite-field elimination. The paths-and-cover enumeration is independent of the producer's permutation-splitting catalogue. The companion controls check nine direct Hamilton permutation expansions on N=8, an exact nonzero Walsh extraction from fixed cuts, and six malformed context classes. All controls give identical output with Python normal and -O. Exterior coefficients are finite rational combinations of fixed-cut evaluations, not individual attachment evaluations. The proof inherits Session145's coefficient formula and all-order upper bound and the independently certified order-seven rank.

Reproduce:

```sh
cat contexts8.part00 contexts8.part01 > contexts8.txt
g++ -O3 -std=c++17 -o verify verify8.cpp
./verify 8 contexts8.txt reps8.txt
./verify 8 contexts8.txt reps8.txt --delete-last
python3 -B controls.py > controls.normal.txt
python3 -B -O controls.py > controls.optimized.txt
cmp controls.normal.txt controls.optimized.txt
sha256sum -c SHA256SUMS
```

If the bank is split for distribution, concatenate numbered chunks in their prescribed order to a temporary file, check its SHA-256 above, and use the temporary file as the second argument of `verify`. The representative list is compared with the independently constructed catalogue, not taken on trust.
