# Exact order-seven attachment dimension and an all-order lower bound

The setting, definitions, all-order carrier upper bound and Hamilton coefficient law are those in the Session 145 `POOLED_PARITY_CARRIER.md`. In particular, attachments fix all arcs with new endpoints independently of the old source. Source measures, updates and dimensions are rational linear; no claim is made for a restricted S/B-only port language, source-dependent attachments, conditioning or shared-source products.

## Theorem

At seven named old vertices, the dimension of the span of all fixed attachment mean responses is exactly **2668**. For normalized laws the constant coordinate is known, leaving exactly **2667** varying coordinates.

The prior pooled-parity carrier has 2668 coordinates and is closed, giving the upper bound. For the lower bound the attached certificate contains 2668 exterior coefficient contexts. Its independent verifier constructs the 2668 representative signed linear forests, computes the coefficient matrix using the exact full Hamilton coefficient law, and reports ranks 2668 modulo 101 and 103, with determinants 48 and 66 respectively. A nonzero determinant modulo either odd prime proves a nonzero rational determinant, since all matrix entries have powers of two as denominators. Boolean multilinear interpolation expresses each exterior coefficient row as a rational linear combination of actual fixed attachment responses. Thus the lower and upper bounds match.

The certificate has 2668 distinct contexts, total enlarged vertex counts between 7 and 23, and 0 through 19 exterior edges. Concatenate `contexts7.part00`, then `contexts7.part01`, without inserting bytes. The SHA-256 of the result is `ee0c332d4dd580eb3073cf95eba3ca098d005cd1f3f26e14c0aaaa9c6671baef`. Each line is `N k u1 v1 ... uk vk`; old vertices are 0 through 6 and all listed edges have at least one new endpoint. `generate_rank7.cpp` is the deterministic discovery program. Independently verifying the saved bank, not replaying the random search, is the certificate check.

## Lifting and lower bounds at every higher order

For every n ≥ 7, let d(n) be the rational response dimension and c(n) the pooled-parity carrier dimension. Then

    1 + binom(n,2) + 3 binom(n,3) + 9 binom(n,4)
      + 40 binom(n,5) + 120 binom(n,6) + 546 binom(n,7)
    ≤ d(n) ≤ c(n).

In particular, 10795 ≤ d(8) ≤ 13028, 36199 ≤ d(9) ≤ 65908, and 103096 ≤ d(10) ≤ 347626. The lower bounds include the constant coordinate; subtract one for normalized laws. No equality is claimed at n ≥ 8.

**Proof.** Choose any seven-element set A of the n old labels, in increasing order. Take an order-seven exterior-coefficient context G and relabel its old vertices onto A, placing its exterior vertices after all n old vertices. For each of the m=n−7 remaining old labels v, add two fresh leaves a,b with a<b, each joined to v, and leave t further fresh exterior vertices isolated. An old forest F containing an edge at any such v has coefficient zero because its union with the fixed leaves has degree at least three there. Every other old F is supported on A.

For an F supported on A with nonzero original coefficient, write c for the number of components including isolates in its union with G. The m added leaf paths are disjoint two-edge components. Relative to the original coefficient, the new coefficient has the multiplier

    (-1)^m (c+1)(c+2)...(c+m) · 2^(-t)
       · (c+m+1)(c+m+2)...(c+m+t).

The sign comes from orienting each new path from its smaller leaf through the smaller-numbered old vertex to its larger leaf; the products are 1 when their index is empty. For a fixed original context, c lies in the finite set {1,...,N}, where N is its enlarged vertex count. As t runs from 0 to N−1, the polynomials `(c+m+1)...(c+m+t)` form a basis for polynomials of degree at most N−1. Rational interpolation therefore chooses a linear combination of the t contexts whose multiplier is 1 at every c in this finite set. This constructs, in the n-vertex response span, the original row extended by zero outside the forests supported on A. Vanishing original coefficients remain zero.

Because the certified seven-vertex coefficient matrix is invertible, these lifted rows span each individual pooled-parity coordinate supported on A. Varying A produces every n-vertex signature whose non-isolated support has at most seven vertices. Distinct signatures have disjoint forest-monomial support and hence give linearly independent functions on the Boolean orientation cube. Their number is the displayed polynomial: the counts with exactly k non-isolated named vertices for k=0,...,7 are 1,0,1,3,9,40,120,546. They follow by binomial inversion from c(0),...,c(7) = 1,1,2,7,28,126,571,2668. The stated upper bound is the previously proved all-order closed carrier. □
