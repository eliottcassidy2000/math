# A three-parameter physical family explaining the ternary clock

## Physical family

Take the six-vertex quotient tournament

`Q = [2,4,41,19,39,11]`

with vertices, in this order, `a,b,U,c,V,R`. Replace `U,V,R` by forward transitive chains of positive lengths `p,q,m`; keep `a,b,c` singleton. Edges between blocks follow Q. Call this actual tournament `G(p,q,m)`. The quotient is strong, so every graph in the family is strong. This is a fixed quotient and literal module substitution, not a formal assignment of response numbers.

The previously audited eight-vertex ternary clock is exactly `G(2,2,1)` up to reversing the labels inside the two length-two chains. Keeping `a,b,c` exposed, one rightward birth beating those three ports and losing to every other old vertex replaces `m` by `m+1`, retaining the same ports.

## Exact all-order formula

For every positive `p,q,m`, put

`A = 4(3^p-2^p)(2^q-1)`,

`B = 2*3^p+3*2^p-6 - (2*3^p-3*2^p)2^q`,

`C = (2^p-3)(2^q-3)-2`.

Then

**H(G(p,q,m)) = A 3^m + B 2^m + C.**

Here is an exact finite-dimensional proof of this identity. The root-skeleton decomposition for a transitive module gives a linear combination of `j^t`, where the possible number of module runs is one plus at most the size of a vertex cover of its out-to-in crossing graph. Between successive module runs, a disjoint internal fragment must cross from an external out-neighbor to an external in-neighbor. These selected crossing edges are disjoint, proving the run bound directly.

The following bounds hold uniformly in all three chain lengths:

- For module U, every crossing edge either starts at the singleton c or ends at the singleton b. Hence U has at most three runs.
- For module V, the only external in-neighbor is singleton c. Hence V has at most two runs.
- For module R, every crossing edge starts at singleton b or singleton c. Hence R has at most three runs.

Consequently H lies in the tensor span of the 18 functions `i^p j^q k^m`, with `i,k in {1,2,3}` and `j in {1,2}`. To justify the tensor statement without assuming simultaneous independence, first recover the p-coefficients as fixed rational linear combinations of the evaluations at p=1,2,3. Each such evaluation separately has the stated q- and m-spans. Repeat for q and m.

The proposed expression also lies in this tensor span. The positive-evaluation matrices for distinct bases are invertible Vandermonde matrices, multiplied by nonzero diagonal matrices. Thus two such expressions agree everywhere if they agree on the full grid `p=1,2,3; q=1,2; m=1,2,3`. The supplied standalone verifier counts Hamiltonian paths of all 18 actual graphs by ordinary subset DP and confirms that grid exactly. Its receipt records every graph parameter and count. This proves the identity for all positive integer parameters. Seven checks outside the interpolation grid are additional controls, not premises for extrapolation.

## Unique cancellation within this entire quotient family

The coefficient of the binary mode is zero if and only if `p=q=2`.

For p=1, B=6. For p>=2, the denominator `D=2*3^p-3*2^p` is positive. Setting B=0 would force

`2^q = (2*3^p+3*2^p-6)/(2*3^p-3*2^p)`.

Because q>=1, this ratio must be at least two, requiring `2*3^p <= 9*2^p-6`. This fails at p=4 and for every larger p: the positive difference at p=4 increases under the next step, as is seen either directly or by induction. For p=2 the ratio is four, giving q=2. For p=3 the ratio is 12/5, not a power of two. Thus the only binary cancellation is p=q=2. At those values A=60 and C=-1, recovering the returning clock.

This supplies an all-order structural explanation of the numerical cancellation. It does not show that this quotient is unique among tournaments, or that no other context can give a returning ternary clock.

## Reproduction and scope

Run `python3 verify.py` and `python3 -O verify.py` in this directory. All checks use explicit exceptions. The code checks the tensor grid, outside controls, graph isomorphism to the audited seed, crossing-edge covers, and finite controls of the cancellation algebra. The uniform cover arguments and the displayed inequality prove the all-order statements; the finite controls alone would not suffice.

The parameters p and q are fixed core choices when using the one-birth returning interface. Allowing target-dependent changes to those parameters is a different operation. No method here chooses parameters for every desired odd count, and no universal ternary grammar or missing-spectrum result is asserted.
