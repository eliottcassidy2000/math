# The same Frobenius carry: cellular kernels, finite logarithms, and a cyclic norm

**Status:** PROVED elementary identities, with classical attribution and exact finite controls. These are a carrier-level synthesis and scoped derivations; no novelty or prime-finiteness claim is made.

The finite logarithm is classical. Relevant primary references are Philippe Elbaz-Vincent and Herbert Gangl, *On Poly(ana)logs I*, Compositio Mathematica 130 (2002), 161–214, [DOI](https://doi.org/10.1023/A:1013757217319), [arXiv](https://arxiv.org/abs/math/0008089), and their author-hosted [finite-polylogarithm account](https://www.maths.dur.ac.uk/users/herbert.gangl/Entropy_and_Polyanalogs_finalV1.pdf). A related primary paper is [*On the Zeros of Fermat Quotients and Mirimanoff Polynomials*](https://arxiv.org/abs/math/0604427). The formulas used here are proved below rather than imported from an uninspected theorem number.

## 1. Exact specialization map

For an odd prime p put

`H_p(X)=((1+X)^p-1-X^p)/p in Z[X]`.

This is exactly the carry polynomial `C_p(1+X)` from the scalar cellular-automaton note. If the homogeneous ternary defect is

`F_p(U,V)=U^p-V^p-(U-V)^p`,

then the exact integral specialization is

`F_p(1+X,1)/p=H_p(X)`.                                   (1.1)

The companion period-polynomial session realizes F_p as a cyclic norm of a specified binary form; that additional realization is not needed for (1.1). The specialization retains the full univariate polynomial, while evaluation at one point retains only a scalar and can hide cancellation.

Modulo p,

`binom(p,j)/p = (-1)^(j-1)/j` for `1<=j<p`.

Consequently

`H_p(X)=-L_(1,p)(-X)`,

where `L_(1,p)(T)=sum_(j=1)^(p-1) T^j/j` is the classical finite logarithm. This identification is an exact polynomial identity, not an approximate similarity.

## 2. Every non-endpoint rational root is exactly double

Differentiation gives the integer identities

`H_p'(X)=(1+X)^(p-1)-X^(p-1)`,

`H_p''(X)=(p-1)[(1+X)^(p-2)-X^(p-2)]`.

For every `r in F_p \ {0,-1}`, Fermat's theorem therefore gives

`H_p'(r)=0,  H_p''(r)=1/[r(r+1)] != 0`.                  (2.1)

Since 2 is invertible, every such r that is a root of H_p has multiplicity exactly two. The endpoints are different:

`H_p(0)=H_p(-1)=0, H_p'(0)=1, H_p'(-1)=-1`,

so both endpoints are simple. No assertion is made that every non-endpoint is a root; all are critical points, which is a weaker predicate.

At X=1,

`H_p(1)=(2^p-2)/p=2q_p(2)`,

where `q_p(2)=(2^(p-1)-1)/p` is the Fermat quotient. Hence

`p^2 | 2^(p-1)-1`

if and only if 1 is an exact double root of H_p modulo p. This is the Wieferich condition translated into a root predicate. It does not imply that H_p is zero as a polynomial: its coefficient of X is always 1. In the cellular interpretation, an aggregate specialization can cancel while the spatial carry persists.

## 3. The intrinsic order-three symmetry

Two exact integral identities are

`H_p(-1-X)=H_p(X)`,

`X^p H_p(1/X)=H_p(X)`.

On roots other than 0 and -1, the transformations `r -> -1-r` and `r -> 1/r` generate the usual S3 action, whose product has order three. Non-endpoint rational roots therefore occur in these orbits. If p>3 is Wieferich at base 2, the orbit of 1 is

`{1,-2,-1/2}`,

and every point in it is an exact double root. This orbit statement concerns one polynomial over F_p, not a new prime-distribution theorem.

There is also an exact characteristic-zero cyclotomic factor. Let omega be a complex primitive cube root of unity. For every prime p>3,

`H_p(omega)=0`,

because `1+omega=-omega^2` and the three cube roots sum to zero. Moreover

- if `p=1 mod3`, then `H_p'(omega)=0` and `H_p''(omega)=p-1!=0`;
- if `p=2 mod3`, then `H_p'(omega)=omega^2-omega!=0`.

Thus `X^2+X+1` divides H_p exactly twice when p=1 mod3 and exactly once when p=2 mod3. The quotient is integral because the cyclotomic polynomial is monic. At split primes p=1 mod3, these two cyclotomic roots lie in F_p and are consistent with the exact-double-root law (2.1). At p=2 mod3 they lie over F_(p^2), outside the rational-point hypothesis of (2.1).

## 4. Connection contract

Source: the integer Frobenius defect of the scalar rule `1+X`. Target 1: a classical finite logarithm, by coefficient reduction. Target 2: the homogeneous ternary defect, by `(U,V)=(1+X,1)`. Target 3: the base-2 Fermat quotient, by X=1. Preserved: exact integer identities, the polynomial-level carry, the root-multiplicity predicate at the stated points, and the explicit S3 symmetry. Lost at the scalar specialization: spatial coefficients and any full polynomial lifting obstruction. Required sidecars: field, prime, derivative normalization, endpoint exclusions, and the chosen integral lift.

In particular, none of the following follows from a single Wieferich specialization:

- vanishing of the Frobenius carry polynomial;
- a Frobenius-compatible lift of the normalized rule `1+X` modulo p^2;
- vanishing of an entire cyclic norm on an arithmetic lattice;
- a finiteness or infinitude theorem for Wieferich or Wilson primes.

The normalized scalar rule `1+X` has two distinct exponent classes modulo p, so the complete lift theorem in `SCALAR_TEMPLATE_DEPTH_CLASSIFICATION.md` proves that no integral Laurent lift can remove its full Frobenius defect modulo p^2, even at a Wieferich prime. The distinction between scalar cancellation and polynomial compatibility is therefore witnessed inside the same example.
