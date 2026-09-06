# Independent audit of the endpoint-33 six-channel return theorem

**Status: PASS — complete core analytic proof and independent polynomial
certificate accepted after review of the finalized producer prose.** This audit concerns the endpoint-33
family in [second_20260906_laurent.md](second_20260906_laurent.md), not the
general trinomial doubled-row problem. The reviewer is the separate Smith
agent. The independent verifier imports no producer or SymPy code.

## The actual first and doubled rows

Write `g=x+16`, where `x>=1` is an integer, and consider support
`(-33,2g-33,3g-33)` with three nonzero complex coefficients `alpha,beta,gamma`.
If `gcd(g,33)=1`, every support return has mass divisible by `g`: all three
exponents are `-33` modulo `g`. The first possible mass `g` is nonempty.
Writing `tau=alpha*gamma^2/beta^3` and `X=alpha^x beta^15 gamma`, the complete
mass-`g` count triples are precisely

```text
(x+j,15-3j,1+2j),       j=0,...,5.
```

Thus the monic first polynomial is exactly

```text
CT(f^g)=X*binom(g,11)*p_x(tau),
[p_x]_j = 11! (x+5)_(5-j)/[(15-3j)!(1+2j)!].
```

At mass `2g`, the full triples are

```text
(2x+e,30-3e,2+2e),       e=-1,...,10.
```

In particular the `e=-1` lower carry is real and is allowed for every `x>=1`.
It cannot be omitted. With `K=(2g)_22>0`,

```text
CT(f^(2g))=X^2*K*Qbar_x(tau),
[Qbar_x]_e=(2x+10)_(10-e)/[(30-3e)!(2+2e)!].
```

These formulas were checked against direct enumeration of **all** nonnegative
count triples at 51 different values of `x`, not by substituting only the
predicted channels. The first-mass gcd condition is independent of the
characteristic certificate: when it fails, earlier masses can occur, even
though the displayed rows and characteristic identity remain valid.

## Polynomiality, the degree bound, and why positivity decides the target

Since `p_x` is monic of degree five and its constant coefficient is nonzero,
`tau` is invertible modulo `p_x`. The only possible nonconstant denominator in
the remainder of `Qbar_x` is removed by the exact identity

```text
[Qbar_x]_-1 / [p_x]_0
 = 64*15!/(33!*11!) * x*product_(j=0)^4 (2x+2j+1).
```

This is a polynomial of degree six. Reducing the inverse of `tau` and all
higher powers gives `R_x(tau)=sum_(j=0)^4 R_j(x)tau^j` with
`deg_x R_j<=10-j`. Indeed `[p_x]_j` has degree `5-j`; monic reduction preserves
total degree when `x` and `tau` each have weight one. The lower-carry quotient
has degree six and multiplies `[p_x]_(j+1)` of degree `4-j`, giving the same
bound.

Let `M_x` be the matrix of multiplication by `R_x` on the basis
`1,tau,...,tau^4`. Its `(i,j)` entry has degree at most `10+j-i`.
Every term of a `k`-by-`k` principal minor therefore has degree at most `10k`:
the selected row and column index sums cancel. Hence

```text
det(zI-M_x)=z^5+c_1(x)z^4+...+c_5(x),
deg c_k <= 10k.
```

This independently establishes the polynomial degree bound used by the
finite interpolation certificate. It is not inferred from sampled values.

The parameter substitution into
[THM-4436 / complete factorial-row simple negative roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md)
is exactly `A=2,B=3,h=5,r=0,z=1`, with the same nonnegative integer `x`.
Thus `p_x` has five distinct negative real roots. Evaluation at those roots
diagonalizes the multiplication algebra over the reals, so the five eigenvalues
of `M_x` are precisely the five real numbers `Qbar_x(tau_i)`.

The independent certificate proves that every coefficient of every polynomial
`c_k(y+1)` is strictly positive. For `x>=1`, each `c_k(x)>0`; consequently
`det(zI-M_x)>0` for every real `z>=0`. Since all five eigenvalues are real,
they must all be strictly negative. Both ingredients matter: positivity of
characteristic coefficients alone would not give real-negative eigenvalues
for an arbitrary matrix.

Therefore at every cancellation of the actual first row,
`Qbar_x(tau_i)<0`, and in particular the actual doubled row is nonzero.
The normalized real sign is
`CT(f^(2g))/(X^2 K)<0`; no sign is assigned to the unnormalized complex
moment. Congruence then proves that the first nonzero moment is exactly
`g` or `2g`. For each allowed `g`, generic phases give `g`, while each of the
five attainable negative roots gives `2g`.

## Independent exact certificate

The [standalone verifier](../../04-computation/second_20260906_laurent_audit.py)
uses only Python's standard library. For `x=1,...,51` it enumerates the complete
first and doubled count fibres, normalizes their exact rational coefficients,
constructs the quotient multiplication matrix, and computes its characteristic
coefficients by sums of principal minors. This is distinct from the producer's
symbolic polynomial arithmetic and Newton trace identities.

For each `k=1,...,5`, Newton forward differences interpolate `c_k(y+1)` from
its first `10k+1` exact values. The analytically proved degree bound makes that
a polynomial identity. Remaining samples independently check the reconstruction
for `k<5`. Every one of the resulting
`11+21+31+41+51=155` monomial coefficients is strictly positive. The full
rational coefficient data are frozen in
[second_20260906_laurent_audit.json](second_20260906_laurent_audit.json).

All five producer `CHAR_CERTIFICATE` rows were separately converted to exact
rational arrays and compared with this interpolation certificate: all 155
coefficients agree. That comparison uses the frozen producer output as data;
the independent computation imports no producer code.

The verifier passes 1,438 always-active exact gates. Normal and optimized
output streams and regenerated certificate bytes agree. Reproduce with

```text
python3 -B 04-computation/second_20260906_laurent_audit.py
python3 -B -O 04-computation/second_20260906_laurent_audit.py
```

Matching output: [second_20260906_laurent_audit.out](second_20260906_laurent_audit.out).

```text
source_sha256:
069897cf1b9ced10fe36f176d07a8b48330d0bed29184f5c9c852b5ce62c4560
output_sha256:
2a01aee2363137e82c3b902b390ca3e869f05fe7d30098959f04b0a90a28363c
certificate_sha256:
52a5fc3095a34639849f540e7c7d4a4819452f5fedd3fa9cad9238c29c233a08
```

The interpolation is a finite exact certificate of universal polynomial
identities under the degree bound. It is not a numerical root sample or an
unsupported extrapolation from finitely many scales.


## Final written review and ancillary scope

The final producer note was read in full. Three prose/type repairs were
requested and applied before acceptance: the symbolic quotient starts over
`Q(x)` before its remainder is shown polynomial; the `g=18` hostile concerns
first admissible support mass; and the unretained symbolic-leading-term claim
was removed from the large-scale multiplier example. The precise `W_x` used
by that finite multiplier test is now displayed. No theorem formula changed.

The derivative-cone appendix is correctly restricted to constant nonnegative
combinations of the phase-function response rays. Its separating linear
functional proves that failure even after reduction at the first roots; it
does not exclude arbitrary phase-dependent weights or all multipliers. The
large-`x` multiplier appendix is likewise a finite exact hostile to one stated
normalized ansatz, not an additional return obstruction. Their formulas and
scope were reviewed against the producer source, but the independent
interpolation verifier certifies the endpoint-33 theorem itself rather than
independently recomputing those two ancillary hostiles. Neither hostile is a
dependency of the positive endpoint-33 proof.
