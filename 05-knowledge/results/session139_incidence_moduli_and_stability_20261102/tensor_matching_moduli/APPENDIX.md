# Tensor compensation in every integral residue coset

Status: **PROVED and independently audited**. This appendix is separate
from the frozen parity/Pfaffian proof and the Eisenstein–Bernstein proof.
It generalizes the elementary tensor lemma only. No stronger tournament
Pfaffian estimate, dimension cutoff, arbitrary-modulus phase proof, or
Riemann-hypothesis statement is asserted.

## Theorem

Fix an integer modulus `M>=2`, an integer residue representative `rho`, and

    d = min_{k in Z} |rho+kM|.

For an integer dimension `s>=0`, take integer coefficients `b_U`, indexed
by subsets of `[s]`, satisfying

    b_empty in rho+MZ,
    b_U in MZ for every nonempty U.

For real `|x_i|<=1`, set

    f_S(x) = sum_{J subset S} b_(S\J) product_{j in J} x_j.

Then the sharp residue-dependent inequality is

    sum_{S subset [s]} f_S(x)^2 >= d^2 product_{i=1}^s (1+x_i^2).  (1)

Empty products are 1. If `rho=0 mod M`, then `d=0` and (1) is the
nonnegativity of a sum of squares; the proof still applies without division
by `d`.

## Proof

Induct on `s`. For `s=0`, the single integer `b_empty` is in the specified
coset and therefore has square at least `d^2`.

For the last coordinate, split the transform into `(a,c+x_s a)`, where
`a` is the transform of coefficients not containing `s` and `c` is the
transform of coefficients containing `s`. The coefficients producing `c`
are all multiples of `M`. Hence the coefficient vectors producing `a`,
`a+c`, and `a-c` each have their constant coefficient in the **same** coset
`rho+MZ` and all other coefficients in `MZ`.

Put

    L=d^2 product_{i<s}(1+x_i^2),
    A=||a||^2,   C=||c||^2,   D=<a,c>.

The induction hypothesis gives `A>=L` and `A+C+2D>=L`, `A+C-2D>=L`.
Consequently `2|D|<=A+C-L`. With `t=|x_s|<=1`,

    ||a||^2+||c+x_s a||^2
      >= (1-t+t^2)A+(1-t)C+tL
      >= (1+t^2)L.

The final comparison uses `C>=0` and the nonnegative coefficients
`1-t+t^2`, `1-t`. This proves (1) in every dimension and for every modulus.
No primality, coprimality, or oddness hypothesis on `M` or `rho` is needed.

## Sharpness and boundary

Choose any minimum-absolute-value representative `b_empty` of the coset,
so `|b_empty|=d`, and set all other coefficients to zero. Then
`f_S=b_empty product_{i in S}x_i`, and equality holds throughout the cube.
This includes the zero-residue case and ties when `M` is even and
`rho=M/2 mod M`.

The cube restriction cannot be omitted. For `M=2`, `rho=1`, one variable,
`b_empty=1`, `b_{1}=-2`, and `x=2`, the left side is `1+(2-2)^2=1`, while
the proposed right side would be `1+2^2=5`.

## Scope relative to tournament Pfaffians

Taking `M=2`, `rho=1` gives the previously audited odd/even tensor lemma.
Ordinary binary tournament Pfaffian grouping establishes exactly that
parity pattern: the base Pfaffian is odd and the other grouped coefficients
are even. It does **not** automatically place those coefficients in the
hypotheses for a larger modulus. Therefore (1) does not, by itself, improve
the matching-support determinant floor or the nonlinear rounding reserve.
The separate Eisenstein–Bernstein argument uses a specific F4 nonvanishing
mechanism; this appendix does not assert that mechanism for every modulus.

`verify_moduli.py` supplies exact integer/Fraction controls across multiple
moduli, residue classes, dimensions, rational points, equality examples,
and explicit invalid-assumption mutations. These controls check the
implementation and edge cases; the induction above, not finite sampling,
is the all-modulus proof. Both normal and optimized runs retain all checks.
