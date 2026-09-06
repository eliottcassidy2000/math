# Independent cubic carried-family proof audit

**PASS: analytic proof plus an independent exact certificate.** The producer's
all-height claim for support `(-21,2g-21,3g-21)`, integral `g>=11` and
`gcd(g,21)=1`, is accepted. It closes a fixed unbounded four-channel family;
general trinomial two-rung separation remains open.

The independently reconstructed source equation is `g(2y+3z)=21m`.
Consequently admissible primitive parameters have no support return before
`m=g` or strictly between `g` and `2g`. Solving in the gamma multiplicity
gives exactly the four and eight displayed fibers. Dividing the first raw
factorial row by its leading coefficient gives

```text
p_j=(g-7)_(3-j) * 7! / [(9-3j)!(1+2j)!],  j=0,1,2,3.
```

Thus `p` is monic with `deg_g p_j=3-j`. Each coefficient and the complete
cubic discriminant has strictly positive coefficients after `g=s+11`.
This independently proves three distinct negative first roots on the full
real parameter interval, without importing the incoming general root theorem.

The raw second row divided by `K=(2g)_14` is `Qbar`. The correct normalized
response is `Qbar/t`; the second source anchor is twice the first minus
`(1,-3,2)`. Dropping that factor reverses the rootwise sign and changes the
predicate. The first monomial factor is nonzero for every allowed coefficient
triple, even when complex; the real sign statement applies to its normalized
doubled value, not to an unordered complex scalar.

The only apparent rational-function denominator in reduction is the constant
term times `t^-1`. The exact factorial cancellation is

```text
Qbar_0/p_0=1152(g-10)(2g-15)(2g-17)(2g-19)/21!.
```

Give `t,g` weight one. All noninverse terms of `Qbar/t` have weight at most
six, and the inverse substitution
`t^-1=-(t^2+p_2 t+p_1)/p_0` has the same bound after the displayed
cancellation. Reduction modulo the monic cubic preserves the bound.
Therefore `R_j`, the reduced coefficients, have degree at most `6-j`.
Multiplication by `R` sends basis monomial `t^j` to coefficients of degree
at most `6+j-i` in row `i`. A principal minor of size `k` consequently has
degree at most `6k`: its row and column indices cancel in every permutation
term. Summing those minors proves degree bounds `6,12,18` for the three
characteristic coefficients. This is the proof that makes interpolation
complete rather than a finite-height extrapolation.

The verifier imports no producer mathematics. At the nineteen integers
`g=11,...,29` it rebuilds both raw fibers from the charge equation and takes
the resultant of `p(t)` and `t*x-Qbar(t)`, normalized to monic in `x`.
This polynomial has precisely the response values `Qbar(lambda)/lambda`
as its roots. Explicit monic normalization removes any resultant orientation
sign. The exact coefficient values match every published shifted polynomial,
and independent rational interpolation reconstructs all three identities.
Nonprimitive sample integers are valid algebraic identity checks; their
first-return interpretation is never used.

All coefficients of the response characteristic polynomial are strictly
positive for `g>=11`, including its constant coefficient. No response value
can be nonnegative. Their reality comes from the already proved first-root
reality. Hence every first cancellation gives a nonzero doubled moment,
and the first detection is exactly `g` or `2g`, with both alternatives
attained. The width is `3g`; no sharp-width equality is claimed.

The audit retains the gcd-dropped `g=12` first mass four and the spectrum
`(-5,1,1)` that defeats trace-plus-norm without the middle coefficient.
It checks the full-domain discriminant, exact denominator cancellation,
nineteen complete source fibers, Sturm root counts, and coefficient identities.
All **218** gates survive optimization, with byte-identical output.

```bash
python3 04-computation/overnight6_20260906_laurent_cubic_carry_audit.py
python3 -O 04-computation/overnight6_20260906_laurent_cubic_carry_audit.py
```

```text
SHA-256, LF bytes:
source edb2a7d3fa43324cd18210a074db76297a753cc4ba1886043c2ad0da838c09e6
output 937a222475587451726438f0b6919a5f12a8d4b4a3316f134965b11c30bd681d
published coefficient record 88d4c05a95fa48066151b72063c5f9fffb0cbecb18f79ca91c975da63f935315
```

The published coefficient entries are positive **rationals**, not all integers;
that presentational correction does not change the certificate. A separate
referee also reviewed the source map and the weighted-degree argument.

**Filing:** root integrated this audited report after `f5f0f7f75`;
portable reproduction paths are shown above. The exact producer and
transcript bytes remain pinned by the sixth manifest.
