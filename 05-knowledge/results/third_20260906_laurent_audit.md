# Independent audit of the forced Laurent factors and endpoint 39

**ACCEPTED: PROVED all-height divisibility; PROVED endpoint-39 consumer from
an exact coefficient certificate and its inherited real-root supplier;
FINITE-EXACT independent coefficient reconstruction.** This audits
[the research note](third_20260906_laurent.md) without importing its producer.
No additional all-height sign theorem, theorem ID, or external priority claim
is made.

## Analytic audit

The normalized first row is monic. The inverse-carry quotient in equation
(2) is polynomial: splitting the doubled falling factorial into its even
and odd factors cancels exactly the first constant coefficient. Consequently
the reduced response and multiplication matrix are polynomial in the
parameter, including at negative integers where the original Laurent
expression cannot be evaluated at zero. This distinction is necessary and
is observed throughout the proof.

At `x=-r`, the first row has exactly an `r`-fold zero root and its constant
coefficient has a simple parameter zero. The actual lower carry has a simple
zero as well. After `delta=epsilon^r`, implicit lifting of the distinct roots
of `a(0)v^r+A` gives `r` small roots `t=epsilon*v(epsilon)`.
The carried response has exact leading order `epsilon^(r-1)`; all regular
terms have higher order. The coprime small and complementary factors of the
first row lift over `C[[delta]]`, so their characteristic coefficients have
integer delta orders. Therefore the `j`th small-block coefficient has order
at least `ceil(j(r-1)/r)=min(j,r-1)`. A coefficient of degree `k` in the full
characteristic polynomial must select at least `k-(h-r)` small eigenvalues.
This gives exactly the claimed exponent

```text
max(0,min(k-h+r,r-1)).
```

Distinct parameter factors are coprime, so their product divides the full
coefficient in `Q[x]`. The equivalent product formulas and degree reductions
follow by summing these exponents. This proves divisibility for every
`h>=1`; it does not assert exact multiplicities or positivity at arbitrary
height. The `r=1` case supplies no forced carried factor. Removing the carry
instead makes every small response have delta order at least one, giving
the different determinant divisor stated in the note.

For endpoint 39, the literal charge equation is `2b+3c=39` at mass `g`
and `2b+3c=78` at mass `2g`. Nonnegative solutions are exactly the seven
first and fourteen doubled channels displayed in the note. Their normalized
multinomial coefficients agree with its two polynomial rows, including the
`e=-1` channel. Reduction of the charge equation modulo `g` forces every
admissible positive mass to be a multiple of `g` when `gcd(g,39)=1`.
The first displayed channel shows that `g` itself is admissible.

The inherited [THM-4436, complete factorial-row roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md)
applies with `A=2,B=3,h=6,r=0,z=1,x=g-19`. Thus all six first roots
are distinct, strictly negative, and real. Positivity of every characteristic
coefficient makes the characteristic polynomial strictly positive on the
nonnegative real axis. Evaluating the real response at any first root gives
a real eigenvalue, which must therefore be negative. The doubled moment
cannot vanish at a first cancellation. Nonzero complex coefficients present
no extra sign assumption: the sign belongs to the normalized response, not
to the complex raw moment. Both asserted detection masses are attained by
the coefficient choices in the note.

## Independent exact reconstruction

The [audit source](../../04-computation/third_20260906_laurent_audit.py) imports
no repository implementation. Its route differs from the producer's symbolic
multiplication matrices and Newton traces:

1. At each integer `x=1,...,58`, enumerate literal nonnegative charge solutions
   for masses `g=x+19` and `2g`, then form their multinomial coefficients.
2. Let `P` be the monic first polynomial and let `Q=t*q` clear the actual
   Laurent carry. Compute
   `Res_t(P,z*t-Q)/P(0)` by exact rational elimination. Because `deg P=6`,
   the resultant equals `P(0)` times the characteristic polynomial of the
   response. Reducing `Q` modulo `P` only shortens the elimination.
3. Divide the six resulting scalar coefficients by their proved parameter
   factors and interpolate with rational forward differences. The proved
   residual degree bound is at most 57, so 58 nodes determine each polynomial.
4. Compare every reconstructed rational monomial coefficient with the frozen
   certificate, then test fresh parameters `x=59,71,83`.

All 208 coefficients agree exactly and are strictly positive. Their degrees
are `11,21,30,38,45,57`. The certificate is identified by its full output hash,
not silently accepted from whichever output happens to occupy its path.
The interpolation nodes alone would not prove positivity; the all-height
factor and degree proof makes the exact polynomial reconstruction decisive.

Separate controls check 45 boundary clusters through `h=9`, including `r=1`;
ten primitive physical parameters' first admissible masses; the nonprimitive
`g=21` first mass seven; and the quadratic retained/deleted-carry distinction.
A direct symbolic Bezout calculation independently recovers the indefinite
coefficient matrix and determinant `-9527204358067041`. This rejects only
the proposed coefficientwise matrix positivity, as scoped in the source.

[The frozen audit output](third_20260906_laurent_audit.out) records 2,436
always-active gates. Reproduce with:

```bash
python3 -B 04-computation/third_20260906_laurent_audit.py
python3 -B -O 04-computation/third_20260906_laurent_audit.py
```

General trinomial two-rung separation and all-height residual positivity
remain OPEN. No repair to the audited mathematical statements was required.

Normal and optimized outputs agree byte for byte. Raw LF SHA256 values:

```text
source 5577addc93b30dda32e1380d086839256482f06a0183186eb1183c7cba3b2e80
output 0a9bb58bc5502cfc0741e45608d4061d60e1ed8d8e82052b90cf80917b811b10
semantic 0bab1c033c32568970ac53a9d8d63ce62f70fbcb8843adc4307af022980382b4
```
