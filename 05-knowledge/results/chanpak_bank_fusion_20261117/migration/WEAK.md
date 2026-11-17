# The exact content boundary for weak root migrations

The [root transfer](PROOF.md) remains exact for every standard
binary word. A weak stage may have genuine physical effects
not captured by the restricted weighted support bound. Its
algebra does, however, give a precise filter for raw prime targets.

Write a standard odd contact address as
u=2^t(v+1)-1, with t>=1 and v even (or v=0).
Let k be the standard bit length of v (k=0 for v=0) and let
c=2^k-rev_k(v), with c=1 at k=0. The root transfer simplifies to

    h'=F(v)h + 2^t(v+1)R,
    R'=c h + 2^(k+t)R,
    S'=h.

Indeed trailing ones leave F unchanged and become leading zeros
in the reverse complement, so 1+RC is still c. Thus the weak
tail scales one column of the two-coordinate read without
the need to change v or any previous module. It preserves
the exact determinant coordinate

    2^k h' - (v+1) R' = (2^k F(v)-(v+1)c)h.

This determinant vanishes only at v=0 by the zero-run proof
in the main note. For a given incoming (h,R,S) and a chosen
even prefix v, a prescribed raw count N therefore chooses its
weak tail exactly: (N-F(v)h)/((v+1)R) must be a power 2^t
with t>=1. The graphs are physical and keep the unit end.

## Integer content: exactly when a weak stage can lock out primes

Define content as gcd(h,R,S), and restrict to the standard
(nonpadded) alphabet and singleton stages. Counts h are odd.
For a nonempty standard **even** word, the new content equals
the old one: a divisor of the new triple divides h=S', then
R by R'−c h=2^k R (the divisor is odd), and then S by the
first row; the reverse divisibility is immediate.
For a nonempty standard **odd** word, the new content is exactly
gcd(h,R), by the same argument and the absence of S from the
first row. A singleton stage is unimodular: from its output
the old coordinates are R'−h'+S', h'−S', R'−h', respectively,
so content is also unchanged.

Thus content can strictly increase only at an odd word,
precisely when the hidden S had been the obstruction to the
larger common divisor of h and R. Once a nonunit content is
created, all subsequent raw counts stay divisible by it; under
strict count growth they cannot be prime. Two consecutive
all-ones contact words necessarily lock this way: the first
makes R'=h' and the second inherits nonunit content h' if
h'>1. This explains why the full arbitrary-weak bank is not
just a rich store of interchangeable prime-target resets.
It does not imply that a downstream two-cover or sevenfold
reader has the same divisibility.
