# Independent audit of beta-step propagation and its model hostiles

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS; FINITE-EXACT replay.**
The [all-h note](long_frontier_sep06_allh.md),
[producer](../../04-computation/long_frontier_sep06_allh.py), and
[frozen output](long_frontier_sep06_allh.out) were read and checked.
The all-h parameter identities, both exact model obstructions, and the
eight genuine fixed-row controls pass. No mathematical correction was
required. A general actual all-h doubled-response sign theorem remains open.

## 1. Full fibres, normalization, and the parameter step

For a=6h+3 and g=x+3h+1, the charge equation at mass g reduces to
`2 n_beta+3 n_gamma=6h+3`. Its complete nonnegative solutions are
`(x+j,3h-3j,1+2j)`, 0<=j<=h. At mass 2g they are exactly
`(2x+e,6h-3e,2+2e)`, -1<=e<=2h. In particular the doubled lower carry
has first count 2x-1>=1 and cannot be discarded.

Dividing the literal factorial coefficients by `binom(g,2h+1)` and
`(2g)_(4h+2)=(2g)!/(2x+2h)!`, respectively, gives the printed monic P
and carried Q. The monomial ratios are precisely tau and tau^e, with
the first anchor X and second anchor X^2. Both normalization factors
are positive and have been retained in the actual controls.

Taking coefficient ratios at x and x+1 yields

```text
P ratio at j: (x+h+1)/(x+j+1),
Q ratio at e: (2x+2h+1)(2x+2h+2)/((2x+e+1)(2x+e+2)).
```

These independently verify both operator identities for every permitted
h,x and every full coefficient index, including e=-1 and both top indices.
The two integral kernels then follow from the elementary integrals of
s^(x+j) and s^(2x+e)(1-s). The carry's exponent is 2x-1>=1, so the second
integral is convergent with no regularization or support truncation.
The source's finite symbolic h checks are controls of this derivation,
not the proof of its unbounded h quantifier.

If gcd(g,a)=1, every return mass must be divisible by g because all support
exponents are congruent to -a modulo g. The exhibited mass-g fibre proves
that g is the first support return. No such interpretation is forced onto
the nonprimitive rows.

## 2. Both model sign reversals are exact and correctly scoped

A fresh standard-library Fraction reconstruction, without importing the
producer, multiplied each three-factor starting cubic, applied the
coefficient step, evaluated the genuine first roots, and used the explicit
cubic discriminant formula. It recovered:

| Step | Before response | After response | After discriminant |
|---|---|---|---|
| h=1, x=1 to 2 | -3/43 | 557/5418 | 1152093393330275/56737436781312 |
| h=1, x=3 to 4 | -874/7575 | 221/21210 | 859047951135198760589/2467080636572160000 |

The shifts from Laurent exponent e to polynomial exponent k=e+1 give
exactly the multipliers 30/((k+2)(k+3)) and 90/((k+6)(k+7)). The same
positive scale 720 matches the actual carry and leading coefficients on
both sides of each transition. The two interior coefficients are different
from the genuine row and are not silently identified with it.

The starting factors give three distinct negative roots. The first
transformed cubic's three strict negative sign brackets exhaust its degree.
For the second transformed cubic, positive discriminant gives three distinct
real roots and strictly positive coefficients exclude every nonnegative
root. Thus both transformed negative-root predicates are justified.
The responses are M(t), obtained by dividing the cubic tM(t) by the same
negative first-root phase; omitting that division would reverse the test.

The first example has designated masses 5 and 6, with gcds 1 and 3 against
9. The second has masses 7 and 8, both coprime to 9. It therefore removes
the missing-first-return-gcd explanation while still remaining a model
with incorrect interior factorial amplitudes. The minimality statement
concerns h=1 and the smallest allowed x-step, not an actual Laurent sign
counterexample or a claim of uniqueness of the model.

## 3. The eight genuine controls use the correct quotient response

The exact universe is h in {7,8} and x in {1,2,10,100}: eight complete
coefficient rows and 4*7+4*8=60 first-root phases. Exactly six rows have a
primitive designated first mass; the x=2 row for each h is a coefficient
control. This count was independently checked.

The auxiliary B_j support is the full interval -x<=j<=h. Its convolution
is formed before the completed-alpha factor suppresses exponents below
-1. The binomial factor and divisor in W agree with Q's positive
normalization `(2g)_(4h+2)`. Hence the printed Q-W is the retained skip
response in the declared convention.

The quotient routine treats the carry correctly. Since p(0)>0, the
identity in Q[t]/(p) is

```text
t^-1=-(p1+p2 t+...+ph t^(h-1))/p0.
```

The producer adds this exact inverse to the ordinary remainder before
forming the multiplication matrix. Its eigenvalues are the response
values at P's roots. The inherited
[THM-4436, complete factorial-row simple negative roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md)
applies with A=2,B=3,r=0,z=1, including the nonprimitive coefficient
controls; it supplies the distinct real first roots. Therefore all matrix
eigenvalues are real. Strictly positive characteristic coefficients exclude
every nonnegative eigenvalue and prove all the reported response values
strictly negative. No sign inference is taken merely from nonnegative
coefficients of Q-W itself.

Literal charge enumeration independently reconstructs every first and
doubled coefficient in all eight named rows. The finite signs are not
promoted to an endpoint-parameter theorem or an all-h result.

## 4. Replay and manifest

Normal and optimized runs each pass **142 always-active exact gates** and
agree byte-for-byte with the frozen output. The full source contains no
disabled assertion, numerical root solver, or imported research producer.
SymPy is used for exact symbolic identities and rational characteristic
polynomials; the additional model reconstruction used Fraction arithmetic
and the explicit cubic discriminant expression.

```bash
python3 04-computation/long_frontier_sep06_allh.py > /tmp/allh-normal.out
python3 -O 04-computation/long_frontier_sep06_allh.py > /tmp/allh-optimized.out
cmp /tmp/allh-normal.out /tmp/allh-optimized.out
cmp /tmp/allh-normal.out 05-knowledge/results/long_frontier_sep06_allh.out
```

Checked SHA-256 manifest:

```text
source  7cc83c7518c27e70db33936809e119a315dc16e1ceb9573ea40bce939d42b1ad
output  bb012dd9cf4b00b904d1cf74746e3455dc753db16fb308bd94d67f4e68bdf07e
semantic a9d05da8f843b6365851e248308082bfbd5d57a7da0379b91c5bddd3726a891f
```

The first failed implication is propagation of the old same-root sign from
the retained qualitative data. The exact beta identities survive; the
missing coordinate is the interior factorial amplitude law. Neither the
actual two-first-channel theorem nor any actual Laurent return claim is
refuted by these models. This audit introduces no new exploration.
