# Separate independent audit: direct cubic/quintic Keller exclusion

20 October 2026. **PASS.**

## Scope and frozen sources

This is a separate audit of `CUBIC_QUINTIC_ADDENDUM.md` and
`verify_cubic_quintic.py`. It proves that no pair with respective y-degrees
3 and 5 and constant nonzero leading coefficients has nonzero constant
Jacobian. Applying this to the generic semigroup `<3,5>` additionally uses
the already audited characteristic classification and constant-leading
terminal reduction.

Reviewed candidate hashes:

```
037ed17075b13f84664eedee1586312be25b507b6c081346ff4f9edd4fa432ce  CUBIC_QUINTIC_ADDENDUM.md
892ae3b17084ef3c3898e31aee8fef11d9769034ff388613735d43000a3745f0  verify_cubic_quintic.py
```

The producer's **31** checks pass normally and with `-O`; the byte-identical
outputs match the staged receipts. The separate reviewer checker
`independent_cubic_quintic.py` passes **32** checks normally and with `-O`,
also byte-identically. It uses symbolic functions of x and actual
differentiation, not the producer's manually assembled coefficient
derivation operator. Its hostile identity is verified for an arbitrary
symbolic function u(x), not only a finite list of powers of x.

The proof below is elementary and universal; the finite examples are
controls, not a substitute for the degree argument.

## 1. Lawful normalization and integration

Constant target scalings make both leading coefficients one and multiply
the nonzero constant Jacobian only by a nonzero constant. Centering the
cubic by `t=y+s(x)` is a determinant-one polynomial source automorphism.
The entire centered pair therefore has precisely the stated shape

```
p=t³+a t+b,
q=t⁵+c t⁴+d t³+e t²+f t+g.
```

I independently differentiated this pair. Every raw coefficient, including
the scalar coefficient `f b'−a g'`, is as displayed in the addendum.
The coefficients of t⁶,t⁵ force c constant and d=5a/3+D with D constant.
The shear `q ↦ q−Dp` is legitimate: D is a complex constant and the
coefficient g is renamed after its change by `−Db`.

With d=5a/3, the remaining positive-degree equations in descending order
give

```
e' = (4c/3)a' + (5/3)b',
f' = (10a/9)a' + (4c/3)b',
g' = (4ca/9+10b/9+2E/3)a' + (10a/9)b'.
```

Each displayed antiderivative in equation (1) differentiates to exactly
this expression. Characteristic zero implies that its integration
constants E,F,G belong to C. No division by a, b, their derivatives, or
any other possibly vanishing polynomial occurs. In particular, the
stationary and degenerate cases have not been lost. The scalar terms
b,g are retained throughout.

## 2. The whole residual Jacobian and first integral

Substituting the integrated coefficients into the entire differentiated
pair gives exactly

```
J=(Ua'+Vb')t+Ub'−(aV/3)a',
U=−5a²/9+4cb/3+F,
V=4ca/3+10b/3+2E.
```

The polynomial

```
H=−5a³/27+4cab/3+Fa+5b²/3+2Eb
```

has partial derivatives H_a=U and H_b=V. Hence its actual x-derivative
is the linear-t coefficient of the whole Jacobian. A constant Jacobian
therefore implies H constant and residual value
`j=Ub'−(aV/3)a'`. Conversely, these identities are exact; no lower
coefficient or integration condition was omitted.

## 3. Exhaustive polynomial-degree argument

If exactly one of a,b is nonconstant, H has a unique highest-degree term:
`−5a³/27` when a is nonconstant, and `5b²/3` when b is nonconstant.
This rules out both cases. If both are constant, the residual Jacobian
is zero.

Suppose both are nonconstant, with degrees r,s>0. The inequality

```
r+s < max(3r,2s)
```

holds for every positive r,s: split into s<2r and s≥2r, as in the note.
The linear terms also have lower degree than this maximum. Thus the
only possible leading cancellation in H is between the cubic in a and
the square in b. Necessarily 3r=2s, so r=2k,s=3k with integer k≥1.
If A,B are the respective nonzero leading coefficients, then:

- U has degree 4k and leading coefficient `−5A²/9`;
- V has degree 3k and leading coefficient `10B/3`;
- the two summands of j both have degree 7k−1.

Their leading coefficients sum to

```
−(5A²/9)(3kB) − (A/3)(10B/3)(2kA)
= −35k A²B/9.
```

This is nonzero in characteristic zero, without any exception at the
relation `B²=A³/9` from H's leading cancellation. Therefore j has
degree 7k−1≥6 and cannot be a nonzero constant.

This covers all polynomial degree possibilities. It does not assume
that the first integral alone makes a,b constant.

## 4. Nonvacuous hostile family

The reviewed corrected family is

```
a=u², b=u³/3, c=E=F=G=0.
```

For arbitrary polynomial u(x), H vanishes identically and the complete
Jacobian is

```
J=−(35/27)u⁶u'.
```

In particular u=x^k gives the displayed `−(35k/27)x^(7k−1)`.
The original pairs still have monic degrees 3 and 5; their function
field is C(x)(t) because its extension degree divides both 3 and 5.
The free C(x)[p]-basis `1,q,q²` has degrees `0,5,10`, distinct
modulo three, and consequently exact codimension `(15−3)/3=4`.
Thus the hostile family really lies in the relevant generic curve
class. Constant H is not enough; the final Jacobian equation is essential.

The producer reported an earlier illustrative typo before this audit began
(`a=3u²,b=u³`). It is absent from both reviewed files and both checkers.
No correction to the reviewed theorem, proof, or formulas was required.

## 5. Reproduction and boundary

```
python3 verify_cubic_quintic.py > review/cubic_quintic_producer_replay.json
python3 -O verify_cubic_quintic.py > review/cubic_quintic_producer_replay_optimized.json
cmp review/cubic_quintic_producer_replay.json review/cubic_quintic_producer_replay_optimized.json
python3 review/independent_cubic_quintic.py > review/cubic_quintic_independent_results.json
python3 -O review/independent_cubic_quintic.py > review/cubic_quintic_independent_results_optimized.json
cmp review/cubic_quintic_independent_results.json review/cubic_quintic_independent_results_optimized.json
```

The separate classification audit remains valid after a status-only
promotion, verified by exact reconstruction of its reviewed hash in
`PROMOTION_RECONCILIATION.json`.

This audit supplies no exclusion of the terminal degree-(6,9) pattern,
no all-degree hyperelliptic theorem, and no resolution of JC(2).
