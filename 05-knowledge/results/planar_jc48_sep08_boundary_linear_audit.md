# Independent audit of the full DG boundary-linear exclusion

**Status: PASS — independent analytic and source audit; normal, optimized,
and frozen exact outputs agree.** This audits
[the boundary-linear proof](planar_jc48_sep08_boundary_linear.md) and its
[standalone producer](../../04-computation/planar_jc48_sep08_boundary_linear.py).
The claim excludes a polynomial source mate of unrestricted degree for
every global function of boundary-chart degree at most one in `r`.
It does not close the general quartic layer or JC(2).

## 1. Full class and rational reduction

I independently reconstructed the argument in both full charts. With
`b=-x^2(1+x^2t)` and `r=1/x`, the only negative source power in
`A(b)r+B(b)` is `A(0)/x`; all other terms are polynomials. Thus `b|A`
is necessary and sufficient for globality, and the stated family is
exhaustive. When `A=0`, both derivatives vanish on `E={x=0}`.

For `A!=0`, `C(r,b)=C(F)(b)` with `r=(F-B)/A`. The literal chart
Jacobian is `r^2`, so a source mate must satisfy

    partial_b G|F = lambda (F-B)^2/A^3.

The two derivatives with respect to the independent coordinate `F`
commute with `partial_b`; they produce a rational primitive of
`2lambda/A^3`. This is a necessary condition on any rational mate,
not a polynomial-degree truncation.

The rational-map argument is valid over an algebraic closure of
`C(F)`. If `deg A=M` and its distinct roots have multiplicities `m_i`,
the primitive has precisely the finite pole orders `3m_i-1`. It has
no additional finite poles or positive-degree polynomial part.
Subtracting its value at infinity gives a nonzero rational function
of degree `3M-s` with a zero of order `3M-1` there, where `s` is the
number of distinct roots. Hence `s<=1`. Together with `b|A`, this
forces `A=a b^m`. The source order of `r b^m` is `2m-1`, so `m>=2`
again gives a critical curve on `E`. The only surviving case is
`A=a b`, with `a!=0`. This recovers the pole mechanism of
[THM-2071 / quadratic-fiber-square-parity-gate, §5](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
with the correct cubic denominator and nonconstant chart Jacobian.

## 2. Complete integrability boundary and all repairs

For `g=F-B(0)` and `h=B-B(0)`, the only finite residue of
`lambda(g-h)^2/(a^3 b^3)` is
`(lambda/a^3)(h_1^2-2g h_2)`. Since `g` is transcendental, its
vanishing is exactly `h_1=h_2=0`. For every polynomial
`h in b^3 C[b]`, not merely the coefficients tested in the producer,
the polynomial integrals

    L = integral h/b^3 db,       M = integral h^2/b^3 db

give the complete particular primitive

    G0 = (lambda/a^3)(-g^2/(2b^2)-2gL+M).

Substitution `g=arb+h` independently gives
`G0=-lambda r^2/(2a)+Q(b,rb)` with polynomial `Q`.
For example the cross term `r h/b` equals `(rb)(h/b^2)` and is
global. Thus this primitive has an exact double pole along `E` and
is regular on the entire other chart. The constants of the nonzero
derivation `(ab/r^2)partial_b` on the pure transcendental field
`C(F)(b)` are exactly `C(F)`. Consequently **every** rational mate is
`G0+H(F)` for a rational one-variable function `H`; no additional
source-dependent repair has been omitted.

The actual fibre `F=B(0)` contains both `E` and
`Gamma={1+x^2t=0}` in the original affine plane. It is reduced along
both: the leading term at `E` is `-a x`, and at the generic point of
`Gamma`, `g=b(ar+h/b)` has a unit second factor. A value pole of
`H` of order two is necessary to remove the double pole of `G0`
on `E`. The same value pole then has order two on `Gamma`, where
`G0` is regular. This contradiction excludes every polynomial
source mate. It needs neither a classification of all components
of this fibre nor any assumption that the mate extends to `W`.

This is the decisive retained datum: two reduced components of the
same actual fibre carry incompatible principal parts. Rational
integrability by itself loses precisely that regularity information.

## 3. Hostile control and scope checks

I separately checked the literal example

    F=rb+b^4=b^4-v,
    H=b^2,                    F-H^2=-v,
    G0=-r^2/2-2rb^3-(4/3)b^6.

It has boundary Jacobian `r^2` and source Jacobian `1`, while its
source primitive has a double pole on `E`. The first coordinate and
its quadratic approximate root are global. The first coordinate has
source `t`-degree four and nonconstant boundary restriction `b^4`.
It has no critical point on the original source: its only critical
point in the other chart is `(r,b)=(0,0)` on the added boundary, and
`F_x=-1` on `E`. It is therefore **source-critical-free**, not
critical-free on `W`. The final producer label was repaired to state
this distinction; its mathematical checks were unchanged.

The unrestricted-degree monomial primitive also follows directly
from the integral formula. For `k>=3` and `F=rb+b^k`, it is

    -1/2 [r^2 + 2k/(k-2) r b^(k-1)
                 + k^2/((k-2)(k-1)) b^(2k-2)].

Thus the finite `k=3..9` bank checks an independently derived
all-parameter identity. The theorem covers arbitrary degrees of
`B,C` in `F=rb C(b)+B(b)` and arbitrary degrees of a proposed mate.
It neither covers all global square-prefix coordinates nor asserts
a result for boundary-chart degree at least two.

## 4. Exact source and replay acceptance

I read the entire source. Its checks are explicit exceptions, so
optimization cannot disable them. The stated finite universe is
`rb^m`, `m=1..9`; `B=b^m`, `m=0..9`; pole-degree arithmetic through
degree nine; symbolic residue and primitive coefficients through
degree five; value-pole orders one through five; monomial primitives
`k=3..9`; and the literal quartic control. The two chart Jacobians
are computed independently. The unbounded proof uses the identities
and valuation argument above, not extrapolation from those ranges.

From the repository root I independently ran:

```sh
python3 -B 04-computation/planar_jc48_sep08_boundary_linear.py > /tmp/boundary_linear_independent_normal.out
python3 -B -O 04-computation/planar_jc48_sep08_boundary_linear.py > /tmp/boundary_linear_independent_optimized.out
```

Both runs pass **149 always-active gates**. Their 525 output bytes
are identical to the current
[frozen output](planar_jc48_sep08_boundary_linear.out), including the
corrected source-critical-free label.

| Artifact | SHA-256 |
| --- | --- |
| Primary source, 5,103 bytes | `aa57669f5b20202adb7e9e39e10b2c7d600a32d188e1303963dbf68e0d068d53` |
| Frozen output and both independent replays | `30a014d4d7850aae83ce1b54c13b026236916e79edb0dd7d667168adb690b193` |
| Semantic record digest | `4e937ecc6e3519b5d405386b568a139bd1db7787f16ba8177a5a7ab369183648` |

No mathematical repair is required. The independent audit accepts
the full stated theorem and the finite controls with these scopes.
