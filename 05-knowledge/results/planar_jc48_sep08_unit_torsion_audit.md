# Independent audit of the boundary-linear unit-torsion supplier

**Status: PASS — independent analytic/source audit and identical normal,
optimized, and frozen exact replays.** This audits
[the unit-torsion theorem](planar_jc48_sep08_unit_torsion.md) and
[its producer](../../04-computation/planar_jc48_sep08_unit_torsion.py).
It is a consequence for an actual fixed smooth polynomial and its
response quotient. It is not a new Jacobian-conjecture case.

## 1. Hypotheses and inheritance

I independently checked the family

    F=a rb+c+h(b),     a!=0,     h in b^3 C[b],
    b=-x^2(1+x^2t),   rb=-x(1+x^2t).

The source is the original `A2_(x,t)`, rather than all of `W`.
For `x!=0` the other chart applies. Its derivatives of `F` are
`ab` and `ar+h'(b)`; their only simultaneous zero is `r=b=0`,
which is outside the source. On `E={x=0}`, `F_x=-a`. Therefore
`F` has unit gradient on the entire affine source. It is not
declared smooth on the added boundary.

The rational field is exactly `C(F)(b)`, with
`r=(F-c-h)/(ab)`. The literal chart Jacobian gives
`D_F=(ab/r^2)partial_b|F`, whose constants are `C(F)` in
characteristic zero. This pays both hypotheses of
[the independently audited torsion/connection theorem](planar_jc48_sep06_torsion.md),
including the rational constants field; it does not assume a
geometric-integrality conclusion as an extra input.

The member `h=0` is already the `r=3` one-arm example
`P=x+lambda x^3t` in Section 5 of that inherited theorem, up to
the displayed scalar and target shift. The present extension is
the complete `h in b^3 C[b]` family with the same labelled
principal-part arm. I verified the canonical map and connection
against Sections 2.3 and 3 of the inherited proof, including signs.

## 2. The complete actual component tuple

Put `g=F-c`. For zero-constant polynomial integrals
`L=integral h/b^3 db` and `M=integral h^2/b^3 db`, the rational
primitive is

    G0=a^(-3)[-g^2/(2b^2)-2gL+M]
       =-1/(2a x^2)+Q(b,rb),

where `Q` is polynomial in global functions. The derivation
identity is `D_F G0=1`. Thus its only affine pole is the double
pole on `E`. The special fibre has the distinct components
`E` and `Gamma={1+x^2t=0}`. Smoothness ensures every component
is reduced and that distinct components are disjoint. Any extra
components of this fibre are retained and have zero principal
part, because `G0` is regular off `E`.

The absence of a simple scalar term is load-bearing and checks
directly for the whole family. Namely

    g=-a x-a x^3t+O(x^6),
    a^2/g^2-x^(-2) is regular at E.

Therefore the scalar principal part, in the actual uniformizer
`g`, is exactly `-a/(2g^2)` on `E` and zero on every other
component. In particular there is no hidden `g^(-1)` term.
The common-diagonal quotient cannot remove this tuple, since
the component `Gamma` has zero principal part.

It follows that `[1]` has exact annihilator `(F-c)^2`.
The upper bound also has an actual polynomial witness:
`g^2G0` has no affine poles and is polynomial by normality, with
`D_F(g^2G0)=g^2`. For the lower bound, multiplying the component
tuple by `g` leaves `-(a/2)[e_E]g^(-1)`, which is nonzero.
Thus generic exactness has not been confused with a polynomial
primitive, and the lower bound is not inferred only from the
pole of one particular choice of primitive.

## 3. Canonical derivatives and exact heights

For any polynomial Bezout field `V(F)=1`, the inherited operator
is `nabla[q]=[V(q)+(div V)q]`. It is independent of the choice
and differentiates the full scalar component tuple. Consequently

    nabla^j[1]  <->
       -(a/2)(-1)^j(j+1)! [e_E] g^(-j-2).

The coefficient is nonzero for every integer `j>=0`. This proves
the exact primary order `j+2`, not merely an upper bound.
Different highest pole orders prove linear independence over
`C`. Multiplication by `g` reaches the first level, and transverse
differentiation reaches every higher level of this same arm.
Thus the smallest connection-stable `C[F]` submodule containing
the unit is the complete arm described in the primary note.
In particular the intrinsic divergence class `nabla[1]` has
exact order three. No identification with a moving-source
Hamiltonian deformation has been smuggled into this statement.

## 4. Exact source and independent replays

The source uses explicit exceptions, with no optimization-disabled
assertions or imports of earlier producers. Its finite universe is
exactly `h=0,b^3,b^4,2b^3-b^5,b^3+b^6`, arbitrary symbolic nonzero
`a`, and derivative levels `0..12`. The source checks the rational
primitive in one chart and the actual polynomial killing witness
by the independent source-chart Jacobian. It also checks the full
scalar coefficient, vanishing simple coefficient, regular
remainder, actual `Gamma`, and the surviving pole after one
value factor. The unbounded result follows from the formulas
and inherited isomorphism, rather than this finite bank.

I independently ran, from the repository root,

```sh
python3 -B 04-computation/planar_jc48_sep08_unit_torsion.py > /tmp/unit_torsion_independent_normal.out
python3 -B -O 04-computation/planar_jc48_sep08_unit_torsion.py > /tmp/unit_torsion_independent_optimized.out
```

Both pass **86 always-active gates** and reproduce all 398 bytes
of [the frozen output](planar_jc48_sep08_unit_torsion.out).

| Artifact | SHA-256 |
| --- | --- |
| Source, 2,469 bytes | `f6f9f99d62216f78e366a30a1452a392ee046966ed914aa1cee88c55268d1ba2` |
| Output and both independent replays | `6aec5ecea6bebbf93d59dbc70b5bcc117db3b41cf7e91320343368eafcec4d2b` |
| Semantic record | `63e99c734749aaaf78f4e6ccf6f48112b62a49c960456d3d1194b2575864fa94` |

No mathematical correction is required. The whole stated family,
the exact component labels, and the unbounded derivative conclusion
are accepted with the preceding scope and inheritance.
