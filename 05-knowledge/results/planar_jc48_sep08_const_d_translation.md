# Moving the finite-six point by an exact intersection of global section spaces

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED.**
This extends the independently audited
[fixed-zero constant-D theorem](planar_jc48_sep08_const_d_quartic.md)
to every finite location of its sixfold boundary zero. The map used
is a symplectic translation of the source plane. It is not an
automorphism of the compactified surface.

## 1. Statement, inheritance, and the precise missing step

Keep the original surface and chart

    W=(P1_x x P1_z) minus {z=x^2},
    s=z-x^2,  t=1/s,  x=1/r,  t=-r^2-r^4*b_D,
    D={r=0} in W,     omega=dx wedge dt.

Let `H=N/s^2` and `L=M/s` be global functions in the proved
filtration spaces `L_2,L_1`, and set `F=H^2+L`.

**Theorem.** Suppose, for a finite `p in C`,

    N|S=a*(x-p)^6,       a!=0,
    M(p,p^2)!=0,         F|D is constant.

Then there is no `G in C(x,t)` with `J_(x,t)(F,G)` a nonzero
constant. The binary octic notation means exactly multiplicity six
at `(p,p^2)` and multiplicity two at the distinguished infinity
point, with no other octic zeros. It does not include an eightfold
zero or allow the nonunit case `M(p,p^2)=0`.

The closest proved mechanism is the fixed-zero theorem: complete
local branch analysis reduces that family to an elliptic curve,
whose full allowed primitive space fails an exact coefficient
identity. The corrected near miss was to move a finite point to
zero without paying globality. The canonical hostile is the pole
in Section 4 below. The underused operation is the exact
intersection of two global function spaces inside one rational
source field.

The live concepts are section-box completeness, the source-plane
symplectic map, changes of boundary values, the auxiliary surface,
and rational exactness. The necessary missing coordinate is the
whole translated numerator, including its lower normal jets. A
boundary multiplicity or a formal local coordinate alone is
insufficient.

## 2. The complete global families at an arbitrary finite point

Set `w=x-p`. The original graph is now `z=(w+p)^2`; importantly,
it is not silently replaced by `z=w^2`. Translation of the first
projective coordinate preserves its degree bound, so the complete
global numerator boxes are still

    deg_w N(w,z)<=4,  deg_z N(w,z)<=2,
    deg_w M(w,z)<=2,  deg_z M(w,z)<=1.

First impose only `N|S=a*w^6`. Write

    N=a*w^6+s*B(w)+s^2*C(w),     s=z-(w+p)^2.

Comparison of the coefficients above degree four in the original
`z` box gives exactly

    C=c0+c1*w+c2*w^2,
    B=beta0+beta1*w+b*w^2
                    +[c1+2p(c2-a)]*w^3+(a+c2)*w^4.     (1)

For completeness, an arbitrary section initially has `deg C<=4`
and `deg B<=6`. The degree-eight and degree-seven boundary
equations, together with the `z`-coefficient bound, successively
force `C4=B6=0` and `C3=B5=0`. The degree-six and degree-five
equations are `B4=a+C2` and `B3=C1+2p(C2-a)`. All remaining
coefficients are free. Thus (1) is exhaustive for every `p`,
including values where a generic-rank argument would need care.

The actual chart `w=1/r-p`, `t=-r^2-r^4*b_D` gives the
coefficient of `b_D` in `H|D` as `a-c2`. Similarly, the complete
numerator of a global `L` before imposing constancy is

    M=m0+m1*w+d*w^2+(n1+2p*n2)*w^3+n2*w^4
                                      +s(n0+n1*w+n2*w^2),

and the coefficient of `b_D` in `L|D` is `-n2`.

The filtration bounds give `deg H|D<=2` and `deg L|D<=1`.
Consequently constant `H|D^2+L|D` forces both restrictions to be
constant. In (1) and the displayed `M` this means

    c2=a,              n2=0.

The entire resulting family is therefore

    N=a*w^6+s[beta0+beta1*w+b*w^2+c1*w^3+2a*w^4]
                              +s^2[c0+c1*w+a*w^2],
    M=m0+m1*w+d*w^2+n1*w^3+s(n0+n1*w).                 (2)

All occurrences of `p` have disappeared from these expressions
in `(w,s)`. This cancellation is the new supplier. It is stronger
than agreement of the boundary octic alone and keeps every
coefficient of the original functions.

Conversely, for every `p`, substitution `s=z-(w+p)^2` in (2)
lies in the original global section boxes and has constant
boundary restrictions. Those actual constants are

    H|D=c0-b+2p*c1+4a*p^2,
    L|D=n0-d+2p*n1.                                  (3)

The exact source also pays all-parameter completeness with
constant nonzero minors: nine independent constraints on the
fifteen-dimensional `N` box leave the displayed six-dimensional
space, and the constant-D `M` space has dimension five. Both
basis minors are one for every complex `p`; the constraint minor
is a constant unit. No finite list of locations is used to
infer this statement.

## 3. Source-plane transport into the auxiliary excluded family

Substitute `s=1/t` in (2). The translated source functions are

    H0(w,t)=a*w^6*t^2
               +(beta0+beta1*w+b*w^2+c1*w^3+2a*w^4)*t
               +c0+c1*w+a*w^2,
    L0(w,t)=(m0+m1*w+d*w^2+n1*w^3)*t+n0+n1*w.

They satisfy the exact identities

    H(w+p,t)=H0(w,t),      L(w+p,t)=L0(w,t).           (4)

These same polynomial expressions are global on the *auxiliary*
surface

    W0=(P1_w x P1_zeta) minus {zeta=w^2},
    zeta=w^2+1/t.

On that surface their boundary octic is `a*w^6`, their finite
M-value is `m0=M(p,p^2)!=0`, and their boundary restrictions are
`c0-b` and `n0-d`. Thus `F0=H0^2+L0` belongs to the complete
fixed-zero family. It has the same nonzero `a,m0` required by
that theorem, including all lower-jet and elliptic subcases.

The source-plane isomorphism is

    T_p:(w,t) -> (x=w+p,t),       determinant dT_p=1.

For an arbitrary rational `G(x,t)` let `G0(w,t)=G(w+p,t)`.
The chain rule in the rational function field gives

    J_(w,t)(F0,G0)=J_(x,t)(F,G)(w+p,t).               (5)

A nonzero constant Jacobian would therefore give the forbidden
rational mate in the proved fixed-zero theorem. This proves the
candidate extension.

The source of the connection is the full original global family.
The target is the full auxiliary global family. The map preserves
the rational field, the original polynomial expression `F`, its
fibre parameter, and the volume form. It changes the boundary
constants (3) and does not preserve the original compactification
as a global morphism. The sidecar supplying what is needed is
the complete section-space calculation, not a claimed surface
symmetry. The target rational mate need not be global on either
surface, exactly as allowed by the fixed-zero theorem.

## 4. Exact hostiles and the failure boundary

The plane map is not generally a surface automorphism. The
original global boundary coordinate is `b_D=-x^2-x^4*t`.
Pulling it back by `x -> x-p`, then using the original boundary
chart, gives

    -(x-p)^2-(x-p)^4*t
      =-2p/r+5p^2-4p^3*r+p^4*r^2+(1-pr)^4*b_D.       (6)

For `p!=0` its simple pole has coefficient `-2p`. The source
checks this coefficient symbolically, and also at nonzero real
and nonreal points. It prevents replacing the family-intersection
proof with an invalid assertion about all global functions.

The inherited actual rational-mate hostile remains

    h=x^2+x^4*t,     F=h^4+h,
    G=1/[3*x^3*(4*h^3+1)],     J(F,G)=1.

Its octic has multiplicity eight, not the six-plus-two pattern.
The new source independently rechecks this identity. Thus the
extension changes the location of a specified sixfold point;
it does not erase the multiplicity or M-unit hypotheses.

Named translated members include the residue-free elliptic tuning
`a=-4/27,b=1,m0=-1,d=8/9`, all other coefficients zero, at
`p=0,1,-3/2,i`. Lower normal jets are also retained and checked.
These examples test (2)--(5); they do not replace the symbolic
completeness proof. The remaining boundaries include nonunit M,
other octic multiplicity partitions, and the infinity-eight case.

## 5. Reproduction and audit boundary

```sh
python3 -B 04-computation/planar_jc48_sep08_const_d_translation.py
python3 -B -O 04-computation/planar_jc48_sep08_const_d_translation.py
```

The standalone source imports no inherited mathematical
implementation. It verifies the general preconstant families,
the complete all-point constraint and basis minors, both actual
boundary constants, the literal rational Jacobian transport, the
nonextendability pole, and named inside- and outside-family
controls. The all-parameter no-mate conclusion is the analytic
consumer of the fixed-zero theorem, not an empirical inference.

Normal and optimized runs agree byte for byte with the frozen
output: **51 always-active exact gates, 373 output bytes**.

- [Source](../../04-computation/planar_jc48_sep08_const_d_translation.py):
  7,273 bytes, SHA256
  `69e58278cb6b84ae791c004b22d59f2b9aebd0e191bb6c08f39d422bcdb53dce`.
- [Frozen output](planar_jc48_sep08_const_d_translation.out): SHA256
  `62a3a4180658aa0fcdc9f0d0b7723fff387c5a87625be0f6f4bea33f50b0a0f5`.
- Semantic digest:
  `396216fbb707e9ba5633c9bbacc67cf1f4e58b2c14a5d4876ddf5780909d1950`.

The source and output are frozen. Independent analytic/source
review is pending; the primary remains RESERVED until that review
is complete. The fixed-zero theorem is PROVED and independently
audited.

The [independent audit](planar_jc48_sep08_const_d_translation_audit.md)
passes the complete all-point section calculation, rational transport
and both51-gate frozen replays. Root independently checked the full
proof, source and both modes as well. Audit7340 bytes, SHA-256
`7eb12d1c99a59f38cd6a81b98aafc0fe8816c9d8159a4463821df19cfcd0a2e3`.
