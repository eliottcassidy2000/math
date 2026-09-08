# First-jet obstructions at shared quartic boundary roots

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED.**
These are bounded extensions of the proved common-root necessity.
They do not exclude all shared roots or settle JC(2).

## 1. Setup and three precise conclusions

Retain the actual surface, sections, and complete boundary from
[the quartic common-root theorem](planar_jc48_sep08_quartic_common_root.md):

    W=(P1_x x P1_z) minus S,       S={z=x^2},
    H=N/s^2 in L_2,              L=M/s in L_1,
    F=H^2+L,                    deg_t H=2.

The binary octic `N|S` is nonzero. At a shared boundary zero let
`s=0` define `S`, and let `w` be its local coordinate. Write

    N=A(w)+s B(w)+s^2 C(s,w),
    m=ord_w A,   j=ord_w B,   n=ord_w M(0,w),            (1)

allowing infinite orders. The original compact generic fibre and
relative differential are

    E=N^2+s^3 M-cs^4=0,       eta=omega/dF.

Away from the unique point at infinity on `S`, the relative form
is a unit times `s^2 dw/E_s`, or equivalently a unit times
`s^2 ds/E_w`. At infinity it has the additional canonical factor
`w^2`. The statements below retain that distinction.

**Local regularity.** At a shared zero, the generic relative form
is regular on every normalized branch if `m=1`, or if `N_s(0,0)!=0`.
The latter conclusion holds for every `m`, with no bound on `n`.

**Finite shared-root obstruction.** At a finite shared zero where
`dN(0,0)=0`, a rational mate is impossible if any of the following
three coefficients is nonzero:

    [w^2]N(0,w),       [w]N_s(0,w),       [w]M(0,w).     (2)

Equivalently, a rational mate requires at every finite shared point
either `m=1`, or `N_s(0,0)!=0`, or all of

    m>=3,             j>=2,              n>=2.          (3)

This is a necessary local condition. No mate is inferred when it holds.

**Global small-multiplicity exclusion.** If every zero of the
complete binary octic `N|S` has multiplicity at most two, then `F`
has no rational mate in `C(x,t)`, for arbitrary `M`. Unlike the
M-avoidance proof, this conclusion does not require `F|D` to be
nonconstant.

The closest proved mechanism is the generic Newton/Morse analysis in
the common-root theorem. The new operation retains the first normal
jet `B` and the tangential jet of `M` at the same actual boundary
point. The source is the common point and these jets; the target is
the normalized branches of the original fibre; the map is weighted
rescaling of its complete equation. It preserves the actual relative
form and exactness. Recording only the pair `(m,n)` discards the
normal derivative and can confuse regular and logarithmic behaviour.

The live concepts are common multiplicities, the first normal jet,
generic tangent quartics, complete infinity charts, and compact
exactness. The nonzero-M rational-mate example in Section 4 is the
hostile to any blanket exclusion of shared roots. The corrected near
miss is to transplant the finite logarithmic obstruction to infinity.

## 2. Regularity without M-avoidance

Suppose first that `m=1`. Use `v=N(s,w)` as the coordinate along
the boundary; this is an invertible local coordinate change because
`N_w(0,0)!=0`. At a shared point, `M(0,0)=0`. In coordinates
`(s,v)`, the generic leading equation is

    v^2+(alpha-c)s^4 + higher weighted terms = 0,

where `wt(s)=1`, `wt(v)=2`; terms such as `s^3 v` have higher
weight. For generic `c`, `alpha-c!=0`. There are two smooth
branches `v=(unit)s^2`, and `E_v` has order two in `s`.
Consequently `s^2 ds/E_v` is regular. The coordinate changes have
unit Jacobian. At infinity the extra canonical zero cannot spoil
regularity.

Now assume `N_s(0,0)!=0`, with arbitrary finite `m`. The equation
`N=0` gives an analytic centre `s=psi(w)`, with
`ord_w psi=m`. At that centre put

    ell=ord_w(M(psi(w),w)-c psi(w)).

For generic `c`, `1<=ell<=m`: if the first summand has smaller
order it remains; otherwise the nonzero leading coefficient of
`psi` makes the order at most `m`, with at most one exceptional
value of `c`. The two Puiseux determinations of the fibre have

    ord_w s=m,       ord_w N=(3m+ell)/2.

The displacement from `psi` has order greater than `m`, so it does
not change the displayed leading order of `M-cs`. The term
`2N N_s` strictly dominates `E_s`; the other terms have orders at
least `2m+ell` and `3m`. Therefore

    eta=(unit) w^((m-ell)/2) dw.

This is regular on every actual normalization, including when the
two Puiseux determinations form one ramified branch. The local
Weierstrass degree in `s` is two, so these exhaust all branches.

This first-normal-derivative condition is load-bearing. For example
`N=z=s+x^2`, `M=x` has a shared zero with `(m,n)=(2,1)` but
`N_s=1`, and is in this regular case. It is not in the logarithmic
class merely because the two boundary forms share a point.

## 3. The generic tangent quartic and its actual logarithmic branches

Suppose `dN=0` at a shared point. In local coordinates write its
quadratic part and the linear part of `M` as

    N2=a w^2+b s w+e s^2,       M1=d w+f s.

Substitute `s=wZ`. The degree-four part of the original equation is

    E=w^4[P_c(Z)+higher powers of w],
    P_c(Z)=(a+bZ+eZ^2)^2+dZ^3+(f-c)Z^4.                (4)

If `a!=0`, this polynomial has four simple nonzero roots for generic
`c`. Indeed a repeated nonzero root would satisfy the independent
polynomial equation `Z P0'(Z)-4P0(Z)=0`, where
`P0=P_c+cZ^4`. Its constant term is `-4a^2`, so it has only
finitely many roots, and each determines at most one exceptional
value `c=P0(Z)/Z^4`. The constant term `a^2` excludes zero roots.

If `a=0` and `b!=0`, factor `Z^2`. The remaining quadratic has
constant term `b^2` and has two simple nonzero roots generically by
the same argument. If `a=b=0` and `d!=0`, factor `Z^3`; the
remaining linear factor has one simple nonzero root generically.
These alternatives do not require classifying any branches with
`Z=0`: one nonzero simple root is sufficient for the obstruction.

At any such root `Z0`, the implicit function theorem gives an
actual smooth branch `s=w Z(w)`, with `Z(0)=Z0!=0`. On it

    E_s=w^3(P_c'(Z0)+O(w)),
    s^2 dw/E_s=(Z0^2/P_c'(Z0)) dw/w + regular terms.     (5)

The residue is nonzero. At a finite boundary point the canonical
multiplier is a unit, so the actual form also has a nonzero residue.
A differential of a meromorphic function cannot have such a simple
pole. Thus no rational mate exists whenever `(a,b,d)!=(0,0,0)`.
This proves (2)--(3). At infinity the extra `w^2` makes these
branches regular instead; no finite-point exclusion is asserted
there.

For the global small-multiplicity statement, a root with `m=1`
is covered by Section 2 or the proved M-unit result. At `m=2`,
a shared point with `N_s!=0` is regular; otherwise `a!=0` in (4),
so all four tangent roots are simple and nonzero generically.
The original equation has Weierstrass degree four, hence those
branches exhaust it. The only possible finite poles are simple;
at infinity the canonical factor makes them regular. Nonshared
points are regular by the prior theorem.

On every compact normalized generic component, a rational mate
would therefore have no poles: a pole of its differential of order
at most one cannot be the derivative of a pole. It would be a
holomorphic function, hence constant, whereas `eta` is nonzero.
Every component meets `W` away from `D`, because neither `S` nor
`Dbar` is a generic fibre component. This proves the global
statement without assuming irreducibility or using a finite point
of intersection with `D`.

## 4. An actual sharp shared-root rational mate and a stopping stratum

Let

    h=x^2+x^4t=-b,   H=h^2,   L=h,   F=h^4+h,
    G=1/[3x^3(4h^3+1)].

All of `h,H,L,F` are global on `W`; `G` is rational. Directly
`J(h,x^(-3))=3`, hence `J(F,G)=1`. The actual numerators are

    N=x^4z^2,   M=x^2z,
    N|S=x^8,   M|S=x^4.

They share only the finite point zero, with `(m,j,n)=(8,6,4)`;
both are nonzero at infinity. This satisfies (3). Moreover
`F|D=b^4-b` is nonconstant and

    G=r^3/[3(1-4b^3)].

On each generic compact component `h=constant`, the primitive has
pole degree three and local degree three at `D`. Thus the prior
degree obstruction has the exact correct threshold. This example
refutes removal of the common-root condition even if `M!=0` and
`F|D` is nonconstant. It is not a polynomial or global regular mate,
and it is not a Keller counterexample.

A separate stopping object explains why the global result above
requires **all** octic roots to have small multiplicity. With
nonzero `beta` and `alpha=-4beta^3/27`, take

    N=beta x^2z-beta x^4+alpha x^2z^2,    M=-1.

There is a finite M-unit multiplicity-six point, while the shared
infinity point has multiplicity two. Here `H|D=-beta`, `L|D=0`,
so the finite-D ramification supplier is absent. Its complete
infinity numerators are

    Ninf=alpha r^2+beta r^2q-beta q^2,   Minf=r^2q.

No rational mate or impossibility is inferred for this object.
It prevents silently extending M-avoidance to arbitrary shared
low-multiplicity points while reusing a now-constant `F|D`.

The remaining precise jet class at a finite singular common point
is `m>=3,j>=2,n>=2`. Its common Newton faces and higher principal
parts remain OPEN here. The rational control shows this remaining
class is not empty, and that a blanket rational-mate exclusion
would be false.

## 5. Frozen exact controls

Reproduce from the repository root with:

```sh
python3 -B 04-computation/planar_jc48_sep08_shared_roots.py
python3 -B -O 04-computation/planar_jc48_sep08_shared_roots.py
```

Both executions pass **454 always-active gates**, reproducing all
418 bytes of [the frozen output](planar_jc48_sep08_shared_roots.out).
The [source](../../04-computation/planar_jc48_sep08_shared_roots.py)
checks normal-unit representatives `m=1..8`, `n=1..9` and infinite
order, the three generic tangent regimes and the first unsupported
face, actual global controls with different normal derivatives,
both full charts of the rational mate, and the constant-D stopping
object. The arbitrary-jet and generic-normalization statements are
the analytic proof above, not an extrapolation from these ranges.

| Artifact | SHA-256 |
| --- | --- |
| Source, 6,255 bytes | `695b0da5d8fb5dec95b85d0d32086432ac6aaa8cb1fd8bc00fa2eef27f4da2ad` |
| Frozen output and both replays | `d156179f74f97998556c7c60f6c2f27166543abe36321cf229a5a38c6a24d085` |
| Semantic record | `5819ccc48388d0db4be9d8a86782c2d2362232652cbf273ffea80111a1a665b8` |

The [independent analytic/source audit](planar_jc48_sep08_shared_roots_audit.md)
accepts all local and componentwise global arguments, the actual
rational-mate and stopping controls, and normal/optimized/frozen
replay. Its SHA-256 is
`92932658bf46f1fa0d441aa9af66664d12cd290d6449e2192b8ef771bfe0a195`.
Root also independently read and accepted the full primary proof.
No source or output changed during promotion.
