# The DG boundary forces primitive pole degree at least three

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED;
standalone controls FINITE-EXACT.** The general gate below is an application
of the classical degree of a meromorphic function. The exact quartic
application excludes a rational Jacobian mate despite residue-free
relative poles. It is not a new polynomial Keller family: the displayed
polynomial already has an affine critical point. JC(2) remains OPEN.

## 1. Inheritance and the precise general gate

Use the specified DG surface and its complete charts:

    X=P1_x x P1_z,  S={z=x^2},  W=X\S,
    U0=A2_(x,t),  Uinf=A2_(r,b),
    x=1/r,  t=-r^2-r^4 b,
    D={r=0} in W,  omega=dx wedge dt=r^2 dr wedge db.

These are the [proved DG filtration and charts](planar_jc48_sep08_dg_quadratic.md).
Let `F in O(W)` have nonconstant restriction `F_D(b)=F(0,b)` to `D`.
For a generic complex value `c`, let `C` be any connected component
of the normalization of the compact closure of `F=c`. Let `eta_C`
be the relative differential determined by

    dF wedge eta = omega.

Define its boundary pole budget by

    B_C=sum_(P in C over S) max(0,-ord_P(eta_C)-1).       (1)

**General gate.** If a rational function `G in C(x,t)` satisfies
`J_(x,t)(F,G)=lambda`, `lambda in C*`, then every `eta_C` has zero
residues. On every component `C` meeting `D`, the restriction of `G`
is a nonconstant meromorphic function of degree exactly `B_C`, and

    B_C>=3.                                             (2)

Consequently a component meeting `D` with `B_C<=2` rules out any
rational Jacobian mate, even if all residues vanish. The statement
is necessary, not sufficient: neither zero residues nor (2) supplies
the missing global primitive.

The degree argument is inherited rather than new. The exact predecessor
[THM-2071 / quadratic-fiber-square-parity-gate, Section 5](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
compares the multiplicity of a zero of a rational primitive with its
pole degree. The [boundary-linear proof, Section 2](planar_jc48_sep08_boundary_linear.md)
uses the same principle for a primitive of `1/A^3`. The new input here
is the actual order-two zero forced by `omega` on the added DG boundary,
with the normalization and all possible poles retained. The
[compact quartic differential criterion](planar_jc48_sep08_quartic_differential.md)
has already separated holomorphic forms, logarithmic poles, and
residue-free poles of higher order. Its last example is the input
control below; no pinned predecessor is changed.

| Live concept | Preserved data and decisive sidecar |
|---|---|
| Added boundary | `F_D` retains actual transverse points with two-form zero order two |
| Compact normalization | Keeps all branches and component labels, even for reducible fibres |
| Pole degree | Restriction of a rational mate preserves exact local pole orders |
| Residues | Necessary for exactness, but a second-kind pole survives this test |
| Affine criticality | Excludes polynomial mates; it does not exclude rational mates |

The corrected near miss is to stop at residue vanishing, or to use an
affine critical point to rule out a rational mate. The least-used
sidecar is a finite point of `D`, which supplies a zero of the proposed
primitive of multiplicity three after subtracting its local value.
The source-to-target map is restriction of the rational identity to
each compact normalized generic component. It preserves exactness
and the local orders. Omitting either `D` or a boundary branch loses
the degree comparison.

## 2. Proof, including rational mates and reducible generic fibres

Choose `c` outside the finite critical-value set of `F on W` and
outside the finite critical-value set of the nonconstant polynomial
`F_D`. The former set is finite because each irreducible component
of the algebraic critical locus has constant `F` in characteristic
zero. Thus the fibre is smooth on `W`, and every one of its points
on `D` satisfies `F_b!=0`.

If a putative rational `G` is being considered, also avoid the finitely
many values whose fibre contains a vertical component of its pole
divisor. Then `G` restricts to a meromorphic function on every compact
normalized component. Horizontal poles cause no problem for defining
the restriction; the differential equation will rule out their poles
on the smooth part of that component. No assumption that `G` is global
on `W` is used.

At `P=(0,b0)` on `D`, the implicit function theorem makes `r` a local
parameter on the fibre. The actual differential identity gives

    eta_C=-r^2/F_b(r,b(r)) dr.                           (3)

Its order is exactly two. It is not the zero differential, and there
is no normalization ramification at this smooth transverse point.
If `dG=lambda eta_C`, then

    G-G(P)=-lambda/(3 F_b(P)) r^3+O(r^4).               (4)

In particular the local degree of `G:C -> P1` at `P` is exactly three.

At every point of the fibre in `W`, `eta_C` is regular: the fibre is
smooth and `omega` is a regular two-form. This holds independently of
whether `G` is regular on the ambient surface. The meromorphic identity
`dG=lambda eta_C` prohibits a pole of the restricted `G` at any such
point. Indeed a pole of order `m>0` gives a differential pole of order
`m+1` in characteristic zero. All its possible poles are therefore
over `S`.

The same Laurent expansion gives zero residues for `dG`, and at a
pole of `eta_C` of order `k>=2` it forces a pole of `G` of order
`k-1`. A simple pole would have nonzero residue and be impossible.
Thus the pole divisor of `G|C` has degree exactly (1). For a
nonconstant meromorphic function on a compact connected curve this
is the degree of its map to `P1`. Every fibre divisor of that map
has the same degree, so the multiplicity three in (4) cannot exceed
`B_C`. This proves (2).

The argument applies separately to each component meeting `D`.
Such a component exists, because the nonconstant polynomial `F_D`
takes every generic value. It does not assume geometric integrality,
that all components meet `D`, or that different points of `D` have
the same `G` value. In particular we do not sum local degrees at
different `D` points without that last additional information.

The numerical bound is sharp as a compact-curve principle:
`G(u)=u^3` on `P1` has differential zero of order two at zero and
one pole of order four, so its primitive pole degree is three.
This is not an assertion that degree three is attained by a DG pair.

Two controls show why the hypotheses matter. First `F=t`, `G=-x`
has Jacobian one, and `F` is global, but `F|D=0` is constant.
Second, the inherited global `F=b^4+rb` has the rational mate

    G0=-r^2/2-2rb^3-(4/3)b^6.

On `F=c`, its restriction is

    G0=-c^2/(2b^2)-cb^2+b^6/6,
    dG0=(c-b^4)^2/b^3 db.

For generic `c` its pole degree is eight, consistent with (2), and
the four transverse `D` points give differential zeros of order two.
The gate does not falsely exclude this known rational mate.

## 3. The exact residue-free quartic and its complete boundary

Take the existing sextuple-root control

    H=x^2t-(4/27)x^6t^2+(4/27)x^2,
    L=-(1+x^4)t-x^2,
    F=H^2+L.                                           (5)

It corresponds to the honest global sections

    N=x^2z-x^4-(8/27)x^4z+(4/27)x^2z^2,
    M=-(1+x^2z),
    H=N/(z-x^2)^2,  L=M/(z-x^2).

The complete second-chart formulas are

    H=-1-r^2b-(8/27)b-(4/27)r^2b^2,
    L=b+r^2+r^4b,
    F|D=(1+8b/27)^2+b.                                (6)

The discriminant of `F|D-c` is `(256c+1593)/729`, so there are two
distinct transverse points of `D` for generic `c`.

On the complete projective boundary the zero divisor of `N|S`
is six times zero plus twice infinity. The form `M|S=-(1+x^4)`
is nonzero at both points. Consequently these are the only points
where the compact fibre meets `S`. We now exhaust the normalized
branches there, rather than retaining just the two known poles.

### The three finite branches

Put `y=x^2` and `s=z-x^2`. The literal compact equation is

    Ebar(s,y)=(sy-4y^3/27+4s^2y/27)^2
                -s^3(1+y^2+sy)-cs^4.                  (7)

Since `Ebar(s,0)=-s^3(1+cs)`, its local Weierstrass degree in `s`
is exactly three. After `s=y^2 h`, the initial divided equation is

    (Ebar/y^6)|_(y=0)=-(h-4/9)^2(h-1/9).               (8)

At `h=1/9` the derivative is `-1/9`, so one branch is
`s=y^2h(y)`, `h(0)=1/9`. Its original derivative satisfies

    Ebar_s=-(1/9)y^4+O(y^5).

For the orientation `dF wedge eta=omega`, the relative form in this
finite chart is `eta=s^2/Ebar_s dx`. It therefore has leading
value `-(1/9)dx` on this third branch and is regular.

The remaining two branches have

    s=(4/9)y^2+e y^3+O(y^4),
    e^2=-64(65+36c)/19683.                             (9)

For `c!=-65/36` the two values of `e` are nonzero. After substitution
`s=y^2(4/9+yh)` and division by `y^8`, the initial equation and
its derivative are

    -h^2/3-64(65+36c)/59049=0,   derivative=-2h/3.

Thus both branches are genuine convergent power series in `y`, by
the implicit function theorem. They are distinct from each other
and from the first branch. Together with Weierstrass degree three,
this proves there is no missing ramified branch. In the original
curve each is parametrized by `x`, with `y=x^2`.

On either branch in (9),

    Ebar_s=-(2e/3)y^5+O(y^6),
    eta=-8/(27e) x^(-2) dx+O(1) dx.                    (10)

Its coefficient is an even Laurent series in `x` to all orders.
Hence these are exactly two double poles, each with zero residue.
The overall sign is opposite to the form `s^2 ds/E_x` displayed
up to sign in the predecessor; the orientation is fixed explicitly
here. The orders, residues, and obstruction are unchanged.

### The two branches at infinity

Use `r=1/x`, `q=1/z`, and `sigma=q-r^2`. The complete local sections
and two-form are

    Ninf=-r^2 sigma-sigma^2-(4/27)r^2-(8/27)sigma,
    Minf=1+r^4+r^2 sigma,
    omega=r^2 sigma^(-2) d sigma wedge dr.

Their fibre equation is `Einf=Ninf^2+sigma^3 Minf-c sigma^4`.
At `r=0` it has order exactly two in `sigma`. Its two branches are

    sigma=-r^2/2+a r^3+O(r^4),   a^2=729/512.           (11)

Indeed substitution and division by `r^6` gives
`(64/729)a^2-1/8`; its derivative is `128a/729!=0`.
Weierstrass degree two and the two distinct implicit roots exhaust
all infinity branches. On them

    Einf_sigma=(128a/729)r^3+O(r^4),
    eta=r^2 sigma^2/Einf_sigma dr
       =[729/(512a)]r^3 dr+O(r^4)dr.                   (12)

Both differential orders are exactly three, so there is no infinity
pole. The extra `r^2` numerator is retained throughout.

## 4. Rational nonintegrability and the critical-point distinction

Choose `c` generic for `F`, and in particular avoid `-65/36` and
`-1593/256`. Across all compact normalized components the boundary
differential orders are exactly

    (-2,-2,0,3,3).

On `W` the differential has no poles. Every residue is zero, but
the total primitive pole-degree budget is only

    sum max(0,-ord(eta)-1)=1+1=2.                       (13)

Each individual component has budget at most two. At least one,
indeed every component containing a point from (6), meets `D` and
would require primitive degree at least three by the general gate.
Thus `eta` is not a rational exact differential on any such
component. Consequently the literal polynomial (5) has **no rational
Jacobian mate in C(x,t)**. No irreducibility proof or assignment of
the two poles to a common global component is needed.

The same polynomial has affine critical points at

    x^2=54/97,   t=7081/1296,
    H=85/36,     F=-77/36.

Both `F_x` and `F_t` vanish there. This elementary observation already
rules out a polynomial mate and must not be advertised as a new
Keller-family exclusion. It does not imply the rational conclusion:
`F=x^2`, `G=t/(2x)` has Jacobian one despite the critical line of
`F`. The degree calculation pays actual rational nonintegrability
that criticality and zero-residue tests do not pay.

For a second hostile to omission of the order-two boundary zero,
`G=u+1/u` on `P1` has degree two and derivative
`(1-u^(-2))du`. It has two residue-free double poles, but only simple
zeros. Thus the same pole pattern alone is compatible with exactness;
the actual `D` zero in (3) is essential.

There is also a direct rational strengthening of the frozen safe
quartic theorem. Under its hypotheses—`N|S` nonzero with all zero
multiplicities at most two, and `M` a unit at those zeros—the same
proof makes `eta` a nonzero holomorphic differential on every compact
generic component. A rational mate restricts meromorphically there
after avoiding its finitely many vertical pole values; holomorphy
forces its restriction to have no poles and hence to be constant.
Therefore that class has no **rational** mate either. This corollary
uses the proved local regularity and the explicit rational-restriction
argument above; it does not modify the earlier frozen statement.

## 5. Reproduction and exact stopping point

The standalone source imports no inherited mathematical implementation.
It checks the literal quartic and complete charts; the exact three
finite and two infinity branches, including their implicit derivatives;
the total pole budget; the boundary discriminant; and the actual affine
critical point. Its positive controls are the sharp classical degree
three example, the exact two-pole degree two example, the constant-D
exception, and the inherited actual DG rational mate of degree eight.
These controls validate the displayed identities, not the general
meromorphic-function degree principle by enumeration.

    python3 -B 04-computation/planar_jc48_sep08_pole_degree.py
    python3 -B -O 04-computation/planar_jc48_sep08_pole_degree.py

Both runs pass **68 always-active gates**, with byte-identical **634-byte**
output. The source and output are frozen for independent audit; their
provisional status marker is intentionally retained.

| Artifact | SHA256 |
|---|---|
| Source, 8,484 bytes | `58743ec60c599028800c560482dae7a7928bd7b97a9818e531e94d47c5bab151` |
| Frozen output and optimized replay | `3ffdc88c37dac410055cff384b84c7c516dbd044201656ecc3a4fb9932ef13bf` |
| Semantic record | `d70d4a1bbd3accf0ae003d9596b661ce6818ed657dee729d3ca49cdf694058ff` |

The remaining general obstruction is exactness when all residues vanish
and each component meeting `D` has pole budget at least three. This
gate does not settle that case or construct a primitive. The finite
component data, not a pooled count alone, must be retained in a future
extension.

The [independent analytic/source audit](planar_jc48_sep08_pole_degree_audit.md)
passes. The frozen source/output retain their pre-audit RESERVED labels;
the promoted primary and this audit govern the proved status.
