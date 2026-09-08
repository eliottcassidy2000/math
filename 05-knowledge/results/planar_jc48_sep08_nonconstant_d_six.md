# Nonconstant-D finite-six: a polynomial weight gate and an exact trace obstruction

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a complete coefficient-stratum argument on the fixed DG surface.
It is not a reduction of arbitrary Keller maps to that surface. JC(2)
remains open.

## 1. Statement, inheritance, and the preserved distinctions

Work over `C` on

    W=(P1_x x P1_z)\{z=x^2}, t=1/(z-x^2),
    x=1/r, t=-r^2-r^4 b_D, omega=dx wedge dt=r^2 dr wedge db_D.

Let `H in L2`, `L in L1` be actual global sections, and put
`F=H^2+L`. Assume the leading binary octic of `H` has a finite
root of multiplicity six and a double root at the original
infinity, and that `F|D` is nonconstant.

**Theorem.** No polynomial `G in C[x,t]` of any degree
satisfies `J(F,G)=1`. After the necessary active first jets, the
subfamily described by `c=nu=0` below has no rational mate either.
No blanket rational exclusion is asserted for the other coefficients.

The [complete constant-D shared-six theorem](planar_jc48_sep08_constant_d_shared_six.md)
and [the translated constant-D M-unit theorem](planar_jc48_sep08_const_d_translation.md)
already exclude polynomial mates for the complementary constant-D
class at every finite location. If the double point is finite,
[leading exactness](planar_jc48_sep08_leading_exactness.md) excludes
rational mates by its nonzero double-root residue. Thus acceptance
of this candidate combines with those proved inputs to close the
polynomial-mate question for the complete binary `(6,2)` class at
all placements. It does not exclude the known constant-D rational
mates, including genuine nonconstant-lower-row submersions.

The closest mechanisms are the full
[M-unit primitive budget](planar_jc48_sep08_quartic_common_root.md),
[shared-root first jets](planar_jc48_sep08_shared_roots.md), and
the [full global section filtration](planar_jc48_sep08_dg_quadratic.md).
The least-used sidecar here is polynomiality in the original source
ring. A birational field alone loses that property. The first step
uses it before taking a field trace; neither step is substituted
for the other. The canonical hostile is the constant-D rational
submersion in Section 7. The corrected near miss is to use the
constant-D coefficient rows in the present class: the parameters
`A` and `lambda` are retained until their actual field cases split.

The live concepts are original polynomial weights, finite primitive
capacity, the weighted infinity, actual finite field degree, and
generic rational residues. The proof map in Section 3 preserves
the original polynomial ring as a subring of Laurent polynomials.
The later maps preserve the actual rational function field and
volume; they do not claim regularity on the source. Their necessary
sidecars are the explicit field inverses, a genuine quadratic
discriminant, and the original generic fibre constant.

## 2. Complete rows and the necessary active entry

Normalize the nonzero leading scalar and write `u=x-p`, retaining
arbitrary original finite position `p`. This is an ordinary
polynomial source coordinate, not an asserted automorphism of the
compactification. The complete pre-active rows for `N=u^6` are

    P=sum_(i=0)^4 P_i u^i,
    Q=(P4-1)u^2+(P3-2p P4+4p)u+Q0,
    M=sum_(i=0)^4 M_i u^i,
    R=M4u^2+(M3-2pM4)u+R0.

They follow from the original numerator rows
`Q,P-2Qx^2,N-Px^2+Qx^4`, each of degree at most four, and
`R,M-Rx^2`, each of degree at most two. Conversely these bounds
force `deg Q<=2`, `deg P<=4` and the displayed coefficients.
The source derives the unique high-coefficient equations and checks
all rows with symbolic `p`.

The original infinity is regular for the relative differential on
every normalized generic branch. Its leading section is

    r^8 N(1/r-p)=r^2(1-pr)^6.

Normal-unit and M-unit cases are covered by the proved local
suppliers. In the remaining shared case the full tangent equation
is

    (a r^2+b r s+h s^2)^2+s^3(lr+es)-f s^4,
    a!=0.

Putting `s=rZ` gives four nonzero simple roots generically; the
repeated-root eliminant has constant `-4a^2` and is independent
of the generic fibre coefficient. Their unweighted relative-form
order is minus one, and the actual multiplier `r^2` makes it one.
This pays the complete infinity fibre without transporting a finite
residue obstruction to it.

At the finite sixfold point, if `M` is a unit, the complete Morse
analysis bounds total primitive pole degree by two. A generic
transverse point of nonconstant `F|D` gives a relative-form zero
of order two, so a primitive on that component would have local
degree three. This contradicts its pole degree at most two. Thus
this entire M-unit branch has no rational mate. If `M` vanishes
but the normal coefficient is a unit, all boundary relative forms
are regular; a primitive would be constant on every compact
generic component, another contradiction. In the remaining case
the shared-root first-jet theorem forces

    P(0)=P'(0)=M(0)=M'(0)=0.

No generic irreducibility is needed for these componentwise degree
arguments. Every hypothetical polynomial mate is rational, so it
must enter these active rows. Rename `A=P4-2` to obtain exactly

    H=u^6t^2+[(A+2)u^4+b u^3+c u^2]t
         +(A+1)u^2+(b-2pA)u+d,
    L=(lambda u^4+mu u^3+nu u^2)t
         +lambda u^2+(mu-2p lambda)u+e.             (1)

The original slopes of `H|D,L|D` are `-A,-lambda`. Hence

    F|D is nonconstant  iff  (A,lambda)!=(0,0).     (2)

The square of a nonconstant affine function has a quadratic term
that cannot be canceled by an affine function; if `A=0`, the
remaining slope is `-lambda`. No constant-D specialization has
been used in (1).

## 3. The original polynomial-ring gate

Use the injective substitution `w=u^2t`, so

    C[u,t] subset C[u,u^-1,w],
    J_(u,t)(F,G)=u^2 J_(u,w)(F,G).

The full expression (1) has no negative `u` powers in this chart.
Its constant row is

    F=f(w)+O(u), f(w)=(cw+d)^2+nu w+e.             (3)

Here is the complete polynomial argument; it does not require a
finite degree cutoff or a rational primitive calculation. Suppose
`G` is polynomial in the original `u,t`. Write its finite Laurent
expansion

    G=sum_(j>=ell) u^j g_j(w), g_ell!=0.

If `ell<0`, each original monomial contributing to `g_ell` has
positive `t` degree. Thus `g_ell` is a polynomial divisible by
`w`. If `f` is nonconstant, the leading bracket term is

    -ell f'(w) g_ell(w) u^(ell+1).

Its coefficient cannot vanish. A constant bracket equal to one
would force `ell=-1` and `f'g_-1=1`, impossible because `g_-1`
is divisible by `w`. If `ell>=0`, the bracket has positive
`u` order (indeed at least two when `ell=0`) and cannot be one.
Therefore polynomial-mate existence forces `f` constant, or
equivalently

    c=nu=0.                                        (4)

This proof is the first part of the separately prepared
[polynomial-weight gate](planar_jc48_sep08_polynomial_weight_gate.md).
Its independent source bundle is pending audit at the time of
this candidate's writing; the complete argument needed here is
included above, so there is no unproved black-box implication.

For comparison, the next row is

    g(w)=2(cw+d)(bw+b-2pA)+mu w+mu-2p lambda.

The stronger second weight gate would force, after (4),
`mu=-2bd` and `-2p(lambda+2dA)!=0`. Those additional restrictions
are not used below. The two exactness obstructions cover every
coefficient in (1) satisfying (4).

## 4. Actual fields and the A-nonzero quadratic trace

Set

    v=u(1+u^2t), h0(v)=v^2+bv+d.

This is an actual birational field change, with
`t=(v-u)/u^3` and

    omega=u^-3 du wedge dv.

After (4), the complete residual is

    H=h0+A u(v-2p),
    L=lambda u(v-2p)+mu v+e.                        (5)

Assume first `A!=0`, and put `h=H`, `delta=lambda/A`. Then

    u=(h-h0)/[A(v-2p)],
    omega=A^2(v-2p)^2/(h-h0)^3 dh wedge dv.

Both expressions are identities in the original rational field;
no points lost by their denominators are declared regular. The
generic constant `f=F` gives the quadratic equation

    h^2+delta h-delta h0+mu v+e-f=0.                (6)

The discriminant of (6) has degree one in the independent formal
constant `f`, with coefficient four. It is not a square in
`C(f)(v)`. Thus (6) is a genuine separable degree-two extension
of that field, and its full trace is available. This does not
require geometric generic integrality after adjoining all
algebraic constants.

Write

    P(v)=h0(v)^2+mu v+e, B=2h0+delta, C=P-f,
    z=h-h0, so z^2+Bz+C=0.

With the convention `dF wedge eta=omega`, the actual relative
differential on (6) is

    eta=A^2(v-2p)^2 dv/[z^3(2z+B)].

For the two roots of the displayed quadratic, the **full** trace
identity is

    Tr[1/(z^3(2z+B))]=1/C^2-B^2/C^3.              (7)

The normalized trace would be half this value. One elementary
proof of (7) takes the negative residue at zero of
`dz/[z^3(z^2+Bz+C)]`; its residue at infinity is zero. The exact
source independently reconstructs the same full trace with the
literal two-by-two companion matrix. Therefore

    Tr eta=-A^2(v-2p)^2[B^2/(P-f)^3-1/(P-f)^2]dv. (8)

If the original rational mate existed, `eta=dG` relative to
`C(f)`. Trace commutes with this derivation, so (8) would be an
exact rational differential in `C(f)(v)`. Only this necessary
direction is used.

## 5. Generic residues force an impossible degree identity

The quartic `P` has leading coefficient one. For the formal
generic `f`, its four inverse roots are simple and avoid every
zero of `P'`. Use the local parameter `w=P(v)` at any one of
them and define rational functions of `v`

    L1=(v-2p)^2(2h0+delta)^2/P',
    L2=(v-2p)^2/P', D_P=(1/P') d/dv.

Apart from the nonzero factor `-A^2`, (8) becomes

    [L1/(w-f)^3-L2/(w-f)^2]dw.

Its residue is

    (1/2)D_P^2 L1-D_P L2.                          (9)

Expression (9) is a rational function of `v` independent of the
generic constant. Vanishing at a generic inverse root therefore
forces it to vanish identically. Integrating once in `C(v)` gives

    L1'-2(v-2p)^2=kappa P'                         (10)

for a constant `kappa`. This is an identity of rational functions,
not a claim about a single chosen level. Constants may depend on
the fixed coefficients but not on `v`.

At infinity, `L1~v^3` and `P'~4v^3`. More precisely the left
side of (10) has leading term

    (3-2)v^2=v^2.

If `kappa!=0`, the right side has degree three; if `kappa=0`,
it is zero. Both alternatives are impossible. This proves that
every residual (5) with `A!=0` has no rational mate, including
`lambda=0` and without any assumption on `b,d,mu,p`.

## 6. The A-zero rational field and its decay obstruction

If `A=0`, condition (2) gives `lambda!=0`. Now `H=h0(v)` and

    F=P(v)+lambda u(v-2p),
    u=(f-P(v))/[lambda(v-2p)].

Thus the original field is already the rational field `C(f)(v)`.
The exact relative differential is

    eta=lambda^2(v-2p)^2 dv/(f-P(v))^3.             (11)

At the same generic simple inverse roots, residue vanishing in
(11) requires `D_P^2 L2=0`. Integrating twice in `C(v)` gives

    L2=alpha P+beta.

But `L2=(v-2p)^2/P'` is nonzero and tends to zero at infinity,
whereas an affine function of the quartic `P` can tend to zero
only if both coefficients vanish. This is a contradiction.
No algebraic trace, irreducibility assumption, or compact genus
estimate is needed in this case.

Sections 4--6 exclude rational mates throughout the residual (4).
Section 3 forces every polynomial mate into that residual. Together
with the complete active-entry argument, this proves the candidate
nonconstant-D polynomial exclusion. The constant-D proved inputs
and the finite-double leading residue give the all-location binary
`(6,2)` polynomial consequence stated in Section 1.

## 7. Hostiles, exact controls, and remaining rational scope

The polynomial-weight gate genuinely needs polynomiality. For

    q=1+u^2t, v=u q,
    H=u^2q^2+q, L=v,
    G0=1/(2v)-H,

the original global rows hold for every finite `p`, `F|D` is
constant, and `J(H^2+L,G0)=1`. This is a nonconstant-`L`
rational submersion from the proved constant-D boundary table.
Its row (3) is `(w+1)^2`, so it illustrates exactly why the
weight gate cannot be applied to an arbitrary rational mate.
Its actual source poles and polynomial nonentry remain as
proved in the constant-D bundle.

The discarded direct route used successively higher inverse
residues and individual source-critical eliminants. For example
`lambda=0,A!=0` forces `mu=0` by a T10 residue; one remaining
family then has an explicit critical sextic. Those computations
are useful scouts but are not needed by this final proof. The
polynomial source weight and actual trace preserve the missing
information and close all parameters without a per-case census.

The [exact companion](../../04-computation/planar_jc48_sep08_nonconstant_d_six.py)
checks the full pre-active row equations, actual original-chart
slopes and multiplicity, complete exceptional rows, explicit field
inverses and volume, the actual quadratic discriminant, full versus
normalized companion-matrix trace, local residue coefficients,
the leading degree mismatch, the rational-field decay, and the
constant-D rational hostile. Small negative-weight monomials are
controls only; the unbounded polynomial proof is Section 3.

The normal and optimized producer replays are complete and match
the frozen output byte-for-byte. Reproduce from the repository root:

```sh
python3 -B 04-computation/planar_jc48_sep08_nonconstant_d_six.py
python3 -B -O 04-computation/planar_jc48_sep08_nonconstant_d_six.py
```

Both pass **64 always-active gates**. The trace control first verifies
the actual centered quadratic coefficients and then uses its universal
two-coefficient companion matrix; this avoids needless expression
expansion without changing the exact identity or specialization.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Source | 7,363 | `946cae045bfd1fd3b4ce8c054cf444196c2a2e99d00d0265ef7e070684bdaa3f` |
| Frozen output and both producer replays | 399 | `16d2315be218c959d98169b153b6f5c2395f3fb12cbf4f8804c45496278693c6` |

The semantic gate digest is
`5821fc3ad3231a3db28a6c9cc3652432ba765959a9c219b76354a2cc578082a8`.
The [independent full audit](planar_jc48_sep08_nonconstant_d_six_audit.md)
accepts every full-coefficient branch, the polynomial gate, both actual
fields and volume forms, the complete trace and both generic-residue
contradictions. All64 gates match in both modes. Combined with the
constant-D theorem, this closes the entire6+2 polynomial boundary class.
