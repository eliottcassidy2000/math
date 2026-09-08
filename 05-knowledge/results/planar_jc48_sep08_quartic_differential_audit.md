# Independent audit of the compact quartic differential criterion

**Status: independent analytic/source audit PASS; normal, optimized, and
frozen-output replay PASS. Audit frozen.** The exact controls are
FINITE-EXACT. The arbitrary-coefficient local argument and its compact
exactness consumer are proved analytically, not inferred from those
controls. No source or prior frozen artifact was edited for this audit.

The audited primary is
[planar_jc48_sep08_quartic_differential.md](planar_jc48_sep08_quartic_differential.md).
The domain is the specified complex surface
`W=(P1_x x P1_z)\{z=x^2}`, with `omega=dx wedge dt`, not an
arbitrary Keller finite envelope. Write `s=z-x^2`,

    H=N/s^2 in L_2,  L=M/s in L_1,  F=H^2+L.

The sufficient exclusion hypotheses are precisely: `N|S` is nonzero,
each of its projective zeros has multiplicity at most two, and `M`
is nonzero at each such zero. The conclusion is absence of a
**global** `G on W` with `dF wedge dG=lambda omega`, `lambda!=0`.
It imposes no degree bound on `G`. It does not exclude a nonglobal
affine-chart mate, close the whole quartic stratum, or settle JC(2).

## 1. Inherited geometry and the full boundary

The [DG filtration](planar_jc48_sep08_dg_quadratic.md) supplies the
actual global sections `N in H0(O(4,2))` and `M in H0(O(2,1))`.
Their restrictions to the graph `S` have degrees eight and four:
the intersection of bidegree `(a,b)` with `(2,1)` is `a+2b`.
In particular a finite polynomial restriction alone does not determine
the zero multiplicities without accounting for infinity.

I independently checked the two-form and transition signs. On the
finite compact chart, `omega=s^(-2) ds wedge dx`. On the corner
chart `r=1/x, q=1/z, sigma=q-r^2`,

    omega=r^2 sigma^(-2) d sigma wedge dr,
    Ninf=r^4 q^2 N(1/r,1/q),
    Minf=-r^2 q M(1/r,1/q).

Thus `H=Ninf/sigma^2` and `L=Minf/sigma` exactly. The sign in
`Minf` is necessary and is retained. These formulas give
`div(omega)=2 Dbar-2S`; the divisor class agrees with
`K_(P1 x P1)=(-2,-2)`. There is no further pole on the finite-x,
infinite-z chart: there `t=q/(1-x^2q)` and its differential is regular
at `q=0`. The extra numerator `r^2` at `S intersect Dbar` can only
improve the orders used in the proof.

The compact fibre equation in any boundary trivialization is

    E=N^2+s^3M-cs^4=0.

Its intersection with `S` consists precisely of the zeros of `N|S`.
Since the latter section is not identically zero, `S` is not a
component of the compact closure, and `F` has a genuine pole of order
four at the generic point of `S`.

## 2. The complete local calculation for multiplicity at most two

At a specified zero of `N|S`, the unit `M-cs` has a local analytic
cube root for every fixed `c`. The map

    v=s(M-cs)^(1/3),  w=w

is a local analytic coordinate change. It transforms the equation
to `Ntilde(v,w)^2+v^3=0`. The normalized base cusp has coordinates
`v=-tau^2, Ntilde=tau^3`, so its pullback is

    Psi(tau,w)=Ntilde(-tau^2,w)-tau^3=0.

This does not replace the original branch by an unrelated finite
cover: away from the boundary the inverse is `tau=-Ntilde/v`.
Consequently the normalizations of the corresponding local branches
agree. No component of the pulled-back curve is contained in
`tau=0`, since `Ntilde(0,w)` is not identically zero.

Fix the convention `dF wedge eta=omega`. On the fibre,
`dF=s^(-4)dE`; changing from `s` to `v` changes the two-form by a
unit. Since `E_w=2 Ntilde Ntilde_w`, pulling the relative form back
gives, up to its harmless orientation sign,

    eta=(holomorphic multiplier) tau^2 d tau/Psi_w.

This derivation uses the actual differential and its Jacobian, not
only a named singularity type. Away from `Dbar` the multiplier is a
unit; at the infinity intersection it is a unit times `r^2`.

Let `k=ord_w Ntilde(0,w)`. The following cases are exhaustive.

* For `k=1`, `Psi_w` is a unit. The implicit function theorem makes
  `tau` a parameter and gives differential order two.
* For `k=2` and `Ntilde_v(0,0)!=0`, the quadratic part is
  `a tau^2+b w^2`, with `ab!=0`. The two distinct tangent lines
  give two smooth branches. Both have `ord(tau)=ord(w)=1` and
  `ord(Psi_w)=1`; the differential order is one on each branch.
* For `k=2` and `Ntilde_v(0,0)=0`, the initial weighted equation is
  `b w^2-tau^3`, with weights `(2,3)`. The terms `tau^2 w`,
  `tau^4`, `w^3`, and all other terms have strictly larger weight.
  Its coprime Newton pair gives one branch with orders
  `(ord(tau),ord(w))=(2,3)`. The leading `Psi_w=2b w+...` has
  order three; `tau^2 d tau` has order five. The differential
  therefore has order two.

These arguments allow all higher coefficients. They prove
holomorphy at every normalized boundary point, including infinity;
there is no exceptional parameter locus within the stated hypotheses.

## 3. Every generic component and the exactness contradiction

The set of critical values of `F on W` is finite. Indeed the algebraic
critical locus has finitely many irreducible components, and in
characteristic zero the restriction of `F` to each positive-dimensional
component has zero differential and hence is constant; isolated
critical points contribute only finitely many more values. Nonconstancy
of `F` excludes a two-dimensional critical component.

Choose `c` outside this set and, if `F|D` is constant, outside that
one additional value. Every component of the fibre is smooth on `W`.
On its compact normalization the relative form is regular there,
and Section 2 pays all points outside `W`. No component can be `S`
or `Dbar`, so every component meets the open subset on which `omega`
is nonzero. Thus `eta` is a nonzero holomorphic differential on
**each** compact connected normalized component. Geometric
irreducibility of the generic fibre is neither assumed nor needed.

If a global `G` satisfied the bracket equation, its restriction to
any such component would be meromorphic and obey `dG=lambda eta`.
A pole of order `m>0` would give a differential pole of order `m+1`
in characteristic zero. Therefore `G` has no poles on that compact
component and is constant, contradicting `eta!=0`. This proves the
typed global exclusion. It does not infer a polynomial primitive
from rational integrability or from local residue conditions.

## 4. Independent reconstruction of both global hostiles

Both examples use `M=-(1+x^2z)`, whose boundary form is
`-(1+x^4)` and whose infinity trivialization has value one on the
boundary at infinity. I reconstructed their compact equations from
the literal sections, rather than importing the producer's local
formulas.

For

    N3=xz-(31/27)x^3+x^3z,

the boundary restriction is `x^3(x^2-4/27)`. Its finite
multiplicities are three at zero and one at each root of
`x^2=4/27`; the infinity multiplicity is three. The bidegree bounds
are valid, and `M` avoids every one of these points.

With `s=tau^2` and `x=3tau/2+d tau^2`, my literal expansion gives

    [tau^8]E3=-(4/3)d^2+351/16-c,
    [tau^6](E3)_x=-(8/3)d,
    d^2=(1053-48c)/64.

The final prose correctly treats this as a parameter substitution,
not as the exact `v` coordinate of Section 2. Divide the substituted
equation, with `d` replaced by `h`, by `tau^8`. Its derivative at
`tau=0,h=d` is `-8d/3`, so for generic `c` the implicit function
theorem supplies the full branches. Since `x'(0)=3/2`, these are
actual smooth normalized branches, not a hidden double cover.
The displayed relative form has residue `-3/(4d)`, up to the chosen
overall orientation sign, and is genuinely logarithmic. Thus `M`
being a unit alone does not imply regularity.

For

    N6=x^2z-x^4-(8/27)x^4z+(4/27)x^2z^2,

the boundary restriction is `-4x^6/27`, with multiplicities six at
zero and two at infinity. It again satisfies the exact global
bidegree bounds. The full second-chart expression

    H6=-1-r^2b-(8/27)b-(4/27)r^2b^2

confirms that no forbidden local coefficient was inserted. Put
`y=x^2`; the original equation becomes

    E6bar=(sy-4y^3/27+4s^2y/27)^2
           -s^3(1+y^2+sy)-cs^4.

My independent expansion at `s=4y^2/9+e y^3` gives

    [y^8]E6bar=-e^2/3-64(65+36c)/59049,
    e^2=-64(65+36c)/19683,
    [y^5](E6bar)_s=-2e/3.

These are the corrected frozen coefficients. The producer detected
an earlier transcription of the split constant on its initial
replay, and repaired it before freeze; my reconstruction independently
confirms the repaired value. The pole coefficient did not change.

For generic `c`, division by `y^8` after
`s=y^2(4/9+yh)` leaves a simple root `h=e`, so `s` is an actual
power series in `y`. With `x` the branch parameter, the relative
form is

    -s^2/(E6bar)_s dx=8/(27e) x^(-2) dx+O(1) dx.

Its coefficient is an even Laurent series in `x` to every order.
Hence its residue is identically zero, while its pole has order
two. This proves the second hostile to the claim that `M`-unit
fibres always have at most logarithmic poles. It does not assert
that this form has a global primitive or that the example has a mate.

## 5. Source audit, complete replay, and scope

I read the complete standalone verifier and independently replayed

    python3 -B 04-computation/planar_jc48_sep08_quartic_differential.py
    python3 -B -O 04-computation/planar_jc48_sep08_quartic_differential.py

Both completed with **78 always-active gates**. Their full outputs
are byte-identical to the frozen primary output, **513 bytes**.
The source has no imported inherited mathematical implementation.
Its universe is three explicit safe global sections, the two
global hostile sections, their full chart and boundary checks, and
the displayed symbolic generic split and polar coefficients. The
safe sections use the global choice `M=1`, and have no boundary
zero at infinity. They realize simple roots and both double-root
cases. The finite parity controls illustrate the all-integer parity
argument above; they are not a finite substitute for it.

| Artifact | SHA256 |
|---|---|
| Source, 6,350 bytes | `865628f671b1eb80984d95fc3d9ba5df5b1ede389a98f11f82841e143c63436a` |
| Frozen output and both independent replays | `931be65a167336b2aca148e02f4e218813d9d30b8566dd940bfb33e7bad14603` |
| Semantic record | `822ffe67992fdc7aec4f0a00d128aebd03824b003972cfbcc73a291e1775c80b` |

The final primary text was accepted after the actual-branch wording
repair. No mathematical correction remains. The necessary unresolved
boundary is a common zero with `M` or multiplicity at least three;
the examples distinguish nonzero residues from higher, residue-free
principal parts. They do not establish sufficiency of that boundary
for a mate. Stronger infinity multiplicity allowances and the general
higher-principal-part problem remain outside this theorem.
