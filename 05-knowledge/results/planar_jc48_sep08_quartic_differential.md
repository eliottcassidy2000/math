# A compact differential exclusion for simple and double quartic boundary roots

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED.**
This concerns global functions on the specified DG surface. It does not
assert that this surface is a Keller finite envelope, or settle JC(2).

## 1. The actual surface, sections, and necessary boundary predicate

Let `X=P1_x x P1_z`, let `S={z=x^2}`, and let `W=X\S`, with the
full charts and two-form of
[the DG filtration](planar_jc48_sep08_dg_quadratic.md):

    U0=A2_(x,t),       t=1/(z-x^2),
    Uinf=A2_(r,b),     x=1/r, t=-r^2-r^4b,
    D={r=0} subset W,  omega=dx wedge dt=r^2 dr wedge db.

On `X`, `div(omega)=2Dbar-2S`, where `Dbar` is the closure of `D`.
In particular the numerator at `Dbar intersect S` has an extra zero,
not an extra pole. In the finite chart put `s=z-x^2`. Write

    H=N/s^2 in L_2,    L=M/s in L_1,    F=H^2+L.

Here `N` and `M` are global sections of `O_X(4,2)` and `O_X(2,1)`,
respectively; the displayed expressions use their finite-chart
trivializations. Their restrictions to `S=P1` are binary forms of
degrees eight and four. All multiplicities below are on this complete
projective boundary, including infinity.

**Theorem.** Suppose `N|S` is not identically zero, every zero
of `N|S` has multiplicity at most two, and `M` is nonzero at every
such zero. Then no global `G on W` satisfies

    dF wedge dG = lambda omega,             lambda != 0.

Thus a hypothetical global pair in this quartic square-prefix class
must have either a common zero of `N|S` and `M|S`, or a zero of
`N|S` of multiplicity at least three. This is a necessary predicate,
not a construction of a mate. The theorem deliberately retains the
same conservative multiplicity bound at infinity.

The closest mechanisms are the global square-prefix reduction in
[the quartic boundary note](planar_jc48_sep08_quartic_boundary.md),
the exact relative-primitive obstruction in
[the boundary-linear note](planar_jc48_sep08_boundary_linear.md),
and the actual compact differential rather than an intermediate
coefficient statistic. The source of this connection is `F=H^2+L`
on `W`; its target is the normalization of a compact generic fibre;
the map is restriction of `omega/dF`. It preserves exactness for an
actual global mate. Forgetting the boundary principal parts destroys
the regularity information needed for that consumer.

The live concepts are: the binary octic and quartic on `S`; cusp
normalization; the divisor of `omega`; logarithmic versus higher
principal parts; and exactness on a compact curve. The two hostiles
below distinguish these concepts. The corrected near miss is that
`M` being a unit would by itself make the differential regular, or
even logarithmic. The least-used sidecar is the full normalization
and its valuations at every boundary point.

## 2. The local proof, including double roots and infinity

The compact fibre over `c` has local equation

    E=N^2+s^3 M-c s^4=0.                                  (1)

Its boundary points are precisely the zeros of `N|S`. At any of them,
`M-cs` is a unit. Choose a local analytic cube root of that unit and
make the invertible coordinate change

    v=s(M-cs)^(1/3).

Equation (1) becomes `Ntilde(v,w)^2+v^3=0`, with `w` a coordinate
along `S`. The model cusp has normalization

    v=-tau^2,       Ntilde=tau^3.

Its pullback is the actual local plane curve

    Psi(tau,w)=Ntilde(-tau^2,w)-tau^3=0.                    (2)

Away from `v=0`, `tau=-Ntilde/v`; hence this pullback is birational
to (1) and does not introduce an unmarked cover of a generic branch.
Normalizations can therefore be computed using (2).

Up to sign and a holomorphic multiplier, the relative form
`eta=omega/dF` is

    s^2 ds/E_w  =  (holomorphic multiplier) tau^2 d tau/Psi_w.   (3)

The multiplier is a unit away from `Dbar`, and may vanish at
`Dbar intersect S`. The change of normal coordinate has invertible
Jacobian. Differentiating the cusp equation proves (3) directly;
no inference from topological cusp type alone is needed.

Let `k=ord_w Ntilde(0,w)`.

* If `k=1`, `Psi_w` is a unit. The pullback curve is smooth, and
  the model differential has order two.
* If `k=2` and `Ntilde_v(0,0)!=0`, the quadratic part of `Psi`
  is `a tau^2+b w^2` with `ab!=0`. It has two distinct tangent
  lines. On either smooth branch `tau` is a parameter and
  `Psi_w` has order one; the differential has order one.
* If `k=2` and `Ntilde_v(0,0)=0`, the weighted leading equation is
  `b w^2-tau^3`, `b!=0`. All other terms have higher weights for
  `wt(tau)=2, wt(w)=3`. This is an ordinary `(2,3)` branch:
  `ord(tau)=2`, `ord(w)=3`, `ord(Psi_w)=3`. The numerator in (3)
  has order five, so the differential has order two.

These exhaust all cases with `k<=2`. The local assertions also
follow directly by the implicit function theorem in the first two
cases and the Newton parametrization in the last. They include
arbitrary higher coefficients, not a truncated family assumption.

For clarity, the complete infinity calculation uses
`r=1/x`, `q=1/z`, `sigma=q-r^2`. In this chart

    omega = r^2 sigma^(-2) d sigma wedge dr,
    Ninf=r^4 q^2 N(1/r,1/q),
    Minf=-r^2 q M(1/r,1/q),
    H=Ninf/sigma^2,       L=Minf/sigma.

Thus (1)--(3) apply there as well, with the additional numerator
`r^2`. No finite-chart calculation silently drops the point at
infinity or the intersection with `Dbar`.

## 3. The compact exactness consumer

Choose `c` outside the finite set of critical values of `F on W`
and outside the value of `F|D` if the latter is constant. Such
values exist: in characteristic zero every irreducible component
of the critical locus has constant `F`, and there are finitely many
components. The assumption `N|S != 0` also makes `F` nonconstant
and prevents `S` from being a component of a fibre closure.

On each compact normalized component of this generic fibre, `eta`
is regular at all boundary points by Section 2. It is regular on
`W` because the fibre is smooth there and `omega` is regular.
Every component meets `W` away from `D`: neither `S` nor `Dbar`
is a generic component. The form is therefore nonzero on every
component, since `omega` is nonzero there and `dF!=0`.

If a global mate existed, its restriction would be a meromorphic
function on each compact normalization and would satisfy
`dG=lambda eta`. A pole of order `m>=1` of a meromorphic function
produces a pole of order `m+1` of its differential in characteristic
zero. Since `eta` is holomorphic, `G` has no poles, and hence is
constant on that compact connected component. This contradicts
the nonzero differential. No irreducibility of the generic fibre,
source quasifiniteness hypothesis, or identification of `W` with
an actual finite envelope is used.

The compact-curve fact used here is elementary; the displayed
argument supplies the full exactness implication. No literature
priority claim is made.

## 4. Two globally admissible failed extensions

Both controls use

    M=-(1+x^2 z),
    L=-[(1+x^4)t+x^2].

This is a global `L_1` function. Its boundary form is `-(1+x^4)`
and is nonzero at infinity. In both controls it avoids every zero
of the corresponding binary octic, not just the displayed finite
point. They are not Keller candidates or counterexamples to the
safe theorem.

### 4.1 Multiplicity three: a genuine logarithmic pole

Take

    N3=xz-(31/27)x^3+x^3z
      =sx-(4/27)x^3+sx^3+x^5,
    H3=(x+x^3)t+(x^5-(4/27)x^3)t^2.

The bidegree is at most `(4,2)`, so `H3` is global. On the complete
boundary, `N3` has multiplicity three at zero and infinity, and
simple roots at `x^2=4/27`. The chosen `M` avoids all four points.

At zero make the substitution `s=tau^2`. Two branches of the
generic fibre have

    x=(3/2)tau+d tau^2+...,
    d^2=(1053-48c)/64.

For generic `c`, `d!=0`. To justify the actual branches, substitute
`x=3tau/2+tau^2 h` in (1) and divide by `tau^8`. At `tau=0`, the
derivative of this equation with respect to `h` at `h=d` is
`-8d/3!=0`. The implicit function theorem gives each actual branch,
and `x'(0)=3/2` makes its map to the original curve a smooth
normalization parameter. The two signs of `d` distinguish the
branches. No unmarked double cover is inferred from the substitution.
Exact leading coefficients in the original equation (1) are

    E3(tau^2,(3/2)tau+d tau^2)
       =[-(4/3)d^2+351/16-c] tau^8+O(tau^9),
    (E3)_x=-(8d/3)tau^6+O(tau^7).

Consequently

    s^2 ds/(E3)_x = -3/(4d) d tau/tau + regular terms.

The nonzero residue refutes the proposed extension from simple
and double roots to arbitrary multiplicity with `M` a unit. In
this control the residue itself also prevents a rational exact
primitive; failure of the holomorphic proof is not existence of
a mate.

### 4.2 Multiplicity six: a residue-free pole of order two

Take

    N6=x^2z-x^4-(8/27)x^4z+(4/27)x^2z^2
      =sx^2-(4/27)x^6+(4/27)s^2x^2,
    H6=x^2t-(4/27)x^6t^2+(4/27)x^2.

Its complete boundary divisor has multiplicities six at zero and
two at infinity; `M` is nonzero at both. A direct full-chart check is

    H6=-1-r^2b-(8/27)b-(4/27)r^2b^2

on `Uinf`. Thus no forbidden local coefficient has been introduced.

Put `y=x^2` and write `E6bar(s,y)` for its original fibre equation.
Two generic branches admit

    s=(4/9)y^2+e y^3+...,
    e^2=-64(65+36c)/19683.

For `c!=-65/36` the coefficient `e` is nonzero and the branches
are distinct. Their expansions are ordinary power series in `y`:
after setting `s=y^2(4/9+y h)`, the divided equation has a simple
root at `h=e`, so the implicit function theorem applies. The exact
leading derivative is

    (E6bar)_s=-(2e/3)y^5+O(y^6).

Using differentiation of the original equation along a branch gives

    s^2 ds/(E6)_x = -s^2/(E6bar)_s dx
                  = 8/(27e) x^(-2) dx + O(1) dx.

The coefficient is nonzero. Moreover the complete Laurent
coefficient is an even function of `x`, because `s` is a power
series in `x^2`. It has **zero residue identically**, not just at a
special value of `c`. This refutes the second proposed extension,
that generic `M`-unit fibres always have at most logarithmic poles.
It also shows why a residue test alone does not settle the remaining
cases. No global primitive or mate is asserted for this example.

## 5. Exact controls, reproduction, and stopping point

The standalone source checks the full chart identities, the
projective boundary multiplicities and avoidance, all displayed
leading coefficients and residues, and safe simple/double controls.
It does not use a finite scan to prove the arbitrary-coefficient
local theorem. Reproduce with

```sh
python3 -B 04-computation/planar_jc48_sep08_quartic_differential.py
python3 -B -O 04-computation/planar_jc48_sep08_quartic_differential.py
```

Both runs pass **78 always-active gates**, with all 513 output bytes
identical to [the frozen output](planar_jc48_sep08_quartic_differential.out).
The initial exact replay caught and repaired a transcription of the
`m=6` split constant before this freeze; the displayed value is the
coefficient reconstructed from the original unsimplified equation.

| Artifact | SHA-256 |
| --- | --- |
| [Source](../../04-computation/planar_jc48_sep08_quartic_differential.py), 6,350 bytes | `865628f671b1eb80984d95fc3d9ba5df5b1ede389a98f11f82841e143c63436a` |
| Frozen output and both replays | `931be65a167336b2aca148e02f4e218813d9d30b8566dd940bfb33e7bad14603` |
| Semantic record | `822ffe67992fdc7aec4f0a00d128aebd03824b003972cfbcc73a291e1775c80b` |

The [independent analytic/source audit](planar_jc48_sep08_quartic_differential_audit.md)
accepts the full local and compact proof, both global hostile controls,
and normal/optimized/frozen replay. Its frozen SHA-256 is
`6f5fc40052b81728ea96037a14f46a7e00a36422fddf6335710147389c5413fb`.
Source and output were unchanged by status promotion.

The precise stopping point is a higher-principal-part problem on
the compact normalization. At multiplicity six, global section
constraints permit a second-kind pole with no residue. A further
consumer must compare all principal parts and global exactness;
it cannot infer either from `M`-unit, local residue vanishing, or
the existence of a global square prefix. A separate cheap question
is whether the additional `r^2` at infinity permits a stronger
local multiplicity allowance there. Both extensions remain OPEN.
