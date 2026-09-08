# All-finite quartic boundary type (4,3,1)

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This concerns the complete declared square-prefix entry on the fixed DG
surface. It does not claim general quartic or Jacobian conjecture closure.
The fourfold point at infinity is not included.

## 1. Statement, actual map, and inheritance

Use the surface and source volume

    W=(P1_x x P1_z) minus S,  S={z=x^2},
    t=1/(z-x^2),  x=1/r,  t=-r^2-r^4 b_D,
    D={r=0},  omega=dx wedge dt=r^2 dr wedge db_D.

Let `H in L_2`, `L in L_1`, `deg_t H=2`, and `F=H^2+L`.
Assume the binary octic given by the leading coefficient of `H`
has three distinct finite zeros of multiplicities `4,3,1`.

**Theorem.** No polynomial `G in C[x,t]` has
`J_(x,t)(F,G)` a nonzero constant, irrespective of its degree.
More precisely, a rational mate forces `L` constant. The final
obstruction when `L` is constant is explicitly polynomial, and
Section 8 gives a rational mate in that same boundary stratum.

The source is the complete pair of global sections. The maps are
restriction to the original compact generic fibre and the actual
coordinate `q=1+(x-p)^2t`. They retain the volume, fibre constant,
and induced normal coefficients. Leading exactness alone loses
the normal jets; residues alone lose the space of possible global
primitives. The final sidecar is that complete primitive space.

Dependencies are [leading exactness](planar_jc48_sep08_leading_exactness.md),
the [degree-at-most-eight classification](planar_jc48_sep08_boundary_exactness.md),
[M-unit Newton/Morse bounds](planar_jc48_sep08_quartic_common_root.md),
the [finite first-jet obstruction](planar_jc48_sep08_shared_roots.md),
and the [complete global filtration](planar_jc48_sep08_dg_quadratic.md).
The compact degree mechanism is already in the
[pole-degree note](planar_jc48_sep08_pole_degree.md).

Live concepts are the marked divisor, inverse-series coefficients,
zeros as well as poles of the relative form, geometric integrality,
and complete primitive spaces. The least-used sidecar is genus-one
Riemann--Roch. The hostile `h=x^2+x^4t`, `F=h^4+h`, with rational
mate `1/[3x^3(4h^3+1)]`, shows that nonconstant `L` does not alone
exclude composition. Its square leading coefficient is the missing
sidecar in that attempted shortcut. Another corrected near miss is
to conflate degree four of `F|D` with local degree three of a
primitive at `D`.

## 2. Normalize the finite points and retain every coefficient

Write `H=N(x)t^2+P(x)t+Q(x)` and `L=M(x)t+R(x)`.
Distinguish the leading coefficients from the section numerators:
with `s=z-x^2`, put `calN=N+sP+s^2Q`, `calM=M+sR`.
Thus `H=calN/s^2`, `L=calM/s`, and the first normal derivative
`calN_s|S=P`; the polynomial `N(x)` itself has no `s` argument.
The full global spaces are

    deg N<=8,
    P=2N8 x^6+2N7 x^5+sum_(i=0)^4 p_i x^i,
    Q=N8 x^4+N7 x^3+(p4-N6)x^2+(p3-N5)x+q0,
    deg M<=4,  R=M4 x^2+M3 x+r0.                      (1)

If `M=0`, then `L` is constant and the nonconstant polynomial
factor `2H` in `J(F,G)=2H J(H,G)` excludes polynomial mates.
Henceforth `M!=0`.

Leading exactness and the complete `(4,3,1)` residue condition
force the points to be `p,p+3d,p-d`, with `d!=0`. The actual
surface scaling `x=d X`, `t=d^-2 T`, and nonzero constant
rescalings of `H,F,G`, normalize `d=1` and the leading scalar to
one. Its volume multiplier is nonzero. Retain the resulting `p`
and put `u=x-p`, so

    N=u^4(u-3)^3(u+1)=u^8-8u^7+18u^6-27u^4.          (2)

The translation to `u` is only an ordinary source coordinate; it
is not claimed to be an automorphism of `W`.

A rational mate also forces `M dx/N` rational-exact. Indeed, in
`F(x,T)=v^-4`, choose `a^2=N` and expand

    T=a^-1 v^-1-P/(2a^2)+(P^2-4a^2Q)v/(8a^3)+T2 v^2+... .

The coefficient of `v^3` in the equation multiplied by `v^4` is
`4a T2+M/a`, giving `T2=-M/(4N)`. The coefficient of `v^6` in
the proposed mate has derivative `T2/2`. Normalized field trace
to `C(x)` commutes with differentiation and yields the claimed
rational primitive. No lower-row term has been omitted.

## 3. The fourfold point must be active

Call a finite point **active** when it is shared by the sections
and `calN_s` vanishes there. Under a rational mate, the first-jet
theorem then requires `m>=3`, `j=ord calN_s|S>=2`, `n=ord M>=2`.
The M-unit and normal-unit models remain separate.

At multiplicity four, an M-unit is unbalanced (`4!=3j`) and
all forms are regular; a normal unit is regular too. The simple
point is regular. A nonactive triple point has no possible
primitive poles: its balanced M-unit case gives regular forms
or forbidden nonzero logarithms.

If only the triple point is active, its equation is

    E=(a w^3+s B(w)+s^2 C(w))^2
         +s^3(M0(w)+(R(w)-zeta)s),
    ord B>=2, ord M0>=2.

Putting `s=w^(3/2)Z` gives the generic leading quartic
`(a+C(0)Z^2)^2+(R(0)-zeta)Z^4`, with four simple nonzero
determinations. They form two branches after `w=tau^2`; each
has differential order `-2`. This exhausts the local degree
four and gives total primitive pole capacity two. Meanwhile
`F|D` has degree four, since `N8=1`. A generic transverse
point on `D` has differential order exactly two, forcing local
degree three of its primitive. This exceeds even the total
capacity two, and hence the capacity on that same component.
No connectedness assumption is needed for this argument.
If there is no active point, the primitive instead has no poles
and is constant on every compact component, also impossible.

Thus the fourfold point is active and

    P(0)=P'(0)=M(0)=M'(0)=0.

The `T2` residue at `u=-1` forces `M(-1)=0`. Write
`M=u^2(Au^2+Bu+C)`. This gives `C=B-A`, and the residue at
zero gives `B=0`. Consequently

    M=lambda u^2(u^2-1),  lambda!=0,
    M/N=lambda(u-1)/(u^2(u-3)^3)
       =d/du [-lambda/(3u(u-3)^2)].                   (3)

In particular `M(3)=72lambda` is a unit. Equivalently, two
active multiple points and the simple root would impose five
counted zeros on `M` of degree at most four. Its complete lower
row, up to a constant, is

    L=lambda Z,  Z=(u^2-1)q+1-2pu,  q=1+u^2t.        (4)

## 4. The moving-root residue fixes the induced normal rows

Before this residue, (1)--(2) and the active jets give exactly

    P=u^2[2u^4+(-4p-16)u^3+A u^2+B u+b],
    Q=u^4-4(p+2)u^3+(A+4p^2-18)u^2
        +(-2Ap+B+32p^2+72p)u+c.                       (5)

Let `Q1` be the coefficient of `u` in `Q`. In `(u,q)`,
`omega=u^-2 du wedge dq`, and

    H0=-27(q-1)^2+b(q-1)+c,
    H1=B(q-1)+Q1,
    F=f(q)+u g(q)+O(u^2),
    f=H0^2-lambda(q-1),
    g=2H0[B(q-1)+Q1]-2p lambda.                       (6)

At an original generic root `f(q0)=zeta`, the actual root has
`q(u)=q0-g(q0)u/f'(q0)+...`. In the convention
`dF wedge eta=omega`, one has `eta=-du/(u^2 F_q)` and residue
`(g'f'-g f'')/(f')^3`. Exactness for every original generic root
forces `(g/f')'=0`, hence `g=C f'` with scalar `C`. Complete
coefficient comparison, using `lambda!=0`, gives

    C=2p,  B=-108p,  Q1=2pb.                          (7)

For `p=0`, only `B=Q1=0` follows and `A` stays free. Both
original source derivatives vanish on `u=0`, excluding a
polynomial mate immediately. We do not divide by `p`; Section 7
also excludes rational mates in the complete anchored family.

For `p!=0`, instead `A=16p-18-b`, and the exact full rows are

    P=u^2[2u^4+(-4p-16)u^3+(16p-18-b)u^2-108pu+b],
    Q=u^4+(-4p-8)u^3+(4p^2+16p-36-b)u^2+2pb u+c.     (8)

## 5. A zero of order four removes the normal-unit triple case

Put

    k=b+20p+54,
    h=(2p-1)b+c-108p^2+108p-27.

Then `P(3)=-72(k+16p)`. At the active fourfold point the
substitution `s=u^2Z` has leading polynomial

    (a+bZ+C0 Z^2)^2-lambda Z^3+(R0-zeta)Z^4,
    a=-27!=0.                                        (9)

It has four simple nonzero roots for generic `zeta`: a repeated
nonzero root satisfies `Z P0'-4P0=0`, whose constant is
`-4a^2`, and each of its finitely many roots determines at most
one exceptional fibre value. These four actual smooth branches
exhaust the local degree four. On each, `E_s` has order six
and `s^2du` order four, so `eta` has order `-2`. The total
possible primitive pole degree is at most four.

If `P(3)!=0`, the M-unit triple point has `j=0`. Its complete
cancellation pair has `s~w^3` and `eta=(unit)w^(3/2)dw`.
It is one actual branch `w=tau^2` with differential order four,
forcing local degree five of a primitive. This is greater than
its own component's pole capacity, which is at most four.
There are no other primitive poles. Hence a rational mate needs

    k=-16p,  b=-36p-54,  P'(3)=864p!=0.                (10)

This zero-degree obstruction replaces an incomplete critical-point
case split; it retains information that residues alone miss.

The same degree argument works before dividing by `p`. When
`p=0`, the condition `P(3)=0` gives `A=30-b/9`, rather than the
specialization of (8). The resulting full anchored family is

    H=D0 q^2+(b+54)(1-u^2/9)q+c-b-27,
    D0=(u-3)^3(u+1),  Z=(u^2-1)q+1.                  (11)

Its triple derivative is `P'(3)=-6(b+54)`. Thus its normal
order is one, or at least two when `b=-54`. Both possibilities
are retained below.

## 6. Geometric integrality and the full genus-one primitive space

Retain `lambda!=0` and either the family (8),(10) for `p!=0`,
or (11) for `p=0`. The geometric generic fibre is integral.
Indeed, if `F=Phi(K)` with a nontrivial outer polynomial, its
degree in `t` forces `deg Phi=2` or `4`. In degree four, `K`
is linear in `t`, so `N^2` is proportional to the fourth power
of its leading coefficient. This would make `N` square up to
scalar, contrary to its odd multiplicities three and one. In
degree two, complete the outer square to `F=K1^2+gamma` and
choose the top sign so `K1+H` has degree two in `t`. Then
`(K1-H)(K1+H)=L-gamma` has degree at most one, forcing
`K1=H` and `L` constant, contrary to `lambda!=0`.

The noncomposite-to-integral route is classical:
[Arzhantsev--Petravchuk, Theorem 1 and Lemma 3](https://arxiv.org/pdf/math/0608157v2)
identify noncomposite polynomials with closed polynomials and
relative algebraic closure of `C(F)` in `C(u,t)`. In
characteristic zero the extension is regular, giving geometric
generic integrality. This is the cited route recovered in
[THM-3827, generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases](../../01-canon/theorems/THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases.md);
it does not assert that every special fibre is irreducible.

At the triple point, `m=3` and `M` is a unit. Normal order
`j>=2` has three simple leading roots and `eta` a unit. For
`j=1` the leading cubic is `(a+bZ)^2+M0 Z^3`, with no triple
root. Its simple roots give a unit. A double root has generic
Morse critical-value order at most two: the derivative in the
fibre parameter is a nonzero unit times `w^2`. Order one gives
one ramified branch with `eta` a unit, and order two gives
nonzero logarithmic poles, already impossible under a mate.
Thus under the hypothetical mate all triple branches have
differential order zero. At the shared simple point, use
`v=calN` as local coordinate: the two branches have `v~s^2`
and `s^2 ds/E_v` a unit, so the order is also exactly zero.

The complete differential divisor now has four poles of order
two at the fourfold point and four zeros of order two on `D`.
It has no others. The infinity point of `S` is absent from the
generic curve since `N8!=0`; on `W\D`, generic smoothness and
the nonvanishing volume make the relative form a unit. Its
canonical degree is zero, and the connected compact normalized
generic curve therefore has genus one.

Let `E0` be the sum of the four points above `u=0`, each once.
A primitive has poles bounded by `E0`, of degree four.
Riemann--Roch gives `dim L(E0)=4`. A complete basis is

    1,  U=1/u,
    J=u(u-3)^2(u+1)t+u^2-(2p+5)u,
    V=(u-3)H/u.                                       (12)

The section `J` is global `L_1` by (1). All three functions
have at most simple poles at the fourfold point. At the triple
point `s~(u-3)^2`, so the factor `(u-3)^2` in `J` cancels
the pole of `t`, and `u-3` cancels the pole of `H`. This holds
both for the simple leading roots and the possible ramified
Morse branch. At the simple point `u+1=O(s)`, making `J`
finite, while `H=calN/s^2` is finite. At `D`, `J,H` are global,
`1/u` vanishes and `(u-3)/u` tends to one. The only possible
additional affine denominator is `u=0`; that line has constant
`F` and is absent from the generic fibre. These are complete
pole checks, not merely formal leading terms at one point.

Independence over the generic constant field follows from
`t` degrees: two uniquely selects `V`, one selects `J`, and
then `1,1/u` are independent. A nonzero polynomial in `t` of
degree below four cannot vanish modulo the geometrically
integral generic equation. The same argument after extending
constants pays independence geometrically; uniqueness descends
the coefficients of a rational mate back to `C(zeta)`.

## 7. Exact coefficient rows exclude every primitive

A mate restricted to the original generic fibre must have

    G=A(zeta)U+B(zeta)J+C(zeta)V+D(zeta).

Use the actual source bracket
`bracket(G)=u^2(F_u G_q-F_q G_u)` and reduce modulo `F-zeta`.
After scaling, the desired bracket is one. The coefficient
functions of `zeta` commute with `F` and do not add derivatives.

For `p!=0`, with (8),(10), the `q^3` row is

    4(u-3)^6(u+1)^2 [A-3(2p+3)B],

so `A=3(2p+3)B`. The resulting `q^2` row is

    (u-3)^4(u+1)(u+3)[4(h+192p^2)B+lambda C],

so `C=-4(h+192p^2)B/lambda`. Finally the coefficient of
`u^4 q` equals

    B [48(h+192p^2)zeta/lambda
        -(32p+48)(h+192p^2)-lambda].                  (13)

The bracket is nonzero in `C(zeta)`: either its fibre-parameter
coefficient is nonzero, or `h+192p^2=0` leaves `-lambda!=0`.
Thus `B=0`, then `A=C=0`, contradicting bracket one.

For `p=0`, use the full anchored family (11), with free `b,c`.
The same basis has `J=(u-3)^2(u+1)q/u-3-9/u`. The `q^3`
and `q^2` rows give, respectively,

    A=-b B/6,
    C=-(b^2+108c)B/(27lambda).

The remaining `q` coefficient is

    B (u-3)^3(u+1)
      [108(b^2+108c)zeta+2b lambda(b^2+108c)-243lambda^2]
      /(243lambda).                                  (14)

Its bracketed factor is nonzero for the same reason: if the
fibre slope vanishes, `-243lambda^2` remains. Again all three
coefficients vanish, which cannot produce bracket one.

This obstruction is in the complete Riemann--Roch space, not a
bounded mate-degree ansatz. Consequently every rational mate
in the all-finite class forces `L` constant. The polynomial
factor from Section 2 finishes the theorem.

## 8. Sharp rational boundary, controls, and stopping scope

There are actual rational mates in the same boundary partition
when `L` is constant, for every normalized finite `p`. Put

    A0=(u-3)(u+1),
    Z0=u^2(u-3)t+u-2p-3,
    H=A0 Z0^2+kappa,
    G_H=(u-2)/(12u(u-3)Z0).

This is a global `L_2` section with exactly the leading (2).
Directly `J(H,G_H)=1`. Hence for `F=H^2+c0`, the rational
function `G_H/(2H)` is a mate. Globality is checked in the
original second chart with `x=u+p`; no translation of the
surface is assumed. The hostile validates both the all-`p`
coefficient scope and the sharp polynomial/rational boundary.

The exact source checks the full symbolic section spaces,
inverse coefficient, all residue and coefficient reductions,
leading local equations, both complete primitive systems, and
the same-partition rational controls. The proof supplies local
branch completeness, geometric integrality and Riemann--Roch;
finite tests do not replace those arguments. Reproduce with

```sh
python3 -B 04-computation/planar_jc48_sep08_four_three_one.py
python3 -B -O 04-computation/planar_jc48_sep08_four_three_one.py
```

Both executions pass **107 gates** and produce the same 426 bytes
as [the frozen output](planar_jc48_sep08_four_three_one.out).
The [source](../../04-computation/planar_jc48_sep08_four_three_one.py)
is 14,646 bytes, SHA256
`2d4051d677c1c1ad72d603d87a8da881cd8bf30a10c38aa3438952e3a0c3ef8d`.
Output SHA256:
`ca712f951b47665fb2314accfa8a231bfc9d03b99cb895823529d7c9f0bffc34`.
Semantic gate digest:
`7919a3c37f7cfaf65c9124c94a1aa86765bcceb13276aa3a34e31d16f6095c14`.

The geometry sibling independently read the noncomposition,
complete local divisor, genus, basis poles and coefficient descent
and reported analytic PASS. It did not independently reconstruct
the three coefficient rows. Root's full proof/source audit has passed, including independently
reconstructed original-source coefficient rows. The primary is promoted.

If the triple or simple root is infinity, the degree-five
`(4,1)` or degree-seven `(4,3)` leading differential is already
excluded by the proved classification. If the fourfold point is
infinity, the remaining finite degree-four `(3,1)` differential
passes the leading gate. Its constant-`D` subclass remains a
separate OPEN target here. No projective-coordinate transport
or all-location closure is inferred from this all-finite result.

The [independent complete audit](planar_jc48_sep08_four_three_one_audit.md)
accepts all analytic steps and both107-gate replays. Direct original-source
Jacobian reductions independently recover both final primitive systems.
The same-partition rational control passes literally for every finite p.
