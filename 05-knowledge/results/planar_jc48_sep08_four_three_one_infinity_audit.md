# Independent audit: the fourfold-infinity (4,3,1) boundary class

**Status: PROVED / INDEPENDENT ANALYTIC AND EXACT-SOURCE AUDIT PASS.**
This is the geometry sibling's independent audit of
[the primary proof](planar_jc48_sep08_four_three_one_infinity.md)
and its [exact producer](../../04-computation/planar_jc48_sep08_four_three_one_infinity.py).
The primary was `RESERVED` when read. This sidecar accepts its theorem;
the root agent owns status promotion and integration. No producer file
was edited, and no mathematical correction was required.

## 1. Accepted statement and precise domain

The source surface is the actual fixed

    W=(P1_x x P1_z) minus {z=x^2},
    t=1/(z-x^2), omega=dx wedge dt.

For complete global sections `H in L2`, `L in L1`, let `F=H^2+L`.
If the binary boundary octic of `H` has multiplicity four at the
original infinity point, and distinct finite roots of multiplicities
three and one, then a rational `G in C(x,t)` with `J(F,G)=1` forces
`L` constant. Consequently there is no polynomial mate of any degree.
The theorem does **not** exclude every rational mate: the primary
gives actual constant-`L` examples in the same binary partition for
every finite triple-point position.

Multiplying a nonzero constant Jacobian to normalize it to one is
harmless. Source/surface scalings normalize the finite separation and
the nonzero leading scalar. Writing `u=x-p` retains the arbitrary
original position `p`; no source translation is asserted to extend to
an automorphism of `W`.

The proof requires the entire global section rows, the original volume,
and generic normalized components. It does not classify arbitrary
quartics in `C[x,t]`, infer a Keller map from rational mates, or prove
JC(2). Its final all-location corollary additionally uses the separate
all-finite and leading-exactness suppliers. The local theorem audited
here does not depend on that integration corollary.

## 2. Full coefficient rows and the global degree drop

I independently recovered the rows from the full global section box,
rather than treating the displayed coefficients as a selected family.
Write

    H=[A(x)z^2+B(x)z+C(x)]/(z-x^2)^2,
    deg A,B,C <=4.

For fixed degree-four leading polynomial `N`, the identity
`N=A x^4+B x^2+C` forces `deg A<=2`. Put

    A=q2 x^2+q1 x+q0,
    B=-q2 x^4-q1 x^3+e2 x^2+e1 x+e0.

Then `C=N-A x^4-B x^2` has degree at most four, and
`P=2A x^2+B`, `Q=A`. This gives all six free coefficients,
with `[x^2]Q=[x^4]P` and `[x]Q=[x^3]P`. After writing `u=x-p`,
these are exactly

    P=sum_(i=0)^4 A_i u^i,
    Q=A4 u^2+(A3-2p A4)u+C0.

The complete `L1` box independently gives the corresponding row
`R=B4 u^2+(B3-2p B4)u+E0` for `M=sum B_i u^i`.
The producer verifies their actual second-chart polynomiality for
symbolic `p` and every free coefficient. My separate coefficient
calculation also checked the shifted relation directly.

For `x=1/r`, `t=-r^2-r^4 b_D`, the retained divisor `D={r=0}`
has `H|D` affine in `b_D` with slope `-A4`, and `L|D` affine
with slope `-B4`. Thus a constant restriction of `F` forces first
`A4^2=0`, then `B4=0`. In particular `deg M<=3`. No assumption
about the genus or number of generic components enters this step.

## 3. Independent local and componentwise degree audit

I read the actual weighted infinity supplier in
[four-four transport, Section 2](planar_jc48_sep08_four_four_transport.md),
the complete [M-unit classification](planar_jc48_sep08_quartic_common_root.md),
and the [shared-root local lemmas](planar_jc48_sep08_shared_roots.md).
Their local statements, rather than a stronger global conclusion from
one former boundary partition, are the dependencies used here.

The actual inversion gives `r=1/x`, `T=-x^2-x^4t` and
`omega=r^2 dr wedge dT`. The transformed leading polynomial is

    r^4(1-pr)^3(1-(p+1)r).

Its coefficient after removing `r^4` is a unit for every `p`, including
`p=0` and `p=-1`. The canonical `r^2` must be retained when interpreting
the normalized local form.

For the shared, normal-nonunit cases at infinity, I recomputed the
following branch balances from the full equation

    Phi=(r^4 a(r)+s B(r)+s^2 C(r,s))^2
         +s^3(D(r)+s E(r,s))-zeta s^4,
    eta=(unit) r^2 s^2 dr/Phi_s.

* If `ord B=1`, the two simple low determinations have `s~r`
  and weighted differential order one. The cancellation centre has
  `s~r^3`, and `N_s` has order one there. The generic order `ell`
  of `calM-zeta s` at that centre is between one and three, since
  the last term has a nonzero coefficient at order three. Thus
  `ord N=(9+ell)/2`, `ord Phi_s=(11+ell)/2`, and the weighted
  form is a unit times `r^((5-ell)/2)dr`. It is regular after the
  actual possible quadratic ramification. The centre displacement
  has strictly higher order and does not invalidate this count.
* If `ord B>=2` and `ord D=1`, there is one simple low root and
  three high determinations with `(ord r,ord s)=(3,7)`. In the
  latter normalization, the orders of the numerator, including
  `dr`, and denominator are respectively `22` and `17`. The
  differential therefore has order five.
* If both first jets vanish, `s=r^2 Z` gives the complete tangent
  quartic in the primary. Its nonzero constant and the nonzero
  constant of `ZP0'-4P0` give four simple nonzero generic roots.
  The weighted form has order zero on all four branches.

Normal-unit and M-unit cases are independently covered by the cited
lemmas; multiplicity four cannot satisfy the balanced M-unit equation
`4=3j`. The local Weierstrass degrees pay exhaustion in every case.
Higher terms of the leading unit have higher face weight. The generic
parameter argument is explicit, not an inference from one specialization.
Therefore infinity contributes no primitive poles.

The finite simple point is regular whether shared or M-unit. A
nonactive finite triple point is regular or has a forbidden nonzero
logarithmic residue; the balanced M-unit contact bound is at most two.
At a shared triple point with normal derivative zero, the finite
first-jet lemma forces

    P(0)=P'(0)=M(0)=M'(0)=0

under a rational mate. Thus either the mate is already impossible,
or this is the only possible active point.

For the active triple, I independently substituted `u=tau^2`,
`s=tau^3 Z`. The leading quartic is

    (a+C(0)Z^2)^2+(R(0)-zeta)Z^4.

Its four generic nonzero simple determinations are paired by the
ramification into two actual normalized branches. On each, `s^2 du`
has order seven and `Phi_s` order nine, giving differential order
minus two. Hence each possible primitive pole has order one and the
total possible primitive pole degree is two. This uses the complete
degree-four local equation; no additional branches are discarded.

Outside the boundary roots the generic compactification adds no
further poles. On the original smooth source the generic fibre is
smooth and the relative form is regular. If no active point exists,
a rational primitive would be constant on each compact normalized
component, although the relative form is nonzero on its source part.

If `F|D` is nonconstant, a generic level has a finite transverse point
on `D`. The original `r^2` volume factor gives the relative form order
exactly two there. Its primitive therefore has local degree three on
that component. This exceeds even the **total** pole capacity two,
so also the capacity on that same component. The argument remains
valid for reducible or geometrically split generic fibres; it never
replaces a componentwise degree by a total degree lower bound.
It follows that `F|D` is constant and the coefficient reduction in
Section 2 applies.

## 4. Formal exactness and the decisive two residues

I checked the full square-prefix inverse coefficients. In the field
`C(x)(a)`, with `a^2=N`, solve `F(x,T)=v^-4`. Including the preceding
coefficients gives the exact `v^3` equation

    4 a T2+M/a=0,    T2=-M/(4N).

The formal identity for a rational mate is
`Gtilde_x=(1/4)v^5 T_v`; its `v^6` coefficient satisfies
`g6'=T2/2`. Normalized field trace commutes with the derivation and
therefore supplies a rational primitive of `M dx/N` in `C(x)`.
No assertion about convergence is needed, and no lower-row term is
lost. This agrees with the independently audited coefficient in the
all-finite companion.

The active and constant-`D` conditions leave

    M=u^2(B3 u+B2),
    M/N=(B3 u+B2)/(u(u-1)).

The residues at `u=0,1` are respectively `-B2` and `B2+B3`.
Solving both equations directly gives `B2=B3=0`; this was also
recomputed separately from the producer. The primary's sequential
version first gets `B2=-B3`, then the remaining differential `B3 du/u`,
and reaches the identical conclusion. Keeping only the first residue
would leave a spurious residual family.

Thus `M=0`, and full globality makes `L` constant. Finally
`J(H^2+constant,G)=2H J(H,G)` cannot equal one for polynomial `G`,
because `H` is nonconstant. This final step does not exclude rational
functions having the necessary poles.

## 5. Independent trace mechanism and sharp controls

The primary explicitly labels the intermediate trace proof as
redundant for the short closure. I checked it independently because
it retains useful source-field information.

For its residual `H,Z`, direct substitution gives
`A u^2+B u+z^2=0`, `y=2Au+B`, `y^2=B^2-4Az^2`, and
`J(H,Z)=u y`. The inverse formula for `t` is rational. The Jacobian
has a nonzero coefficient independent of the free parameters, so
`h,z` are algebraically independent. The monic quadratic in `h`
has discriminant `16z^2(z-1)(z+a+b-1)`, a nonzero polynomial for
every `a,b,c`; hence the ambient field degree is exactly two, not
one. Even coinciding roots of that displayed discriminant in `z`
do not turn the quadratic in `h` into a square.

Replacing `z` by `(zeta-h^2)/lambda` is a birational base-field
change for `lambda!=0`. It preserves the degree-two field extension.
The relative form follows from the original Jacobian, with sign

    eta=-dh/(lambda u y)
       =(1+B/y)dh/(2lambda z^2).

Tracing under `y -> -y` gives `lambda dh/(zeta-h^2)^2`.
After adjoining only the generic constant `rho^2=zeta`, its residue
at `h=rho` is `-lambda/(4rho^3)`, nonzero. Trace of an exact
differential is exact, so this supplies a second obstruction for
that intermediate nonconstant-row family without a genus or
geometric-integrality hypothesis.

For the sharp control, `H=u^3(u-1)t^2+c` is global in the original
shifted chart and `G_H=-1/(u^2 t)` satisfies `J(H,G_H)=1`.
The primary's second-chart polynomial verifies the boundary partition
for every `p`. Dividing this mate by `2H` gives a rational mate for
`H^2+c0`. These examples pay the strict distinction between the
constant-`L` necessary condition and the polynomial-only exclusion.

## 6. Exact-source audit, reproduction, and freeze

I read the complete producer, including all conditions in the local
case table and all shifted coefficients. Its **87 always-active
gates** are exact symbolic identity and boundary controls, not a
finite parameter census used to infer local completeness. The proof
supplies branch exhaustion, genericity, and the componentwise degree
argument analytically. The source retains the normal-unit, M-unit,
first-normal-jet, first-M-jet and higher-jet alternatives; both
residues; actual quadratic-field dominance and degree; and rational
controls for arbitrary original `p`.

Independent normal and optimized replays completed successfully:

```sh
python3 -B 04-computation/planar_jc48_sep08_four_three_one_infinity.py
python3 -B -O 04-computation/planar_jc48_sep08_four_three_one_infinity.py
```

Both reproduced the frozen output byte for byte. I additionally
derived the global coefficient constraints from the full numerator
box and solved the two unconditioned residue equations independently.
A targeted correction search in the current mistakes ledger found
no applicable correction to the inherited local or inverse-series
suppliers.

Frozen files at acceptance:

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Producer source | 10,605 | `b4bf39bad5a819d2435379795d4bee97fbbd6c79e271b507192bd4d3217df77c` |
| Frozen output | 399 | `92a0617b4594bc84c293a743ac79d6921310400924f06e8341c7894363e05268` |
| Primary proof, before promotion | 12,754 | `1a612853e78919661ff8da9f5dac816eeb03cb75cc3d787dcbee7455e2d9d911` |

The output's semantic gate digest is
`a747eb90da21add105fe24156322e989d0f302a74337e88862db4aee6ab5f057`.
No producer edit or Git mutation was made during this audit.

**Final finding: PASS.** The fourfold-infinity theorem and its sharp
rational boundary are accepted. No remaining mathematical or replay
repair is required.
