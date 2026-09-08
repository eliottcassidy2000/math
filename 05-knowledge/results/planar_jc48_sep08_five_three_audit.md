# Independent audit: the all-finite (5,3) square-prefix class

**Status: PROVED / INDEPENDENT ANALYTIC, MATRIX, AND SOURCE AUDIT PASS.**
The geometry sibling accepts the complete
[five-three proof](planar_jc48_sep08_five_three.md) and
[producer](../../04-computation/planar_jc48_sep08_five_three.py).
The primary was `RESERVED` at acceptance. Root owns its promotion.
No producer file was edited, and the final repaired proof requires
no further mathematical or source correction.

## 1. Statement, full coefficients, and preserved coordinates

The accepted theorem concerns `H in L2`, `L in L1`, and `F=H^2+L`
on the fixed surface `(P1_x x P1_z) minus {z=x^2}`, with original
volume `dx wedge dt`, `t=1/(z-x^2)`. Its boundary octic has two
distinct **finite** roots of multiplicities five and three. A
rational Jacobian mate forces `L` constant. Thus no polynomial
mate of any degree exists. The primary's all-position constant-`L`
rational examples prevent a blanket rational exclusion. Neither
placement at infinity nor an arbitrary quartic is included.

Normalizing a nonzero leading scalar and the nonzero separation uses
actual scalings. The arbitrary original location `p` is retained in
`u=x-p`; a translation of the surface is not asserted. The complete
global rows, not a restricted prefix, are used in every bracket.

I independently recovered their coefficient structure from the full
numerator box `N=A(x)x^4+B(x)x^2+C(x)`, with `deg A,B,C<=4`.
For fixed `N` of degree eight and `P=2A x^2+B`, `Q=A`, one has

    P6=2N8, P5=2N7,
    Q4=N8, Q3=N7, Q2=P4-N6, Q1=P3-N5.

The remaining five coefficients of `P` and the constant of `Q`
are free. Conversely these relations reconstruct the full box.
The six-dimensional fixed-leading kernel is therefore complete.
The analogous `L1` row has `R2=M4`, `R1=M3`, with its constant
free. The source checks every relation in the original `x`
coordinates for symbolic `p`, including the two pre-active jets,
and verifies actual polynomiality on the second chart.

The coefficient `N8=1` has two important consequences beyond those
rows: the original `S` point at infinity is not a boundary zero,
and `H|D` has quadratic leading coefficient one. Hence `F|D`
is a monic quartic. Its four generic roots on the retained `D`
are simple; the source volume's `r^2` gives differential order
exactly two at each. These are four actual zeros of the relative
form, not four missing points or a new source-chart volume.

## 2. Necessary active entry and the complete moving residue

The local dependencies are the proved shared-root first-jet and
normal-unit lemmas, the M-unit Newton/Morse classification, and the
same-component primitive-degree argument. At multiplicity five an
M-unit is never balanced with `3j`, so its forms are regular.
A shared normal unit is also regular. A nonactive triple point is
regular or has forbidden nonzero logarithms. An active triple has
only total primitive pole capacity two.

If the fivefold point were nonactive, a generic `D` zero would force
primitive local degree three, larger than even the total available
capacity two. This also exceeds the capacity of that same component.
If there were no primitive poles, the primitive would be constant
on every compact component while its relative differential is
nonzero on the source part. Neither alternative needs an assumption
of generic connectedness. Thus the fivefold point must be active.
The finite first-jet theorem then forces `P(0)=P'(0)=M(0)=M'(0)=0`.

The complete inverse coefficient `T2=-M/(4N)` and normalized trace
make `M du/N` rational-exact. Its residue at zero is the negative
of `mu4+3mu3+6mu2`; vanishing gives the stated relation. This is
used before asserting that the triple point is M-unit.

I independently recovered the first two complete rows after setting
`w=u^2t`, with volume `u^-2 du wedge dw`. They agree with the
primary's `f(w)` and `g(w)`. Exactness along all simple roots of
`f(w)=zeta` requires

    g'f'-g f''=0.

For `b!=0`, `f'` is linear and the polynomial identity implies
`g=C f'`; the cubic coefficient of `g` is `-2b`, so this is
impossible. Therefore `b=0`. If `mu2!=0`, `f` is linear and
`g` must be constant, forcing first `c=0`, then `mu3=0`.
Together with the `T2` relation this gives precisely case I.
If `mu2=0`, the same relation gives case II unless `M=0`,
which is the separately handled constant-`L` case. This argument
does not invent outer roots when `f` is constant.

Only after this reduction do the triple values become the units
`-5lambda` and `-2lambda`. The two active lower rows are therefore
the complete residue survivors. Absorbing the additive `r0` into
the original fibre value is harmless; the source keeps `r0` in
the higher inverse-residue checks as an additional control.

## 3. Independent branch, zero, and genus inventories

I reconstructed the local rescalings with the complete original
equation `calN^2+s^3 calM-zeta s^4`.

In **case I**, `m=5`, normal order at least three, and lower order
two. Because the moving residue has already forced `Q(0)=0`,
the outer face `s=u^2 Z` is
`Z^3(lambda-zeta Z)` after removal of the higher-weight terms.
It supplies one simple outer branch. Its numerator `s^2 du`
and denominator have orders four and six, so the form has a
pole of order two. The high face uses
`u=tau^3`, `s=tau^8 Z` and has polynomial `a^2+lambda Z^3`.
Its three determinations form one actual branch, since the
ramification orders are coprime. Numerator and denominator orders
are `16+2=18` and `22`, giving a pole of order four. One plus
three exhausts the Weierstrass degree four. The primitive divisor
has orders one and three, total degree four.

If the M-unit triple has nonzero normal value `P(1)`, its
cancellation pair is one normalized branch with `u-1=tau^2`
and `s~tau^6`. The actual form has order four there. A primitive
would have local degree five, exceeding the entire case-I pole
capacity four. Thus case I necessarily has `P(1)=0`.

In **case II**, lower order three and normal order at least three
give `u=tau^2`, `s=tau^5 Z`. Before the later inverse coefficient
sets `Q(0)=0`, the full face is
`(a+Q(0)Z^2)^2+(R(0)-zeta)Z^4`; it has four simple nonzero
generic determinations. After that coefficient reduction the source
checks the corresponding simplified face. The four determinations
form two normalized branches. Their differential order is
`10+1-15=-4`, so the primitive divisor has orders three and three.

At an M-unit triple with `P(1)=0`, normal order at least two
gives unit forms. Normal order one is the balanced cubic case:
simple leading roots give units, while a double root has generic
Morse contact at most two. Contact one is a unit after quadratic
normalization; contact two gives nonzero logarithms and is excluded
under the proposed mate. If `P(1)!=0`, the sole boundary branch
instead has the zero of order four described above.

All remaining zeros are the four order-two points on `D`. On
the smooth generic fibre in `W\D`, the source volume is a unit
and the relative form has neither zeros nor poles. A nonzero
boundary value at the remaining `S` points prevents the generic
fibre from approaching them. Thus the inventories are complete:

| Case | Pole orders | Additional triple zero | Degree of canonical divisor | Genus |
|---|---|---|---:|---:|
| I | 2,4 | none | `8-2-4=2` | 2 |
| II, `P(1)=0` | 4,4 | none | `8-4-4=0` | 1 |
| II, `P(1)!=0` | 4,4 | order 4 | `8+4-4-4=4` | 3 |

For converting these degrees to genera, geometric integrality is
paid separately. A nontrivial polynomial composition of the
nonconstant-`L` quartic has outer degree two or four, since its
`t` degree is four. Outer degree four would make the polynomial
`N` square up to a scalar, contradicted by its odd root
multiplicities. After completing the outer square in degree two,
choose the leading sign so that the inner quadratic and `H`
have the same `t^2` coefficient. The factor `K+H` has `t`
degree two, whereas `(K-H)(K+H)=L-constant` has degree at most
one. Consequently `K=H`, forcing `L` constant, a contradiction.

I checked this reduction against the already recovered precise
noncomposite/closed-polynomial route: Arzhantsev–Petravchuk,
Theorem 1 and Lemma 3, as routed in
[THM-3827](../../01-canon/theorems/THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases.md).
Relative algebraic closure in characteristic zero gives geometric
generic integrality. No claim that every special fibre is integral
is needed.

The three primitive divisors consequently have Riemann–Roch dimensions
`4+1-2=3`, `6+1-1=6`, and `6+1-3=4`. In each case its degree
is strictly greater than `2g-2`, so no special-divisor correction
is omitted.

## 4. Complete primitive spaces and the affine-safe correction

The functions `J` and `J0` in the primary are genuine global `L1`
sections. Their induced linear terms are respectively
`-(2+2p)u` and `-(3+2p)u`, as required by their leading
polynomials `u^2(u-1)^2` and `u(u-1)^3`. This pays regularity
at `D` for every original `p`.

In **case I**, `U=1/u` has pole orders one and three at the outer
and high points. The function `J` is regular at the outer point
and has pole order two at the high point. At the triple, its
double factor cancels the possible order-two `t` pole. Thus
`1,U,J` lie in the required divisor space. The only affine
denominator is `u`; its line has constant `F` and is absent from
the generic fibre. Independence follows from `t` degree one,
then independence of `1,1/u`, in the degree-four integral generic
field. Their number equals the Riemann–Roch dimension.

In **case II with `P(1)=0`**, `H` is finite at both active
points. The functions `U,J,UJ,V,VJ`, with `V=(u-1)H/u`, have
pole bounds `2,1,3,2,3`. At the triple, `J` is regular and
the factor `u-1` makes the possible pole of `H` harmless.
All are regular on `D` and on the ordinary affine points of
`u=1`. The only remaining denominator is again `u`, absent
from the generic fibre. Their `t` degrees first isolate `VJ`
at degree three, then `V` at degree two, then `J,UJ`, whose
coefficients differ by the nonconstant factor `u`. Finally
`1,U` are independent. This proves that all six functions are
a basis, rather than just a low-degree ansatz.

For **case II with `P(1)!=0`**, I explicitly checked the
rejection of the earlier candidate `(J+d/(u-1))/u`.
At ordinary finite source points over `u=1`, `J` is polynomial
and `d=P(1)!=0`; the added term has a genuine pole. The generic
fibre has two such affine points, since its equation there is
a quadratic in `t` with nonzero leading coefficient `d^2`.
That candidate cannot belong to the primitive space, even though
it cancels a boundary principal part. It is not used in any final
matrix.

The replacement `W0=H*Adj` is valid. On the normal-unit boundary
branch, with `w=u-1`, one has

    t=-P/N-(H-Q)/P+O(N(H-Q)^2/P^3),
    H=(unit)w^(-3/2).

Thus `J0=K(u)+O(w^(3/2))`, with the primary's rational `K`.
I independently computed both its jets and the exact subtraction:

    K-k0+k1(1/u-1)=(u-1)^2(2p-u+1)/u.               (1)

Hence `Adj` has order at least `3/2` on that branch and `H*Adj`
is regular. The only affine denominator in `Adj` is `u`, so
this repair is also regular at every ordinary affine point over
`u=1`. It avoids the failed candidate's lost source addresses.

At the active points `J0` has pole order three, `U` order two,
and `H` is finite, so `W0` has at most order three. At `D`,
`H,J0` are global and `U` tends to zero. Consequently
`1,U,J0,W0` all belong to the required divisor space. Their
independence follows from the unique `t` degree three of `W0`,
then degree one of `J0`, and the final pair `1,U`. This is a
complete four-dimensional basis.

The coefficient arguments remain valid after extension of the generic
constant field: the generic polynomial is geometrically integral and
the basis functions have `t` degree below four. The determinants
below stay nonzero over that extension. No unproved rational descent
of individual branch labels or primitive coefficients is required.

## 5. Independent higher inverse coefficients and matrix reconstruction

I independently derived the formal pair-sum equation. For the centered
quartic `(Ny^2+D0)^2+My+E0=v^-4`, factor into
`(y^2+sigma y+qa)(y^2-sigma y+qb)`. Matching coefficients gives

    qa+qb=sigma^2+2D0/N,
    qb-qa=M/(N^2 sigma),
    qa*qb=(D0^2+E0-v^-4)/N^2.

Elimination, followed by `sigma=M v^2 V(v^4)/(2N)`, gives exactly
the primary's implicit equation for `V`. The solution with constant
one has the stated first two coefficients. The two opposite formal
inverse branches sum to `-sigma`; their even terms therefore yield

    T6=-ME0/(8N),
    T10=-M(3E0^2-D0 M^2/N)/(32N).

These coefficients are rational in the original base coordinate.
The formal mate identity makes each `T_k du`, for `k!=0`, exact
in the radical field; normalized trace descends the rational rows.
I independently recomputed the two residues in case II:

    Res T6 du=lambda^2(B+12p+2)/16,
    Res T10 du=lambda^3 c/32 after B=-12p-2.

Thus the coefficient constraints used before the final matrices are
necessary for every rational mate. They retain all original induced
lower coefficients.

For a separate matrix path, I differentiated the original `(u,t)`
polynomials **before** replacing `t=w/u^2`. I then reduced each
Jacobian by explicit leading-coefficient long division in `w`.
This verifier imported no producer code and did not call its
`bracket` or `sympy.rem` routine. Every resulting remainder was
independently polynomial in `u,w`; no denominator clearing changed
the displayed monomial addresses.

The independent path reproduced:

* Both case-I rows: the cubic row `4a u^2(u-1)^6`, then the
  linear row `-2b lambda u(u-1)^2` after its first coefficient
  vanishes.
* The genus-one determinant
  `192 lambda zeta(-lambda^2+32p^2 zeta)`.
* The complete genus-three determinant (16) and its separate
  anchored determinant (17), with exactly the primary's signs
  and factors.

All selected rows have zero target coefficient for the constant
bracket one. In case I the successive rows force both primitive
coefficients to zero. In genus one the determinant is nonzero in
`C(zeta)` for all `lambda!=0`, including `p=0`. In genus three
with `p!=0`, a vanishing fibre slope forces `A=8p+6`; substitution
leaves `4096 lambda p^3` in the determinant, nonzero. At `p=0`,
the alternative determinant has nonzero fibre slope `-240lambda`.
Thus all nonconstant primitive coefficients vanish in each case,
contradicting bracket one. This is a complete primitive-space
argument, not a finite degree bound on a proposed mate.

The independent matrix verifier passed 18 exact controls. Its
semantic check-label digest is
`23b1acb71e13fc5c4808375e6e1a85565b2496c67b5d602fda29ddf282401b08`.
Its separate ordinary Jacobian/long-division core is reproduced below;
the inputs are precisely the source polynomials and bases stated in
the primary, so all matrices can be reconstructed without importing
the producer.

```python
def rem_original(F, G):
    T = S.cancel(S.diff(F,u)*S.diff(G,t)-S.diff(F,t)*S.diff(G,u))
    T = S.Poly(S.expand(T.subs(t,w/u**2)),w)
    D = S.Poly(S.expand((F-zeta).subs(t,w/u**2)),w)
    coeff = {j:S.cancel(T.nth(j)) for j in range(T.degree()+1)}
    dcoeff = {j:S.cancel(D.nth(j)) for j in range(5)}
    for j in range(T.degree(),3,-1):
        ratio = S.cancel(coeff.get(j,0)/dcoeff[4])
        if ratio == 0:
            continue
        for k in range(5):
            coeff[j-4+k] = S.cancel(coeff.get(j-4+k,0)-ratio*dcoeff[k])
        if coeff[j] != 0:
            raise RuntimeError('long division top')
    ans = S.cancel(sum(coeff.get(j,0)*w**j for j in range(4)))
    if S.fraction(ans)[1] != 1:
        raise RuntimeError('nonpolynomial remainder')
    return S.Poly(ans,u,w)
```

## 6. Sharp rational boundary and final exact-source audit

I checked the constant-`L` sharp family independently in the
coordinate `Z0`. If `A0=u(u-1)` and
`S(u)=1/(3u^2)+5/(3u)+1/(u-1)`, its original Jacobian reduces to

    u^2(u-1)[-A0' S-2A0 S']=1.

This proves the displayed rational mate, with the all-position
globality verified by the complete second-chart identity. The leading
polynomial remains precisely `u^5(u-1)^3`. Dividing the mate of `H`
by `2H` gives the rational mate of `H^2+constant`; multiplication
by the nonunit `2H` still excludes every polynomial mate.

I read every final source gate. The full pre-active coefficient rows,
actual blow-up volume, residue cases, local faces and normalizations,
formal factor resolvent, higher trace residues, original affine-line
constancy, corrected jets, polynomial remainder checks, complete
matrices, and sharp examples agree with the proof. The source's
order-bookkeeping checks supplement, rather than establish, the
analytic branch exhaustion and Riemann–Roch arguments audited above.

Independent normal and optimized replays both completed:

```sh
python3 -B 04-computation/planar_jc48_sep08_five_three.py
python3 -B -O 04-computation/planar_jc48_sep08_five_three.py
```

Each passes **117 always-active gates** and reproduces the frozen
output byte for byte. The independent matrix reconstruction and inverse
residue calculations additionally agree.

| Artifact at acceptance | Bytes | SHA256 |
|---|---:|---|
| Producer source | 13,167 | `0ceac3b9198d66241883179a4fc0b7fc598f2e838af87561f42d1d966951dba4` |
| Frozen output and both replays | 437 | `3013be961d8e1c884adb298f858ef01e4d721d8b9521e08ffa4faf2d92e40c9b` |
| Primary before promotion | 16,578 | `5c27b0d355c8c547e4dde5b119d411d5520eb70142ad1a1c5854214f99f8421d` |

Producer semantic digest:
`a1fc2ea5ca623ebc16b2c697ee095b8c019938f119491a01e669a98e99a95bcf`.

**Final audit: PASS.** The all-finite rational conclusion, polynomial
exclusion, affine-safe genus-three repair, and exact rational boundary
are accepted. Infinity placements remain outside this theorem.
Source/output were preserved; root owns promotion and integration.
