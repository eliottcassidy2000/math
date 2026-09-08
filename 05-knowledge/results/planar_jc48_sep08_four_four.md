# The distinguished finite/infinity 4+4 quartic boundary class is excluded

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This concerns the actual DG global section class with boundary octic
`alpha*x^4`, `alpha!=0`. Its two roots are the distinguished finite
point zero and infinity. No translation or arbitrary-pair normalization
on the surface is asserted.

## 1. Full entry and the source/target contract

Keep the fixed surface and source chart

    W=(P1_x x P1_z) minus S,   S={z=x^2},   t=1/(z-x^2),
    omega=dx wedge dt.

Let `H in L2`, `L in L1` be global, `F=H^2+L`, and suppose
its complete boundary octic is `N|S=alpha*x^4`, `alpha!=0`.

**Theorem.** No polynomial `G in C[x,t]`, of any degree,
has `J(F,G)` a nonzero constant.

The full entry, paid by the rank-nine octic restriction map and the
six-dimensional `L1` kernel, is

    w=x^2*t,   v=x+x^3*t,   h=x^2+x^4*t,
    H=alpha*w^2+k0+kt*t+kxt*x*t+k2*w+k3*v+k4*h,
    L=l0+lt*t+lxt*x*t+l2*w+l3*v+l4*h.                 (1)

The closest mechanisms are the complete local Newton classification in
[common-root necessity](planar_jc48_sep08_quartic_common_root.md), the
finite logarithmic first-jet obstruction in
[shared roots](planar_jc48_sep08_shared_roots.md), and the exact global
section spaces of [the DG filtration](planar_jc48_sep08_dg_quadratic.md).
The [infinity-octuple proof](planar_jc48_sep08_infinity_octuple.md)
suggests retaining the actual volume weight; none of its
conclusions is a dependency of the proof below.

The live concepts are weighted infinity branches, residues at the finite
point, the full `L1` coefficient space, critical points in the source
chart, and the distinction between polynomial and rational mates.
The map to a compact generic fibre preserves a proposed mate's exact
relative differential. A genus or pole-order table alone loses the
residues. The second map, `(x,t)->(x,q=1+x^2*t)`, is used only at
`x!=0`, where its Jacobian is `x^2`; this restriction is essential.

If `L` is constant, `J(F,G)=2H*J(H,G)` excludes polynomial mates
immediately. Henceforth assume `L` nonconstant, and normalize a
hypothetical polynomial mate to `J(F,G)=1`.

## 2. Complete weighted regularity at the infinity root of order four

The actual surface involution `(x,z)->(1/x,1/z)` gives

    u=1/x,  T=-x^2-x^4*t,   omega=u^2 du wedge dT.

It preserves the `4+4` boundary condition and the entire affine family
(1), with changed coefficients. For example the six basis elements
`(1,t,xt,w,v,h)` become `(1,-h_new,-v_new,-1-w_new,-uT,-T)`.
The base `alpha*w^2` becomes `alpha*w_new^2+2alpha*w_new+alpha`.
These are literal identities, not a source-plane polynomial automorphism.
In the new second chart the weighted form is `dr wedge db`.

It suffices to analyze the new finite point `(u,s)=(0,0)`, where
`s=zeta-u^2`. Write the new coefficients again as in (1):

    N=alpha*u^4+s(kt+kxt*u+k2*u^2+k3*u^3+k4*u^4)
      +s^2(k0+k3*u+k4*u^2),
    M=lt+lxt*u+l2*u^2+l3*u^3+l4*u^4
      +s(l0+l3*u+l4*u^2),
    E=N^2+s^3*M-c*s^4,
    eta=u^2*s^2 du/E_s.                              (2)

All generic local normalized branches have regular `eta`, for every
coefficient value with `alpha!=0`. Here is the full case split.

If `kt!=0`, the normal-unit local lemma gives two branches and
unweighted regularity for arbitrary `M`; the weight improves it. If
`kt=0,lt!=0`, the `M`-unit lemma has `m=4`. The balanced case
`4=3j` is impossible for an integer normal order. Its one low plus
two cancellation determinations (`j=1`), or three dominant cubic
determinations (`j>=2`, including infinite order), are regular.
These are only local uses of the proved lemmas.

Now `kt=lt=0`, so

    E(0,s)=(k0^2+l0-c)*s^4.                          (3)

Thus the exact local Weierstrass degree is four.

If `kxt!=0`, there are two low branches `s=u*Z+...` from

    Z^2[(kxt+k0*Z)^2+lxt*Z+(l0-c)*Z^2].

Their two nonzero roots are simple for generic `c`; their unweighted
forms are logarithmic, so the weight makes them regular of order one.
The other two determinations cancel `N` at
`s=-(alpha/kxt)u^3+...`. At that analytic centre put
`ell=ord_u(M-cs)`. Generically `1<=ell<=3`: lower terms remain,
and otherwise the coefficient of `c` supplies a nonzero term at order
three. Since `N_s` has order one, the weighted relative form has
exponent `(5-ell)/2` in `u`. Its actual normalized order is positive,
even when the two determinations form one ramified branch. All four
sheets in (3) have now been accounted for.

If `kxt=0,lxt!=0`, there is one low branch from
`Z^3[lxt+(k0^2+l0-c)Z]` at `s=uZ`, again weighted-regular.
The other three determinations have `s` of order `7/3`, with the
simple cubic face `alpha^2+lxt*Z^3`. On the actual normalization
`u=tau^3`, `s=tau^7` times a unit. The leading order of `E_s` is
17, and the numerator `u^2*s^2*du` has order 22. Hence `eta` has
order five. One plus three exhausts (3).

Finally let `kxt=lxt=0`. At `s=u^2Z` the complete face is

    P(Z)=(alpha+k2Z+k0Z^2)^2+l2Z^3+(l0-c)Z^4.        (4)

It has four simple nonzero roots generically. Indeed `P=P0-cZ^4`,
`P0(0)=alpha^2`, and any repeated nonzero root satisfies
`Z*P0'-4P0=0`, a nonzero polynomial with constant term `-4alpha^2`.
There are only finitely many exceptional `c` values. The generic
leading coefficient of `P` is nonzero as well. Every branch has
`E_s` of order six, exactly the order of `u^2*s^2`; its form is
regular. This completes the local infinity proof without restricting
`M` or dropping higher coefficients.

## 3. The finite point forces a full coefficient residue identity

There are only two boundary points. If at the finite point `lt!=0`,
the same `m=4`, `M`-unit lemma makes the unweighted form regular.
If `kt!=0`, the normal-unit lemma does so. Section 2 pays the other
boundary point. On every compact generic component the form would be
holomorphic and nonzero, which cannot be a rational derivative.
Indeed no generic component equals `S` or `D`; each meets the original
source where the volume form is nonzero. Generic fibres are smooth there
in characteristic zero, so their relative form is nonzero. A rational
primitive of a holomorphic differential on a compact component has no
poles and is constant. Thus `lt=kt=0`.

If `kxt!=0`, either of the two simple nonzero low roots in the
quadratic of Section 2 gives an unweighted simple pole with nonzero
residue. If `kxt=0,lxt!=0`, the low simple root there gives such a
pole. Equivalently, this is the already proved finite first-jet
obstruction. Hence

    lt=kt=kxt=lxt=0.                                (5)

In the birational chart `(x,w=x^2t)`, the original form is
`x^-2 dx wedge dw`. The entire function has expansion

    F=f(w)+x*g(w)+O(x^2),
    A(w)=alpha*w^2+k2*w+k0,
    f=A^2+l2*w+l0,
    g=(1+w)(2k3*A+l3).                              (6)

For each generic root `f(w0)=c`, `f'(w0)!=0`, the actual branch is
`w=w0-g(w0)x/f'(w0)+O(x^2)`. On this branch

    eta= -dx/(x^2 F_w),
    Res eta = (g'f'-f''g)(w0)/f'(w0)^3.              (7)

As `c` varies, vanishing of every such residue implies the polynomial
identity `g'f'-f''g=0`, equivalently

    g=C*f'                                         (8)

for a scalar `C`. This uses the original fibre and all of its first
coefficient corrections. The higher `O(x^2)` rows do not affect (7).
The face (4), now unweighted, also shows these are all four generic
local branches; no discarded branch is needed to obtain (8).

The full coefficient comparison in (8) is short but restrictive.
Since `alpha!=0`, `f'` has degree three.

- If `k3=0`, leading coefficients force `C=0`, then `l3=0`.
  Under (5), the original polynomial `F` has both partial derivatives
  zero along the source line `x=0`. This excludes polynomial mates.
- If `k3!=0`, the degree-three coefficient gives
  `C=k3/(2alpha)`. The degree-two, degree-one, and constant
  coefficients successively give

      k2=2alpha,    l3=0,    l2=0.                  (9)

  Since `L` is nonconstant, its remaining coefficient `l4` is nonzero.

The condition (8) is necessary, not sufficient for exactness. The
remaining cases are excluded by actual source critical points next.

## 4. Every residue-compatible nonconstant-L remainder is critical

Put `q=1+w`, `c0=k0-alpha`. Equations (5) and (9) give

    B=alpha*q+k3*x+k4*x^2,
    H=c0+q*B,
    F=H^2+l4*x^2*q+l0,       alpha*k3*l4!=0.         (10)

We work only at `x!=0`, where `(x,q)` is an actual regular source
coordinate chart, with inverse `t=(q-1)/x^2` and Jacobian `x^2`.
A critical point of (10) in this chart is therefore a critical point
of the original source polynomial.

On `q=0`, the criticality condition is

    2c0*k3+(2c0*k4+l4)*x=0.                         (11)

If `c0!=0` and `2c0*k4+l4!=0`, this has a nonzero solution, giving
an immediate actual critical point. Two exceptional coefficient cases
remain, and both can be solved without a parameter scan.

In both cases set

    q=-k3*x/(4alpha),   Cx=k3+2k4*x.

Then `F_x` and `F_q` vanish together precisely when
`H*Cx+l4*x=0`, for `x*q!=0`.

If `c0=0`, choose any complex root of

    k3*x*(3k3+4k4*x)*(k3+2k4*x)=16alpha*l4.          (12)

This is a nonconstant polynomial equation even when `k4=0`. Its
nonzero right-hand side ensures `x!=0`, `q!=0`, and both remaining
linear factors are nonzero. Substitution gives `H*Cx+l4*x=0`.
Thus it supplies an actual source critical point for every `k4`.

If `c0!=0` and `l4=-2c0*k4`, then `k4!=0`. Choose a root of

    x^2*(3k3+4k4*x)*(k3+2k4*x)=16alpha*c0.           (13)

Again the equation is nonconstant and its nonzero right-hand side
forces all necessary factors nonzero. Direct substitution gives
`H*Cx+l4*x=0`, hence an actual source critical point. These two
cases and (11) exhaust all coefficients in (10).

Every critical point contradicts `J(F,G)=1`. This proves the
polynomial-mate exclusion for the full distinguished `4+4` class.

## 5. Boundaries and verification

The local residue-free condition alone survives a named hostile:
`alpha=k3=l4=1`, `k2=2`, `k0=1`, `k4=0`, all other nonconstant
coefficients zero. It satisfies (8), yet (12) gives the explicit source
critical point

    x=16/3,    q=-4/3,    t=-21/256.

Thus residue cancellation does not supply a Keller mate. There is also a genuine rational-mate hostile with nonconstant `L`:

    H=w^2,    L=w,    F=w^4+w,
    G=1/[x(4w^3+1)],    J(F,G)=1.

Indeed `J(w,1/x)=1`. This object satisfies the residue identity with
`g=0` and is source-critical along `x=0`, while its mate has poles.
Thus the final conclusion must remain polynomial even with nonconstant
`L`; no rational-mate exclusion is claimed for the whole class.

No actual global affine-pair realization, general quartic closure, or
normalization of arbitrary two boundary points is inferred. The weighted
infinity lemma does apply to every coefficient in this fixed section
space. Other root partitions remain a separate question.

## 6. Exact controls and reproduction

The new standalone [source](../../04-computation/planar_jc48_sep08_four_four.py)
and [frozen output](planar_jc48_sep08_four_four.out) pass **96 always-active
exact gates**, with ordinary and optimized Python outputs byte-identical.
No inherited mathematical implementation is imported.

The exact universe is the entire six-dimensional global correction space
for both `H` and `L`, all unit/first-jet strata of the local multiplicity-four
model, the complete generic tangent quartic, and every coefficient case
of the residue-compatible criticality equations. The source checks actual
inversion and volume, section-space rank, the local degree-two/three/four
suppliers, all displayed leading faces and derivative orders, the moving
root's residue including its displacement, every coefficient of (8), and
both exceptional critical-point identities. It includes named generic
quartic faces, the explicit residue-compatible critical point, and the
nonconstant-`L` rational-mate hostile.

The analytic arguments prove branch exhaustion and existence of a complex
root in (12) and (13); finite tests do not replace those steps. Likewise,
the local pole/residue argument retains the original fibre throughout.

Reproduce from the worktree root:

```bash
python3 04-computation/planar_jc48_sep08_four_four.py
python3 -O 04-computation/planar_jc48_sep08_four_four.py
```

Frozen SHA256 pins:

- Source, 10,853 bytes:
  `f03a0843dd2521a35b0df92481a05db8af0c8130e210a5a5fb315afb9d5e8089`.
- Output, 405 bytes:
  `a3c7c6823ea672b09a459803e75111f1392e810d13d8897d10882724e92ccb90`.
- Semantic digest:
  `e6de99bdcde1e30ca8f3d177c22ec9b24f73de73ade3eaff6a5e47f1f7e49eae`.

The [independent full analytic/source audit](planar_jc48_sep08_four_four_audit.md)
accepts all coefficient strata and all96 normal/optimized gates. It separately
reconstructs the critical-line elimination and both exceptional equations.
No canonical ID is assigned; the scope remains the distinguished pair.
