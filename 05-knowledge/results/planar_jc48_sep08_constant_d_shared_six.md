# Constant-D shared finite-six: complete polynomial exclusion and rational boundary

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a theorem about the complete global square-prefix class on the
fixed DG surface. It retains rational mates and an actual rational
submersion. No unrestricted quartic or JC(2) conclusion is claimed.

## 1. Statement, inheritance, and the complete global entry

Work over `C` on

    W=(P1_x x P1_z) minus {z=x^2},
    t=1/(z-x^2), x=1/r, t=-r^2-r^4 b_D,
    D={r=0}, omega=dx wedge dt=r^2 dr wedge db_D.

Let `H in L2`, `L in L1`, and `F=H^2+L`. Suppose the boundary
octic of `H` has a finite root of multiplicity six and the original
infinity root of multiplicity two, `F|D` is constant, and the boundary
section of `L` vanishes at that finite root.

**Theorem.** Such an `F` has no polynomial constant-Jacobian
mate, of any degree. After the necessary local entry reduction, the
table below completely classifies its rational mates' existence.

Actual source/surface scalings and a scalar normalization of `H` put
the leading polynomial in the form `N=u^6`, where `u=x-p` and the
original finite position `p` remains arbitrary. The notation `u`
is a polynomial coordinate change for calculations, **not** a claimed
automorphism of the compactification.

The complete constant-`D` coefficient rows are

    H=u^6 t^2+(beta0+beta1 u+b u^2+c u^3+2u^4)t
                       +(b+k)+cu+u^2,
    L=(m0+m1 u+d u^2+n u^3)t+(e+d)+nu.                 (1)

The shared-root hypothesis is `m0=0`. If a rational mate exists,
the inherited local lemmas force

    beta0=beta1=m1=0.                                 (2)

For the resulting active class put

    q=1+u^2 t,     v=uq,
    h=H=v^2+cv+bq+k,     L=dq+nv+e.                   (3)

Thus `k=c0-b`, `e=n0-d` in the original constant-row notation.
No parameter relation among `b,c,k,d,n,e` has been imposed.

The complete rational-existence table for (3) is:

| Disjoint coefficient case | A rational mate exists exactly when |
|---|---|
| `b*d!=0` | Never |
| `b=0,d!=0` | `2kc+n=0` |
| `b!=0,d=0,n!=0` | `c=k=0` |
| `b!=0,d=n=0` | `c=0` |
| `b=d=0` | Always |

Every row excludes polynomial mates. Outside (2), there is no rational
mate in the stated shared class.

The closest proved suppliers are the
[full global filtration](planar_jc48_sep08_dg_genus.md),
[shared-root regularity and first-jet obstruction](planar_jc48_sep08_shared_roots.md),
and the [constant-D finite-six M-unit theorem](planar_jc48_sep08_const_d_quartic.md).
That last theorem is kept separate: it assumes the lower section is
a unit at the finite root and is stated at its specified original
finite position. Its conclusion is not transported by translation.

For the final regularity obstruction, the direct antecedents are
[THM-3978 — linear-seam submersion and rational-mate pole obstruction](../../01-canon/theorems/THM-3978-linear-seam-submersion-rational-mate-pole-obstruction.md)
and [octuple transport, Section 5](planar_jc48_sep08_octuple_transport.md).
Their reusable mechanism is the incompatibility of one invariant
integration constant with two divisors of the same fibre. The
present polynomial and its two divisors are derived explicitly below;
no theorem is transported merely because a rational chart resembles
the old chart.

The live concept board is complete global rows, the actual volume,
residue tests in rational coordinates, source criticality, and
same-fibre pole repair. The canonical hostile is the rational
submersion in Section 6. The corrected near miss is to keep analyzing
outer residues after an inner conic-infinity logarithm already
excludes the whole `b*d!=0` case. The least-used sidecar is the
original retained source divisor where a rational correction fails.

## 2. Full coefficient completeness and the local entry

Here is an independent derivation of the global rows. In the complete
quadratic section box write

    H=[A(x)z^2+B(x)z+C(x)]/(z-x^2)^2,
    deg A,B,C<=4,
    N=A x^4+B x^2+C=(x-p)^6.

The degree constraints force `deg A<=2`. If
`A=q2 x^2+q1 x+q0`, the top terms of `B` are
`(1-q2)x^4+(-6p-q1)x^3`; its three remaining coefficients
are arbitrary. This gives all six parameters of the fixed-leading
section space. In `u` coordinates, before imposing constant `D`,

    P=sum_(i=0)^4 A_i u^i,
    Q=(A4-1)u^2+(A3-2p A4+4p)u+C0.                   (4)

The complete linear section box similarly gives

    M=sum_(i=0)^4 B_i u^i,
    R=B4 u^2+(B3-2p B4)u+E0.                         (5)

On the actual `D`, the slopes of `H|D` and `L|D` in `b_D`
are `2-A4` and `-B4`. The square of a nonconstant affine function
cannot be canceled by an affine function. Thus constancy of `F|D`
forces `A4=2`, `B4=0`, giving exactly (1). In particular,

    H|D=4p^2+2pc+k,      L|D=2pn+e.                  (6)

These generally depend on `p`. Their explicit values preserve the
original compactification while the active polynomial formulas (3)
happen to be independent of `p`.

Infinity contributes no primitive poles. To see this for the full
rows, use the actual surface inversion `r=1/x`,
`T=-x^2-x^4t`, and `s=1/T`. Its numerator has boundary value
`r^2(1-pr)^6`. Constant `D` implies that the quadratic part of
the numerator of `H` and the linear part of that of `L` have the
form

    N2=r^2+a r s+h_D s^2,    M1=l r+l_D s.

For the original generic fibre `F=zeta`, putting `s=r Z` gives

    (1+aZ+h_D Z^2)^2+l Z^3+(l_D-zeta)Z^4.             (7)

The constant is one and the leading coefficient is `F|D-zeta`.
The eliminant `ZP0'-4P0`, with `P0` obtained by removing the
`-zeta Z^4` term, has constant minus four. Hence all four roots
are nonzero and simple for generic `zeta`. They exhaust the
local Weierstrass degree four. The original relative form has
orders `2+2-3=1`, including the canonical `r^2`, on every
normalized branch. Higher terms and the two nonactive finite
jets do not affect this conclusion.

At the finite shared root, if `beta0!=0`, the normal-unit lemma
makes every branch regular. There would then be no primitive
poles anywhere on the compact generic components. This is
impossible: their nonzero relative form cannot be the differential
of a constant. If `beta0=0`, the finite shared first-jet lemma
forces `beta1=m1=0` under a rational mate. Otherwise it produces
a nonzero logarithmic residue. These two alternatives prove (2).
They use the actual finite point, where the canonical multiplier
is a unit, rather than the infinity version of the lemma.

Generic components meet the original source, and the relative
form is nonzero there. Poles of a rational mate's denominator
containing an entire special fibre may be discarded by taking
the generic fibre. No irreducibility or genus assumption is
needed for this entry argument.

## 3. The actual birational chart and its retained volume

The first chart in (3) has rational inverse

    u=v/q,       t=(q-1)q^2/v^2.

Moreover

    J_(u,t)(v,q)=u^2q=v^2/q,
    omega=(q/v^2) dv wedge dq.                       (8)

If `b!=0`, the pair `(h,v)` is also a rational field coordinate
system, because

    q=(h-v^2-cv-k)/b,
    omega=q/(b v^2) dv wedge dh.                     (9)

These equations preserve the source field and the exact relative
differential. They do not transfer polynomial regularity through
the denominators of the inverse map. Whenever regularity is used
below, it is checked in the original `(u,t)` plane.

## 4. Two nonzero coefficients: universal conic logarithms

Assume `b*d!=0`. Put `delta=d/b`, `alpha=n-delta*c`. Eliminating
`q` gives the complete original first coordinate

    F=h^2+delta h-delta v^2+alpha v+e-delta k.          (10)

On its generic fibre, (9) gives

    eta=-q dv/[b(2h+delta)v^2].                       (11)

Let `w=1/v`, `xi=h/v`. The complete equation at infinity is

    xi^2-delta+w(delta xi+alpha)
                    +w^2(e-delta k-zeta)=0.           (12)

At `w=0` there are two simple points `xi=epsilon*sqrt(delta)`,
`epsilon=+-1`, since `delta!=0`. The implicit derivative is
`2xi!=0`, so they give actual branches with local parameter `w`.
The square root is a nonzero complex constant, not an extension
of the moving fibre parameter. The form (11) is exactly

    [(xi-c)w-1-k w^2] dw
      /[b^2 w(2xi+delta w)].                         (13)

Its residues at these branches are

    -epsilon/(2b^2 sqrt(delta)),

both nonzero. Thus there is no rational mate, for any values of
`c,k,n,e`. No outer-residue condition, sample coefficient test,
or possible genus label is needed. The branches belong to the
actual field (9), so an exact rational differential would have
zero residues there.

## 5. Vanishing b: a residue becomes an actual critical line

First let `b=0`, `d!=0`. Define

    P(v)=(v^2+cv+k)^2+nv.

Now `F=P(v)+dq+e`, so `(F,v)` generates the full rational source
field. Keeping `F=f` fixed, (8) gives

    eta=-(f-P(v)-e)dv/(d^2 v^2).                     (14)

Its only possible finite logarithmic coefficient is
`(2kc+n)/d^2`. Rational exactness is therefore equivalent to

    2kc+n=0.                                        (15)

Sufficiency is explicit: under (15), the rational function

    G=[v^3/3+c v^2+(c^2+2k)v+(F-e-k^2)/v]/d^2        (16)

has Jacobian one. Thus this row is not a rational exclusion.

On the original source line `u=0`, however, `q=1`, `v=0`,
and direct differentiation gives

    F_t=0,       F_u=2kc+n.

The necessary condition (15) makes the entire source line critical.
There can be no polynomial mate. This argument checks the source
gradient itself rather than transferring regularity from `(F,v)`.

If `b=d=0`, then `F=P(v)+e` is a polynomial composition of
outer degree four. The nonconstant factor `P'(v)` in any
polynomial Jacobian excludes a polynomial mate. Nevertheless
`J(v,1/(2u^2))=1`, and hence

    G=1/[2u^2 P'(v)]                                (17)

is a rational mate for every choice of the coefficients. This
also handles `n=0` without any division by `n`.

## 6. Vanishing d: the rational submersion and the original-source repair failure

Suppose `b!=0`, `d=0`, `n!=0`. Then

    F=h^2+nv+e,     v=(F-h^2-e)/n,

so `(F,h)` generates the entire rational source field. The exact
relative form at fixed `F=f` is

    eta=[h-v^2-cv-k]dh/(n b^2 v^2),
    v=(f-h^2-e)/n.                                  (18)

For `rho^2=f-e`, its residues at `h=epsilon*rho` are

    epsilon*(nk+2c rho^2)/(4b^2 rho^3).               (19)

Since `f` is the generic parameter and `n!=0`, both residues
vanish exactly when `c=k=0`. In this case

    H=v^2+bq,       F=H^2+nv+e,
    G0=1/(2b^2v)-H/(n b^2)                          (20)

satisfies `J(F,G0)=1` identically in the original source field.
Every rational mate is `G0+R(F)` with `R in C(T)`: in the
coordinate field `C(F)(h)`, the Hamiltonian derivation is a
nonzero scalar multiple of `partial_h`, whose kernel is exactly
`C(F)`.

This is an actual polynomial submersion. Away from `v=0`, it
follows from the regular Jacobian identity with `G0`. On `u=0`
the derivative `F_u=n` is nonzero. On `q=0` (where `u!=0`),
`F_t=n u^3` is nonzero. These exhaust `v=uq=0`.

Its polynomial-mate obstruction is global. In original coordinates
the special fibre factors as

    F-(b^2+e)=u K(u,t),
    K=u(q^2+bt)(H+b)+nq.                            (21)

One has

    K(0,t)=n,
    K(u,0)=u^3+2bu+n,
    K(u,-1/u^2)=-b^2/u.                             (22)

Thus `K` is nonconstant and has an irreducible curve factor;
its zero set avoids both `u=0` and `q=0`. Along `u=0`, `G0`
has a simple pole with coefficient `1/(2b^2)`, and `F-(b^2+e)`
has order one. To cancel this pole, `R(F)` would need a pole
at `F=b^2+e`. Along any component of `K=0`, however, `G0` is
regular because `v=uq` is a unit, whereas that same `R(F)`
would have a pole. It cannot be canceled. Therefore no rational
mate is polynomial.

This is the same-fibre mechanism recovered from THM-3978 and
octuple transport, now with the actual new field kernel, fibre
factorization, and retained divisors paid explicitly. It also
explains why absence of source critical points would not suffice
for this class.

## 7. Constant L and the complete rational table

The remaining case with `b!=0` is `d=n=0`, so `F=H^2+e`.
A rational mate of `F` is equivalent to a rational mate of `H`:
multiply by `2H` in one direction and divide by `2H` in the
other. At fixed `H=h`, (9) gives

    eta_H=[1+c/v+(k-h)/v^2]dv/b^2.                  (23)

Its residue is `c/b^2`. Thus a rational mate exists exactly
when `c=0`; then

    G_H=[v+(H-k)/v]/b^2,
    G_F=G_H/(2H)                                    (24)

are explicit mates. For a polynomial mate of `F`, the nonunit
factor `2H` excludes every coefficient choice, including those
passing (23). Together with Sections 4–6, this proves every
row of the rational-existence table and the complete polynomial
exclusion.

The five cases are disjoint and exhaustive. No limiting argument
is used at `b=0`, `d=0`, or `n=0`; those denominators are handled
by separate exact coordinate fields. The lower constants and
the arbitrary original position `p` are retained throughout.

## 8. Exact controls, status, and reproduction

The producer checks the full shifted global box and constant-`D`
coefficient constraints, both actual rational chart inverses and
Jacobians, the full generic infinity face before the finite entry
reduction, all three residue mechanisms, explicit rational mates
on every accepted row, the polynomial factor obstructions, and
the complete original-source fibre and transverse-gradient
controls for the rational submersion.

All checks are symbolic in the free coefficients. They do not
replace the analytic case exhaustion, generic branch argument,
or field-kernel and divisor proof with a finite scan. The source
uses always-active exceptions rather than Python assertions.

```sh
python3 -B 04-computation/planar_jc48_sep08_constant_d_shared_six.py
python3 -B -O 04-computation/planar_jc48_sep08_constant_d_shared_six.py
```

Both normal and optimized replays pass **68 always-active gates** and
produce the same 392 bytes as
[the frozen output](planar_jc48_sep08_constant_d_shared_six.out).
The [source](../../04-computation/planar_jc48_sep08_constant_d_shared_six.py)
has 9,140 bytes and SHA256
`798b7824988fa6cbdcaaccd2689997d29ad57473e64206ca8b7f42a1e31d4aa6`.
Output SHA256:
`ecd0a552a4bde15e90bebc8dd5764f6105dc24d05463107e7d30404b048c410c`.
Semantic gate digest:
`27e6a0636da93ea56848bc4aeda7c8dd4871c78f04c6bde3d0ae2137b88f8f3d`.

Root supplied the complete active coefficient target and the explicit
nonconstant-lower-row rational submersion. The geometry sibling
independently derived the full case split, the direct conic logarithms,
the residue/critical-line identification, and the rational boundary
formulas, and reconstructed the same-fibre obstruction on the original
source. The exact source/output are frozen for independent audit.
The [independent complete audit](planar_jc48_sep08_constant_d_shared_six_audit.md)
accepts every coefficient case, all original-source divisors, the full
rational table and both68-gate replays. Combined with the independently
proved [M-unit transport](planar_jc48_sep08_const_d_translation.md), this
closes the entire constant-D finite6/infinity2 polynomial-mate class
at every finite location. The nonconstant-D class remains separate.
