# A finite/infinity 4+4 boundary pair is excluded at every finite point

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The new argument handles a nonzero finite point. The distinguished zero
point is the separate [4+4 theorem](planar_jc48_sep08_four_four.md).
That distinguished theorem is proved and independently audited, so the
combined polynomial-mate exclusion holds for every finite p. No translation is asserted to be an automorphism of the surface.

## 1. Full shifted global entry

Fix `W=(P1_x x P1_z) minus {z=x^2}`, with source
`t=1/(z-x^2)`. Let `H in L2`, `L in L1`, and `F=H^2+L`.
Suppose the complete boundary octic is

    N|S=alpha*(x-p)^4,      alpha!=0,  p!=0.

**Theorem.** No `G in C[x,t]`, of any degree, has
`J(F,G)` a nonzero constant.

Put `u=x-p`. The full corrected global basis is

    w=u^2*t,   v=u(1+w)-2p,   h=u^2(1+w)-2pu=u*v,
    K=k0+kt*t+kxt*u*t+k2*w+k3*v+k4*h,
    L=l0+lt*t+lxt*u*t+l2*w+l3*v+l4*h,
    H=alpha*w^2+K.                                   (1)

The six basis functions are global on the original `W`. Their section
coefficient determinant is a nonzero constant for every `p`; the octic
restriction has rank nine, with exactly this six-dimensional kernel.
These are the same complete section-space mechanisms as in the
[all-finite octuple transport](planar_jc48_sep08_octuple_transport.md).
The present base section is `alpha*w^2`, not `alpha*h^2`.

The source-to-target map first takes the actual relative differential
on the compact generic fibre, retaining the moving-root residue and all
corrected section coefficients. It then passes to a particular critical
curve of a polynomial map `(P,z)`. That second map need not be invertible;
its zero Jacobian is the supplier of simultaneous criticality. The only
chart inversion used to get a source point is `t=(q-1)/u^2` at `u!=0`.
Three forbidden values in a normalized parameter encode exactly the
lost chart addresses. A polynomial root-count argument pays those values.

Closest mechanisms are the [distinguished 4+4 residue/criticality proof](planar_jc48_sep08_four_four.md),
full shifted global sections, and the generic boundary differential
lemmas in [shared roots](planar_jc48_sep08_shared_roots.md).
The corrected near miss is to translate the surface itself. The live
concepts are actual section completeness, weighted infinity, residues,
a critical-curve parameter, and the integer spectrum of blocked roots.

As before, constant `L` is excluded directly by the nonunit factor `2H`
in `J(H^2+L,G)`. Hence assume `L` nonconstant.

## 2. The infinity supplier survives the actual unit factor

Under the actual surface inversion `(x,z)->(1/x,1/z)`, put
`r=1/x`, `T=-x^2-x^4t`. The original volume is `r^2 dr wedge dT`,
and the boundary octic becomes

    alpha*r^4*(1-pr)^4.

Near `r=0` this is an order-four root times a unit. The weighted
regularity proof for the distinguished case still applies, for the
following explicit reason. In local coordinates `(r,s)`, write

    N=r^4 a(r)+s B(r)+s^2 C(r,s),    a(0)=alpha!=0,
    M=D(r)+s E(r,s),
    eta=r^2 s^2 dr/[N^2+s^3M-cs^4]_s.                (2)

Normal-unit and `M`-unit cases are covered by the proved local lemmas.
Otherwise the local degree is four. If `ord B=1`, there are two low
branches of order `s~r`, and two cancellation determinations `s~r^3`.
The actual generic order of `M-cs` at the cancellation centre is at most
three; their weighted exponents are `(5-ell)/2>0`.
If `ord B>=2` but `ord D=1`, there is one low branch and three
`(ord r,ord s)=(3,7)` determinations, whose weighted form has order
five. Finally when both first jets vanish, the entire tangent face is

    (alpha+B2 Z+C0 Z^2)^2+D2 Z^3+(E0-c)Z^4,
    s=r^2Z.

It has four simple nonzero roots generically, since the repeated-root
supplier `ZP0'-4P0` has constant term `-4alpha^2`. Its weighted
forms have order zero. The higher terms of `a(r)` enter at strictly
higher weight in every displayed face. They do not change the generic
centre bound, where the term `-c*s` supplies the order. Thus every
normalized branch at the original infinity point is regular.

At the finite point `u=0`, the exact boundary value is `alpha*u^4`.
The same unweighted regularity and finite first-jet residue argument
therefore forces

    lt=kt=kxt=lxt=0.                                 (3)

This conclusion uses the two actual boundary points together in the
regular cases; a compact generic component meets the original source,
where the relative form is nonzero. In the first-jet cases the local
nonzero residue alone prevents exactness.

## 3. All shifted residue coefficients and the remaining polynomial map

In the chart `(u,w=u^2t)`, put

    A=alpha*w^2+k2*w+k0-2p*k3,
    f=A^2+l2*w+l0-2p*l3,
    g=2A[k3(1+w)-2p*k4]+l3(1+w)-2p*l4.              (4)

The complete function is `F=f(w)+u*g(w)+O(u^2)`, and the actual
volume is `u^-2 du wedge dw`. Consequently the original moving-root
residue is `(g'f'-f''g)/f'^3`. A proposed mate forces

    g=C*f'                                          (5)

with constant `C`, since generic values of `f(w)` cover a dense set.
All induced first-row terms in (4) are retained.

If `k3=0`, the degree comparison makes `C=0`, hence `g=0`.
Then the original `F_u` and `F_t` both vanish at every source point
`u=0`, excluding a polynomial mate.

If `k3!=0`, the complete coefficient comparison gives

    C=k3/(2alpha),
    k2=2alpha-4alpha*p*k4/k3,
    l3=0,
    l2=-4alpha*p*l4/k3.                              (6)

Nonconstant `L` now means `l4!=0`. Define

    a=4alpha/k3,    q=1+w,    D=p*a,
    z=(u^2-D)q-2pu,
    P=alpha*q^2+k3*u*q,
    c0=k0-alpha-2p*k3+p*a*k4.

After (6), the entire function, with no row omitted, is

    F=(P+k4*z+c0)^2+l4*z+(l0-l2).                   (7)

The exact Jacobian in `(u,q)` is

    J(P,z)=u[2p*k3-q(k3*u+4alpha*q)].                 (8)

We use the nonzero-`q` branch of this critical curve:

    u=2p/q-a*q,
    P=2p*k3-3alpha*q^2,
    z=a^2*q^3-3p*a*q.                               (9)

Write `H=P+k4*z+c0`. Along (9), set

    R(q)=H[k3*q+2k4(p-a*q^2)]+l4(p-a*q^2).          (10)

Then the literal two partial derivatives are

    F_u=2R,
    F_q=(4p/q^2-a)R.                                (11)

Thus every root of `R` with `q!=0` and `u!=0` supplies an actual
critical point of the original source polynomial. The coordinate
change `(u,t)->(u,q=1+u^2t)` has Jacobian `u^2`; no conclusion is
transported through its zero locus.

## 4. The critical polynomial cannot have only forbidden roots

Choose either nonzero square root `s` of `2p/a`, and put `q=sQ`.
The forbidden values `q=0` or `u=0` are exactly

    Q=0, 1, -1.

Normalize parameters by

    B=2a*s*k4/k3,
    C=4+2c0/(p*k3),
    L=2l4/(k3^2*s)!=0.

Up to the nonzero factor `p*k3^2*s/2`, the complete polynomial (10) is

    Rbar(Q)=[C-3Q^2+B(2Q^3-3Q)]
            [Q+(B/2)(1-2Q^2)]+L(1-2Q^2).           (12)

It always has a root outside `{0,1,-1}`. The proof is elementary and
retains all multiplicities.

If `B!=0`, its leading three coefficients are

    [Q^5]Rbar=-2B^2,   [Q^4]Rbar=5B,
    [Q^3]Rbar=4B^2-3.

If every root were forbidden, there would be nonnegative integers
`e,r,s0` with `e+r+s0=5` such that

    Rbar=-2B^2 Q^e(Q-1)^r(Q+1)^s0.

Put `d=r-s0` and `v=r+s0<=5`. The degree-four coefficient gives
`B=5/(2d)`, so `d!=0`. The degree-three coefficient then gives

    13d^2=25v-100.                                  (13)

If `v<4`, the right side is negative. If `v=4`, it forces `d=0`.
If `v=5`, it says `13d^2=25`, impossible for an integer `d`.
Thus none of the 21 complete multiplicity patterns can occur. No
constraint on the lower coefficients or on `C,L` was needed here.

If `B=0`, then

    Rbar=-3Q^3-2LQ^2+CQ+L.

Since `L!=0`, zero is not a root. If every root were `1` or `-1`,
write `Rbar=-3(Q-1)^r(Q+1)^s0`, `r+s0=3`. The constant term and
degree-two coefficient respectively give

    L=-3(-1)^r,    r-s0=2(-1)^r.

The left side of the second equality is odd and the right side even.
This contradiction handles all four remaining multiplicity patterns.

Select the guaranteed allowed root. Equations (9)--(11) give a critical
point with `u!=0`, and `t=(q-1)/u^2` is finite. This contradicts
`J(F,G)=1`, proving the theorem for every `p!=0`.

## 5. Controls, scope, and verification boundary

The forbidden-root analysis is an exact finite algebraic universe, not
an assumption that one sampled critical polynomial is generic. It covers
21 quintic patterns when `B!=0` and four cubic patterns when `B=0`.
The top-coefficient contradiction pays every complex parameter value.

The simple case `B=0,C=0,L=1` supplies a named allowed-root control.
The case `B=0,L=0` is a necessary hostile to dropping nonconstant `L`:
for `C=3`, `Rbar=-3Q(Q-1)(Q+1)` has only forbidden roots. Constant
`L` was excluded by the original polynomial factor before (12), so it
cannot be silently admitted to the cubic argument.

This is a changed-global-family proof, not an automorphism argument for
`x->x-p` on `W`. The map `(P,z)` also need not be a coordinate system.
Its critical locus is deliberately used, and the source-coordinate loss
is paid by the forbidden-root theorem. The conclusion remains about
polynomial mates; no generic rational-mate exclusion is inferred.

## 6. Exact verification and reproduction

The new standalone [source](../../04-computation/planar_jc48_sep08_four_four_transport.py)
and [frozen output](planar_jc48_sep08_four_four_transport.out) pass **151
always-active exact gates**. Normal and optimized Python outputs are
byte-identical; no inherited mathematical implementation is imported.

The exact universe includes the entire all-`p` global section space, the
higher-unit terms in the infinity multiplicity-four local model, every
coefficient of the residue identity, the complete critical-curve map and
normalization, all 21 quintic forbidden-root patterns, and all four cubic
patterns. A triangular rank-nine minor is identically one for all `p`,
so section completeness is not inferred from a generic specialization.
The test of every multiplicity pattern uses the displayed exact coefficient
and parity contradictions, not numerical roots or a numerical scan.

Named controls retain an allowed-root cubic and the `L=0` hostile whose
three roots are all forbidden. The analytic proof supplies compact
exactness and the existence of a complex critical-polynomial root; the
finite controls do not replace those arguments.

Reproduce from the worktree root:

```bash
python3 04-computation/planar_jc48_sep08_four_four_transport.py
python3 -O 04-computation/planar_jc48_sep08_four_four_transport.py
```

Frozen SHA256 pins:

- Source, 10,516 bytes:
  `4c284b826ed72bcb0e06775488f08c89a1fed1f492bbb7f118d113918d080851`.
- Output, 420 bytes:
  `a420c162add140c10de64ff4290bf3929102b4aae72ad47d519ab5a1832a0897`.
- Semantic digest:
  `ad47627c90d50ba6afc0079324c4d069634efab507670e15816c229022bdc0f8`.

The [independent complete analytic/source audit](planar_jc48_sep08_four_four_transport_audit.md)
accepts all151 normal/optimized gates, full shifted sections, higher-unit
infinity branches and the complete critical-polynomial argument. A separate
exact original-source critical-point control also passes. The conclusion
remains polynomial-only and treats one finite root paired with infinity.
