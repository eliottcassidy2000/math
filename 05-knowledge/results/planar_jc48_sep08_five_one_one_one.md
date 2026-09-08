# The entire quartic boundary type 5+1+1+1 has no polynomial mate

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The formal leading-coefficient dependency is now proved and independently
audited, and this consumer has passed its full independent audit. This is
an unrestricted-degree **polynomial**
mate exclusion on the declared DG source. In the all-finite case the
rational conclusion stops at `L` constant. Section 6 supplies an actual
rational-mate family in this exact boundary partition, making that
distinction sharp.

## 1. Actual family, complete sections, and all placements

Retain the fixed surface and actual source coordinates

    W=(P1_x x P1_z) minus S,      S={z=x^2},
    t=1/(z-x^2),                 omega=dx wedge dt,
    x=1/r,  t=-r^2-r^4 b,        D={r=0},
    omega=r^2 dr wedge db.

Let `H in L2`, `L in L1` be global functions, with `deg_t H=2`, and put
`F=H^2+L`. Suppose the **complete binary-octic** leading boundary section
of `H` has exactly four distinct zeros of multiplicities `5,1,1,1`.

**Candidate theorem.** There is no `G in C[x,t]`, of any degree, with
`J_(x,t)(F,G)=lambda!=0`. This does not require `G` to extend to `W`.
There is no normalization of arbitrary boundary points by an automorphism
of `W` in the argument.

Write the entire source polynomials as

    H=N(x)t^2+P(x)t+Q(x),       L=M(x)t+R(x).              (1)

The [proved complete DG filtration, Section 3](planar_jc48_sep08_dg_quadratic.md)
gives the full coefficient spaces:

    deg N<=8,
    P=2N8*x^6+2N7*x^5+sum_{i=0}^4 p_i*x^i,
    Q=N8*x^4+N7*x^3+(p4-N6)*x^2+(p3-N5)*x+q0;          (2)

    M=sum_{i=0}^4 m_i*x^i,
    R=m4*x^2+m3*x+r0.                                  (3)

Here `N_i` are all nine coefficients of `N`, and the five `p_i`, `q0`
are independently free. In (3) the five `m_i`, `r0` are independently
free. In particular `deg M<=4`, and **M=0 forces L to be constant**.
That last implication uses globality; an arbitrary source polynomial
of degree zero in `t` need not be constant.

For clarity, (3) also follows directly by substituting the second chart:
the coefficients of `r^-2` and `r^-1` in `Mt+R` are `R2-m4` and
`R1-m3`; after their removal the entire expression is polynomial in
`r,b`. Completeness of the degree bounds comes from the whole section
box `H^0(O(2,1))`, not from a chosen selection of generators.

The degree of `N` equals eight minus the boundary multiplicity at
infinity. Therefore the distinct-point hypothesis leaves exactly three
placement cases:

1. All four points are finite: `deg N=8`, partition `5+1+1+1`.
2. The fivefold point is infinity: `deg N=3`, three simple finite roots.
3. One simple point is infinity: `deg N=7`, finite partition `5+1+1`.

The final section pays the last two cases in the original affine
coordinate, without transporting the source volume.

## 2. A second inverse coefficient with a genuine rational primitive

The [separate leading-exactness argument](planar_jc48_sep08_leading_exactness.md)
constructs the formal inverse and embeds every rational mate in its
Laurent series field. We repeat the needed coefficient calculation
explicitly so that neither a frozen prefix nor an unproved hierarchy
is substituted for it.

Choose `a=sqrt(N)` and put `K=C(x)(a)`. Solve `F(x,T)=v^-4`, using
`T=a^-1 v^-1+O(1)`. Write `U=vT`. The exact equation is

    (N U^2+P U v+Q v^2)^2+M U v^3+R v^4=1.            (4)

The derivative with respect to `U` at `v=0`, `U=a^-1`, is `4a!=0`.
The unique formal coefficients through order three give

    T_-1=1/a,
    T_0=-P/(2N),
    T_1=(P^2-4NQ)/(8N^(3/2)),
    T_2=-M/(4N).                                       (5)

Thus all induced coefficients in `P,Q` have been retained in the
calculation; their cancellation from `T2` is an identity, not a choice.
The constant row `R` first enters (4) at order four and does not change
any coefficient in (5).

Suppose now that `J(F,G)=1` for an arbitrary `G in C(x,t)` and let
`Gtilde=G(x,T)`. Holding `v` fixed, the exact formal chain rule gives

    partial_x Gtilde=(1/4)v^5 partial_v T.              (6)

The source field embeds because any nonzero polynomial denominator in
`t` has a unique lowest Laurent term after substitution. Differentiation
in `x` extends uniquely to `K` and acts coefficientwise, as in the
leading-exactness proof. No degree or pole bound on `G` is involved.

If `g6=[v^6]Gtilde`, (6) yields

    g6'=T2/2,
    d(2g6)=-M dx/(4N),          2g6 in K.               (7)

It is essential that this initially gives a primitive in `K`, not
necessarily in `C(x)`. But the differential on the right is rational
in `x`. Its primitive descends by the normalized trace:

    B=(1/[K:C(x)]) Tr_{K/C(x)}(2g6),
    dB=-M dx/(4N),              B in C(x).              (8)

Trace commutes with the uniquely extended characteristic-zero
derivation. This formula also covers a square `N`, when the field
degree is one. Alternatively, at a simple root of nonsquare `N`, the
ramification `x-q=s^2` doubles a rational residue and cannot remove it.

At each simple root `q` of `N`, (8) forces

    Res_q[-M dx/(4N)]=-M(q)/(4N'(q))=0,
    hence M(q)=0.                                      (9)

Since `N'(q)!=0`, this is an actual evaluation constraint. On the
radical normalization the residue would be `-M(q)/(2N'(q))`, with the
same consequence. A nonzero constant Jacobian is reduced to one by
scaling `G`.

## 3. The fivefold point must be a singular shared point with two zero jets

Assume for this section that all four boundary points are finite. Let
`p` be the fivefold point and `q1,q2,q3` the three distinct simple
points, all different from `p`.

Put `s=z-x^2` locally. The full numerators, not only their boundary
values, are

    Ncal=N(x)+sP(x)+s^2 Q(x),
    Mcal=M(x)+sR(x),
    E=Ncal^2+s^3 Mcal-c s^4=0.                         (10)

Here `c` is a generic fibre value. At a finite point the relative form
is, up to an everywhere nonzero coordinate multiplier,

    eta=s^2 dx/E_s,           dF wedge eta=omega.        (11)

We use the audited local conclusions in
[quartic common-root necessity, Sections 2.1–2.2](planar_jc48_sep08_quartic_common_root.md)
and [shared roots, Sections 1–3](planar_jc48_sep08_shared_roots.md),
with their actual normal derivatives retained.

First, if `M(p)!=0`, every normalized branch at the fivefold point
has regular `eta`. To check the entire local scope, write

    Ncal=A(u)+sB(u)+s^2 C(s,u),
    u=x-p,  ord A=5,  j=ord B,

allowing `j=infinity`. The M-unit cases are all unbalanced because
`5` is not divisible by three:

* For `j=0`, the local Weierstrass degree in `s` is two. Both
  cancellation determinations have `s` of order five and
  `eta=(unit)u^(5/2)du`, regular on their actual normalization.
  The apparent third Newton root has nonzero `s` at `u=0` and does
  not pass through the point.
* For `j=1`, the Weierstrass degree is three. One branch has `s`
  of order two and regular differential. The other two have `s`
  of order four and `eta=(unit)u du`. These exhaust that degree.
* For every `j>=2`, including infinity, the leading balance is
  `A(u)^2+M(p)s^3=0`, with three simple nonzero Puiseux
  determinations and `s` of order `10/3`. The term `3M(p)s^2`
  dominates `E_s`, so `eta` is a unit times `du` and is regular
  after the possible cubic normalization. These exhaust the
  Weierstrass degree three.

Second, if `M(p)=0` but `P(p)=Ncal_s(p,0)!=0`, all branches are
again regular. The analytic centre `Ncal=0` has `s=psi(u)` of
order five. For generic `c`, the order `ell` of
`Mcal(psi(u),u)-c psi(u)` satisfies `1<=ell<=5`. Its two
determinations give

    eta=(unit)u^((5-ell)/2)du,

regular after normalization. The local Weierstrass degree is two,
so no additional branch is missing.

Every other boundary point is simple, and is regular for arbitrary
`Mcal`: if shared, use `Ncal` itself as the tangential coordinate;
the generic leading equation is `Ncal^2+unit*s^4=0`, and (11)
is regular on its two branches. If not shared, the M-unit local
calculation applies. The point at infinity is not a zero of `N|S`
in the all-finite case, so it contributes no boundary branch.

Consequently, in either of the two regular alternatives at `p`,
`eta` is holomorphic on every compact normalized generic component.
This contradicts a rational mate. Indeed, choose a generic value so
the fibre is smooth in `W` and avoids the finitely many exceptional
values of the rational function's denominators. The relative form
is regular in `W`, including `D`, and nonzero on each component
meeting `W` away from `D`. Every generic component does so: it
cannot equal the fixed `D`, and `N|S` is nonzero, so no component
of the closure equals `S`. The restriction of a rational mate
satisfies `dG=eta`. A pole of that restriction would give a pole
of its derivative; thus a holomorphic `eta` would have a global
holomorphic, hence constant, primitive. This is impossible for
the stated nonzero form. No geometrical irreducibility or
nonconstant `F|D` assumption is needed here.

Therefore a rational mate in this all-finite case requires

    M(p)=0,      P(p)=0.                                (12)

This makes the shared point singular for `Ncal`: its tangential
derivative already vanishes because its multiplicity is five.
The proved finite first-jet obstruction now applies at this
actual point. It forces

    P'(p)=0,     M'(p)=0.                               (13)

For completeness, the quadratic part of `Ncal` and linear part
of `Mcal` at the point give, on setting `s=uZ`,

    E=u^4[(a+bZ+eZ^2)^2+dZ^3+(f-c)Z^4+...],
    a=[u^2]N=0,   b=P'(p),   d=M'(p).

If `b!=0`, the nonzero part of the tangent polynomial is a
quadratic with two simple nonzero roots for generic `c`; if
`b=0,d!=0`, it has one such root. Each gives a nonzero logarithmic
residue `Z0^2/P_c'(Z0)` in (11). This is an obstruction to a
rational primitive, not a count of selected branches. Thus (13)
is valid without imposing higher-jet conditions.

## 4. Three simple-root values plus two jets force the whole linear row to vanish

Let `C(x)=prod_i(x-qi)`. Equations (9) give `C|M`; (12)–(13)
give `(x-p)^2|M`. These factors are coprime. Since `deg M<=4`
by the complete global space (3), it follows that

    M=0,      L=R=r0 in C.                              (14)

This is an exact degree-five counted-zero argument, valid at all
distinct complex positions. It is not a numerical interpolation.
Equivalently write `M=C(Ax+B)`. The two remaining jet equations
have determinant `-C(p)^2!=0` in the unknowns `A,B`. The source
also checks the independent confluent Vandermonde determinant
for evaluations at `p,p,qi`, with one derivative at `p`.

If `G` were polynomial, (14) would yield the polynomial identity

    J(F,G)=2H J(H,G)=lambda!=0,

which is impossible because `H` is nonconstant (`deg_t H=2`).
This proves the all-finite polynomial-mate exclusion, of every
degree. The intermediate result is precise: **any rational mate
would force L constant**. That result does not itself exclude
rational mates. For example, outside this boundary partition,
`F=t^4`, `G=-x/(4t^3)` has Jacobian one with `L=0`.

No assumption from the elliptic exactness locus is used. In
particular the argument covers both `N=x^5(x^3+1)`, whose leading
differential is exact, and the nonexact members of the same
multiplicity type.

Two literal global controls show why both parts of the argument
are needed. Take

    N=x^5(x^3+1),  P=2x^6,  Q=x^4-x,
    H=N t^2+P t+Q.

Its second-chart expression is `b^2+2br+b^2r^3`, so all induced
global corrections are present. The functions

    L1=(x^3+1)t+x,
    L2=x(x^3+1)t+x^2

are both global and their leading coefficients vanish at every
simple root of `N`. The first still has `M(0)=1`; the second has
`M(0)=0`, `M'(0)=1`. Neither three simple-root evaluations alone
nor four counted zeros force `M=0`. In the second control the
literal tangent polynomial is `Z^3-cZ^4`; its nonzero root has
residue `-1`. This supplies an actual logarithmic boundary
obstruction instead of silently deleting the remaining jet.

## 5. Infinity placements are excluded by the leading differential

The leading-exactness theorem applied to `F`, whose leading
coefficient is `N^2`, requires

    dx/sqrt(N) exact in C(x)(sqrt(N)).                  (15)

If the fivefold point lies at infinity, `N` has degree three
and three distinct simple finite roots. On its normalized
quadratic curve, the differential in (15) is nonzero holomorphic:
it has order zero at each finite branch point and at its sole
infinity point. A rational exact differential cannot have this
property, so in this case there is no rational mate at all.

If a simple point lies at infinity, `N` has degree seven with
finite multiplicities `5,1,1`. Its differential has just one
possible primitive pole, of order three, over the fivefold
root. At infinity the differential has order four, so a
nonconstant primitive would have local degree five there.
Its total map degree is at most three, a contradiction.
Again this excludes rational mates, not only polynomial ones.

These direct calculations are also the corresponding two
entries of the [complete boundary exactness classification](planar_jc48_sep08_boundary_exactness.md). The present proof
does not use its more delicate elliptic coefficient condition.
Together with the all-finite case, they exhaust the specified
binary partition in all positions. They do not identify the
finite and infinite volume forms or assert a projective source
symmetry.

## 6. A sharp rational-mate family inside this same partition

Let `delta!=0` and let `beta,gamma,c0` be arbitrary complex constants.
Put

    h=x^2+x^4t,
    H=(x^4+delta*x)(1+x^2t)^2+beta*h+gamma,
    F=H^2+c0.                                          (16)

This is an actual global family: its second-chart expression is

    H=(1+delta*r^3)b^2-beta*b+gamma.

Its leading coefficient is `N=x^5(x^3+delta)`. Thus its fivefold point
is zero and its three other roots are distinct, nonzero and simple,
exactly the boundary partition of the theorem. Direct calculation gives

    J(H, 1/(3delta*h))=1,
    J(F, 1/(6delta*H*h))=1.                            (17)

Both mates are rational. The second one is an actual counterexample to
any attempted rational-mate upgrade of the theorem, without leaving its
stated global coefficient class or its boundary multiplicity type.
Neither mate is polynomial: its displayed reciprocal denominator is a
nonconstant source polynomial.

The mechanism is transparent after a field map with its degree retained.
On `x!=0`, set `u=x^-3` and keep `h`. Then

    du wedge dh=-3 dx wedge dt,
    H=(1+delta*u)h^2+beta*h+gamma.                      (18)

The equation is linear in `u`; differentiating `1/(3delta*h)` pays the
volume factor `-1/3` in (18), proving the first identity in (17). The
second follows by dividing that primitive by `2H`.

This is **not a birational coordinate change**. The original field is
`C(x,h)` because `t=(h-x^2)/x^4`, while the subfield is `C(u,h)` with
`x^3=u^-1`. Its degree is three: `u^-1` is not a cube in `C(u,h)`, as
its valuation at `u=0` is `-1`. No source translation, inversion of the
DG surface, or plane polynomial automorphism is asserted. This explicit
retention of the field degree is what permits the simple primitive to
lift to the original rational source without creating a polynomial mate.

For the normalized `delta=1` instance, this family has
`p0=p1=p2=0`, `p3=2`, and `p4=beta` in (2). Thus it preserves all
induced coefficient rows. Its leading differential is the exact
elliptic survivor already paid by the boundary classification, and now
it also has a full rational mate. The earlier generic example
`t^4` remains a simpler control, but (16) pays the stronger within-class
sharpness assertion.

## 7. Inheritance, lost data, and exact verification

The closest proved mechanisms are the full DG section space,
the shared-root first-jet obstruction, and the formal leading
inverse. The new operation extracts an additional **actual**
inverse coefficient and combines its three separate residues
with the two necessary local jets. The source is the entire
polynomial `F=H^2+L`, the target is a rational differential in
`x`, and the map sends it to `-M dx/(4N)`. It retains all lower
rows until their cancellation has been proved, then discards
everything except `M,N`. The missing information is supplied
by normalized branches at the original fivefold point.

The live concepts are formal coefficient primitives, trace
descent, labelled simple-root residues, complete global
coefficient degree bounds, and finite normal jets. The controls
above distinguish the separate proof obligations. The
same-partition rational family (16) prevents upgrading (14)
to a rational-mate exclusion. The corrected
near miss is to assume an algebraic primitive was already
rational without paying trace or the doubled residues.

The source is
[planar_jc48_sep08_five_one_one_one.py](../../04-computation/planar_jc48_sep08_five_one_one_one.py).
Its universe is the symbolic complete global coefficient
space (2)–(3), arbitrary four distinct complex boundary
points, and all three infinity placements. It verifies the
inverse through the exact required order, the entire second
chart cancellation, labelled residues and both independent
five-zero determinants, local order thresholds, and the
literal global controls, including both exact identities in (17)
and the actual second chart. No finite bound on a mate or sample
coefficient search is substituted for the proof.

```bash
python3 04-computation/planar_jc48_sep08_five_one_one_one.py
python3 -O 04-computation/planar_jc48_sep08_five_one_one_one.py
```

Every gate uses an explicit exception and remains active in
the optimized run. Normal and optimized replays are byte-identical:
**62 gates**, 394 output bytes. The frozen output is
[planar_jc48_sep08_five_one_one_one.out](planar_jc48_sep08_five_one_one_one.out).

SHA256 pins:

* Source, 8,859 bytes:
  `029193e4fc766ebbfc7bf49a8af1bbc928da6cd88ea04f77dd1067d784e9a6fd`.
* Output, 394 bytes:
  `c6ef3ec938cffd949ad4a755ef8732c52c0a5bafce736415eafb3992ae90f776`.
* Semantic gate manifest:
  `de0cbef63edcb5a98cf214274d92dd91417c5d4bef999d43548d5dd5cd36fecc`.

All three artifacts are frozen for independent review. The formal
leading-exactness and boundary classification dependencies are now
proved and independently audited. The [complete independent audit](planar_jc48_sep08_five_one_one_one_audit.md)
accepts every local regime, full section rank, T2/trace descent, the five-zero
argument, all boundary locations and the genuine same-partition rational
hostile. Source and output are frozen and independently accepted.
