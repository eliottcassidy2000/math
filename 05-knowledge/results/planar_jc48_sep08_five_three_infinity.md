# The two infinity placements of the (5,3) quartic boundary class

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This note concerns the complete global square-prefix entry on the fixed DG
surface. It does not assert a general quartic or Jacobian-conjecture theorem.

## 1. Statement and exact coefficient contract

Let

    W=(P1_x x P1_z)\{z=x^2},  t=1/(z-x^2),
    omega=dx wedge dt=r^2 dr wedge db_D,
    x=1/r,  t=-r^2-r^4 b_D.

Take every global `H in L2`, `L in L1`, and `F=H^2+L`, with
`deg_t H=2`. Suppose its leading binary octic has partition `(5,3)`
and one of its two zeros is infinity.

**Theorem.** No polynomial `G in C[x,t]` of any degree satisfies
`J(F,G)=1`. This statement includes both placements and every position of
the finite point. Nonconstant `L` can admit rational mates: Sections 5 and 7
supply two actual global families and prove why their rational mates cannot
be repaired polynomially.

Together with the independently audited
[all-finite theorem](planar_jc48_sep08_five_three.md), this would close the
polynomial-mate question for the complete binary `(5,3)` class at all
locations. That combination is only as strong as the two proved inputs.

Normalize the nonzero leading scalar by scaling `H` and the output, and
write `u=x-p`. This is an ordinary polynomial source coordinate, with
`du wedge dt=dx wedge dt`. It is **not** asserted to be an automorphism
of `W`. The full section equations are imposed in the original `x`.
For `d=3` or `5`, put `N=u^d` and write

    P=A4 u^4+A3 u^3+A2 u^2+A1 u+A0,
    Q=A4 u^2+(A3-2p A4-e_d)u+C0,
    M=B4 u^4+B3 u^3+B2 u^2+B1 u+B0,
    R=B4 u^2+(B3-2p B4)u+E0,
    H=N t^2+P t+Q,  L=M t+R,                       (1)

where `e_3=0`, `e_5=1`. These are complete coefficient spaces, not
selected prefixes. They follow by expanding the proved original-coordinate
rules `Q2=P4`, `Q1=P3-N5`, `R2=M4`, `R1=M3` in
[the full global filtration](planar_jc48_sep08_dg_quadratic.md).
An additive constant of `L` can be absorbed in the original generic value
`zeta`.

The inherited mechanisms are the
[finite first-jet obstruction](planar_jc48_sep08_shared_roots.md),
[unit-M Newton/Morse analysis](planar_jc48_sep08_quartic_common_root.md),
and [complete inverse coefficients](planar_jc48_sep08_five_three.md).
The canonical hostile is a rational mate with incompatible principal parts
on two actual components of one source fibre. The corrected near miss is
transporting the finite local differential without its infinity factor.
The live concepts are full section rows, weighted infinity, compact pole
degree, higher inverse residues, and actual same-fibre pole repair.

The source/target map in the local arguments is normalization of the
original generic fibre, retaining both its constant and its volume. In the
field arguments it is an explicitly given finite or birational map, with
its degree and transformed volume stated. Boundary multiplicities alone
lose the normal jets; fibrewise poles alone lose ordinary affine points.
Both are retained below.

## 2. Complete weighted infinity regularity for multiplicities three and five

Use the actual surface involution

    r=1/x,  T=-x^2-x^4t,  omega=r^2 dr wedge dT.

Near its boundary point, with local boundary equation `s=0`, write

    calN=r^m a(r)+s B(r)+s^2 C(r,s),  a(0)!=0,
    calM=D(r)+s E(r,s),
    Phi=calN^2+s^3 calM-zeta s^4,
    eta=(unit) r^2 s^2 dr/Phi_s,                    (2)

with `m=3` or `5`. The unit accounts for harmless chart trivializations.
The factor `r^2` is the actual canonical numerator; omitting it changes the
answer. Write `j=ord B`, `n=ord D`, allowing infinity.

**Weighted local lemma.** Every normalized generic branch in (2) has
regular `eta`, for arbitrary higher coefficients.

If `B(0)!=0`, the proved analytic-centre argument gives two already
regular unweighted branches for every `m,n`. If `D(0)!=0`, use the
unit-M analysis: for `m=5` all cases are unbalanced (`5!=3j`) and regular;
for `m=3` the only balanced case `j=1` has at most logarithmic unweighted
poles. The extra `r^2` makes those regular as well.

It remains to treat `j,n>=1`, where the local Weierstrass degree in `s`
is four. The following list exhausts all four determinations in each
case; ramified determinations are counted as actual branches when giving
orders.

* For `j=1`, there are two low branches `s~r`, each with weighted order
  one. The other two determinations cancel `calN` around an analytic
  centre of order `m-1`. If `ell` is the order of `calM-zeta s` there,
  genericity gives `1<=ell<=m-1`; the `-zeta s` term pays the upper bound.
  Their form is a unit times
  `r^((m-3-ell)/2+2) dr`, nonnegative for both values of `m`.

* For `m=3,j>=2,n=1`, there is one low branch `s~r` and three high
  determinations `s~r^(5/3)`. The latter form one branch with
  `(ord r,ord s)=(3,5)` and weighted differential order five. For
  `m=3,j>=2,n>=2`, put `r=tau^2,s=tau^3 Z`. The quartic in `Z` has four
  simple nonzero roots generically, paired into two actual branches;
  each weighted differential has order two.

* For `m=5,j=2,n>=2`, there are two low branches `s~r^2`, weighted order
  zero, and two cancelling determinations centred at order three.
  Their generic contact has `1<=ell<=3` and their form is a unit times
  `r^((3-ell)/2) dr`, regular on either normalization.

* For `m=5,j>=3,n=2`, one low branch `s~r^2` has weighted order zero.
  The other three determinations have `(ord r,ord s)=(3,8)` and
  weighted order two. For `m=5,j>=3,n>=3`, four determinations
  `s~r^(5/2)` form two actual branches, each of weighted order zero.

* Finally take `m=5,n=1,j>=2`. There is one low branch `s~r` of
  weighted order one. The three high determinations use `s=r^3 Z`.
  For `j>=3` their face is `a0^2+D1 Z^3`, with simple roots and
  weighted order one. For `j=2` the full face is

      (a0+B2 Z)^2+D1 Z^3.                           (3)

  Its roots are nonzero and it has no triple root. Its only possible
  repeated root is a double root; explicitly it requires
  `Z=-3a0/B2`, `D1=4B2^3/(27a0)`, and the second derivative is nonzero.
  Parameterized Morse reduction gives a local equation
  `v^2=R(r,zeta)`. At its critical section,
  `R_zeta` is a unit times `r^2`, because `-zeta s^4` is two weights
  above the order-ten initial face. Hence generic contact is at most
  two. The weighted form is a unit times `r dr/v`. At contact one its
  actual ramified order is two; at contact two the two smooth branch
  orders are zero. The remaining simple high root has order one.

The nonzero low faces are generically simple because their repeated-root
eliminants have a nonzero constant term, exactly as in the proved tangent
quartic argument. In the fractional faces the displayed monomials have
simple nonzero roots. Higher terms have strictly higher weight or enter
the paid analytic centre/Morse contact. Thus the list includes all
coefficient degenerations and every branch, including the balanced (3).

On a smooth generic fibre in `W\D`, `eta` is a unit. At an intersection
with `D` it is regular. Accordingly all possible primitive poles in these
two placements occur over the single finite boundary point. This is a
statement about the complete compactification, not just a source chart.

## 3. Quintuple infinity: finite pole degree forces the complete constant-D stratum

Here `N=u^3`. Suppose a rational mate exists and `L` is nonconstant.
At the finite triple, an M-unit is regular or has forbidden nonzero
logarithms; a normal unit is regular. Together with Section 2, absence of
an active finite point would leave no poles of a rational primitive on
any compact generic component, although `eta` is nonzero there.
The finite first-jet obstruction therefore requires

    P(0)=P'(0)=M(0)=M'(0)=0.                         (4)

At this active triple, `u=tau^2,s=tau^3 Z` produces four simple nonzero
determinations, giving exactly two actual differential poles of order two.
Their total primitive pole capacity is two. If `F|D` were nonconstant,
a generic transverse point of `D` would make a primitive have local
degree three, exceeding this capacity even on the same component.
Consequently `F|D` is constant. In (1), the coefficient of `b_D` in
`H|D` is `-A4`, and in `L|D` it is `-B4`. A square of a nonconstant
linear polynomial cannot be cancelled by a linear one; hence
`A4=B4=0`.

The exact inverse coefficient `T2=-M/(4N)` implies that `M du/N`
is rational-exact. Indeed the corresponding Laurent coefficient of a
rational mate has derivative proportional to `T2`, and normalized trace
from `C(u)(sqrt N)` commutes with differentiation. Its residue sets
`B2=0`. Under the standing nonconstant-L assumption, write `B3=lambda!=0`.
The complete residual, for every `p`, is therefore

    H=u^3 t^2+u^2(a u+b)t+a u+c,
    Z=u^3 t+u,  L=lambda Z,  F=H^2+lambda Z.          (5)

This same expression is obtained from the full global coefficient
intersection. No surface translation has been used.

## 4. The a-nonzero residual: incompatible degree-two and degree-three maps

Assume `a!=0` in (5). At infinity the first normal order is one and
`calM` has order one. More explicitly,

    calN(r,0)=r^5(1-pr)^3,
    [s]calN=-a r+O(r^2),
    calM(r,0)=-lambda r+O(r^2).

The two low branches `s~r` each give a zero of order one. The two high
cancelling determinations have centre order four and contact one; they
form one ramified branch `r=tau^2` with a differential zero of order six.
There are no points on `D` in the generic fibre, since its restriction
is constant. The two finite poles have orders two. Thus the divisor of
`eta` has degree `1+1+6-2-2=4`.

The geometric generic curve is integral. A nontrivial polynomial
composition must have outer degree two or four because `deg_t F=4`.
An outer quartic forces its leading `N` to be a square up to scalar,
contradicted by `N=u^3`. An outer quadratic, after completing the square,
would give `(K-H)(K+H)=L-constant`; comparison of t-degrees forces
`K=H` and `L` constant. The classical closed-polynomial criterion then
pays geometric generic integrality in characteristic zero; see
[Arzhantsev--Petravchuk, Theorem 1 and Lemma 3](https://arxiv.org/pdf/math/0608157v2)
and [THM-3827, generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases](../../01-canon/theorems/THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases.md).
The normalized compact generic curve consequently has genus three.

There is an actual degree-three map to the `h=H` line. Put `U=1/u`
and `z=Z`; this is a birational field change and

    h=z^2 U^3-2z U^2+(1+bz)U+a z+c-b,
    omega=-U dU wedge dz,
    z=(zeta-h^2)/lambda.                             (6)

As a polynomial in `U`, (6) gives degree three over `C(h,z)`;
geometric generic integrality preserves this degree after extending the
generic constant field. A rational primitive has total pole degree at
most two. Degree one would make the curve rational. If its degree is two,
the function fields generated by it and by `h` generate the entire
curve field: the remaining degree divides both two and three.
The resulting map to `P1 x P1` is birational onto a curve of bidegree
`(2,3)`. By adjunction that curve has arithmetic genus
`(2-1)(3-1)=2`; normalization cannot increase genus. This contradicts
genus three. Thus this subcase has no rational mate.

The ordinary cubic trace of `U/h_U` is zero. It would give no obstruction;
the complete degree and divisor data in (6) are the needed supplier.

## 5. The a-zero residual: all critical alternatives and a sharp conic exception

Set `v=ut`, so the actual source form is `du wedge dv/u` on `u!=0`.
For (5) with `a=0`,

    H=u v(v+b)+c,
    F=H^2+lambda(u^2 v+u).

For a critical point with nonzero denominators, elimination gives

    u=-(2v+b)/[v(3v+b)],
    H=lambda/[2v(3v+b)],
    C(v)=-4v^3+6(c-b)v^2+2b(c-b)v-lambda=0.           (7)

Every root of (7) outside `v=0,-b/3,-b/2` supplies an actual point
with `u!=0`, hence a source critical point. The constant `-lambda`
excludes zero. If `b=0` all forbidden addresses are zero, so a usable
root exists. If `b!=0`, compare the four possible factorizations

    C(v)=-4(v+b/3)^j(v+b/2)^(3-j),  0<=j<=3.

Coefficient comparison leaves exactly

    c=b/3,  lambda=4b^3/27,  j=3.                    (8)

Thus all cases except (8) have a source critical point and no polynomial
mate. This finite list exhausts possible root concentrations of a cubic;
it is not a parameter sample.

Normalize (8) by the invertible polynomial scaling
`u=3U/b`, `t=b^2 T/9` and by the constant output factor. This reduces to
`b=3,c=1,lambda=4`, preserving polynomial-mate existence. The resulting
literal functions, written again in `(u,t)`, are

    v=ut,  H=u^3t^2+3u^2t+1,
    A=1+u(v+1)(v+4),  B=1+uv(v+1),
    F=H^2+4(u^3t+u)=AB,
    G0=2/(A^2-F)=1/[2A u(v+1)],  J(F,G0)=1.          (9)

The map to `(A,B)` is birational:

    v=4(B-1)/(A-B),
    u=(A-B)^2/[4(A+3B-4)].

Hence the constants of the fibre derivation are exactly `C(F)`, and every
rational mate is `G0+K(F)`. On the special fibre `F=1`, the two disjoint
source divisors `E={u=0}` and `Gamma={1+ut=0}` have

    [(F-1)G0]|E=2,  [(F-1)G0]|Gamma=1.              (10)

Both are simple poles and cannot be cancelled by a common rational
function of `F`. Higher poles of `K` would themselves survive. No
polynomial mate exists. This is an actual nonconstant-L rational mate
in the exact infinity-five partition, and it validates the distinction
between rational and polynomial conclusions. All original global rows
are retained for every `p`; the scaling used only a plane polynomial
change after that reduction.

## 6. Triple infinity: exact higher coefficients and actual critical points

Here `N=u^5`. By Section 2 the entire infinity fibre is regular.
At the finite quintuple an M-unit is unbalanced (`5!=3j`) and regular;
a normal unit is regular. Thus a rational mate again forces the active
jets (4). The residue of `M/u^5` then sets `B4=0` in (1). Rename
coefficients to obtain the complete rows

    H=u^5t^2+u^2(a u^2+b u+k)t
        +a u^2+(b-2pa-1)u+c,
    L=n u^2t+m(u^3t+u).                             (11)

In the actual chart `w=u^2t`, the first two coefficients are

    f(w)=(kw+c)^2+n w,
    g(w)=2(kw+c)(w^2+bw+b-2pa-1)+m(w+1).

For every simple root of `f(w)=zeta`, the moving-root residue forces
`g=C f'`. If `k!=0`, its cubic coefficient `2k` cannot match the linear
`f'`; hence `k=0`. If `n!=0`, `f` is linear and the quadratic/linear
coefficients of `g` force `c=m=0`. Then in original `(u,t)` the function
`F` is divisible by `u^2`, so every point on `u=0` is source critical.
This already excludes polynomial mates, for all `a,b,p`.

Suppose therefore `n=0,m!=0`. The higher inverse coefficients are
load-bearing. With

    D0=Q-P^2/(4N),  E0=R-MP/(2N),
    T6=-M E0/(8N),
    T10=-M(3E0^2-D0 M^2/N)/(32N),                    (12)

their rational exactness follows from the full formal inverse, by
normalized trace as for `T2`. Formula (12) can be recovered by centering
`t=y-P/(2N)` and summing the two opposite inverse roots; the exact
pair-sum identity and all factors are proved in Section 5 of
[the all-finite (5,3) theorem](planar_jc48_sep08_five_three.md).
It concerns the original `F=zeta`, not a modified first row. In (11),

    Res_0 T6=m(2an+bm-2m)/16,
    Res_0 T10 |(n=0,b=2)=c m^3/32.

Consequently `b=2,c=0`. Put `q=1+u^2t`, `v=u q`; the full residual is

    H=v^2/u+a u(v-2p),  L=m v,  F=H^2+m v.          (13)

For `a!=0`, use the genuine chart `(u,v)` with `u!=0`. Its coordinate
Jacobian is `u^3`, so vanishing of both coordinate derivatives is
exactly source criticality. The fold is

    a(v-2p)u^2=v^2.

On it `H=2v^2/u` and

    F_v=12a v^2-16ap v+m,  F_u=0.                   (14)

The quadratic in (14) has a root outside `{0,2p}`: zero is excluded by
`m!=0`, and concentration of both roots at `2p` would require `p=0`
from its linear coefficient and then `m=0` from its constant. For an
allowed root choose either square root
`u^2=v^2/[a(v-2p)]`. It is nonzero, and (14) gives an actual source
critical point. This pays every value of `p` and every nonzero `a,m`.

## 7. The remaining rational family and its incompatible second principal parts

In (13) with `a=0`, the literal functions are

    q=1+u^2t,  v=u q,  H=u q^2,
    F=H^2+m v,  G0=1/(6v^3),  J(F,G0)=1.            (15)

Both `H` and `L=m v` satisfy the full global constraints for every finite
`p`, in the triple-infinity partition. The field change `(u,t)->(h,v)`
is birational with `u=v^2/h`; since `F=h^2+m v`, the generic field is
`C(F)(h)`. Every rational mate therefore differs from `G0` by `K(F)`.

Consider the same special fibre `F=0`. On `E={u=0}`, the source coordinate
`t` is regular and `H^2=v^2+O(v^4)`. On
`Gamma={1+u^2t=0}`, with `u` a unit, `H^2=v^4/u^2`. These two disjoint
irreducible divisors give respectively

    G0 = m^3/(6F^3)+m/(2F^2)+O(F^-1)       on E,
    G0 = m^3/(6F^3)          +O(F^-1)       on Gamma. (16)

The nonzero second coefficient cannot be simultaneously cancelled by
`K(F)`, while higher poles would survive on both divisors. Thus no
polynomial mate exists even in this exact rational family.

If `L` was constant at the start of either placement, the factor
`2H` in `J(H^2+constant,G)` already excludes polynomial mates. This
completes the two-placement polynomial proof. Statements (9) and (15)
show why its rational scope cannot be enlarged. No assertion about
other binary partitions or all quartic polynomials follows.

## 8. Reproducible exact controls and review boundary

The source is
[planar_jc48_sep08_five_three_infinity.py](../../04-computation/planar_jc48_sep08_five_three_infinity.py).
Run it normally and with `-O`. It uses always-active checks, not removable
assertions. The named exact universe consists of the full symbolic
coefficient spaces (1), all local face regimes of Section 2, the four
possible cubic root-concentration patterns, both symbolic residue
identities, both actual critical-point eliminations, and the two
rational families with distinct same-fibre principal parts. The source
checks original-coordinate globality and field inverses independently
of the displayed local arguments. It does not infer a theorem from a
bounded coefficient census.

The normalization/completeness, generic-component and genus arguments
are analytic proof obligations, not consequences of the gate count.
The [independent complete audit](planar_jc48_sep08_five_three_infinity_audit.md)
accepts both complete infinity placements, every local degeneration, the
genus/degree argument, all actual critical points and both incompatible
same-fibre principal parts. Both107-gate replays match the frozen output.
Together with the all-finite companion this closes every polynomial
placement of binary5+3, while retaining the rational exceptions below.


Reproduce from the repository root:

```sh
python3 -B 04-computation/planar_jc48_sep08_five_three_infinity.py
python3 -B -O 04-computation/planar_jc48_sep08_five_three_infinity.py
```

Both executions pass **107 always-active exact gates** and agree with all
349 bytes of [the frozen output](planar_jc48_sep08_five_three_infinity.out).

| Artifact | SHA-256 |
| --- | --- |
| Source, 8,514 bytes | `a693f7ad98a082d6897fa247316904f89056a146a771fe04862eccd3067c15a1` |
| Frozen output and both replays, 349 bytes | `05ee85e0d98edf3f36d14b8cfe7116aff7824b814c5e7b0abfe0262d3088cceb` |
| Semantic record | `d2c7bd4f6fead41fa623987705af5b830f098bf14243d9a1f53d9b44d70fcf73` |
