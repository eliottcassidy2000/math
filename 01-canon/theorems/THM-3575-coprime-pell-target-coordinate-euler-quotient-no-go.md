---
id: THM-3575
title: "Coprime Pell target-coordinate Euler and quotient no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For coprime nonzero
  P,Q in C[b], the target coordinate G_{P,Q} has an exact linear-times-
  quadratic pullback.  When Delta=(bP-Q)^2+3Q^2 and E=bP-2Q are nonzero,
  the residual component has compactly supported Euler characteristic
  -3+3n_P+3n_Q+2n_Delta+n_E-2[Q(0)=0].  Collision compatibility never
  gives Euler characteristic one.  Conversely, Euler characteristic one
  forces, up to common scale, P=t and Q=b.  The resulting residual S_t is
  not A2 for any irreducible parameter: its C*-quotient is a non-line
  cubic curve.  The split parameters t=1+-i*sqrt(3) remain a separately
  typed three-factor boundary.
source: kps-s188 + subagent-zeta
contributor: /root/resonant_graph_escape
audit: >
  PASS.  The target-coordinate identity, strict vertical fibres, Euler
  formula, collision obstruction, exhaustive Euler-sharp classification,
  torus quotient, parameter boundaries, and finite hostile atlas were
  independently rederived.  Normal, optimized, and stored transcripts are
  byte-identical.
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-1345-jc2-equivariant-category-poisson-reframing-dc1-shadow
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3560-jelonek-euler-gate-monomial-target-shear-no-go
  - THM-3570-universal-pell-conic-target-graph-factor-compiler
related:
  - THM-3565-resonant-linear-a-target-graph-factor-classification
  - THM-3568-reducible-target-graph-component-euler-no-go
  - THM-3571-quadratic-target-graph-euler-no-go
  - THM-3574-universal-reducible-target-graph-component-unit-no-go
companion: 04-computation/jacobian_coprime_pell_target_coordinate_euler_quotient_no_go_resonant_thm3575.py
output: 05-knowledge/results/jacobian_coprime_pell_target_coordinate_euler_quotient_no_go_resonant_thm3575.out
script_sha256: 72f727c29f2ce52a8a13c5d06a912af9c5cc3734540fcd4f1943f66f52485c57
output_sha256: ec5f89558074247aaf0722111e4ba96ed3599b413031b58e21a291b90052116c
hash_basis: LF-normalized bytes
---

# THM-3575 -- coprime Pell target-coordinate Euler and quotient no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
compiles every coprime polynomial numerator/denominator in THM-3570's Pell
factor parameter into a target coordinate.  The surprising boundary is
sharp: collision-compatible rows all fail the Euler test, while every
Euler-sharp row belongs to one one-parameter family which misses the known
triple collision.  Euler characteristic alone does not finish that family;
its inherited torus quotient does.

All varieties are over `C`, and `chi` means compactly supported Euler
characteristic.  Write the fixed THM-1300 map as

```text
F=(A,B,C),                         u=1+xy.                 (1)
```

## 1. A coprime Pell graph is a target coordinate

Let `0!=P,Q in C[b]` be coprime and define

```text
G_(P,Q)=Q^3 c-2P^3 a+bP^2Q-2PQ^2.                        (2)
```

Choose `r,s in C[b]` with

```text
(-2P^3)r+Q^3s=1
```

and put `K=-sa+rc`.  Derivatives of `P,Q,r,s` occur only in the middle
`b`-column, so direct expansion gives

```text
Jac_(a,b,c)(G_(P,Q),b,K)=(-2P^3)r+Q^3s=1.                (3)
```

Thus `(G_(P,Q),b,K)` is a Keller triangularization and `G_(P,Q)` is a
target coordinate.  In fact the inverse linear system in `(a,c)` has
determinant one, so it is a polynomial automorphism.  Consequently

```text
T_(P,Q)=V(G_(P,Q)) ~= A2.                                (4)
```

The pair `(P,Q)` is defined only up to common nonzero scale for purposes of
the zero set in `(2)`.

## 2. Exact linear-times-quadratic source factorization

Put

```text
p=P(B),                         q=Q(B).
```

Then direct expansion, without localization, gives

```text
G_(P,Q)(F)=(qx-pu) S_(P,Q),                               (5)

S_(P,Q)=
 z(2p^2u^2-pqxu-q^2x^2)
 +2p^2y^2(3u+1)-pqy(3u-2)+q^2(5-3u).                    (6)
```

On `Q!=0`, this is exactly THM-3565's rational parameter `h=P/Q`.
Unlike that earlier graph, `(2)` clears the denominator while retaining a
target coordinate because the `a`- and `c`-coefficients are coprime.

The hypersurface `V(G_(P,Q)(F))` is smooth: `F` is etale and `G_(P,Q)` is a
coordinate.  Hence distinct irreducible factors in `(5)` have disjoint zero
sets.  This justifies Euler subtraction between the linear component and
the residual component; it is not an inclusion-exclusion approximation.

## 3. The reduced Jelonek section with vertical corrections

Set

```text
D=3aP^2-bPQ+Q^2,
H=12a^2P^2-4abPQ+(16a-b^2)Q^2,
Delta=(bP-Q)^2+3Q^2,                 E=bP-2Q.             (7)
```

Modulo the target equation `(2)`, the Jelonek polynomial `L` satisfies

```text
Q^6L=D^2H.                                                (8)
```

The identity is polynomial before restriction: `Q^6L-D^2H` is a multiple
of `G_(P,Q)`.  It must not be divided by `Q` on a vertical fibre.  Two
further identities control the cleared strict transforms:

```text
disc_a(H)=64Q^2 Delta,
P^2H+E^2Q^2=4(aP^2+Q^2)D.                               (9)
```

Let `n_R` denote the number of distinct complex roots of a nonzero
polynomial `R in C[b]`, and let

```text
delta_0=[Q(0)=0].                                        (10)
```

Assume for the moment that `Delta` and `E` are nonzero polynomials.  The
projection of the two reduced cleared curves `V_T(D)` and `V_T(H)` to the
`b`-line has the following exact fibres.

* `V_T(D)` has one point away from `V(PQ)`, no point over `V(P)`, and a
  whole affine line over each root of `Q`.  Therefore
  `chi(V_T(D))=1-n_P`.
* `V_T(H)` has two points generically, one over each root of `P` or
  `Delta`, and a whole affine line over each root of `Q`.  The only
  `Q`--`Delta` overlap is `b=0` when `delta_0=1`.  Hence

```text
chi(V_T(H))=2-n_P-n_Q-n_Delta+delta_0.                   (11)
```

* The reduced intersection consists of those `Q`-vertical affine lines
  together with one point over each root of `E` off `b=0`.  Thus

```text
chi(V_T(D) intersect V_T(H))=n_Q+n_E-delta_0.            (12)
```

Their reduced union therefore has the Euler characteristic in `(13)`
below.  This union equals the actual Jelonek section on `Q!=0`, but it has
an extraneous vertical affine line over every nonzero root `b_0` of `Q`.
The actual target equations there are

```text
a=0,                         L=b_0^2(b_0c-1),             (13a)
```

so the true Jelonek fibre is one point.  Replacing the extraneous `A1` by
that point does not change compactly supported Euler characteristic.  At
the possible root `b_0=0`, both the cleared and actual fibres are the whole
`c`-line.  Hence the actual reduced Jelonek section has

```text
chi(T_(P,Q) intersect V(L))
 =3-2n_P-2n_Q-n_Delta-n_E+2delta_0.                      (13)
```

The correction terms and the strict-transform replacement in `(11)--(13)`
are load-bearing.  Treating `P/Q` as a polynomial parameter would erase
precisely the vertical behaviour over `Q=0`; treating `(8)` as an equality
of reduced zero sets would add false vertical components.

## 4. Omitted curve, the linear component, and the Euler formula

Parameterize THM-2473's omitted curve using its nonzero `b`-coordinate:

```text
a=b^2/12,                         c=4/(3b),               (14)
```

Then the exact cube identity is

```text
6b G_(P,Q)|_(14)=-E^3.                                   (15)
```

Thus the reduced omitted-curve intersection has

```text
chi(T_(P,Q) intersect E_omit)=n_E-delta_0.               (16)
```

The linear component `R_(P,Q)=V(qx-pu)` maps isomorphically to the
complement of `V_T(D)`, with the nonzero-`b` vertical lines over roots of
`Q` restored.  If `Q(0)=0`, its zero vertical line has no point because
`u=0` and `B=0` are incompatible.  Explicitly, on `Q D!=0`, put
`h=P/Q`, `D_h=D/Q^2`; the inverse is

```text
x=h/D_h,                  y=b-3ah,
z=aD_h^3-y^2D_h(D_h+3).                                  (17a)
```

At a nonzero root `b_0` of `Q`, the restored vertical inverse is

```text
x=2/b_0,             y=-b_0/2,
z=b_0^2(10-b_0c)/8.                                      (17b)
```

These formulas also show injectivity on every displayed stratum.  Therefore

```text
chi(R_(P,Q))=n_P+n_Q-delta_0.                            (17)
```

Apply THM-3560's exact `3/1/0` Euler integration to the coordinate plane
`T_(P,Q)`, then subtract the disjoint linear component using `(13)`, `(16)`,
and `(17)`.  The residual surface in `(6)` has

```text
chi(S_(P,Q))
 =-3+3n_P+3n_Q+2n_Delta+n_E-2delta_0.                   (18)
```

The residual quadratic over `C(T_(P,Q))` has nonsquare part `H`.  By `(9)`,
it is irreducible whenever `Delta!=0`; hence `(18)` is the Euler
characteristic of one smooth irreducible component in the nondegenerate
scope.

## 5. Collision-compatible rows never pass Euler

The known triple collision value is `(-1/4,0,0)`.  Substitution in `(2)`
gives

```text
G_(P,Q)(-1/4,0,0)
 =P(0)(P(0)-2Q(0))(P(0)+2Q(0))/2.                       (19)
```

Thus collision compatibility is equivalent to

```text
P(0)=0                         or P(0)=+-2Q(0).           (20)
```

At the three collision sources, the factor `qx-pu` contains exactly one
point: the fixed source when `P(0)=0`, the `x=1` source when
`P(0)=-2Q(0)`, and the `x=-1` source when `P(0)=2Q(0)`.  The residual
component contains the other two.  Coprimality and `(20)` also force
`delta_0=0`.

If `(18)` were one, then

```text
3(n_P+n_Q)+2n_Delta+n_E=4.                              (21)
```

When `P(0)=0`, equation `(21)` forces
`(n_P,n_Q,n_Delta,n_E)=(1,0,0,1)`.  But then `Q` is constant and
`P=lambda b^m`, so `Delta=(bP-Q)^2+3Q^2` is nonconstant, a contradiction.
When `P(0)=+-2Q(0)`, either `P,Q` are both constant, in which case
`(n_Delta,n_E)=(2,1)` and `chi(S)=2`, or `n_P+n_Q>=1`.  Equality in the
latter case again requires `Delta` to be constant.  Its factorization

```text
Delta=[bP-(1+i sqrt(3))Q][bP-(1-i sqrt(3))Q]             (22)
```

would make both factors units, hence make both `Q` and `bP` constant,
which is impossible.  Therefore

```text
collision compatible  ==>  chi(S_(P,Q))!=1.             (23)
```

This excludes the smooth irreducible residual component from being `A2`.
It also explains the failure mechanism: the linear factor can isolate only
one collision source, while the vertical/root census charges the component
carrying the collision pair.

## 6. Classification of every Euler-sharp pair

Now drop collision compatibility but retain `Delta E!=0`.  Solving `(18)=1`
gives

```text
3(n_P+n_Q)+2n_Delta+n_E=4+2delta_0.                      (24)
```

There is no solution with `delta_0=0`: the constant pair has
`(n_Delta,n_E)=(2,1)`, while a nonconstant pair would force the impossible
constant `Delta` from `(22)`.  Hence `delta_0=1`.  Equation `(24)` then
has `n_Q>=1`; moreover `(22)` makes `Delta` nonconstant.  Thus
`n_P+n_Q>=2` would make the left side at least eight, and `(24)` forces
`n_P+n_Q=1`.  Hence `P` is constant and `Q=lambda b^m`.  After common
scaling write `P=t`, `Q=b^m`.  Direct root counts give

```text
n_Delta=2m-1,                         n_E=m.              (25)
```

Substitution in `(24)` forces `m=1`.  Consequently every nondegenerate
Euler-sharp pair is, up to common scale,

```text
P=t,                         Q=b,
t!=0,2,                         t^2-2t+4!=0.              (26)
```

In particular, the unique Euler escape abandons the collision: `Q(0)=0`
and `P(0)!=0`, so `(19)` is nonzero.

The two degeneracies omitted while deriving `(24)` introduce no additional
shape.  If `E=0`, then `bP=2Q`; coprimality forces `P` to be a unit and
`Q` to be proportional to `b`, which is the parameter `t=2` below.  If
`Delta=0`, one factor in `(22)` vanishes identically, so the same argument
forces `P` to be a unit and `Q` proportional to `b`; after common scaling
these are exactly `t=1+-i sqrt(3)`.  Thus the next family also exhausts
the degenerate Euler-sharp boundaries.

## 7. The sharp family and why Euler alone is insufficient

Extend `(26)` to every `t!=0` and put

```text
G_t=b^3c-2t^3a+t(t-2)b^2,
G_t(F)=(Bx-tu)S_t.                                       (27)
```

The linear component is exactly a two-torus.  In the invariant coordinates

```text
v=xy,                         w=x^2z,                    (28)
```

the equation `Bx-tu=0` forces `x!=0` and `u=v+1!=0`; it solves uniquely for
`w`.  Conversely `(x,u) in (C*)^2` reconstructs `y,z`.  Hence

```text
V(Bx-tu) ~= (C*)^2.                                      (29)
```

Eliminate `a` from `G_t=0`.  The Jelonek restriction factors as

```text
L|_(T_t)=
 b^2(3bc+t^2-4t)^2
 (3b^2c^2+4t(t-1)bc-4t^2)/(4t^6).                       (30)
```

Its reduced zero set has Euler characteristic one for every `t!=0`: the
line `b=0` supplies one, while every nonvertical component is a `C*`
hyperbola.  At `t=4`, the middle factor becomes `bc`; its extra `c=0` line
meets `b=0` once and leaves the same total.  By `(15)`, the omitted curve is
disjoint from `T_t` unless `t=2`; at `t=2` the whole omitted `C*` lies in
`T_t` and still has Euler characteristic zero.  THM-3560 and `(29)` give

```text
chi(S_t)=1                         for every t!=0.         (31)
```

Thus `t=1`, for example, defeats every Euler-only no-go despite being a
smooth irreducible residual component.

## 8. The torus quotient supplies the missing coordinate

The inherited action has weights `(1,-1,-2)` on `(x,y,z)`.  The polynomial
`S_t` has weight `-2`, so `J_t=x^2S_t` is invariant.  Since every invariant
monomial is generated by `(v,w)` in `(28)`, write the graded source ring as
`A=sum_j A_j`.  Every weight-two monomial is divisible by `x^2`, hence
`A_2=x^2A_0`, and therefore

```text
(S_t) intersect A_0=S_t A_2=(x^2S_t) in A_0.            (32a)
```

Exactness of `C*`-invariants gives

```text
C[S_t]^(C*)=C[v,w]/(J_t).                                (32)
```

Put `U=v+1`.  As a cubic in `w`, the quotient equation has

```text
lc_w(J_t)=-9U^4,
disc_w(J_t)=108t^2U^6 H_t(v),                            (33)

H_t=
 3(t-2)^2(t^2-2t+4)U^4
 +2(7t^3-12t^2+24t+8)U^3
 -12(t^2-2t+3)U^2-24U-4.                                (34)
```

At `v=-1`, the cubic drops to one simple root and
`partial_w J_t=-4`.  Moreover `H_t(-1)=-4`, while the fixed `-24U` term
makes `H_t` nonconstant.  It therefore has a root away from `v=-1`; over
such a root the cubic has at most two distinct points.  If `m>=1` is the
number of distinct roots of `H_t` and the corresponding reduced fibres have
sizes `k_i<=2`, Euler integration over the `v`-line gives

```text
chi(V(J_t))=3(-m)+1+sum_i k_i <=1-m<=0.                  (35)
```

In particular the quotient curve is not `A1`.  When
`t^2-2t+4!=0`, the residual quadratic is irreducible, hence `S_t` is smooth
and irreducible.  If it were `A2`, linearization of its nontrivial
`C*`-action (the equivariant `A2` mechanism used in THM-1345) would make
its one-dimensional categorical quotient `A1`, contradicting `(35)`.
Therefore

```text
t!=0 and t^2-2t+4!=0  ==>  S_t is not A2.                (36)
```

This is the missing sidecar in the Euler-sharp row: the quotient remembers
the weighted orbit coordinate which Euler characteristic destroys.

## 9. Exact parameter boundaries and the `t=1` hostile

The boundary ledger is as follows.

* `t=0`: `G_0=b^3c` is not a coordinate; it lies outside the theorem.
* `t=1+-i sqrt(3)`: `t^3=-8` and `t(t-2)=-4`, so both parameters give the
  same target coordinate

```text
G_t=b^3c+16a-4b^2.                                       (37)
```

  Here the residual quadratic splits.  Thus `S_t` itself is reducible and
  cannot be `A2`, but the present theorem does **not** classify each of
  its two additional irreducible factors as a possible source coordinate.
  This is the exact second-denominator boundary from THM-3565, not evidence
  that a generic residual quadratic splits.
* `t=2`: the omitted curve is contained in `T_t`; the residual is still
  irreducible, and its quotient curve has `chi=-2`.
* Away from `t=2`, `t^2-2t+4=0`, and

```text
9t^3-18t^2+36t-8=0,                                     (38)
```

  the quartic `H_t` is squarefree and the quotient has `chi=-3`.  At the
  three roots of `(38)`, `H_t` acquires one double root, but the associated
  `w`-fibre becomes a triple root; the two corrections cancel and the
  quotient still has `chi=-3`.

The clean hostile is `t=1`:

```text
G_1=b^3c-2a-b^2,
H_1(v)=9v^4+90v^3+192v^2+126v+11.                       (39)
```

The target is a coordinate plane; `S_1` is smooth and irreducible, maps
generically etale of degree two to that plane, and has `chi(S_1)=1`.
Nevertheless `H_1` is squarefree and the quotient has `chi=-3`, so `S_1`
is not `A2`.  This is a hostile to the assumption that smoothness plus the
correct Euler characteristic makes a residual Keller component a plausible
coordinate plane.

## 10. Exact verification and scope

Run

```bash
python3 04-computation/jacobian_coprime_pell_target_coordinate_euler_quotient_no_go_resonant_thm3575.py
python3 -O 04-computation/jacobian_coprime_pell_target_coordinate_euler_quotient_no_go_resonant_thm3575.py
```

The companion verifies `(2)--(9)`, `(15)`, `(27)`, `(30)`, and
`(32)--(39)` symbolically.  It checks the reduced-root formula on explicit
vertical, collision, repeated-root, Euler-sharp, and nonsymmetric controls.
An independent finite hostile enumerates all `532` coprime nondegenerate
ordered pairs of nonzero degree-at-most-two polynomials with coefficients
in `{-1,0,1}`: no collision row is Euler-sharp, and the only four sharp
ordered pairs are the signed representatives `P=+-1,Q=+-b`.  This finite
atlas is a control, not the proof of the all-degree classification.

No planar Jacobian counterexample is constructed.  The theorem closes the
irreducible quadratic residual inside the coprime Pell-graph compiler.  It
does not classify the individual extra factors at `(37)`, a non-Pell target
coordinate, or an irreducible pullback with no rational Pell factor.

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
