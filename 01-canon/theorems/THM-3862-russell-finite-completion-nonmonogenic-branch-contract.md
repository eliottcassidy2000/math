---
id: THM-3862
title: "Russell finite-completion nonmonogenic branch contract"
status: >
  PROVED + CITED + INDEPENDENTLY PROOF-AUDITED.  This is a proved conditional
  contract, not an existence theorem: conditional on a hypothetical Darboux
  morphism from the THM-3785 Russell pseudo-plane, its finite normalization
  over the target plane is finite flat, contains the Russell surface as an
  open etale locus, and is not globally monogenic.  Every codimension-one
  branch component has affine normalization A1; over C all branch components
  and the Russell-arm Jelonek component share the unique target point at
  infinity.  If the arm image itself is a branch component, the completion
  degree is at least four.  JC(2) remains OPEN.
source: root / iut_structure_audit / Russell completion and branch lane, 2026-08-23
audit: >
  INDEPENDENT ROOT PROOF AUDIT PASSED after the theorem was written by the
  iut_structure_audit lane.  The audit separately reconstructed finiteness of
  normalization, S subset B, finite type and quasi-finiteness in Zariski Main,
  two-dimensional miracle flatness at every prime, the global monogenic
  derivative-unit contradiction and its D(u) hostile, tame DVR ramification
  over each generic branch component, deleted-divisor valuation transfer,
  the normalization consequence of polynomial uniruledness, Nguyen's
  common-infinity conclusion, and the L+D+E sheet count.  It also checked the
  exact hypotheses against the cited primary sources and retained the
  conditional/existence boundary.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3843-russell-arm-birational-immersion-and-forced-self-identification
  - THM-3849-russell-arm-conductor-polynomial-and-residual-contact-graph
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3847-two-place-cubic-deformation-monogenic-unit-debt
  - THM-3850-nonconstant-cubic-profile-irreducible-branch-puncture-formula
  - THM-3851-tricuspidal-quartic-rank-two-two-place-tradeoff
citation:
  - "Z. Jelonek and M. Lason, Quantitative properties of the non-properness set of a polynomial map, arXiv:1411.5011v2, Theorem 3.2."
  - "Nguyen Van Chau, Note on the Jacobian condition and the non-proper value set, Ann. Polon. Math. 84 (2004), 203--210, DOI 10.4064/ap84-3-2; arXiv:math/0305088."
---

# THM-3862 -- the Russell finite completion must be nonmonogenic with polynomial branches

**PROVED + CITED + INDEPENDENTLY PROOF-AUDITED (conditional contract; JC(2)
OPEN).**  Work over `C`.  Fix `c in C*` and let

```text
B=C[r,z,e]/(r^2e-z^3+c^3r),                 Y=Spec B,             (1)
L=V(r,z) ~= A1_e.                                                   (2)
```

THM-3785 supplies the surjective etale cubic atlas

```text
phi:A2_(x,y) -> Y.                                                   (3)
```

Suppose, conditionally, that there is a Darboux pair `P,Q in B`,

```text
{P,Q}=lambda in C*,                                                  (4)
```

and write

```text
Psi=(P,Q):Y -> A2_(U,V),
R=C[P,Q] ~= C[U,V],
K=Frac(B),
d=[K:Frac(R)].                                                       (5)
```

The bracket in `(4)` is the bracket induced by `(3)`.  Thus `dP,dQ` are
independent everywhere, `Psi` is etale, and

```text
F=Psi phi:A2_(x,y) -> A2_(U,V)                                     (6)
```

is a polynomial Keller map with Jacobian `lambda`.  In particular `P,Q`
are algebraically independent and the identification of `R` in `(5)` is
literal.  Dominance in equal dimension makes `K/Frac(R)` finite.  THM-3785
also gives the all-degree floor

```text
d>=3,                    [C(x,y):C(P,Q)]=3d>=9.                     (7)
```

Let

```text
S=the integral closure of R in K,
Xbar=Spec S,
pi:Xbar -> Spec R=A2_(U,V).                                        (8)
```

Then the following conditional conclusions hold.

1. `S` is finite flat of rank `d` over `R`, and the inclusion `S subset B`
   induces an open immersion

   ```text
   j:Y -> Xbar,                       Psi=pi j.                      (9)
   ```

2. The finite `R`-algebra `S` is not globally monogenic.  That is, there is
   no `theta in S` with `S=R[theta]`.
3. Every irreducible codimension-one component `Gamma` of the branch locus
   of `pi`, equipped with its reduced curve structure, has affine
   normalization

   ```text
   Gamma^nu ~= A1.                                                  (10)
   ```

4. The projective closures of all such branch components, together with the
   arm image

   ```text
   Gamma_L=Psi(L),                                                   (11)
   ```

   meet the target line at infinity at one common point.
5. If `Gamma_L` is itself a branch component of `pi`, then

   ```text
   d>=4.                                                            (12)
   ```

In particular a degree-three finite completion cannot use the arm image as
one of its branch components, even though THM-3843 forces that image to be a
component of the nonproperness set of the composite Keller map `(6)`.

## 1. The finite normalization contains `Y` as an open set

The ring `R` is a polynomial ring over `C`, hence excellent.  Since
`K/Frac(R)` is finite, excellence makes its integral closure `S` finite over
`R`.  Its fraction field is `K`: after clearing the coefficients in an
algebraic equation for any `a in K`, a nonzero `h in R` can be chosen so that
`ha` is integral over `R`, and then

```text
a=(ha)/h in Frac(S).                                                (13)
```

There is an inclusion in the direction needed for `(9)`.  If `s in S`, then
`s` satisfies a monic equation with coefficients in `R subset B`.  It is
therefore integral over `B`.  Since

```text
s in K=Frac(B),                    B is normal,                     (14)
```

one has `s in B`.  Hence

```text
S subset B.                                                          (15)
```

The resulting morphism `j:Spec B -> Spec S` is birational by `(13)`.  It is
also of finite type, a hypothesis which must not be left implicit in the use
of Zariski Main.  Indeed, choose finite `C`-algebra generators
`b_1,...,b_m` of `B`.  Since `C subset R subset S subset B`,

```text
B=S[b_1,...,b_m].                                                    (16)
```

It is affine and therefore separated.  Finally, `j` is quasi-finite.  The
map `Psi` is etale and hence quasi-finite, while `pi` is finite; for every
`xi in Xbar`,

```text
j^(-1)(xi) subset Psi^(-1)(pi(xi)),                                 (17)
```

and the set on the right is finite.  Thus `j` is a separated, finite-type,
quasi-finite birational morphism to the normal scheme `Xbar`.  The
birational form of Zariski Main gives the open immersion `(9)`.

It remains to check flatness rather than silently infer it from finiteness.
The domain `S` is a normal ring of dimension two.  Serre's `S2` condition
therefore makes it Cohen--Macaulay at every prime.  The base `R` is regular
of dimension two, and the finite dominant map `pi` has zero-dimensional
fibres.  Equivalently, after localization at a base prime, all source local
dimensions equal the base local dimension; this also follows from
going-down for the integral extension of the normal domain `R`.  The local
miracle-flatness criterion now applies at every source prime.  Hence `S` is
finite flat over `R`.  Its rank is locally constant on the connected scheme
`Spec R` and equals its generic rank, namely

```text
rank_R(S)=[K:Frac(R)]=d.                                            (18)
```

This proves assertion 1.

## 2. The derivative of a global generator would be a forbidden unit

Suppose for contradiction that `S=R[theta]`.  The kernel of
`R[T] -> S`, `T |-> theta`, is a height-one prime in the UFD `R[T]`.
It is principal.  Since `theta` is integral, the kernel contains a monic
polynomial; its prime generator can therefore be scaled to a monic
polynomial `f`.  Tensoring with `Frac(R)` identifies `f` with the minimal
polynomial of `theta`, so `deg(f)=d`.  Hence there is an isomorphism

```text
S=R[T]/(f),                     T |-> theta.                        (19)
```

In particular `(1,theta,...,theta^(d-1))` is an `R`-basis of `S`.  The
relative differential module is

```text
Omega_(S/R)=S dtheta/(f'(theta)dtheta).                             (20)
```

On the open set `j(Y)`, the finite map `pi` restricts to the etale map
`Psi`.  Therefore

```text
delta=f'(theta)                                                     (21)
```

has no zero on `Y`.  Its image in the affine coordinate ring `B` is a unit.
THM-3785 computes the complete unit group,

```text
B*=C*.                                                              (22)
```

Consequently `delta=a` in `B` for some `a in C*`.  The injection `(15)`
makes the same identity hold in `S`.  But `f'(theta)-a` has degree at most
`d-1` in the power basis `(19)`, and its coefficient of `theta^(d-1)` is
the nonzero scalar `d`.  This is impossible over characteristic zero.
Assertion 2 follows.

Both adjectives in *globally monogenic* are load-bearing.  Finite
extensions remain locally monogenic at many primes, and localization can
create precisely the unit used in `(21)`.  The sharp elementary hostile is

```text
R_0=C[A,C],
S_0=C[u,v],                       A=u^n, C=v,       n>=2,
Y_0=D(u)=Spec C[u,u^(-1),v].                                      (23)
```

The finite completion `S_0/R_0` is globally monogenic, while its restriction
`Y_0 -> D(A)` is etale and has the nonconstant derivative-unit
`n u^(n-1)`.  Thus the contradiction above would fail without `(22)`, and
the argument does not forbid local power-basis presentations of `S`.

## 3. A generic branch divisor supplies an omitted ramification valuation

Let `Gamma` be an irreducible codimension-one component of the branch locus
of `pi`, and let `p` be its height-one prime in `R`.  The localization
`R_p` is a DVR.  The height-one localizations `S_(q_i)` over `p` are the
DVRs in the integral closure of `R_p` in the separable field extension
`K/Frac(R)`.  Write their ramification indices and residue degrees as

```text
e_i=e(q_i/p),                    f_i=[kappa(q_i):kappa(p)].          (24)
```

Finite flatness gives the fundamental equality

```text
d=sum_(q_i over p) e_i f_i.                                       (25)
```

At least one `e_i` exceeds one.  Indeed, if every `e_i=1`, then every
residue extension in `(24)` is separable because its characteristic is
zero.  Each finite DVR map is then unramified, and finite flatness makes
`pi` etale over the generic point of `Gamma`, contrary to the definition of
a branch component.  Choose such a prime and denote its divisor by `E`.
Then

```text
e_E>=2,                         pi(E)=Gamma.                        (26)
```

The divisor `E` is omitted from `Y`.  Its generic point cannot belong to
`j(Y)`, because `pi|_(j(Y))=Psi` is etale.  Since `j(Y)` is open, a nonempty
intersection `E intersect j(Y)` would be a nonempty open subset of the
irreducible divisor `E` and would contain its generic point.  Therefore

```text
E intersect j(Y)=empty.                                             (27)
```

Characteristic zero, finite flatness, and the fact that `Gamma` is a
generic codimension-one branch component are essential here.  In positive
characteristic, inseparable residue extensions can invalidate the inference
from branch to `e_i>1`.

## 4. The deleted valuation is transferred exactly to the plane Jelonek set

Put

```text
U_E=Xbar minus E.                                                    (28)
```

Equation `(27)` gives a morphism

```text
j phi:A2 -> U_E.                                                     (29)
```

It is dominant because `phi` is dominant and `j(Y)` is dense in `Xbar`.
Thus the reusable lemma of THM-3841 applies literally to `(29)`; one does
not have to identify `Y` with the whole complement `(28)`.

For clarity, its valuation mechanism is as follows.  Extend the divisorial
valuation `v_E` of `K` to a valuation `w` of `C(x,y)`.  If `w` had a center
on the affine source of `(29)`, functoriality would give `v_E` a center in
`U_E`, while its center on the separated surface `Xbar` is the generic point
of `E`.  This is impossible.  The restriction of `w` to `Frac(R)` is
centered at the generic point of `Gamma` by `(26)`.  Properness of the
composite near that generic point would give `w` a source center by the
valuative criterion.  Consequently

```text
Gamma subset S_F,                 F=pi j phi=Psi phi.               (30)
```

The nonproperness set of a generically finite polynomial plane map is pure
of dimension one when nonempty.  Hence the closed irreducible curve
`Gamma` is an irreducible component of `S_F`.

Jelonek--Lason, Theorem 3.2, says that `S_F` is covered by polynomially
parametrized curves.  A general point of the component `Gamma`, away from
the other components, therefore lies on the image of a nonconstant
polynomial map

```text
A1 -> Gamma,                                                         (31)
```

which is necessarily dominant.  This forces `(10)`.  Here is the short
curve argument, included to fix the exact conclusion.  The normal source
`A1` lifts `(31)` to the affine normalization `Gamma^nu`.  After completing
both curves, the lift extends to a nonconstant morphism

```text
P1 -> overline(Gamma^nu).                                           (32)
```

The target of `(32)` has genus zero.  Moreover, every projective point
deleted from `Gamma^nu` must have a preimage under the surjective map `(32)`,
but no such preimage lies in the source `A1`.  All deleted points would
therefore have to be images of the single source point at infinity.  There
is at most one; there is at least one because `Gamma^nu` is an affine curve.
Thus its projective completion is `P1` with one point removed, proving
`Gamma^nu ~= A1`.

The complex field is used again for the common target infinity statement.
The map `F` in `(6)` is nonsingular.  Nguyen Van Chau proves that the
projective closure of its nonempty nonproperness curve has exactly one point
on the target line at infinity.  THM-3843 already makes the arm curve
`Gamma_L` an irreducible component of `S_F`; equation `(30)` does the same
for every branch component of `pi`.  Therefore all their projective closures
pass through the same unique target point at infinity.  This proves
assertions 3 and 4.

## 5. Two visible arm sheets plus ramification force degree at least four

Restrict `P,Q` to the arm:

```text
p(e)=P|_L,                         q(e)=Q|_L,
gamma=(p,q):L -> Gamma_L.                                           (33)
```

THM-3843 proves that `gamma` is the finite birational normalization of its
image and has nowhere-zero differential.  Let `Delta(U,V)` be the
irreducible implicit equation of `Gamma_L`.  THM-3849 proves the reduced
pullback decomposition

```text
div_Y(Delta(P,Q))=L+D,
D!=0,                              [D]=2[L] in Pic(Y),               (34)
```

with no component of `D` equal to `L`.

Now localize the target at the DVR prime `(Delta)`.  Because `j` is an open
immersion, the generic points of `L` and of every component `D_i` of `D`
give distinct height-one primes of `S` over `(Delta)`.  The reduced
coefficient of `L` in `(34)` gives

```text
e_L=ord_L(Delta(P,Q))=1,           f_L=1,                            (35)
```

where `f_L=1` is exactly the birationality in `(33)`.  Each `D_i` maps
dominantly to `Gamma_L`: it lies over that curve, and the quasi-finite map
`Psi` cannot contract a divisor to a point.  Reducedness in `(34)` gives

```text
e_(D_i)=1,                         f_(D_i)>=1.                        (36)
```

Suppose `Gamma_L` is also a branch component of `pi`.  Section 3 then
supplies a further prime `E` over `(Delta)` which is absent from `Y` and
satisfies

```text
e_E f_E>=2.                                                         (37)
```

It is distinct from the visible primes in `(35)--(36)`, which are
unramified.  Using one component `D_i` in the fundamental equality `(25)`
gives

```text
d=sum_(Z over (Delta)) e_Z f_Z
 >=e_Lf_L+e_(D_i)f_(D_i)+e_Ef_E
 >=1+1+2=4.                                                         (38)
```

This proves assertion 5.  In degree three, the forced occurrence of
`Gamma_L` in `S_F` is supplied by the nonproper cubic source atlas `phi` and
cannot be reinterpreted as evidence for a ramification divisor of `pi` over
the same curve.  The argument does not classify a possible omitted
*unramified* boundary prime over `Gamma_L`.

## 6. What the arm conductor preserves and what it forgets

THM-3849 defines the unique polynomial `kappa in C[e]` by

```text
Delta_U(p,q)=kappa q',              Delta_V(p,q)=-kappa p'.         (39)
```

It proves

```text
D intersect L=V_L(kappa),
deg(kappa)=2 sum_y delta_y.                                        (40)
```

Thus `kappa` preserves the finite collision addresses of the immersed arm,
the total weighted contact degree incident to each address, and the
aggregate scheme intersection `D intersect L`.  The individual pairwise
edge weights require the completed target-branch factorization in addition
to `kappa`.  The scalar polynomial alone does **not** retain:

1. the allocation of incident weight among distinct edges when three or
   more branches collide;
2. whether `Gamma_L` is a branch component of the finite map `pi`;
3. any omitted ramification prime `E`, or its sheet data `(e_E,f_E)`;
4. any omitted unramified boundary prime over `Gamma_L`;
5. components of `D` disjoint from `L`;
6. the common target point or branch directions at infinity; or
7. the global index form, different ideal, or monogenicity class of `S/R`.

Accordingly, the sheet decomposition

```text
{(Z,e_Z,f_Z,visible/deleted):Z over Gamma_L}                       (41)
```

and a global index-form/different sidecar are both necessary; neither can be
reconstructed from `kappa`.  The exact THM-3849 node, tacnode, and ordinary
triple-point controls make the loss visible: their arm normalizations are
all `A1`, while `deg(kappa)` equals two, four, and six.  Nontrivial conductor
therefore remains compatible with polynomial uniruledness.  Independently,
THM-3844 gives a one-place branch with two cusps and a node whose cubic
completion is nevertheless globally monogenic.  Curve-normalization
conductor and cover monogenicity are different invariants.

## 7. Scope and audit boundary

This is a conditional completion theorem.  It does not construct a Darboux
pair on `Y`, a finite completion with the stated properties, or a planar
Jacobian counterexample.  It gives necessary conditions on any such pair.
Nonmonogenicity is necessary here but is not sufficient; affine-line branch
normalizations and a common infinity point are also necessary but not
sufficient.  The word *branch* always means the reduced irreducible curve
underlying a generic codimension-one component of the non-etale locus of the
finite flat map `pi`.  No claim is made about an isolated non-etale point.

Sections 1--3 and 5 are algebraic in characteristic zero once the THM-3785,
THM-3843, and THM-3849 inputs are available.  Section 4 deliberately works
over `C` in order to use the cited Jelonek--Lason and Nguyen theorems without
an unstated field-transfer argument.  The nonmonogenic conclusion is global
over `R`; local monogenic presentations are explicitly not excluded.

The independently written proof audit recorded in the metadata reconstructed
all five implications and the hostile boundary.  QED.
