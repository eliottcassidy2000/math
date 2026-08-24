---
id: THM-3972
title: "Simple-collision affine-P graph debt has an exact blowup normalization"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. When the THM-3969
  collision polynomial is
  nonconstant and squarefree, every collision is one transverse basepoint
  of the relative-P1 graph. Blowing up those points and deleting the strict
  cubic pole multisection gives the smooth finite normalization. If the pole
  multisection is irreducible and there are n collisions, its class group is
  Z^n and both the canonical and ramification classes equal the primitive
  fibre-section class H. Thus simple collisions exactly repair the old
  Z/3 debt and can pass the canonical-vector gate. Euler rigidity forces
  any possible irreducible pole to have the one-support Kummer form
  a=lambda(t-alpha)^d with 3 not dividing d. Every squarefree row with
  a=t and constant c,r is excluded by its ramification normalization; the
  first r=0 row has a genuine zero-class ramification prime. General
  canonical-compatible simple-collision rows remain open; the pulled-back
  target volume is exact, so de Rham exactness alone does not close them.
source: jc-extra-debt-local / post-THM-3969 simple-resultant collision packet, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-cohn3709, 2026-08-24). The audit checked
  the resultant and infinity root count, transversality, identification of
  the graph closure with the blowup, quasi-finiteness and proper finiteness,
  the preimage of target infinity, exceptional ramification charts, and
  normal-field bridge, the Picard,
  canonical, and ramification classes, the homogeneous incidence identity,
  Euler/Kummer rigidity, every constant-c,r puncture computation and double
  collision seam, and the exact-volume scope. The characteristic-zero Euler
  base-change sentence and maximal-base-ideal graph argument were made
  explicit before promotion. A second all-frontiers hostile pass extended
  the squarefree/irreducible-pole assumption ranges through the Euler gate,
  supplied both local exceptional-ramification Jacobians, and replaced the
  tautological infinity gate by the actual four derivatives. Normal and
  optimized 51-gate replays match the refreshed frozen output and hashes.
depends_on:
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
  - THM-3968-canonical-vector-different-affine-plane-boundary-obstruction
  - THM-3969-affine-p-graph-debt-relative-p1-normalization
  - THM-3966-euler-rigid-split-a1-fibre-affine-plane-boundary
related:
  - THM-3964-polynomial-graph-hidden-double-root-normalization
  - THM-3967-quadratic-p-depth-natural-cubic-conductor-debt-closure
script: 04-computation/jc2_simple_collision_graph_blowup_thm3972.py
output: 05-knowledge/results/jc2_simple_collision_graph_blowup_thm3972.out
script_sha256: 88961bea53a000a39378aed8683ffbe6571119bd354286235954d293d5c879b4
output_sha256: f238e984b667f3973f29a09524db7c86220ffce248184d4e937d7a9860605a97
semantic_sha256: 947829aee929640e6e59e5c050819d81de73f1f972a75330ede8dcb16bef22ea
hash_basis: raw LF bytes
---

# THM-3972 -- a simple collision is one blowup, not an untyped conductor

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Let `a,c,r in k[t]`
and put

```text
y=P-r^2,
q=3rP-r^3+y^2(c+ay),
F=T^3-3PT-q.                                             (1)
```

Use the THM-3969 collision polynomial

```text
Xi=c^3+27a^2r^3-27acr+27a.                              (2)
```

Assume throughout Sections 1--6 that `Xi` has positive degree and is
squarefree. The theorem computes the exact finite normalization at all of
its simple collision fibres. When the cubic pole multisection is
irreducible it also computes the complete class and canonical ledger. This
does not exclude every row in that ledger; Section 7 closes the first sharp
control and Section 9 states the residual honestly.

## 1. The graph on the relative projective line

Put `x=T+r`. Then `(1)` becomes

```text
x^3-3rx^2-3xy-cy^2-ay^3=0.                              (3)
```

On the line `y=vx`, define

```text
D=1-av^3,                    N=3r+3v+cv^2.              (4)
```

Away from `D=0`, the residual intersection is

```text
x=N/D,                       y=vN/D.                    (5)
```

On `S=A1_t times P1_[V:W]` homogenize `(4)` as

```text
D_h=W^3-aV^3,
N_h=3rW^2+3VW+cV^2,
P_h=r^2D_h+VN_h.                                        (6)
```

The target rational map is

```text
phi:S --dashrightarrow P1_P,             phi=[P_h:D_h]. (7)
```

Its base scheme is exactly `D_h=N_h=0`: on `D_h=0` one has
`V!=0`, so `P_h=0` is equivalent to `N_h=0`. Direct elimination gives

```text
Res_v(1-av^3,3r+3v+cv^2)=Xi.                            (8)
```

Thus the roots of `Xi` are precisely the collision fibres, including the
point `[V:W]=[1:0]` when `a=c=0`.

## 2. Squarefree Xi means one transverse point per collision

Consider first a finite basepoint `v`. There

```text
av^3=1,                   r=-v-cv^2/3.                  (9)
```

Let `J_b=det partial_(t,v)(D,N)`. Writing primes for `t` derivatives,
exact substitution in `(2)` gives

```text
Xi'=[(c^2v^2+3cv+9)/v^3] J_b.                           (10)
```

The prefactor is not an accidental denominator. If `omega` is a primitive
cube root of unity, then under `(9)`

```text
N(omega v)N(omega^2v)
 =3v^2(c^2v^2+3cv+9).                                   (11)
```

Consequently its vanishing means that a second pole root is also a
basepoint in the same fibre. Equations `(10)--(11)` show that a simple root
of `Xi` has exactly one finite basepoint and that `D,N` are local parameters
there.

At the infinity point use `w=W/V`. The two rows are

```text
w^3-a,                    3rw^2+3w+c.                   (12)
```

At `a=c=w=0`, equation `(2)` has derivative `Xi'=27a'`.  For the two rows
in `(12)`, the four derivatives are

```text
D_t=-a',  D_w=3w^2,  N_t=3r'w^2+c',  N_w=6rw+3,
```

so their Jacobian at the basepoint is `-3a'`. A simple infinity collision
is therefore transverse as well. Hence if `n=deg Xi`, the base scheme of
`(7)` consists of exactly `n` distinct reduced points.

## 3. Blow up the basepoints and retain the exceptional affine lines

Let

```text
beta:S_tilde=Bl_Z(S) -> S                                (13)
```

be the blowup of the reduced base scheme. Section 2 makes this the blowup
of a smooth surface at `n` distinct points, so `S_tilde` is smooth. It is
also the graph closure of `(7)`: at every basepoint `D_h,N_h` are regular
parameters, so the base ideal `(P_h,D_h)=(N_h,D_h)` is the maximal ideal
and its Rees blowup is precisely the closure of the graph.  Thus `(7)`
resolves to a morphism

```text
Phi:S_tilde -> P1_P times A1_t.                         (14)
```

This morphism is finite. Away from the exceptional curves, a target point
`(P,t)` has preimage cut on `P1_[V:W]` by the binary cubic

```text
(P-r^2)D_h-VN_h=0.                                      (15)
```

The coefficient of `V^2W` in `(15)` is the constant `-3`, so this cubic
never vanishes identically. On an exceptional curve, the two independent
linear parts of `D_h,N_h` give a degree-one map to the target `P1` fibre.
Thus no exceptional curve is contracted. The proper morphism `(14)` is
quasi-finite and therefore finite, of generic degree three.

Let `E=V(D_h)` and let `E_tilde` be its strict transform. A local blowup
chart shows that the inverse image of target infinity in `(14)` is exactly
`E_tilde`: an exceptional curve meets it in one point but is not itself
removed. Hence

```text
X=S_tilde minus E_tilde = Phi^(-1)(A2_(P,t))             (16)
```

is affine, smooth, and finite of degree three over the target. Over each
collision it contains the exceptional curve minus its infinity point, an
affine line mapping isomorphically to the vertical target line.

For completeness, `(16)` maps birationally to the old order `(1)`. On the
finite-`v` chart put

```text
y=P-r^2,                    G=yD-vN,
x=N+av^2y.                                               (17)
```

Then the cancellation-free identities are

```text
vx-y=-G,                    Dx-N=av^2G.                  (18)
```

Thus `T=x-r` is regular even at `D=N=0`. On the infinity chart it is
`T=w(P-r^2)-r`, and the two formulas glue. Generically `v=y/x`, so

```text
k(X)=k(P,T,t),                [k(X):k(P,t)]=3.           (19)
```

Therefore `F` is irreducible, `Spec A` is a domain, and the smooth finite
birational model `X` is its full finite normalization.

## 4. An irreducible pole gives the exact free class lattice

Assume from now through Section 6 that `D_h` is irreducible in
`k[t,V,W]`. Let `H` be the class of a section of `S -> A1_t`, and denote
the exceptional curves by `E_1,...,E_n`. Then

```text
Pic(S_tilde)=Z H direct-sum direct-sum_i Z E_i,
[E_tilde]=3H-sum_i E_i.                                 (20)
```

The coefficient `-1` at every exceptional generator records the transverse
basepoint multiplicity. Weil localization along the single removed prime
`E_tilde` gives

```text
Cl(X)
 =(Z H direct-sum direct-sum_i Z E_i)/(3H-sum_iE_i)
 isomorphic to Z^n.                                     (21)
```

The relation is primitive because `n>=1`. This is the exact topological
transition missed by the collision-free model: before a collision an
irreducible cubic pole leaves `Z/3`; the first exceptional curve turns its
coefficient vector `(3,-1)` primitive and repairs the torsion.

Since `K_S=-2H`, the blowup formula and `(20)` give

```text
K_(S_tilde)=-2H+sum_iE_i,
[K_X]=H in Cl(X).                                       (22)
```

Thus the simple-collision normalization is neither a hypersurface nor
canonically trivial. It pays exactly the nontrivial canonical debt required
by THM-3968.

## 5. The ramification divisor pays the same class

On the chart `(17)`, differentiation at fixed `(P,t)` gives the saturated
ramification row

```text
rho=av^2y+cv^2+2v+r.                                    (23)
```

Its generic projective equation is the quartic

```text
R_h=aV^4+2arV^3W+cV^2W^2+2VW^3+rW^4.                   (24)
```

The exact identity

```text
D rho-[av^4+2arv^3+cv^2+2v+r]=av^2G                  (25)
```

is the required saturation check.  The assertion that an exceptional curve
is not a ramification component needs more than its degree-one restriction,
so fix a finite basepoint `(t0,v0)` and write
`J_b=det partial(D,N)/partial(t,v) != 0`. In the blowup chart
`D=e, N=ez`, the exceptional divisor is `e=0` and
`P=r(t)^2+vz`. Inverting the Jacobian of `(D,N)` gives on `e=0`

```text
partial t/partial e=(N_v-D_vz)/J_b,
partial t/partial z=0,             partial P/partial z=v0,
det partial(t,P)/partial(e,z)=v0(N_v-D_vz)/J_b.         (25a)
```

This determinant is not identically zero because `v0!=0`. Its unique zero is

```text
z=N_v/D_v=v0+2r0,                                      (25b)
```

using `av0^3=1` and `r0=-v0-cv0^2/3`. Identity `(26a)` below pulls back as
`R_h=e(z-v-2r)`, so `(25b)` is exactly the point where the strict
ramification transform meets the exceptional line.

At an infinity basepoint use the chart `N=e, D=ez`. There
`J_inf=-3a'!=0`, `P=r^2+1/z`, and on the finite-target part `z!=0`,

```text
det partial(t,P)/partial(e,z)=-3/(J_inf z) !=0.          (25c)
```

Moreover `R_h/e=w^2-(1+2rw)z`, whose exceptional intersection is `z=0`,
precisely the target-infinity point removed with the strict pole. Thus no
exceptional divisor is a ramification component. The relative ramification
Cartier divisor on `X` is the strict transform of `(24)`, with multiplicities
retained, and

```text
[Ram(Phi)]=4H-sum_iE_i=H=[K_X].                         (26)
```

Equation `(26)` independently recovers finite duality and the canonical
vector of THM-3968. It also shows why the whole simple-collision packet
cannot be rejected by canonical class alone: `H` is primitive in `(21)` and
can be one coefficient-one entry of a boundary basis.

There is a useful compactification identity that records the infinity
exception rather than hiding it:

```text
R_h=W^2N_h-(V+2rW)D_h.                                  (26a)
```

Thus finite intersections of the ramification and pole multisections are
the collision basepoints, where their strict transforms separate.  The
factor `W^2` is load-bearing: over a zero of `a`, ramification can still meet
the pole at `[V:W]=[1:0]` without producing a basepoint.  No global
disjointness is inferred from the finite-`v` identity.

## 6. Euler rigidity forces a one-support Kummer pole

The geometry above supplies a stronger necessary gate. After base change
and descent to `C` (equivalently, using compact-support etale cohomology),
compactly supported Euler characteristic gives

```text
chi(X)=chi(S)+n-chi(E)=2+n-chi(E).                       (27)
```

Indeed each point blowup raises Euler characteristic by one, and blowing up
a smooth point of the pole does not change its strict transform.  If a
Keller `A2` were open in `X`, `(21)` and THM-3922 would give exactly `n`
prime boundary curves.  The Hartogs/affine-intersection purity step in
THM-3966 rules out isolated boundary points.  THM-3920 makes each boundary
curve rational and unibranch, so each has Euler characteristic at most one;
intersections can only lower the Euler characteristic of their union.
Additivity with `chi(A2)=1` yields

```text
2+n-chi(E)=chi(X)<=1+n,
chi(E)>=1.                                               (28)
```

An irreducible affine curve has Euler characteristic at most one.  Equality
in `(28)` therefore forces

```text
chi(E)=1,
normalization(E)=A1,
E has no multibranch identifications.                    (29)
```

This has an exact Kummer translation.  On `V!=0`, the pole is

```text
E: w^3=a(t).                                             (30)
```

For the Euler argument one may first descend the finite coefficient set to
a finitely generated characteristic-zero subfield and embed it in `C`;
geometric irreducibility, normalization, puncture counts, and the displayed
Euler identities are unchanged by algebraically closed base extension.
Equivalently, the entire argument may be read with compactly supported
`l`-adic Euler characteristic.  Thus no unstated assumption `k=C` enters
`(27)--(30)`.

For the projective normalization of its function field, let `B` be the
number of places, including infinity, where the valuation of `a` is not
divisible by three.  Kummer Riemann--Hurwitz gives

```text
genus(Ebar)=B-2.                                         (31)
```

Condition `(29)` forces `B=2` and makes infinity one of the two branch
places, so there is exactly one finite cube-free support point.  A cube
factor at a different finite point would create three normalization
addresses identified on `E`, strictly lowering `chi(E)`.  Consequently any
simple-collision row that can contain a Keller plane must satisfy

```text
a=lambda(t-alpha)^d,          lambda in k*,
d>=1,                         3 does not divide d.       (32)
```

This conclusion is necessary, not sufficient.  It reduces the irreducible
pole search from an arbitrary polynomial to one Kummer tower.

## 7. Every constant-(c,r) row at the first Kummer height is closed

Translate and scale the first case of `(32)` to

```text
a=t,                         c,r in k.                   (33)
```

If `r=0`, then

```text
Xi=c^3+27t,
R_h=V[tV^3+cVW^2+2W^3].                                 (34)
```

There is one collision, so `Cl(X)=Z H`.  The section `V=0` is a
ramification prime of class `H`.  The other factor is primitive and linear
in `t`; on its affine normalization

```text
tv^3+cv+2=0,
v(tv^2+c)=-2.                                           (35)
```

Thus it is a genuine second prime with normalization `G_m` and class
`3H-E_1=0`.  THM-3922 forbids a zero-class boundary prime.  This excludes
all constants `c` in the `r=0` row, including the control `c=1`.

For that control the unique basepoint and the saturated exceptional check
are especially transparent:

```text
(t,v)=(-1/27,-3),             det partial(D,N)=-81,
(1-tv^3)T=v(v+3),             P=vT,
J=T(1+2tv^3)+2v^2+3v.                                  (36)
```

On the exceptional line, `J=3(T+3)` selects one point.  Generically

```text
J=3v(tv^3+v+2)/(1-tv^3),
[K_X]=[Ram]=H+0.                                        (37)
```

Now suppose `r!=0`.  The collision polynomial and its discriminant are

```text
Xi=27r^3t^2+27(1-cr)t+c^3,
disc_t(Xi)=-27(cr-3)^2(4cr-3).                          (38)
```

Hence the simple-collision hypothesis excludes both `cr=3` and
`4cr=3`.  The ramification quartic is primitive and linear in `t`:

```text
t v^3(v+2r)+cv^2+2v+r=0.                                (39)
```

Its two coefficients are coprime because their values at the possible
common roots are

```text
(cv^2+2v+r)|_(v=0)=r,
(cv^2+2v+r)|_(v=-2r)=r(4cr-3).                          (40)
```

Thus `(39)` is irreducible, and its normalization in `X` has at least the
two distinct missing finite parameters `v=0,-2r` (as well as the projective
infinity address).  It is not `A1`.  Here `E:w^3=t` has `chi(E)=1`, so
`(27)` is the equality case of the Euler-rigid boundary theorem: every
Keller boundary prime would have `A1` normalization.  The ramification
prime `(39)` contradicts that requirement.

The two excluded discriminant seams explain the failure boundary.  At
`4cr=3`, `(39)` cancels `v+2r` and

```text
Xi=27(8r^3t+1)^2/(64r^3).                               (41)
```

At `cr=3`, two collision addresses lie over the same fibre and

```text
Xi=27(r^3t-1)^2/r^3.                                    (42)
```

Both are multiple-resultant packets, not counterexamples to the
simple-collision theorem.

## 8. The de Rham exact-volume gate is silent here

For the named `a=t,c=1,r=0` surface, second algebraic de Rham cohomology is
nonzero.  The blowup formula gives

```text
H^2_dR(S_tilde)=k H direct-sum k E_1.
```

The strict pole is `A1`, and its Gysin class is `3H-E_1`; hence

```text
H^2_dR(X)= (kH direct-sum kE_1)/k(3H-E_1) isomorphic to k. (43)
```

Nevertheless the volume form dictated by the fixed target pair is globally
exact on every such finite completion:

```text
dP wedge dt=d(P dt).                                    (44)
```

On a hypothetical Keller source chart `(u,z)`, one would have
`dP wedge dt=lambda du wedge dz` for a scalar `lambda`; equation `(44)`
would provide an exact global extension of that source volume.  Thus the
exact-volume necessary gate passes.  A different nonzero class in `(43)`
cannot be substituted for the actual Keller volume without first proving
that it restricts to `du wedge dz`.  De Rham rank alone does not close the
collision packet.

## 9. Exact residual

For a nonconstant squarefree `Xi` and an irreducible pole multisection, the
normalization and coarse invoices are now fixed:

```text
X=Bl_(n points)(A1 times P1) minus strict pole,
Cl(X)=Z^n,
[K_X]=[Ram]=H,
exceptional collision curves are A1.                    (45)
```

Any Keller boundary would have exactly `n` prime components.  Its pole must
have the one-support form `(32)`.  All distinct ramification primes must
occur among the boundary components, their classes must be primitive and
independent members of that basis, their normalizations must meet the Euler
invoice, and their tame coefficients must sum to `H`.  Section 7 closes the
entire first-height constant-`(c,r)` family.

What remains is the canonical-compatible packet with `(32)` at higher
height or moving `c,r`, at most `n` suitably shaped ramification primes, and
enough possible unramified boundary primes to complete the basis.  Multiple
roots of `Xi` and reducible pole multisections are separate.  None is closed
by `(45)` or by the exact form `(44)`.  This theorem does not prove `JC(2)`.
**QED candidate.**

## Reproduction

```bash
python3 04-computation/jc2_simple_collision_graph_blowup_thm3972.py
python3 -O 04-computation/jc2_simple_collision_graph_blowup_thm3972.py
sha256sum 04-computation/jc2_simple_collision_graph_blowup_thm3972.py \
  05-knowledge/results/jc2_simple_collision_graph_blowup_thm3972.out
python3 agents/check_docs.py
```
