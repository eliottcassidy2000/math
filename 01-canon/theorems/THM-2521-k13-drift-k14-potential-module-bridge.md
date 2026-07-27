---
id: THM-2521
title: "K13 collision drift and the K14 potential-module bridge"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The centred
  thirteen-predecessor profile of THM-2519 embeds, by the signless-incidence
  map p -> (p_x+p_y), into the universal degree-visible potential submodule
  of THM-2514's factor-balanced K_14 edge space.  This submodule has
  dimension 12 on the marked physical slice; ordinary degree is exactly
  12p, its edge norm is 156 times the collision drift, and its degree norm
  is 1,872 times that drift.  In every affine cut chart one aligned Radon
  marginal is injective on the submodule.  At every rational nonuniform
  predecessor fibre the pulled-back defect has all 72 mixed modes and all
  5,184 vector-cut modes.  The six Hamilton columns are exactly the six
  antipodal collision-gap energies; together with the star they give a
  nonconstant rational seven-vector and hence all six additive septimal
  cut characters.  The construction is signed and pointwise/Hilbert-valued,
  the marked vertex is barycentric, and the affine torsor identification is
  chosen rather than semantically supplied.  It does not orient the
  self-cospan, force the multiplicative chi_7 probe, identify source,
  arrival, owner, or deep sheets, exclude a row, or prove LRC(14).
source: codex-2026-07-27-k13-k14-potential-bridge
depends_on:
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2514-cyclic-k14-factor-chart-and-six-phase-ordinary-degree-reconstruction
  - THM-2519-last-digit-collision-drift-and-k13-dirichlet-boundary
related:
  - THM-2509-antipodal-radon-cospan-and-lossless-septimal-chart
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
  - THM-2520-rational-jump-crt-dichotomy-and-delayed-owner-forcing
script: 04-computation/lrc14_k13_k14_potential_bridge_thm2521.py
output: 05-knowledge/results/lrc14_k13_k14_potential_bridge_thm2521.out
script_sha256: 41b6f8b0f9757981bf077661108664fc8d9c8a30d0929d8ca2426fbca560aefc
output_sha256: e41bdc283ba88c3da185784ac846e3bef089494f8da6e89879b870a0dce316a9
hash_basis: working-tree bytes (LF)
---

# THM-2521 -- collision drift enters the degree-visible `K_14` potential module

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2519 produces a genuine physical graph: over a marked future point, the
thirteen last predecessors carry masses whose `K_13` Dirichlet energy is the
last-digit collision drift.  THM-2514 produces a different graph: a signed
`13 x 7` strip is the cyclic one-factor chart of `K_14`.  Equality of the
numbers `13` and `14` did not identify those objects.

There is nevertheless an exact map once one retains the uncontracted
predecessor profile and chooses an affine identification of its `F_13` torsor
with the finite vertices of the static chart.  Centre the thirteen masses,
put zero at the marked vertex `Omega`, and assign to an edge the sum of its
endpoint potentials:

```text
W({x,y})=p_x+p_y.                                             (1)
```

Every perfect matching then has total weight zero.  Thus (1) is already in
THM-2514's factor-balanced space.  More strongly, ordinary degree recovers
the potential by the scalar `12`.  This identifies the precise submodule in
which collision drift lives and proves that THM-2514's large marginal kernel
is not the obstruction on that submodule.

The map is signed.  It uses predecessor masses before the quadratic
self-correlation forgets their phase, and it is naturally valued in a
weighted `L^2` space over the future base.  Those qualifications are the
remaining mathematics, not presentational caveats.

## 1. The universal signless-potential submodule

Let

```text
V=F_13 disjoint_union {Omega},                 |V|=14,
E=E(K_14).                                                    (2)
```

First work over a characteristic-zero scalar field.  Put

```text
V_0={p in K^V:sum_(v in V)p_v=0}.                             (3)
```

Define the signless-incidence map

```text
S:V_0 -> K^E,
(S p)_{ {x,y} }=p_x+p_y.                                     (4)
```

Let `{M_h:h in F_13}` be the cyclic one-factorization of THM-2514.  Since
each `M_h` is a perfect matching, every vertex occurs once and

```text
sum_(e in M_h)(S p)_e=sum_(v in V)p_v=0.                      (5)
```

Hence `S(V_0)` lies in the factor-balanced edge space

```text
mathcal B
 ={W in K^E:sum_(e in M_h)W_e=0 for every h}.                 (6)
```

The thirteen equations in (6) have disjoint edge supports.  They are
independent, so

```text
dim mathcal B=91-13=78.                                       (7)
```

For an edge weighting define ordinary degree by

```text
(Deg W)_v=sum_(e incident to v)W_e.                           (8)
```

If `p in V_0`, then

```text
(Deg S p)_v
 =sum_(u!=v)(p_v+p_u)
 =13p_v+(sum_u p_u-p_v)
 =12p_v.                                                      (9)
```

Thus `S` is injective and degree is a scalar isomorphism on its image.
There is an exact splitting.  Summing (6) over the thirteen matchings gives
`sum_e W_e=0`, hence

```text
sum_v(Deg W)_v=2sum_e W_e=0.                                  (10)
```

Therefore `Deg W in V_0`, and

```text
W=S((Deg W)/12)+[W-S((Deg W)/12)].                            (11)
```

The bracketed term has degree zero.  Consequently

```text
mathcal B
 =S(V_0) direct_sum (mathcal B intersection ker Deg),

dim S(V_0)=13,
dim(mathcal B intersection ker Deg)=65.                       (12)
```

This explains THM-2514's one-phase degree rank `13` and supplies an explicit
right inverse.  It is independent of the cyclic factorization: equation
(5) holds for every perfect matching of `K_14`.

For real or complex Hilbert-valued potentials, all identities hold
coefficientwise.  With the unnormalized vertex and edge norms,

```text
||S p||_E^2
 =sum_(x<y)||p_x+p_y||^2
 =12sum_v||p_v||^2,                                          (13)

||Deg S p||_V^2
 =144sum_v||p_v||^2.                                         (14)
```

Indeed every `||p_v||^2` occurs thirteen times in (13), while

```text
2sum_(x<y)<p_x,p_y>
 =||sum_vp_v||^2-sum_v||p_v||^2
 =-sum_v||p_v||^2.                                           (15)
```

## 2. The physical twelve-dimensional slice

Use the THM-2519 notation at one collision depth.  Its predecessor masses
are

```text
q_r(z)=h((z+r)/13),                         r in F_13,

mu(z)=1/13 sum_r q_r(z),

p_r(z)=q_r(z)-mu(z),
p_Omega(z)=0.                                                (16)
```

Pointwise,

```text
sum_(v in V)p_v(z)=0.                                        (17)
```

The potentials lie in the marked slice

```text
U={p in V_0:p_Omega=0},                  dim U=12.             (18)
```

Let `G>=0` be the future weight and put

```text
mathcal H=L^2(G(z)dz).
```

Then the thirteen coordinate functions `p_r` define one element of
`mathcal H^13`.  THM-2519 equation (17) is

```text
D_13
 =integral G(z)[1/13 sum_r p_r(z)^2]dz.                       (19)
```

Apply Section 1 in `mathcal H`.  If

```text
W(z)=S(p(z)),                                                 (20)
```

then the exact invoices are

```text
integral G||p||^2=13D_13,

integral G||W||^2=156D_13,

integral G||Deg W||^2=1872D_13.                              (21)
```

Equivalently,

```text
D_13=||W||_(L^2(G;E))^2/156
    =||Deg W||_(L^2(G;V))^2/1872.                            (22)
```

Thus positive last-digit drift is exactly a nonzero Hilbert-valued ordinary
degree current on `S(U)`.  No integration of a signed root phase has been
taken in (22), so cancellation over the future base cannot erase this
statement.

The degree vector is especially simple:

```text
(Deg W)_r=12(q_r-mu),                    r in F_13,
(Deg W)_Omega=0.                                            (23)
```

The zero at `Omega` is barycentric.  It is not the value of an owner
indicator.

## 3. Pullback to every affine cut chart

Choose one affine identification of the physical predecessor torsor with
the finite vertex set `F_13`.  Let

```text
tau in F_13^*,
a in F_7^*,
c in F_7                                                     (24)
```

and use THM-2514's chart.  Write

```text
s=ar+c,
rho(s) in {0,1,...,6}                                        (25)
```

for the ordered representative.  Pull (20) back through the edge bijection:

```text
d_p(h,r)=W(Phi_(tau,a,c)(h,r)).                              (26)
```

Explicitly,

```text
d_p(h,r)=p_h,                                      if s=0,

d_p(h,r)=p_(h+tau rho(s))+p_(h-tau rho(s)),         if s!=0. (27)
```

The row sum is zero because the seven entries are the edge weights of the
perfect matching `M_h`.  The column sums vanish too:

```text
sum_h d_p(h,r)=0,                                             (28)
```

since a spoke column sums `p_h`, while a finite column sums two translates
of `p`.  Hence `d_p` is doubly centred and has only mixed finite characters.

The Radon marginal aligned with the same chart has a closed form.  Directly
from (27),

```text
R_(tau,a,c)d_p(v)
 =7p_v+sum_(s=1)^6 p_(v-2tau rho(s)).                         (29)
```

The cut parameters `a,c` disappear because `s=ar+c` merely reindexes the
seven columns.  With

```text
p_hat(alpha)=sum_v p_v zeta_13^(-alpha v),                   (30)
```

equation (29) has multiplier

```text
m_tau(alpha)
 =7+sum_(s=1)^6 zeta_13^(-2alpha tau rho(s)).                 (31)
```

For every `alpha!=0`, the triangle inequality gives

```text
|m_tau(alpha)|>=7-6=1.                                       (32)
```

Thus a single aligned oriented Radon marginal is injective on `S(U)`, for
every one of the `504` affine cut charts.  In normalized `ell^2`, it even
obeys

```text
||R_(tau,a,c)d_p||_2>=||p||_2.                               (33)
```

In particular, the potential submodule intersects THM-2514's
`54`-dimensional paired-marginal kernel trivially.  Replacing `tau` by
`-tau` conjugates (31).  The chart supplies that orientation; the
self-correlation drift does not choose between the two signs.

## 4. Exact mixed spectrum of the pulled-back potential

The bridge is stronger than degree detection.  Let `xi=zeta_7`, let
`a_bar=a^(-1) mod 7`, and use the unnormalized transform

```text
d_tilde(alpha,beta)
 =sum_(h,r)d_p(h,r)
    zeta_13^(-alpha h)xi^(-beta r).                           (34)
```

Put

```text
gamma=beta a_bar,

L(lambda,gamma)
 =1+sum_(s=1)^6 xi^(-gamma s)
       [zeta_13^(-lambda rho(s))+zeta_13^(lambda rho(s))].    (35)
```

Reindexing `r` by `s=ar+c` in (27) gives the exact factorization

```text
d_tilde(alpha,beta)
 =xi^(beta a_bar c)
   L(alpha tau,gamma)p_hat(alpha).                            (36)
```

THM-2514 proves

```text
L(lambda,gamma)!=0                    for lambda,gamma!=0.    (37)
```

Now assume the THM-2519 step data are rational-valued.  At every refinement
cell, `p(z) in Q^13`.  If `p(z)` is nonzero and
`p_hat(alpha,z)=0` for some `alpha!=0`, then

```text
P_z(X)=sum_r p_r(z)X^r                                      (38)
```

vanishes at a primitive thirteenth root.  Since `deg P_z<=12`, rational
cyclotomic irreducibility gives

```text
P_z=c Phi_13.
```

But `P_z(1)=sum_r p_r(z)=0`, so `c=0`, a contradiction.  Therefore, on
every nonuniform rational predecessor fibre,

```text
p_hat(alpha,z)!=0                         for all alpha!=0.   (39)
```

Equations (36)--(39) imply, simultaneously on that same fibre,

```text
d_tilde(alpha,beta,z)!=0
                   for all alpha,beta!=0.                    (40)
```

Thus the one potential defect has all

```text
12*6=72
```

mixed modes pointwise.  Positive `D_13` means the nonuniform locus has
positive `G`-measure, so every mode in (40) is a nonzero element of the
complexification of `mathcal H`.

Applying the THM-2508 cut factor coefficientwise to `d_p` gives all

```text
12*6*12*6=5,184                                             (41)
```

primitive cut coefficients as `mathcal H`-valued functions.  This is an
exact static algebraic lift.  It is not an integrated scalar current: a
signed integral of one function in (40) may still cancel, while its
`L^2(G)` norm cannot.

The distinction between pointwise rational and merely Hilbert-valued data is
load-bearing.  Equations (21)--(36) hold for arbitrary real Hilbert-valued
potentials.  Without pointwise rationality, a nonzero potential can have
some zero root characters, so `D_13>0` alone does not imply (39)--(41).

## 5. The six Hamilton cycles are the six collision gaps

There is also a nonnegative quadratic shadow.  At a fixed future point put

```text
E_0(z)=sum_r p_r(z)^2                                       (42)
```

on the star column.  For `s!=0`, put the gradient energy on the corresponding
THM-2514 Hamilton cycle

```text
E_s(z)
 =sum_h [p_(h-tau rho(s))(z)-p_(h+tau rho(s))(z)]^2.          (43)
```

The six cycles partition all finite edges of `K_13`.  The complete-graph
identity for a mean-zero thirteen-vector is

```text
sum_(r<t)(p_r-p_t)^2
 =13sum_rp_r^2-(sum_rp_r)^2
 =13E_0.                                                      (44)
```

Therefore

```text
sum_(s=1)^6 E_s=13E_0.                                       (45)
```

After integrating with the future weight, write

```text
mathcal E_s=integral G E_s.
```

Equations (19), (42), and (45) give

```text
mathcal E_0=13D_13,

sum_(s=1)^6 mathcal E_s=169D_13,

sum_(s=0)^6 mathcal E_s=182D_13.                             (46)
```

More precisely, THM-2519's needle gaps satisfy

```text
mathcal E_s
 =26[B_0-B_(2tau rho(s))],                   s!=0.            (47)
```

Indeed the right side is the same cyclic squared difference divided by
`2*13`.  Thus THM-2514's six Hamilton columns are literally the six
antipodal last-digit toothpick energies, not only isomorphic graphs.

If `D_13>0`, the rational seven-vector

```text
mathcal E=(mathcal E_0,...,mathcal E_6)                       (48)
```

is nonconstant.  Otherwise every entry would equal `mathcal E_0`, and the
six-cycle sum would be `6mathcal E_0`, contradicting (46), where it is
`13mathcal E_0`.  If the step data are rational, one vanishing nontrivial
additive `F_7` transform of (48) would make its coefficient polynomial a
multiple of `Phi_7`, hence make all seven entries equal.  Consequently

```text
Ehat(beta)!=0                         for every beta in F_7^*. (49)
```

This supplies all six additive septimal cut characters in a nonnegative
quadratic package.  It does not supply an oriented root phase: (43) is even
under converse.

Nor does it force the multiplicative quadratic-character probe.  The sharp
centred-delta control

```text
p_0=12/13,
p_r=-1/13                         for r!=0                     (50)
```

has

```text
mathcal E=(12/13,2,2,2,2,2,2)                                (51)
```

before any common positive scalar weight.  Every normalized nontrivial
additive `F_7` coefficient is `-2/13`, while

```text
sum_(s!=0)chi_7(s)mathcal E_s=0.                              (52)
```

Thus the still-open multiplicative `chi_7`/Fano orientation probe is not
silently solved by (49).

## 6. Representation and gauge ledger

The map has a precise representation-theoretic location.

- `S(V_0)` is the signless-incidence copy of the thirteen-dimensional
  standard representation of `S_14`.
- The physical `U` in (18) is the augmentation representation of the
  two-transitive action

  ```text
  AGL_1(F_13) acts on F_13.                                  (53)
  ```

  It is irreducible over `C`: the translation subgroup supplies the twelve
  nontrivial characters and `F_13^*` permutes them transitively.  It is also
  irreducible over `Q`, because a proper rational submodule would complexify
  to a proper complex submodule.  Under translations alone, the rational
  irreducibility is equivalently the irreducibility of `Phi_13`; over `C`
  that restriction splits into twelve lines.
- `U` is not invariant under all of `S_14`, because an arbitrary vertex
  permutation can move `Omega`.  The containing module `V_0` is.
- An affine change of the physical predecessor coordinate merely relabels
  the finite vertices.  Equations (4), (9), and (21) are
  `AGL_1(F_13)`-equivariant and their norms are gauge-independent.
- Pullback to the strip requires a chart `(tau,a,c)`.  Changing the
  septimal cut permutes which semantic column is the star and which are
  cycles.  Equations (27), (36), and (48) remain covariant only when that
  cut torsor is retained.

The affine identification of the physical predecessor torsor with
THM-2514's static finite vertex set is a choice.  Existing LRC theorems do
not prove that it agrees with a target label, a runner label, a source root,
an arrival root, or a deep root.  Equivariance says different affine choices
are controlled gauges; it does not manufacture the missing semantic map.

## 7. Exact gain and sharp stopping boundary

The proved connection is now:

```text
positive rational last-digit collision drift
  -> nonzero centred predecessor profile in L^2(G)
  -> 12-dimensional signless-potential submodule of K_14
  -> one nonzero ordinary degree and one injective aligned Radon marginal
  -> all 72 pointwise mixed modes and all 5,184 H-valued cut modes
  -> all six nontrivial additive septimal column-energy characters.       (54)
```

This materially removes one false candidate obstruction.  The
`54`-dimensional paired-marginal kernel and the need for six degree phases in
the full `78`-dimensional chart do not obstruct this physical potential
component.  It lies in the easiest degree-visible summand, and even one
aligned phase recovers it.

The theorem does **not** close the live LRC seam.

1. Positive `D_13` on an actual lawful response is not proved here.
   THM-2520 turns that question into an endpoint-current dichotomy and a
   delayed-owner forcing test; its charged branch still has to be verified
   for the live packet.
2. The map uses the uncontracted vector `q_r`.  The scalar drift and the
   even autocorrelation vector `B_u` do not determine `p`; reconstructing it
   from them would be a phase-retrieval error.
3. `p`, `W`, and `d_p` are signed centred masses.  They are finite rational
   linear combinations of Boolean ancestry factors, not measures of one
   Boolean owner/arrival event.
4. The nonnegative energies in Section 5 are squares of conditional masses,
   not idempotent Boolean intersections.
5. `p_Omega=0` is a barycentric convention.  The future owner event `G`
   weights the Hilbert norm but does not become a signed owner-loop charge.
6. A chart orientation `tau` makes (29) oriented; the self-cospan remains
   antipodally even and does not select `tau` against `-tau`.
7. Source, arrival, target, old deep, and future-owner sheets remain distinct.
   No role-wise ancestry carry or temporal atom intertwiner is supplied.
8. An `L^2`-nonzero mixed mode need not have a nonzero signed scalar integral.
   Taking its norm restores positivity but erases the phase needed for an
   owner-loop current.

Accordingly no scalar row is excluded and LRC(14) remains open.  The next
decisive bridge is not another spectral census.  It is a Boolean or polarized
pairing which couples one of the signed potential modes in (40) to the actual
owner/source/deep ancestry current without collapsing its cut phase.

## 8. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_k13_k14_potential_bridge_thm2521.py
python3 -O 04-computation/lrc14_k13_k14_potential_bridge_thm2521.py
```

Both executions must reproduce

```text
05-knowledge/results/lrc14_k13_k14_potential_bridge_thm2521.out
```

byte-for-byte.  The companion uses exact integer and rational arithmetic and
the finite field `F_547`, which contains primitive seventh and thirteenth
roots.  It verifies:

- matching-constraint rank `13`, factor-balanced dimension `78`, degree rank
  `13`, and degree-kernel dimension `65`;
- the full and physical potential dimensions `13` and `12`, together with
  the edge and degree Gram scalars `12` and `144`;
- all `504` affine charts on all twelve physical basis vectors, including
  `78,624` aligned-Radon entries;
- `435,456` direct mixed-multiplier identities, with every geometric factor
  and mixed mode nonzero;
- all `5,184` THM-2508 cut-lift factors on a physical potential control;
- `1,728` coefficientwise six-cycle quadratic-form identities; and
- the exact centred-delta boundary (51)--(52).

Normal and optimized runs reproduce the stored transcript exactly.  The
finite-field checks referee the algebraic identities; nonvanishing in the
theorem is proved over the cyclotomic fields by (32), (37), and rational
irreducibility, not inferred from reduction modulo one prime.

An independent line audit rederived the factor-balanced dimension and direct
sum, both Hilbert norm constants, the predecessor centring, the chart
pullback, the aligned-Radon multiplier and its lower bound, the direct mixed
cyclotomic factorization, pointwise rational all-or-all law, the constants
`26`, `169`, and `182` in the six-cycle energy ledger, the additive
`F_7` dichotomy, and the rational/complex representation claims.  It also
confirmed that the affine torsor choice, barycentric `Omega`, signed
linearization, Hilbert-versus-scalar distinction, antipodal orientation loss,
and multiplicative `chi_7` hostile are all explicit stopping boundaries.
**QED.**
