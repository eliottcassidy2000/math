# Planar JC(2): squeeze the counterexample through target arithmetic, intrinsic width, and principal edge kernels

**Research synthesis, 2026-08-25.** The planar Jacobian conjecture remains
**OPEN**. This session proves four narrower results: THM-4120 on the smooth
maximum-weight-eight `2:3` response, THM-4122 on intrinsic nonproperness
width, THM-4124 on integral-degree target shears, and THM-4126 on vertical
nonproperness inside that smooth seam. No statement below promotes a formal
scaffold to a global Keller map.

## 1. Inheritance and portfolio

The inheritance pass was:

- **closest proved mechanism:** THM-4103's boundary divisor and five-valued
  response ledger on the first smooth theta seam;
- **canonical hostile:** the quadratic base change
  `q=a^3/2+r^2`, where the target elliptic curve acquires `(a/2,r)` and a
  Mordell--Weil conclusion over `C(q)` ceases to control the quadratic
  punctures;
- **corrected near miss:** THM-3997's all-row Hasse system is necessary at
  collisions but cannot be replaced by a smooth simple-root jet;
- **least-used sidecar:** the normalization degree hidden inside Nguyen's
  polynomial parametrization of a nonproperness component.

The portfolio was:

- **Anchor:** collapse THM-4103's five smooth responses by computing the
  target Mordell--Weil group;
- **Niche:** divide polynomial-cover inflation out of the nonproperness
  parametrization and compare its intrinsic degree with source-direction
  width;
- **Wildcard:** ask whether Lang's polygon similarity synchronizes
  coefficients strongly enough to produce an actual target shear.

All three lanes produced exact statements. The Anchor gives the largest local
squeeze; the Niche and Wildcard change the global normal form in which future
counterexamples should be searched.

## 2. Live concept board

```text
E  target elliptic surface:   sections over C(q), singular fibres, base change;
B  boundary divisor:          punctures, residue fields, ramification indices;
K  edge kernel:               evaluation w->T and (w-T)-adic re-entry cost;
W  intrinsic width:           normalized target pole orders, not cover degree;
S  target shear:              an invertible triangular response operation;
H  collision system:          Hasse rows and normalized wall branches.        (1)
```

The main new interaction is

```text
B --maps to--> E --collapses rational punctures--> exact degree 21
                                      |
                                      v
                         K converts pole order into degree floors.             (2)
```

The global interaction is different:

```text
source escape width --bounds--> W --types--> target nonproperness curves,
Lang polygon scaling --plus coefficients--> S --reduces--> target orbit.       (3)
```

THM-4126 now connects the first arrow to the smooth local seam, but only by
using complete boundary exhaustion. The target-shear arrow still does not
enter that nonintegral `2:3` gauge automatically.

## 3. The smooth maximum-weight-eight seam is much thinner

In THM-4103's smooth nonresonant theta-only branch, with `b=d=0`, the target
generic fibre is

```text
E_q: V^2=U^3-(3a^2/4)U+q-a^3/4.                          (4)
```

Its discriminant is `-432q(q-a^3/2)`. The two finite fibres are `I1`, while
the minimal fibre at infinity is `II*`. The rational elliptic surface has
Neron--Severi rank ten; the section/fibre block together with the `E8` root
lattice already has full rank and discriminant one. Shioda--Tate therefore
gives

```text
E_q(C(q))={O}.                                             (5)
```

Miranda--Persson's `X_211` row independently identifies the same
`II*+I1+I1` surface and records section group order one.

Five boundary punctures are rational over `C(q)` and have total index

```text
1+3+7+3+3=17.                                             (6)
```

The remaining two form one quadratic closed point and contribute either zero
or four sheets. Hence only degrees `17,21` survive target rationality. The
inherited Eisenstein-norm gate excludes `17`, forcing

```text
degree=21,                    every puncture maps to O.    (7)
```

This turns a five-valued response ledger into one value and fixes target pole
divisors of degrees `42,63`.

At the index-seven `DE` puncture, fixed weight layers are polynomials in

```text
w=x^6t^7=s^6t,                  w|_DE=T=-gamma/theta.      (8)
```

Restriction is evaluation at `T`, so the vanishing ideal is principal:
`(w-T)`. The exact source equation gives re-entry order one when `phi!=0` and
order two when `phi=0`; in the latter case the inherited theta-row identity
forcing `epsilon+kappa!=0` is load-bearing. The Keller bracket couples the
two top-layer multiplicities rather than allowing separate cancellation. The
baseline square/cube response and the higher-weight alternative together
raise the necessary total-degree floors to

```text
deg A>=28,                       deg C>=31.                (9)
```

This does not empty the seam. It says any survivor must pay a much larger
global degree bill than its first local jet suggests.

Full boundary exhaustion supplies a second, global consequence. Removing the
target origin and its complete seven-point inverse image makes the generic
affine-fibre map finite of degree `21`. Clearing denominators in monic
integrality equations spreads finiteness away from finitely many values of
`q=E(U,V)`. Hence every target nonproperness component is vertical. THM-4122
excludes every smooth elliptic fibre, leaving only the two nodal `I1` fibres

```text
E(U,V)=0,                       E(U,V)=a^3/2.             (10)
```

On their normalizations `(U,V)` has pole pair `(2,3)`, so `rho=1` and the
global polynomial-degree pair is `(2G,3G)`. The local floor `deg A>=28` now
forces

```text
G>=14,                    (deg A,deg C)>=(28,42),         (11)
```

or `(30,45)` in the higher-`DE`-weight branch. Thus `21` is also the global
function-field degree on this seam, and the possible nonproperness set is
exactly one or both nodal fibres.

## 4. Global counterexample normal form

Let a nonproper Keller pair have degrees `(m,n)=G(d,e)` with `gcd(d,e)=1`,
and let `E=max(d,e)`, `D=max(m,n)`. For each irreducible nonproperness
component, normalize it and write its intrinsic target pole pair as

```text
(rho_C d,rho_C e).                                        (12)
```

Nguyen's raw parametrization multiplier is `M=(cover degree)rho_C`; it is not
intrinsic. Jelonek--Lason's total-degree and fixed-source-direction
parametrizations then give

```text
rho_C E <= min(D-1,w_X(F)).                               (13)
```

Consequences include `G>=2` and, in every affine-linear source chart,

```text
w_X(F)>=max(4,E).                                         (14)
```

For the global reduced pair `2:3`, Lang similarity makes every directional
degree pair `(2r_j,3r_j)`, so `w_X=3r_min`; THM-3544 gives `r_min>=2` and

```text
rho_C<=min(35,r_min)              at degrees (72,108).    (15)
```

A genuinely attained chart with `w_X=6` would therefore force only pole pairs
`(2,3)` and `(4,6)`. This is **CONDITIONAL**. The fact that `(4,6)` is the
first unexcluded width in THM-3586 is not proof that such a chart is attained.

Independently, if one component degree is an integer multiple of the other,
Lang similarity plus leading-form dependence upgrades polygon scaling to a
single coefficient-synchronized shear:

```text
deg Q=r deg P  =>  Q -> Q-cP^r lowers degree
                    and cancels every nonzero outer vertex.                 (16)
```

Therefore a degree-sum-minimal representative under two-sided triangular
target moves has neither degree dividing the other, so its reduced pair has
`d,e>=2`. This is a normalization of the target orbit, not a claim about every
presentation.

The resonant `2:3` expression does not extend `(16)`:

```text
Jac(P,Q^2-cP^3)=2Q Jac(P,Q).                              (17)
```

It is an observable, not an invertible target move. The current extremal
`2:3` seam therefore remains genuinely nonintegral and cannot be shortened by
this shear theorem.

## 5. Exact counterexample-space partition

| region | current squeeze | status | cheapest next decisive test |
|---|---|---|---|
| arbitrary planar Keller pair | choose a target-orbit representative with reduced `d,e>=2`; every nonproper component has `(rho d,rho e)` and obeys `(13)` | PROVED/CITED | search for a source chart with width close to `E`, recording componentwise `rho` |
| reduced `(72,108)` | each chart has `w_X=3r_min>=6` and `rho<=min(35,r_min)`; `rho in {1,2}` only if `r_min=2` is attained | PROVED + CONDITIONAL | compute the smallest `r_min` allowed by the reduced support constraints |
| oriented `2:3`, `b=d=0`, maximum residual weight `8`, smooth nonresonant theta-only | global degree `21`; `S_F` is one/both nodal fibres; `rho=1`; `(deg A,deg C)>=(28,42)` | PROVED RELATIVE + CITED + VERIFIED-EXACT | enumerate rank-21 sheet defects and node monodromy for the two possible components |
| same weight-eight face on three collision walls | smooth response proof unavailable | OPEN | normalize each repeated-root boundary and run the full THM-3997 Hasse system |
| same reduced depth with maximum residual weight `>=9` | no finite support closure | OPEN | enumerate the first new edge-kernel generators and test whether pole-degree 21 still forces a coupled floor |
| other depth cells or chart entry | not reached by the local reduction | OPEN | prove an invariant routing theorem into the oriented B2 chart, or produce a hostile cell |

This table is the honest squeeze. THM-4126 identifies the local and global
`2:3` descriptions only inside the smooth seam by the full-boundary spreading
argument. No such identification is proved on a collision wall, another
depth cell, or before chart entry.

The numerical carriers remain distinct even after that bridge: `21` is the
function-field/generic-fibre mapping degree; `(14,21)` is the pole pair at the
single `DE` source puncture; `(42,63)` is the total source-curve pole-divisor
pair; and `(2,3)` is the pole pair on a target nonproperness normalization.
The bridged value is `rho=1`, not `7` or `21`.

## 6. Connection contracts

### 6.1 Boundary punctures to Mordell--Weil sections

- **source:** boundary closed points of the generic source fibre;
- **target:** `E_q(C(q))` or the Mordell--Weil group over a residue extension;
- **map:** evaluate the generic-fibre morphism at each puncture;
- **preserved:** residue field and Galois orbit, ramification index, image
  section;
- **destroyed:** local Laurent coefficient data after passing only to the
  section class;
- **sidecar:** field of definition of every puncture;
- **cheapest hostile:** quadratic base change `(15)` in THM-4120.

### 6.2 Pole response to an edge-kernel filtration

- **source:** maximum `DE`-weight polynomial layers;
- **target:** the local Laurent expansion at the index-seven puncture;
- **map:** `f(w)` maps to `f(T)` and then to its `(w-T)`-adic re-entry;
- **preserved:** weight, evaluation, vanishing multiplicity, first re-entry
  order;
- **destroyed:** cross-layer cancellation unless the full Keller bracket is
  retained;
- **sidecar:** the exact source equation through the first two `z` terms;
- **cheapest hostile:** a top layer killed by `(w-T)` whose next weight layer
  restores the same Laurent pole.

### 6.3 Escaping source directions to target intrinsic width

- **source:** bounded-degree parametric curves produced from escape sequences;
- **target:** normalized components of `S_F`;
- **map:** factor the parametrization through normalization;
- **preserved:** coordinate pole ratio and degree after multiplication by the
  cover degree;
- **destroyed:** intrinsic multiplier if the cover degree is forgotten;
- **sidecar:** `deg h` in `M=(deg h)rho_C`;
- **cheapest hostile:** precompose a fixed normalization by `t->t^h`.

### 6.4 Newton similarity to target descent

- **source:** `N_0(Q)=rN_0(P)` for integral `r`;
- **target:** the response fibre `Q+C[P]`;
- **map:** use the common top scalar and subtract `cP^r`;
- **preserved:** constant Jacobian and every exposed-vertex coefficient
  relation;
- **destroyed:** lower radial layers and the amount of degree drop;
- **sidecar:** the coefficient identity at every scaled vertex;
- **cheapest hostile:** the scaled-polygon non-Keller pair with vertex ratios
  `1,2` in the THM-4124 certificate.

### 6.5 Complete boundary to vertical nonproperness

- **source:** the finite projective generic-fibre map and its full inverse
  image of the target origin;
- **target:** the global nonproperness locus stratified by `q=E(U,V)`;
- **map:** restrict away from the boundary, then clear denominators in monic
  integrality equations;
- **preserved:** mapping degree and finiteness over a dense pencil open set;
- **destroyed:** which exceptional `q` values are nonproper before fibre
  geometry is reintroduced;
- **sidecar:** equality of the seven punctures with the complete boundary;
- **cheapest hostile:** omit one puncture mapping to a finite target point,
  which can restore a horizontal nonproperness component.

## 7. Highest-value next session

1. **Rank-21 nodal census.** For the three possibilities
   `S_F=N_0,N_1,N_0 union N_1`, enumerate generic sheet loss, node monodromy,
   and THM-3996 address packets in the finite-flat completion. The hostile is
   a disconnected extra address packet invisible to total degree.
2. **Collision-wall normalization.** For each of THM-4053's three walls,
   compute the normalized boundary points, their residue fields, and the
   target elliptic surface after stable reduction. Then ask whether a
   Mordell--Weil or component-group gate replaces the smooth argument.
3. **Full `(w-T)`-adic compiler.** Order all support layers by `DE` weight and
   principal-kernel multiplicity, feed their first two re-entry coefficients
   into THM-3997's all-row Hasse system, and search for a decreasing response
   budget. The positive control is THM-4120's square/cube baseline; the
   hostile is cross-layer pole restoration.
4. **Target-orbit peak reduction.** Decide whether a triangular target word
   can first increase degrees from `(2G,3G)` and later lower their sum. The
   hostile `C^2-cA^3` is resonant but not an automorphism; any lawful word must
   retain the explicit pencil and nonproperness sidecars.

The strategic reframing is now precise: do not search for arbitrary
coefficients of a vast Keller pair. Search for compatibility among three
small objects--a target-orbit-minimal degree pair, a bounded intrinsic
nonproperness curve, and a local edge-kernel response--while keeping residue
fields and collision Hasse rows as mandatory sidecars.
