---
id: HYP-6815
title: LRC14 four-far cone and affine-slope threshold suspension
status: EXACT REPRESENTATION LEMMAS + EXACT AUDITS; THM-761 high-support sheet descent proved; raw fragmentation splice refuted inside f=4; peel-relative compactification open
source: codex-2026-07-14-S2
script: 04-computation/lrc14_affine_slope_suspension_codex_S2.py
result: 05-knowledge/results/lrc14_affine_slope_suspension_codex_S2.out
related:
  - HYP-6780
  - HYP-6785
  - HYP-6820
  - HYP-6830
  - HYP-6755
  - HYP-3106
  - HYP-3072
  - HYP-3025
  - HYP-3007
  - THM-509
  - THM-555   # tiling-cycle-moment hierarchy (not the colliding insertion-DP file)
  - THM-668
  - THM-738
  - THM-741
  - THM-742
  - THM-755
  - THM-760
  - THM-761
  - THM-765
  - THM-764
  - THM-771
---

# HYP-6815: LRC14 Four-Far Cone And Affine-Slope Threshold Suspension

## Verdict

There is one literal four-dimensional object and a second, exact
four-coordinate description that must not be confused with it.

1. The first unresolved far-count chart is literally four-dimensional.  After
   THM-738 closes `f<=3`, fix a nine-speed core `C subset {1,...,14}` and let
   the four ordered far speeds vary.  Their gap coordinates form a rank-four
   lattice chart in an ordered orthant.
2. A broader affine family has a mixed four-coordinate suspension
   `(u,t,c,lambda)`: two torus phases, a discrete integer slope, and a
   clearance level, subject to `u=ct`.  For fixed `c` it has only the two
   continuous coordinates `(t,lambda)`.  It is not another four-dimensional
   space, and a general affine family need not stay in one four-far chart.

There is now also a first *proved recursive chart change*, not a third
dimension claim.  THM-761 sends a scaled-core packet `cD union W` to the
finite sheet cycle `Z/cZ`: the core margin is exactly constant on every
sheet, while the exceptions become an inhomogeneous residue-cover problem.
This closes high-support structured scale packets with at most six exceptions.
It does not close the literal four-far chart: on `C union kF`, the four far
speeds are the scaled core and the nine fixed small speeds are exceptions, so
the sheet packet has `r=9`. HYP-6830's proposed raw-fragmentation bridge to
the scale-free/small-scale regime is now refuted inside this same four-far
chart; a peel-relative replacement remains open.

The global moduli space of 13 speeds modulo dilation is not four-dimensional.
Calling all of LRC14 a four-manifold would discard the higher far-count strata.
A true global four-dimensional reduction still requires a peel or descent
theorem.

The exact object on the first chart is therefore not a naked point of `Z^4`.
It is a semilinear rank-four cone chart carrying, fiber by fiber, an
owner-colored threshold loop, rational gap metric, embedded marked relation
lattice, and the action of the next permitted proof operation.  Uniform wall
stratification and transport of those fibers over the chart remain open.

## 1. The outer four-far cone

Fix `C subset {1,...,14}`, `|C|=9`.  In gap coordinates write

```text
w1 = a
w2 = a+d1
w3 = a+d1+d2
w4 = a+d1+d2+d3,
```

with `a>=15` and `d1,d2,d3>=1`.  Let `X_C` be the points
`g=(a,d1,d2,d3)` for which

`S(C,g)=C union {w1,w2,w3,w4}`

is covering.  Every nine-element `C subset {1,...,14}` already has gcd one:
no integer `d>=2` has nine positive multiples at most 14.  Primitivity is
therefore automatic on this exactly-`f=4` chart.  The union of these charts
over the `binom(14,9)=2002` possible cores is exactly the covering,
exactly-nine-small `f=4` stratum, where the core is unique.  THM-741 quantifies
over rows with *at least* nine small speeds, so its 2002 bodies overlap on the
`f<=3` strata; those bodies are not a disjoint model of its full statement.

Let `Q=lcm(2,...,14)=360360`.  Covering is a Boolean combination of the
congruences `q|wj`, `2<=q<=14`.  Hence:

> **Semilinear cone lemma.** `X_C` is a finite union of residue classes
> modulo `Q`, intersected with the positive ordered gap cone.

The quotient has finitely many residue addresses modulo `Q`; `X_C` itself is
infinite, and this is not a finite LRC classification.  The archimedean
heights still move endpoint walls and pair-sum moduli inside a fixed residue
class.

> **Far-dilation monoid lemma.** If `g in X_C` and `k>=1` is an integer, then
> `k*g in X_C`.  Indeed `w_j(k*g)=k*w_j(g)`, so every divisibility witness
> survives.

Thus "cone chart" does not mean an ordinary real convex cone: `X_C` need not
be addition-closed, though it has this exact positive-integer action.  The
action scales only the four far coordinates while fixing `C`, so it is not
global speed dilation and does not preserve exact clearance.  The exact
canary

```text
C={1,...,9}, far=(22,26,28,60)       M=6/61
C={1,...,9}, far=(44,52,56,120)      M=12/121
```

is covering in both rows.  The action on the fiber is therefore data for a
future descent, not a quotient that may be erased.

For `g in X_C`, define the cocharacter

`phi_{C,g}(t)=(s*t mod 1)_{s in S(C,g)}`

and the closed safe cube `K=[1/14,13/14]^13`.  The universal incidence

`Z_C={(g,t): phi_{C,g}(t) in K}`

projects to `X_C`.  Surjectivity of `Z_C -> X_C` for every core `C` is exactly
the remaining primitive-covering, exactly-`f=4` claim after the dilation and
noncovering reductions (including THM-366).  It is not by itself a statement
about the separately dispatched noncovering rows or the overlapping `f<=3`
part of THM-741.

## 2. The inner affine-slope suspension

Let `A=(a_i)` and `R=(r_i)` be integer 13-vectors and let

`V(c)=(c a_i+r_i)`,

for every integer scale `c` at which the coordinates are positive and
distinct.  Primitivity is a separate hypothesis when a covering reduction
needs it; the identity below does not require it.  On the two-torus set

`Phi_{A,R}(u,t)=min_i ||a_i u+r_i t||`.

For the slope-`c` closed geodesic `L_c={(ct,t):t in R/Z}`, direct substitution
gives the exact identity

`M(V(c)) = max_t Phi_{A,R}(ct,t)`.

Thus the mixed four-coordinate incidence object

`X_{A,R}={(u,t,c,lambda): u=ct, Phi_{A,R}(u,t)>=lambda}`

in `(R/Z)^2 x N x [0,1/2]` has the LRC14 assertion as nonemptiness of its
integer-slope fiber at `lambda=1/14`.  It is a stratified continuous/discrete
object, not a four-dimensional manifold.  Runner owners label its strip
strata.  This affine model is broader than a fixed `X_C`: a ray stays in one
core chart only under additional zero-slope/fixed-core conditions, and the AP
dilation and shear audits below intentionally move outside that chart.

When `R=0`, `Phi` is independent of `t`; the suspension is cylindrical and
recovers HYP-6780's dilation law.  Nonzero `R` is transverse shear/holonomy,
not a small perturbation that may be discarded.

This generalizes THM-742's exact slope-geodesic formulation for
`B union (W+J)`.  THM-742 is one polygonal chart; HYP-6815 retains arbitrary
owned offsets and the clearance filtration.

### The first proved recursion: the sheet chart

For a scaled-core packet

`V=cD union W`,

where no exception `w in W` is divisible by `c`, choose a core time `t0` and
write

`t_k=(t0+k)/c`, for `k in Z/cZ`.

Then for every `d in D` and every sheet `k`,

`||(cd)t_k||=||d t0||`,

while each exception has phase

`w t_k=(w t0)/c + wk/c mod 1`.

Thus the core endpoint movie may be *annihilated legally* for this operation:
its required margin is identical on every sheet.  The exception carrier must
retain the common sheet index, owner-labelled residues `w mod c`, offsets
`(w t0)/c`, gcd multiplicities `g_w=gcd(w,c)`, and the closed threshold.  It
cannot split the exceptions into independent clocks, because their bad sets
must be unioned on the same `k`.

THM-761 proves that a free sheet exists whenever

`sum_{w in W} g_w*(floor(c/(7g_w))+1) <= c-1`.

For LRC14 this closes every one-exception packet at every scale and, for at
most six coprime exceptions, every `c>=43`.  The former `c=26` / `q<=25`
canary is closed immediately on its native 26-sheet clock even though its
first ordinary rational witness has denominator 27.  The surviving sheet
debts are exact: seven exceptions can tile the sheet cycle, excessive gcd
multiplicity calls for `c -> c/g` descent, and a peel-relative small-scale/
scale-free splice is HYP-6830's corrected open content.

For the fixed-core far action in Section 1, `C union kF` can be written with
scaled core `D=F` and exception set `W=C`.  Thus the exact sheet
reparameterization still applies, but `|W|=9`: the free-sheet union bound is
beyond its sharp `r<=6` range.  This locates the literal four-far chart inside
the open tiling/overlap residue rather than claiming THM-761 has solved it.

THM-771 makes a second, transverse route to infinity inside the literal chart exact.
For prime `N>110`, set

```text
V_N={1,...,9,15,110,N,1092N}.
```

These rows are primitive, covering, and exactly `f=4`. No divisor is shared
by seven speeds, but the good-set component count of the top-peeled core
`{1,...,9,15,110,N}` is unbounded: the fixed safe interval
`[1/14,111/1540]` is cut by at least `N/1540-8/7` full `N`-teeth. Exact counts are
`66,104,174,310` at `N=211,503,1009,2003`. Thus high-support coherent dilation and
high-frequency tooth insertion are independent directions of noncompactness.
The scalar `c*` sees the first only when at least seven speeds share the scale and
is blind to the second; it can also miss radial dilation of the four far coordinates.

This does not threaten LRC14: the audited `V_N` rows fire the capped-envelope
peel. It identifies a missing coordinate. A transport state for the splice
must retain at least the normalized cap load
`r_P/(v|G'_P|)=pi*v*(P)/v`, rather than `r_P` alone, together with the
divisor-support profile and endpoint/tooth ownership.

There is also an exact compressed local coordinate.  At a fixed time put
`x_i={s_i t}`, `kappa_i=floor(14x_i)`, and `rho_i={14x_i}`.  If
`gamma_ij=floor(14{x_j-x_i})`, then

`gamma_ij=[kappa_j-kappa_i-1_{rho_j<rho_i}] mod 14`.

Thus the nominal edge-sector table is a vertex potential vector plus the
inversion matrix of one global weak order.  This is a better local carrier than
a binary tournament, but it still needs wall flags: the safe endpoint
`x_i=13/14` shares its one-sided sector with dangerous points immediately
above it.  The sector potential and weak order do not replace rational gaps,
scale/residue transport, or marked relations.

## 3. The exact dual resonance picture

For any coefficient vector `z in Z^13`, define

```text
m(z)=z dot A,
n(z)=z dot R.
```

The torus character attached to `z` is

`chi_z(u,t)=exp(2*pi*i*(m(z)u+n(z)t))`.

On the slope-`c` fiber it becomes

`chi_z(ct,t)=exp(2*pi*i*(c*m(z)+n(z))*t)`.

Therefore:

> **Resonance-line lemma.** `z` is an integer relation of `V(c)` exactly
> when its projected point `(m(z),n(z))` lies on `c*m+n=0`.

For a finite trigonometric polynomial, averaging on the slope fiber selects
exactly the full marked coefficient vectors `z` whose projections lie on that
line.  Bare points `(m,n)` are not enough: the projection can have nontrivial
fibers, and Fourier mass depends on their multiplicities and coefficient
weights.  With a justified Abel/Fejer limit, this can recover the
one-dimensional Haar measure of the pulled-back safe indicator.  It does not
recover closed-threshold nonemptiness or an isolated equality witness:
AP/Goddyn-Wong equality cases can be lonely with zero safe measure.  No
absolute-convergence or pointwise-boundary claim is made here.

This identifies the repository's relation-lattice, Fourier, and scale-residue
languages.  Shape keeps only `m`; residue keeps only `n`; the LRC slice uses
their pairing with `c`, while measure needs the marked weighted preimage over
the selected line.  The exact audit finds that the AP shape has 36
support-three relations `e_i+e_j-e_k`, while a single owned offset kills from
6 to 11 of them, depending on its owner.

After fixing the owner labels and target anchor, the full embedded marked
lattice

`ker(Z^13 -> Z, z |-> z dot S)`

determines a primitive positive speed vector up to sign: it is the labelled
rank-12 kernel inside the fixed `Z^13` of a primitive integer functional,
whose normal is unique up to sign; positivity selects the sign.  For a
nonprimitive vector the kernel determines only its primitive normalization,
which preserves `M` but not component count under dilation.  An abstract
isomorphism class, rank, support enumerator, or successive minimum does not
retain that labelled normal or its coordinate-facet pairing with the safe
cube.

## 4. What information is actually sufficient?

No nontrivial task-independent minimal quotient has been demonstrated.  The
identity representation is universally sufficient, but a useful compressed
carrier depends on the next question.

### Fixed-threshold truth

Let

`beta_S(t)_i = 1_{||s_i t|| < 1/14}`.

Starting at `t=0`, record the cyclic sequence of endpoint blocks, with each
block retaining entering owners, exiting owners, and exact coincidences.  If
`B_-` is the danger-owner set immediately before a tie block, with exiting
owners `E` and entering owners `A_in`, then exactly

`B_tie=B_- \ E` and `B_+=(B_- \ E) union A_in`.

Thus the owner-colored event word reconstructs both the piecewise-constant
open-cell loop `beta_S` and every zero-dimensional equality state.  It decides
whether an open cell or boundary tie state hits `0^13`.  Exact gap lengths are
unnecessary for this Boolean predicate.

An aggregate signed current is insufficient: an owner entering exactly when
another exits can cancel in the scalar current even when the tie state carries
the only equality witness.

### Margin, measure, and autocorrelation

To recover exact `M`, safe measure, component count, or THM-731/732, add the
rational endpoint phases and gap lengths, or equivalently retain the full
threshold filtration `lambda |-> G_S(lambda)`.  The signed owner endpoint
divisor reconstructs the Bernoulli `B2` discrepancy.  An exact peak
certificate `(t*,M)` recovers `M` alone only when it includes a
maximality/upper-bound certificate, not merely the lower-bound witness `t*`.
It still does not recover the safe measure or topology.

The exact 552-row endpoint audit supplies sharp canaries:

```text
AP                         M=1/14
{1,...,12,26}              M=2/27
```

They have the same endpoint tournament, divisor mask, and cap sign. Thus
order, covering, and the THM-755 bit do not preserve exact clearance or
metric data.

### Covering and capped-envelope status

Covering needs the divisor mask.  For peeled speed `v_peel` and peeled core
`B`, the THM-755 decision needs the projective comparison
`v_peel*|G'_B|/r_B > 1/pi`, equivalently `pi*v_peel*|G'_B| > r_B`, not endpoint
orientation.  A lift `21 -> 22` in the audit changes the cap decision with zero
endpoint-tournament edge flips.

### Transport across rows

Deletion and peel require owner labels.  Dilation and affine transport require
the scale/residue action on the fiber, not merely a normalized point.  Fourier
or interaction arguments require the embedded relation coefficients and
heights.  Proof assembly additionally needs the next observer, available
certificate, legal discharge, or named residual debt.

THM-765 supplies a sharp transport model: when one runner covers every
super-`L` component of a core, each component's length and midpoint phase must
fit one tooth, and all gcd-deck translations must be retained.  For `L<1/4`,
a genuine margin drop forces the core gcd to divide the added runner.  Neither
component length alone nor a residue label alone preserves that conclusion.
The global hereditary-primitivity corollary is automatic on the literal
four-far chart: deleting a far speed leaves the gcd-one nine-core, while
deleting a small speed leaves eight distinct integers in `[1,14]`, whose gcd
is already one. THM-765 is therefore a carrier model here, not an `f=4`
pruning theorem.

The resulting hierarchy is:

```text
truth carrier    = initial state + signed owner event blocks
metric carrier   = truth carrier + rational phases/gap filtration
functor carrier  = metric carrier + scale/residue action + marked relations
proof carrier    = functor carrier + next observer/certificate/discharge
```

## 5. Geometry-side and obstruction-side finite presentations

The owner-colored cubical loop is a geometry-side presentation.  HYP-6785's
endogenous pair-sum blocker complex is an obstruction-side presentation:
pair-sum peak candidates carry blocker sets, and LRC holds exactly when some
candidate has an empty blocker edge.  Both are independently equivalent to
fixed-row LRC truth.  No pairing, chain map, or transport theorem between them
has been constructed, so calling them a proved primal/dual pair would be too
strong.  The blocker presentation is finite at fixed `S`, but a scale-normal
pressure or descent theorem over `X_C` remains open.

This suggests a bicomplex rather than a flat ledger:

```text
d_geom  = cross an endpoint/contact wall
d_arith = change a residue, valuation, or scale depth
```

Building these operations and a comparison map is a proposed program.  If a
future exact quotient makes them fail to commute, that commutator would locate
a missing owner, height, relation, or certificate sidecar; the displayed
bicomplex is not yet a theorem.

## 6. Exact executable evidence

`lrc14_affine_slope_suspension_codex_S2.py` verifies:

1. the exact slice identity for pure AP dilations;
2. `|G_{cA}|=|G_A|` and `r_{cA}=c*r_A` for `A={1,...,12}`;
3. fixed-core far-coordinate dilation preserves covering but changes exact
   `M` from `6/61` to `12/121` on the displayed core-nine canary;
4. the transverse ray `c{1,...,13}+e_13` at covering scales, with exact
   `M=1/13`;
5. the projected support-three relation loss under owned shear;
6. a metric collision at equal shape, owned offset, and `c mod 14`:
   `c=2` and `c=16` have the same `M=2/15` but different exact safe measures
   `5/84` and `115/1904`, and component counts `4` and `30`;
7. a representation tournament under predicate-first and compression-first
   gauges.  Both are transitive, with 11 edge flips between gauges.

The larger endpoint/tournament audit in
`lrc14_endpoint_tournament_sidecar_audit_structure_S1.py` uses 552 exact rows,
checks 552/552 Bernoulli reconstructions, and supplies the collision counts
quoted above.

THM-764 makes the concurrent HYP-6820 denominator canary exact.  For every
`15<=q<=28`, the task-specific rational-witness carrier is precisely a
zero-owner set plus the owner-resolved signed-unit-pair deck; a witness exists
exactly when there is no zero owner and some signed unit pair is missing.  The
`c=26` transverse ray blocks every denominator `15<=q<=25` and first witnesses
at `q=27`; a gcd-incoherent uncapped residual blocks the same window and first
witnesses at `q=26`.  The deck preserves this bounded rational predicate, not
full LRC truth or metric data.

THM-761 then resolves the apparent conflict: the coherent `c=26` row is
dispatched by its scale-native sheet deck, without any bounded ordinary
denominator claim.  The proof carrier is therefore an atlas of exact charts
(signed-pair decks, scale sheets, endpoint words, blocker complexes), with an
explicit legal transition or residual debt, not the Boolean field "some
`q<=25` works."

## 7. Transfers from seemingly unrelated threads

- Tournament switching and tiling: preserve the group action or gauge
  section, not merely a bijection of underlying sets.
- Coding and matroids: a weight enumerator preserves support counts but loses
  support connectivity, coefficient height, and realizability.  LRC likewise
  needs the embedded marked relation code.
- Ising and Lee-Yang: a first moment or `p0` loses compatible many-body
  packing.  Root surfaces are candidate phase classifiers in the finite banks
  audited in that thread, not a proved global classifier or discharge.
- Resolvent folds: the symmetric inner page is useful only with center,
  antisymmetric leakage, and boundary sidecars.  The proposed LRC dictionary
  treats shape as the inner page and residue shear as leakage; no exact
  resolvent-to-LRC map is claimed.
- CRT and p-adic trees: residue skeleton and valuation/height flex are separate
  coordinates.  Fixed residue does not freeze endpoint geometry.
- Canonical metagraph addresses (HYP-6825): a finite rooted tiling address can
  organize the combinatorial base, but threshold walls, metric gaps, owners,
  and certificate obligations must remain in a stalk. The `V_N` ray makes one
  such stalk visibly noncompact: its divisor-scale type stays bounded while
  new owner-labelled teeth subdivide the safe fiber.
- Safe-component tooth decks (THM-765): component width, midpoint phase, and
  gcd translation sheets form one exact containment condition.  This is a
  proved instance of a metric carrier plus arithmetic action being smaller
  than the full row but larger than either of its projections.
- Observer categories: an anchored cyclic order, metric gap widths, and the
  observer placement are distinct fields.  The exact LRC predicate is marked;
  an observer-blind order class can mix safe and unsafe placements.
- By analogy, THM-509's tournament baby-Hodge integrality gaps warn that a
  continuous relaxation can contain integer-unrealizable holes.  No analogous
  hole has been exhibited for the LRC torus suspension; the methodological
  consequence is to audit a realizability or carry syndrome before using a
  continuous certificate to discharge an integer row.
- The tiling-cycle-moment THM-555's cut/cycle boundary warns that aggregate
  scores and low cycle counts can survive while higher overlap and packing
  data are lost.  Only the local invariant views tested there justify this
  transfer, not a universal failure of local-to-global reasoning.
- Observer-cut and proof-circuit ledgers supply an audit discipline: relative
  to a named next operation, every changed predicate should be reconstructed,
  annihilated, descended, boundary-stopped, or routed to named debt.  This is
  not a universal quotient theorem.

## 8. Open theorem targets

1. **Peel-relative splice.** The former complementarity
   `r_D<=B(c*)` is refuted by the exactly-`f=4` family above. Replace it with a
   scale-normal state containing `r_D/(v|G'_D|)`, the divisor-support profile,
   and endpoint owners; then compose THM-761's high-support sheet regime with
   the capped-envelope band or a strictly smaller descendant.
2. **Seven-exception sheet residue.** Classify the `1/7` cyclic tilings on
   `Z/cZ`, and formalize the terminating `c -> c/g` bookkeeping for excessive
   gcd multiplicity.  These are method walls, not known LRC counterexamples.
3. **Finite jet outside sheet structure.** For affine rays not dispatched by
   THM-761, test whether leading shape, first owned offset, residue address,
   and an adaptive named certificate eventually decide the ray.  THM-764
   rules out a uniform `q<=25` terminal window.
4. **Four-circuit localization.** Outside coherent pack/cluster faces, show
   that blocker completeness forces a bounded-height marked relation circuit
   meeting all four far coordinates.
5. **Far-dilation event cocycle.** For the exact action `g |-> k*g`, classify
   which owner endpoint blocks split, merge, or reorder, and determine the
   smallest sidecar that transports fixed-threshold truth and metric data.
   The `6/61` versus `12/121` canary rules out an `M`-invariant orbit quotient.
6. **Cone descent.** On every residue chamber of every `X_C`, produce a
   threshold point or a well-founded move to a smaller normalized state.
7. **Cubical/blocker comparison in Lean.** Connect zero-state reachability in
   the colored endpoint loop to an empty endogenous blocker edge and state
   exactly what data transports between the two fixed-row presentations.
8. **Information-minimal assembly.** Prove which fields of the truth, metric,
   functor, and proof carriers can be reconstructed or discharged on each
   current theorem lane.

## Tournament Analysis and assumption challenge

The affine audit uses representations as vertices.  Its observable is the
number of proof-critical row pairs separated; its switches are
predicate-first and compression-first; its tie path is

```text
full affine slope packet
> exact peak certificate
> owned offset plus c mod 14
> dual relation signature
> owner-resolved offset
> offset multiset
> shape only.
```

Both gauges have score histogram `0,1,...,6`, zero directed 3-cycles,
singleton SCCs, and one Hamiltonian path; 11 edges flip between gauges.

Alternate vertex sets explicitly considered were runners, gaps, fixed circle
sections, section boundaries, endpoint events, residues, cover arcs, Fourier
modes, marked relation circuits, blocker edges, and proof obligations.  Raw
runners and unweighted endpoints are rejected as sufficient vertices.  The
selected carrier is the representation/quotient itself, because the question
being audited is controlled forgetting.

The status is exact reparameterization plus exact finite evidence.  No finite
or well-founded classification of all integer-slope fibers is proved here.
