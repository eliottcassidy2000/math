---
id: HYP-6815
title: LRC14 four-far cone and affine-slope threshold suspension
status: EXACT REPRESENTATION LEMMAS + EXACT AUDIT; finite descent open
source: codex-2026-07-14-S2
script: 04-computation/lrc14_affine_slope_suspension_codex_S2.py
result: 05-knowledge/results/lrc14_affine_slope_suspension_codex_S2.out
related:
  - HYP-6780
  - HYP-6785
  - HYP-6755
  - HYP-3106
  - HYP-3072
  - HYP-3025
  - THM-668
  - THM-738
  - THM-741
  - THM-742
  - THM-755
---

# HYP-6815: LRC14 Four-Far Cone And Affine-Slope Threshold Suspension

## Verdict

There are two precise meanings of the requested four-dimensional object.

1. The first unresolved far-count chart is literally four-dimensional.  After
   THM-738 closes `f<=3`, fix a nine-speed core `P subset {1,...,14}` and let
   the four ordered far speeds vary.  They form an arithmetic lattice cone in
   four variables.
2. An affine ray in that cone has a mixed four-coordinate suspension
   `(u,t,c,lambda)`: two torus phases, an integer slope, and a clearance level.

The global moduli space of 13 speeds modulo dilation is not four-dimensional.
Calling all of LRC14 a four-manifold would discard the higher far-count strata.
A true global four-dimensional reduction still requires a peel or descent
theorem.

The exact object on the first chart is therefore not a naked point of `Z^4`.
It is an arithmetic four-cone decorated by an owner-colored threshold loop,
rational gap metric, embedded marked relation lattice, and the action of the
next permitted proof operation.

## 1. The outer four-far cone

Fix `P subset {1,...,14}`, `|P|=9`.  In gap coordinates write

```text
w1 = a
w2 = a+d1
w3 = a+d1+d2
w4 = a+d1+d2+d3,
```

with `a>=15` and `d1,d2,d3>=1`.  Let `X_P` be the points
`g=(a,d1,d2,d3)` for which

`S(P,g)=P union {w1,w2,w3,w4}`

is primitive and covering.  The union of these charts over the
`binom(14,9)=2002` possible cores is exactly the `f=4` stratum behind the
THM-741 computation.

Let `Q=lcm(2,...,14)=360360`.  Covering is a Boolean combination of the
congruences `q|wj`, `2<=q<=14`.  Primitivity depends only on the primes
dividing `gcd(P)` and the `wj`.  Hence:

> **Semilinear cone lemma.** `X_P` is a finite union of residue classes
> modulo `Q`, intersected with the positive ordered gap cone.

This is an exact finite arithmetic base, but not a finite LRC classification.
The archimedean heights still move endpoint walls and pair-sum moduli inside a
fixed residue class.

For `g in X_P`, define the cocharacter

`phi_g(t)=(s*t mod 1)_{s in S(P,g)}`

and the closed safe cube `K=[1/14,13/14]^13`.  The universal incidence

`Z_P={(g,t): phi_g(t) in K}`

projects to `X_P`.  The exact `f=4` claim is simply surjectivity of
`Z_P -> X_P` for every core `P`.

## 2. The inner affine-slope suspension

Let `P=(p_i)` and `R=(r_i)` be integer 13-vectors and let

`V(c)=(c p_i+r_i)`,

on the positive, distinct, primitive integer scales `c` under study.  On the
two-torus set

`Phi_{P,R}(u,t)=min_i ||p_i u+r_i t||`.

For the slope-`c` closed geodesic `L_c={(ct,t):t in R/Z}`, direct substitution
gives the exact identity

`M(V(c)) = max_t Phi_{P,R}(ct,t)`.

Thus the mixed four-coordinate incidence object

`X_{P,R}={(u,t,c,lambda): u=ct, Phi_{P,R}(u,t)>=lambda}`

in `(R/Z)^2 x N x [0,1/2]` has the LRC14 assertion as nonemptiness of its
integer-slope fiber at `lambda=1/14`.  It is a stratified continuous/discrete
object, not a smooth four-manifold.  Runner owners label its strip strata.

When `R=0`, `Phi` is independent of `t`; the suspension is cylindrical and
recovers HYP-6780's dilation law.  Nonzero `R` is transverse shear/holonomy,
not a small perturbation that may be discarded.

This generalizes THM-742's exact slope-geodesic formulation for
`B union (W+J)`.  THM-742 is one polygonal chart; HYP-6815 retains arbitrary
owned offsets and the clearance filtration.

## 3. The exact dual resonance picture

For any coefficient vector `z in Z^13`, define

```text
m(z)=z dot P,
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
exactly the Fourier modes on that line.  For discontinuous safe indicators,
the same statement must be taken through a justified Abel/Fejer limit; no
absolute-convergence claim is made here.

This identifies the repository's relation-lattice, Fourier, and scale-residue
languages.  Shape keeps only `m`; residue keeps only `n`; the LRC slice uses
their pairing with `c`.  The exact audit finds that the AP shape has 36
support-three relations `e_i+e_j-e_k`, while a single owned offset kills from
6 to 11 of them, depending on its owner.

The full embedded marked lattice

`ker(Z^13 -> Z, z |-> z dot S)`

determines a primitive positive speed vector up to sign: it is the labelled
rank-12 kernel of a primitive integer functional, whose normal is unique up to
sign.  An abstract isomorphism class, rank, support enumerator, or successive
minimum does not retain that labelled normal or its coordinate-facet pairing
with the safe cube.

## 4. What information is actually sufficient?

There is no task-independent minimal statistic.  The correct carrier depends
on the next question.

### Fixed-threshold truth

Let

`beta_S(t)_i = 1_{||s_i t|| < 1/14}`.

Starting at `t=0`, record the cyclic sequence of endpoint blocks, with each
block retaining entering owners, exiting owners, and exact coincidences.  This
owner-colored event word reconstructs the piecewise-constant cubical loop
`beta_S`.  It decides whether the loop or a boundary tie state hits `0^13`.
Exact gap lengths are unnecessary for this Boolean predicate.

An aggregate signed current is insufficient: an owner entering exactly when
another exits can cancel in the scalar current even when the tie state carries
the only equality witness.

### Margin, measure, and autocorrelation

To recover exact `M`, safe measure, component count, or THM-731/732, add the
rational endpoint phases and gap lengths, or equivalently retain the full
threshold filtration `lambda |-> G_S(lambda)`.  The signed owner endpoint
divisor reconstructs the Bernoulli `B2` discrepancy.  An exact peak
certificate `(t*,M)` recovers `M` alone, not the safe measure or topology.

The exact 552-row endpoint audit supplies sharp canaries:

```text
AP                         M=1/14
{1,...,12,26}              M=2/27
```

They have the same endpoint tournament, divisor mask, and cap sign.  Thus
order, covering, and the THM-755 bit do not preserve the metric predicate.

### Covering and capped-envelope status

Covering needs the divisor mask.  The THM-755 decision needs the projective
ratio `v*|G|/r` (or its exact sign against the threshold), not endpoint
orientation.  A lift `21 -> 22` in the audit changes the cap decision with zero
endpoint-tournament edge flips.

### Transport across rows

Deletion and peel require owner labels.  Dilation and affine transport require
the scale/residue action on the fiber, not merely a normalized point.  Fourier
or interaction arguments require the embedded relation coefficients and
heights.  Proof assembly additionally needs the next observer, available
certificate, legal discharge, or named residual debt.

The resulting hierarchy is:

```text
truth carrier    = initial state + signed owner event blocks
metric carrier   = truth carrier + rational phases/gap filtration
functor carrier  = metric carrier + scale/residue action + marked relations
proof carrier    = functor carrier + next observer/certificate/discharge
```

## 5. Primal and dual finite presentations

The owner-colored cubical loop is the primal presentation.  HYP-6785's
endogenous pair-sum blocker complex is the dual presentation: pair-sum peak
candidates carry blocker sets, and LRC holds exactly when some candidate has
an empty blocker edge.  The dual is finite at fixed `S`, but a scale-normal
pressure or descent theorem over `X_P` remains open.

This suggests a bicomplex rather than a flat ledger:

```text
d_geom  = cross an endpoint/contact wall
d_arith = change a residue, valuation, or scale depth
```

If quotienting makes these operations fail to commute, the commutator is not
noise.  It is the missing owner, height, relation, or certificate sidecar.

## 6. Exact executable evidence

`lrc14_affine_slope_suspension_codex_S2.py` verifies:

1. the exact slice identity for pure AP dilations;
2. `|G_{cP}|=|G_P|` and `r_{cP}=c*r_P` for `P={1,...,12}`;
3. the transverse ray `c{1,...,13}+e_13` at covering scales, with exact
   `M=1/13`;
4. the projected support-three relation loss under owned shear;
5. a metric collision at equal shape, owned offset, and `c mod 14`:
   `c=2` and `c=16` have the same `M=2/15` but different exact safe measures
   `5/84` and `115/1904`, and component counts `4` and `30`;
6. a representation tournament under predicate-first and compression-first
   gauges.  Both are transitive, with 11 edge flips between gauges.

The larger endpoint/tournament audit in
`lrc14_endpoint_tournament_sidecar_audit_structure_S1.py` uses 552 exact rows,
checks 552/552 Bernoulli reconstructions, and supplies the collision counts
quoted above.

## 7. Transfers from seemingly unrelated threads

- Tournament switching and tiling: preserve the group action or gauge
  section, not merely a bijection of underlying sets.
- Coding and matroids: a weight enumerator preserves support counts but loses
  support connectivity, coefficient height, and realizability.  LRC likewise
  needs the embedded marked relation code.
- Ising and Lee-Yang: a first moment or `p0` loses compatible many-body packing;
  the whole threshold/root surface is a classifier, not by itself a discharge.
- Resolvent folds: the symmetric inner page is useful only with center,
  antisymmetric leakage, and boundary sidecars.  Shape is the analogous inner
  page here; residue shear is the leakage.
- CRT and p-adic trees: residue skeleton and valuation/height flex are separate
  coordinates.  Fixed residue does not freeze endpoint geometry.
- Cut/cycle topology: locally consistent views need not glue globally.  The
  owner event loop and blocker complex must agree as global presentations.
- Observer-cut and proof circuits: a quotient is legal only relative to a
  named next operation, with every changed predicate reconstructed,
  annihilated, descended, boundary-stopped, or routed to named debt.

## 8. Open theorem targets

1. **Finite jet at infinity.** Prove that every affine ray `cP+R` is
   eventually dispatched by `P`, its first nonzero owned offset jet, finitely
   many residue classes, and a named certificate.
2. **Four-circuit localization.** Outside coherent pack/cluster faces, show
   that blocker completeness forces a bounded-height marked relation circuit
   meeting all four far coordinates.
3. **Cone descent.** On every residue chamber of every `X_P`, produce a
   threshold point or a well-founded move to a smaller normalized state.
4. **Cubical/blocker duality in Lean.** Connect zero-state reachability in the
   colored endpoint loop to an empty endogenous blocker edge.
5. **Information-minimal assembly.** Prove which fields of the truth, metric,
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
