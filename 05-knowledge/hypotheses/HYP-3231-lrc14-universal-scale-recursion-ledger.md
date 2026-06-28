---
id: HYP-3231
title: Universal scale invariance turns the LRC14 route into a scale-normal recursion ledger
status: SYNTHESIS / route-sharpening ledger; not an LRC14 proof
source: codex-2026-06-28
tangent: T1330
technique: LTI-330
tournament_technique: LTT-230
reflection: 07-reflections/lrc14-universal-scale-invariance-recursion-ledger-codex-20260628.md
related:
  - HYP-3245
  - HYP-3237
  - HYP-3236
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3230
  - HYP-3229
  - HYP-3228
  - HYP-3227
  - HYP-3215
  - HYP-3214
  - HYP-3213
  - HYP-3212
  - HYP-3205
  - HYP-3162
  - HYP-2963
  - HYP-1969
  - THM-573
  - THM-532
  - THM-407
  - OPEN-Q-108
---

# HYP-3231: Universal Scale Invariance As The LRC14 Recursion Ledger

## Claim

Assume the scale invariance used throughout the repo holds universally:

```text
M(cS) = M(S) for every nonzero integer c,
Lonely(c*v, t/c) <=> Lonely(v, t).
```

Then the current LRC14 route should be read as a projective, scale-normal
recursion rather than as a search for a final scalar.  Every successful
sharpening in the repo has had the same form:

```text
normalize scale
-> split the first residue / route / topology coordinate not killed by scale
-> attach the sidecar that makes the quotient legal
-> discharge an easy branch
-> recurse on the surviving primitive packet.
```

The hidden pattern is that the route keeps returning to the same object:

```text
primitive projective shape
+ finite address
+ destroyed-coordinate sidecar
+ named discharge or residual debt.
```

This is why nonprimitive ties such as doubled AP are not counter-signals.
They are gauge copies.  The real signal is always the first coordinate that
does *not* descend through scale: nonunit residue stratum, endpoint owner,
relation height, branch route, order-3 overlap, or certificate image.

Push-time integration: mac-mini S75b concurrently claimed HYP-3230 for the
three-gap/Stern-Brocot recursion in the cap kernel.  This HYP-3231 ledger is
the route-level companion: HYP-3230 makes the order-2 cap-kernel recursion
explicit, while HYP-3231 asks every quotient in the larger LRC14 route to be
scale-normal, sidecar-restored, or named as residual debt.

## How The Route Sharpened

### 1. Proof currencies, not one scalar

HYP-1969 split the proof landscape into finite-sieve, kernel, zonotope, and
endpoint-pressure currencies.  That was the first durable reframing: LRC14
progress is not a scalar contest, but a comparison of proof carriers.

The retained lesson:

```text
route value = predicate retained + cost of forgotten coordinates.
```

### 2. Scale becomes a theorem-level quotient

THM-407 made the scale quotient explicit on the additive residual:
`M(cS)=M(S)` and time reversal generate a shell action.  For n=14,
`2n-1=27=3^3`, so the pair-sum shell ledger collapses from `13` shells to
three scale-stable strata:

```text
gcd(a,27) = 1, 3, 9.
```

The recursive pattern appears here first: the unit stratum descends by the
scale action, while the nonunit strata survive as named residual coordinates.

### 3. Scale-invariant sector mass isolates relation height

THM-532 proves the seven-sector cover mass is scale-invariant and decomposes
it as

```text
meas(S7(E)) = M7(k) + corr(E),
```

where the correction is a relation-lattice tail.  High relation height is the
scale-free generic branch; the low-height AP-rich branch is finite modulo
scale and translation.  The same recursion appears again:

```text
generic high-height branch -> certificate margin
low-height scale-normal branch -> finite residual.
```

### 4. The covering proof tree moves scale to the boundary

The THM-573 level-7 lift sieve and the covering-bound correction sharpen the
global LRC14 perimeter:

```text
WLOG primitive
14-free                         -> t = 1/14
covering, >=7 multiples of 7     -> level-7 lift sieve
covering, <=6 multiples of 7     -> residual
```

The important correction was scale-theoretic: covering rows are not bounded
away from `1/14`, because dilated AP/GW rows sit exactly on the same floor.
Thus the proof is not allowed to use a margin that scale invariance forbids.

### 5. Packet sheaves replace unsafe quotients

HYP-2963 and the HYP-2990..HYP-3046 stack sharpened the route from scalar
features into labelled packet sidecars: route labels, exact `M`, endpoint
owners, topology, primitive decks, Haar zeta, residual capacitors, and
cycle-class images.  The operative law became:

```text
a quotient is legal only if the forgotten coordinate is
fiber-constant, reconstructed, dual-annihilated, descended,
boundary-stopped, or routed to named residual debt.
```

This is the same scale-normal recursion with better bookkeeping.

### 6. Median and branch closures make the proof interface finite

HYP-3070/HYP-3074/HYP-3077/HYP-3082/HYP-3083 sharpened the interface from
route labels to legal proof-state medians and protected branch graphs.  Raw
route labels had empty centers; sidecar-completed route states had legal
centers.  The branch-kernel audit found that scalar contraction creates
naked bridges, while protected branch graphs do not.

So the recursion sharpened again:

```text
raw route quotient
-> sidecar-completed route state
-> median center or named bridge debt.
```

### 7. The k=8 bounded node becomes a certificate-vector problem

HYP-3160 through HYP-3229 turned the bounded core into a compatibility problem
among several scale-stable certificate faces:

```text
Fejer / de-Moivre sector kernel
S75 comb-overlap Gram kernel
Johnson J(14,2) cap kernel
Toeplitz / covariance / Green trap discharge
Joukowski / Hermite-Biehler / cyclotomic cubic
Gamma0(7) coefficient sidecar
```

HYP-3229 is the latest sharpened form: Gamma0(7) is a coefficient engine for
finite LP/Toeplitz rows, while Stark, Beraha, Mahler, and subshift data remain
sidecars until they preserve a named LRC predicate.

Post-rebase companion HYP-3245 adds a lag-space test for that preservation:
after scale normalization, an autocorrelation sidecar must say whether
short-lag deficit and outward-lag surplus descend through the same primitive
packet or whether they create a named scale/fiber debt.

Incoming HYP-3232 is the separate interlocking-recursion companion: it locates
where scale covariance breaks at the apex fold, so HYP-3231's scale-normal
ledger should mark whether a signal survives below-apex three-gap scaling or
falls into the antipode-half deviation sector.

The newest hidden recursion is now order-theoretic:

```text
order-2 PSD / Gram data closes pair overlap,
remaining debt moves to order-3 triple overlap / non-associative cocycle.
```

## New Recursive Pattern To Name

The route has been describing a scale-normal renormalization tower without
quite naming it.  Each level has the same grammar:

| Level | Normalized object | First surviving coordinate | Discharge | Residual |
|---|---|---|---|---|
| Scale | primitive projective shape | nonunit shell / dilation fixed row | quotient by `c` and sign | gcd strata |
| Apex | covering profile | multiples of 7 | `t=1/14` or level-7 lift | `<=6` multiples of 7 |
| Relation | sector-cover orbit | short relation height | high-height certificate | AP-rich finite set |
| Packet | labelled row | route/topology/owner | q-witness, AP/GW, petal, K33, covering | F7/THM-572 debt |
| Sidecar | proof-state median | first missing sidecar | legal center | bridge/cocycle debt |
| Analytic | Fejer/Gram/Toeplitz vector | triple-overlap constant | order-2 PSD certificate | order-3 debt |

This is not merely a history.  It suggests a proof schema:

```text
Every surviving packet must either
  descend one level in the scale-normal tower,
  acquire the sidecar that makes the descent legal,
  or expose the first nonzero residual coordinate.
```

## Proof Targets

1. **Scale-normal packet theorem.**  Add `primitive_scale_gcd`,
   `scale_orbit_representative`, `dilation_fixed_boundary`,
   `nonunit_residue_stratum`, `scale_destroyed_payload`, and
   `renormalization_depth` to the active packet ledger.  Prove every
   post-THM-573 residual has a finite scale-normal address before it enters
   Fejer/Gram/Toeplitz certificates.

2. **Scale-fiber exactness theorem.**  For every quotient used on HYP-2963
   packets, define the emitted cocycle `omega_Q`.  Prove it is exact,
   dual-annihilated, descended, boundary-stopped, or a named residual.

3. **Order-3 overlap theorem.**  HYP-3229 and S75 reduce the new certificate
   debt to triple-overlap constants after order-2 Gram positivity.  The next
   finite LP should make these constants explicit and test them against
   HYP-3227 Green trap weights.

4. **Wide scale-separated theorem.**  Import the HYP-3215/Rosenfeld
   exponential-sum route into the same scale-normal fields: high relation
   height and far-cluster scale separation should be one theorem, not two
   incompatible estimates.

5. **Base flag isolation.**  HYP-3215's induction-base warning is not a
   bounded-core certificate issue.  Keep it in the global route ledger so the
   Fejer/Gram/Gamma0 certificate does not accidentally certify an unverified
   base.

## Tournament Analysis

Vertices are proof-route states and scale fibers, not runners or raw arcs.
Alternate vertex sets considered: runners, gaps, fixed circle sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
Haar rectangles, matroid circuits, endpoint owners, triple overlaps, and proof
obligations.

Chosen vertices:

```text
lean_lonely_scale
primitive_projective_shape
twisted_shell_gcd_strata
level7_lift_sieve
relation_height_split
labelled_packet_sheaf
route_state_medianization
fejer_gram_magic_certificate
green_toeplitz_trap_discharge
gamma0_7_coefficient_engine
raw_scalar_or_sequence_shadow
```

Pairwise observable:

```text
scale-invariant LRC predicate retained,
then first surviving destroyed coordinate exposed,
then certificate checkability,
then residual debt named.
```

Switch/gauge: orient `A -> B` when `A` preserves more of the scale-normal
proof payload with no increase in unnamed residual debt.  Tie path:

```text
primitive_projective_shape
> level7_lift_sieve
> relation_height_split
> labelled_packet_sheaf
> route_state_medianization
> fejer_gram_magic_certificate
> green_toeplitz_trap_discharge
> gamma0_7_coefficient_engine
> raw_scalar_or_sequence_shadow
```

Challenged assumption: the useful tournament vertices are not runners by
default.  For this session the preserved LRC predicate is `M(S)>=1/14` modulo
nonzero scale.  The quotient destroys absolute speed size, reset-calendar
complexity, finite product-bound size, and some endpoint event order; those
must be restored by packet sidecars before any scalar or analogy can be used.

## Status

This is a route-sharpening synthesis, not a proof.  It makes the current route
more exact by naming the recursion already present in THM-407, THM-532,
THM-573, HYP-2963, HYP-3083, and HYP-3229.  The open mathematical work is
unchanged but cleaner:

```text
close the bounded certificate vector,
discharge order-3 overlap debt,
prove the wide scale-separated branch,
and verify or replace the global induction base.
```

-> HYP-3245, HYP-3243, HYP-3242, HYP-3241, HYP-3240, HYP-3239, HYP-3238, HYP-3237, HYP-3236, HYP-3235, HYP-3234, HYP-3233, HYP-3232, HYP-3230, HYP-3229, HYP-3228, HYP-3227, HYP-3219, HYP-3215, HYP-3214,
HYP-3213, HYP-3212, HYP-3205, HYP-3162, HYP-2963, HYP-1969, THM-573,
THM-532, THM-407, T1330, LTI-330, LTT-230, OPEN-Q-108.
