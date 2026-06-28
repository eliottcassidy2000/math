---
id: HYP-3244
title: Tiling and half-tiling interlocking recursions as a controlled-forgetting span
status: EVIDENCE / exact small tournament recursion audit plus LRC proof-interface proposal; not a proof
source: codex-2026-06-28
tangent: T1344
technique: LTI-344
tournament_technique: LTT-244
script: 04-computation/tournament_tiling_half_tiling_interlock_codex_20260628.py
result: 05-knowledge/results/tournament_tiling_half_tiling_interlock_codex_20260628.out
reflection: 07-reflections/tiling-half-tiling-interlocking-recursions-codex-20260628.md
related:
  - HYP-3243
  - HYP-3242
  - HYP-3241
  - HYP-3240
  - HYP-3239
  - HYP-3238
  - HYP-3237
  - HYP-3236
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3231
  - HYP-3230
  - HYP-3229
  - HYP-3228
  - HYP-3227
  - HYP-3220
  - HYP-3219
  - HYP-3218
  - HYP-3216
  - HYP-3214
  - HYP-3053
  - HYP-3052
  - HYP-3051
  - HYP-3149
  - HYP-3143
  - HYP-3199
  - OPEN-Q-108
---

# HYP-3244: Tiling and Half-Tiling Interlocking Recursions

## Claim

The fixed-path tiling model and the half-tiling/orbit model should not be
identified.  They are two interlocking recursions connected by a
controlled-forgetting span.

```text
C_n  --tiling add/free diagonal-->  C_{n+1}
 |                                  |
 v  quotient with fiber sidecar     v
U_n --Aut(parent)-word orbit-----> rooted U_{n+1} -> U_{n+1}
```

Here `C_n` is the fixed-Hamiltonian-path tiling cube and `U_n` is the A000568
orbit set of unlabelled tournaments.  The tiling recursion keeps explicit flip
witnesses.  The half-tiling recursion compresses by parent automorphisms,
incident-word orbits, and `GF(2)` coboundary laws.

The square commutes only after adding sidecars:

```text
path-presentation fiber H(T)/|Aut(T)|
parent automorphism incident-word orbit
rectangle/hourglass coboundary residue
tail/tip deletion signature
```

This is the tournament-side analog of the incoming LRC recursion story.  HYP-3216
names the LRC-side cyclotomic ladder and 2-adic fold, HYP-3230 names the
three-gap/Stern-Brocot recursion inside the cap kernel, HYP-3231 names the
scale-normal route ledger, and HYP-3232 names the modulus-covariance recursion
that breaks at the apex half.  HYP-3233 names the cyclotomic-factor layer,
HYP-3234 names the signed-address chart-change sheaf, HYP-3235 names the
totally-real cyclotomic field packet, HYP-3236 names the Green conductance
face, HYP-3237 names the Vitali-wall/core-construction split, HYP-3238 names
the crossed even-positive/odd-negative compression bridge, HYP-3239 names the
`D_7`/Borsuk-Ulam sign-isotypic refinement plus the bimodal-discrepancy
diagonal, HYP-3220 names the even-odd/positive-negative parity wall, and
HYP-3219 names the Brouwer/sign factorization sidecar.  Incoming HYP-3240 names
the hard-core dilation witness guardrail, where covering-tight dilations move
the witness to a finer `Phi_{14d}` grid.  Incoming HYP-3241 names the
equioscillation saddle-index count `(p-1)/2` and the shared AP/Goddyn-Wong
`Phi_14` witness core, giving the topological sidecar that decides when the
quotient should be read as Borsuk-Ulam forcing rather than a symmetric
Brouwer/SOS fold.  Incoming HYP-3242 names the Euler/Cech view: the cap is the
measured Euler characteristic of the danger-cover nerve and the lonely point is
the hole.  Incoming HYP-3243 turns that geometry into a typed graph proof-route
atlas.  HYP-3244 supplies the
missing tournament chart: lift to the witness-rich tiling cover, compress
through the half-tiling quotient only when sidecars certify that the LRC
predicate descends.

## Exact Readout

The `n=4` fixed-path tiling chart uses the base path `0->1->2->3` and free
flips

```text
a=(0,2), b=(1,3), c=(0,3).
```

The exact fixed-path fibers are

```text
T: ['E']
+: ['a']
-: ['b']
S: ['ab', 'c', 'ac', 'bc', 'abc']
```

So the prompt table is exactly XOR on the named generators followed by the
isomorphism-class quotient:

```text
* | E  a  b  c
E | T  +  -  S
a | +  T  S  S
b | -  S  T  S
c | S  S  S  T
```

The half-tiling chart can be realized with fixed arcs

```text
2->0, 3->0, 1->2, 3->1
```

whose partial outdegree vector is `[0,1,1,2]`, and free bits

```text
x=(2,3), y=(0,1).
```

It gives

```text
* | E  x  y
E | T  +  -
x | +  T  S
y | -  S  T
```

Thus the half chart is a four-state section, while the fixed-path tiling chart
has an `S`-fiber of size `5`.  The missing payload is the smallest visible
canary/filler sidecar.

The fixed-path cover law through `n=6` is:

```text
n=3: fixed_cube=2    U(n)=2   fiber_hist={1: 2}
n=4: fixed_cube=8    U(n)=4   fiber_hist={1: 3, 5: 1}
n=5: fixed_cube=64   U(n)=12  fiber_hist={1: 4, 3: 1, 5: 3, 9: 2, 11: 1, 13: 1}
n=6: fixed_cube=1024 U(n)=56  H(T)/|Aut(T)| check=True
```

The half-tiling incident-word orbit recursion gives:

```text
1->2: word_orbits=2    rooted=2    U(2)=1
2->3: word_orbits=4    rooted=4    U(3)=2
3->4: word_orbits=12   rooted=12   U(4)=4
4->5: word_orbits=48   rooted=48   U(5)=12
5->6: word_orbits=296  rooted=296  U(6)=56
```

At the first large bookkeeping surface:

```text
fixed-path cover at n=6 = 1024
rooted half-tiling/orbit count = 296
unrooted classes = 56
```

These are not contradictory counts; they are the two sides of the span.

The coboundary compression side agrees with HYP-3053:

```text
K_{k,k+1}: lines = k(k+1), rank = 2k, redundancy = k(k-1)
```

The rectangle law is

```text
L(a,b)+L(a,b')+L(a',b)+L(a',b') = 0.
```

A nonzero rectangle or hourglass residue means the quotient has forgotten real
payload: endpoint owner, route status, active support, barcode bar, or
proof-obligation state.

## LRC Translation

Use the two recursions alternately:

1. Lift to the fixed-path tiling cover to get an explicit witness, exchange
   move, tail/tip packet, or finite obstruction.
2. Compress through the half-tiling/orbit recursion only after the fiber,
   automorphism, coboundary, and deletion sidecars certify descent.
3. If the descent fails, keep the failed sidecar as a named proof obligation
   rather than replacing it by a scalar count.

This reframes the LRC14 finite frontier.  HYP-3227/HYP-3229/HYP-3230 are
building Fejer/Gram/three-gap certificate rows, HYP-3231 records the
scale-normal recursion discipline, and HYP-3232 isolates the modulus fold
where the antipode half begins to deviate.  HYP-3244 supplies the corresponding
tournament interface: explicit tiling witnesses should discharge finite traps,
while half-tiling compression should certify that those discharges survive
quotienting.

Incoming HYP-3238 turns this from a naming discipline into a finite audit on
the same anchored k=8 bank.  It finds `18` primitive non-AP rows with zero
negative covariance leakage, so absence of a negative edge is a false descent
certificate.  It also finds `2879` primitive rows with positive `q3` debt but
`0` exchange-margin violations, and the `11` non-AP HYP-3202 traps split into
`8` negative-leakage-plus-odd-`q3` debts and `3` odd-`q3` debts with no
negative covariance leakage.  For HYP-3244, that means a half-tiling quotient
cannot certify an even/positive scalar unless the odd/negative sidecar is
explicitly retained or priced by the endpoint-bimodality exchange.

## Proof Target

The next theorem-facing form is a descent criterion:

```text
tiling witness descends through the half-tiling quotient
  if fiber sidecar is constant or accounted for,
  parent-aut word orbit is named,
  rectangle/hourglass residues vanish or are dual-annihilated,
  and tail/tip deletion packets agree with the target LRC predicate.
```

If proved for the finite HYP-3202/HYP-3227 trap boundary, this would turn the
current trap-discharge computations into a small controlled-quotient theorem:
construct the witness in the tiling cover, then certify its half-tiling descent.

## Tournament Analysis

The tournament vertices are proof carriers, not runners or arcs:

```text
tail_tip_deletion_sidecar
tiling_witness_lift
half_tiling_orbit_compression
parent_aut_incident_word_orbit
coboundary_rectangle_hourglass_law
fiber_debt_H_over_Aut
n4_canary_filler_section
lrc_fejer_gram_recursion_bridge
raw_fixed_path_cube_count
raw_A000568_class_count
```

Pairwise observable:

```text
LRC predicate retained through lift/compress recursion minus quotient risk.
```

Switch:

```text
orient toward the carrier that keeps a witness and certifies descent.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3_cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

## Assumption Challenge

Vertices considered:

```text
runners, arcs, gaps, fixed circle sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits,
flip bits, deletion roots, rectangle defects, and proof obligations.
```

Chosen vertices are proof carriers for the two recursions.

Preserved predicate:

```text
explicit witnesses, orbit legality, fiber multiplicity, coboundary
consistency, and LRC sidecar descent.
```

Destroyed if scalarized:

```text
runner labels, Hamiltonian-path presentation, canary/filler S-fiber mass,
endpoint ownership, and active certificate coordinates.
```

Challenged assumption:

```text
The half-tiling chart is not a replacement for the tiling cube.  It is a
compressed section whose legality depends on sidecars.
```

## Next Pulls

1. Attach this span to the HYP-3227 trap-discharge graph: each trap edge should
   name its tiling lift, half-tiling descent certificate, and failed sidecar if
   descent is not legal.
2. Combine HYP-3230's three-gap cap-kernel recursion, HYP-3232's
   modulus-covariance apex break, HYP-3234's signed chart-change debts,
   HYP-3235's totally-real Fejer-square packet, HYP-3236's Green-resistance
   slack, HYP-3237's Vitali core/bulk wall, HYP-3238's exact false-terminal
   and odd-`q3` payload audit, HYP-3239's `D_7`
   sign-isotypic/Borsuk-Ulam audit, HYP-3220's positive/negative parity
   sidecar, and HYP-3219's sign/degree sidecar
   with the `K_{k,k+1}` rectangle laws: test whether continued-fraction kinks,
   speed-`>n/2` deviations, local-slot cancellations, cyclotomic conductor
   defects, core/bulk witness loss, parity-sign flips, or topological sign
   defects correspond to nonzero rectangle/hourglass residues after
   quotienting.
3. Extend the `n=4` canary/filler chart to `n=5` by naming which fixed-path
   fibers survive under tail deletion, tip deletion, and converse/reflection.
4. Turn the descent criterion into a small finite schema for the HYP-3202
   exchange traps, HYP-3236 Green-conductance traps, HYP-3237 Vitali-core
   witnesses, HYP-3238 crossed-duality sidecars (`18` zero-negative false
   terminals and the trap `8+3` split), HYP-3239 `D_7`/bimodality sidecars,
   HYP-3240 dilation-witness sidecars, HYP-3241 saddle-index/GW-core
   sidecars, HYP-3242 Euler/Cech-hole sidecars, HYP-3243 topology-graph route
   sidecars, HYP-3220 parity sidecars, HYP-3219 sign sidecars, and the
   HYP-3218 Fejer/equidistribution proof-push before attempting a global LRC14
   proof.

-> HYP-3244, HYP-3243, HYP-3242, HYP-3241, HYP-3240, HYP-3239, HYP-3238, HYP-3237, HYP-3236, HYP-3235, HYP-3234,
HYP-3233, HYP-3232, HYP-3231, HYP-3230, HYP-3229, HYP-3228, HYP-3227,
HYP-3220, HYP-3219, HYP-3218, HYP-3216, HYP-3214, HYP-3053, HYP-3052,
HYP-3051, HYP-3149, HYP-3143, HYP-3199, T1344, LTI-344, LTT-244,
OPEN-Q-108.
