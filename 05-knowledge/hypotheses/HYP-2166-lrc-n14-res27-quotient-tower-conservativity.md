---
id: HYP-2166
status: SUPPORTED by S610 quotient-tower synthesis; lift/CRT conservativity open
source: user-2026-06-03; codex-2026-06-03-S610
related:
  - HYP-2175
  - HYP-2165
  - HYP-2164
  - HYP-2163
  - HYP-2162
  - HYP-2161
  - HYP-2160
  - HYP-2142
  - HYP-2107
  - HYP-2105
  - THM-398
  - THM-407
  - THM-406
  - THM-401
---

# HYP-2166: the n=14 `Res_27` improvements form a quotient tower with one lift seam

## Claim

The recent n=14 improvements are not separate tricks.  They form one
coimage/opfibration tower:

```text
clock exit
  -> <2,-1> shell-orbit quotient
  -> least-positive pair-sum/pinch lower-bound quotient
  -> owner-fibre certificate reattachment.
```

The left-moving maps forget data and prove lower bounds cheaply.  The
right-moving owner fibre reattaches exactly the labels needed to discharge
floor rows by cheap pairs or positive measure.

The remaining n=14 proof should therefore be a conservativity theorem for the
lift/CRT seam:

```text
arbitrary integer lifts over the Res_27 coimage
  land in strict pinch rows, AP/V*, normalized 2*AP,
  or the two positive owner-control rows.
```

## S610 Evidence

`04-computation/lrc_n14_res27_quotient_tower_s610.py` composes HYP-2162,
HYP-2163, HYP-2164, and HYP-2165.

It recomputes:

```text
shell fold:                 13 raw shells -> 3 gcd strata {1,3,9}
pinch D/U/N survivors:      27733
pinch proof-obligation types: 148
pinch strict rows:          27730
pinch floor rows:               3
pinch below rows:               0
owner slack rows through 81: 17550
owner full covers:           9506
owner open residuals:           0
```

The floor atoms are:

```text
AP                 pinch floor, owner floor, cheap pair (1,13)
V*                 pinch floor, owner floor, cheap pair (1,13)
2*AP               nonprimitive pinch floor, normalizes to AP
```

The owner-only block-all controls are positive-measure, and the pinch gauge
already sees them strictly above the floor:

```text
slack=(3,9,12,42):   pinch=2/23 at t=4/23
slack=(3,12,27,42):  pinch=3/32 at t=7/32
```

Thus there is no reason to enumerate the `64` fixed classes, `27733` pinch
survivors, or `9506` owner covers again.  The proof frontier is the functorial
step from arbitrary lifts to this finite atom list.

## Structural Reading

Each recent improvement preserves a different predicate:

```text
HYP-2163 clock exit:
  preserves "some runner is lonely at t=1/n" and removes the no-multiple bulk.

HYP-2162 shell fold:
  preserves tightness under scaling and time reversal; destroys shell identity
  inside the gcd strata {1,3,9}.

HYP-2164 pinch quotient:
  preserves least-positive Res_27 lower-bound status; destroys lift choices.

HYP-2165 owner fibre:
  preserves cheap-pair / positive-measure discharge; destroys arbitrary unit
  representative data outside the bounded canonical fibre.
```

So the tower is a pair of adjoint moves:

```text
quotient/forget  -> lower-bound certificate
reattach fibre   -> owner discharge certificate.
```

This is the concrete n=14 form of HYP-2161's coimage/Yoneda slogan.

## Tournament Analysis

S610 uses proof atoms as tournament vertices, not runners.  The observable is:

```text
(open?, exact-floor?, owner-missing?, normalized?, label)
```

The switch orients from lower proof burden to higher proof burden.  The
fingerprint is transitive:

```text
vertices: 11
directed 3-cycles: 0
SCCs: [1]
Hamiltonian path:
  C1 no-multiple clock exit
  owner no-cheap positive controls
  THM-407 gcd=1 shell orbit
  pinch strict rows
  owner cheap-pair floor AP/V*
  pinch floor AP
  pinch floor V*
  pinch floor 2*AP
  lift/CRT conservativity seam
  THM-407 gcd=3 shell orbit
  THM-407 gcd=9 shell orbit
```

This is deliberately not a cyclic residual.  A nontrivial SCC should appear
only if arbitrary lift/CRT states create mutually incompatible owner obligations.

## Assumption Challenge

Candidate vertices considered:

```text
runners,
residues,
gaps,
pair-sum denominators,
fixed round classes,
slack rows,
shell orbits,
lift CRT states,
proof atoms.
```

S610 chooses proof atoms because the current target is not another search
space; it is a proof obligation.  This quotient preserves the LRC discharge
predicate and deliberately destroys raw runner identity, phase order, and
individual lift choices.

The challenged assumption is that the next improvement must be a bigger finite
enumeration.  The evidence says the useful next theorem is smaller and more
structural: prove lift/CRT conservativity for the atom list.

## Honest Status

This is not a proof of LRC n=14.  It is a proof-program compression.

Open theorem:

```text
For every arbitrary integer lift of a Res_27 candidate that survives the known
clock, dominance, shell, pinch, and owner gates, the lift data either preserves
a strict pinch/owner certificate or normalizes to AP/V*/2*AP.
```

That theorem should be attacked with the integer bounded-CRT automaton and
two-block determinant machinery from HYP-2142/HYP-2107/HYP-2105, not by
reopening the least-positive or canonical owner scans.

## See

`04-computation/lrc_n14_res27_quotient_tower_s610.py`,
`05-knowledge/results/lrc_n14_res27_quotient_tower_s610.out`,
`07-reflections/lrc-n14-res27-quotient-tower-s610.md`,
`05-knowledge/hypotheses/HYP-2175-lrc-dimension-descent-salience.md`,
`05-knowledge/hypotheses/HYP-2165-lrc-n14-res27-fixed-class-bridge.md`,
`05-knowledge/hypotheses/HYP-2164-lrc-n14-res27-pinch-certificate.md`,
`01-canon/theorems/THM-407-twisted-involution-shell-reduction-of-the-LRC-additive-residual.md`,
`05-knowledge/hypotheses/HYP-2161-coimage-yoneda-2nm1-resonance-cancellation.md`.
