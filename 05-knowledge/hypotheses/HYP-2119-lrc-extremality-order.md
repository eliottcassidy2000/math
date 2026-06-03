---
id: HYP-2119
status: OPEN proof-program; S586 separates certificate-bearing extremality from decoy extremality
source: codex-2026-06-03-S586
related: [THM-400, HYP-2118, HYP-2117, HYP-2116, HYP-2115, HYP-2114, HYP-2113, HYP-2112, HYP-2108, HYP-2107, HYP-2101, HYP-2100, HYP-2096]
---

# HYP-2119: LRC extremality is ordered by certificate power, not scalar size

## Claim

In the additive/sleeve branch of the Lonely Runner Conjecture, "extremal" is
not a single scalar condition.  The proof-relevant order is:

```text
critical sleeve saturation
> Phi / endpoint-gap terminal gates
> visible low 3-term delta-clock folds
> circuit-free margin floors
> hidden virtual-fold pressure
> dyadic binder profiles
> raw 4-term additive energy.
```

The top axes can certify or nearly certify the LRC boundary.  The lower axes
are diagnostics unless they become labelled input to a terminal gate.

Equivalently: an LRC proof should not ask which row maximizes a familiar
statistic.  It should ask which extremal statistic still has certificate power
after AP, shifted AP, `V*`, and hidden-fold controls are compared in the same
clock.

## Evidence From S586

S586 computes exact maximin `M`, exact safe measure, saturation curves, 3-term
fold counts, ordered additive energy, hidden virtual-fold pressure, clock
binders, and a Tournament Analysis over extremality axes.

At `n=14`:

- AP is truly extremal: `M=1/14`, `mu_safe=0`, `36` visible 3-folds, ordered
  energy `1469`, critical saturation true, clock binders `[1,13]`.
- Far-shift AP is raw-energy extremal but safe: same ordered energy `1469`,
  no visible folds, `M=5/14`, `M-delta=+0.28571`, `mu_safe=0.11734`.
- Unit-shift AP is fold-rich but not critical: same ordered energy `1469`,
  `30` visible folds, but `M=0.125`, `M-delta=+0.05357`,
  `mu_safe=0.06122`.
- `V*=(1,2,3,4,5,6,7,8,9,10,11,13,24)` is the non-AP calibration floor row:
  `M=1/14`, `mu_safe=0`, `31` visible folds, critical saturation true, and
  the same binder stratum `[1,13]`.
- The doubled-apex stress row `(1..12,26)` is nearly critical but not tight:
  `M=0.07407`, `M-delta=+0.00265`, `mu_safe=0.00549`.
- The S584 hidden row `(6,11,14,15,16,18,19,23,28)` has no visible folds and
  a hidden sum `34` of multiplicity `4`, but remains safe:
  `M-delta=+0.13529`, `mu_safe=0.11225`; adding the hidden sum gives
  augmented margin `+0.06909`, so the virtual pressure is real but not a
  certificate by itself.

The deterministic sample audit supports the same split.  At `n=8` and `n=10`,
the correlation of hardness `-margin` with 3-term count is about `+0.55`, while
the correlation with ordered additive energy is negative (`-0.26`, `-0.25`).
Covering hardness `-mu_safe` behaves similarly: about `+0.55` with 3-term
count and near zero/negative with energy.  This is not a proof, but it is a
good warning label: raw energy points in the wrong direction in these controls.

After rebasing over THM-400, this order has a sharper invariant language.
Balanced relations (`epsilon=0`) are translation-invariant and observer-blind;
unbalanced relations (`epsilon != 0`) are observer-coupled.  Thus raw 4-term
energy is a decoy because it is balanced, while visible fold/doubling axes are
hard only when they carry nonzero augmentation into the observer.

## Tournament Analysis

S586 uses extremality/certificate axes as tournament vertices, not runners.

The pair observable is:

```text
(certificate_power, hardness_correlation, extremal_specificity,
 decoy_resistance, maturity).
```

The switch sends the arc to the lexicographically larger observable, with
declaration order as the tie Hamiltonian path.  The resulting tournament is
transitive:

```text
critical_sleeve_saturation
> Phi_endpoint_gap_terminal
> visible_low_3fold_delta_clock
> circuit_free_margin_floor
> hidden_virtual_fold_pressure
> dyadic_fractal_binder_profile
> raw_4term_energy.
```

Score histogram is `{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, directed 3-cycles `0`,
and Hamiltonian path count `1`.

This is the transitive form of the hidden no-shortcut fact: raw energy cannot
jump above hidden pressure or visible folds without losing certificate power.

## Proof Route

1. Prove a critical-saturation lemma: if a row reaches the LRC floor or has
   zero safe measure, its sleeves must be critical in the S578 sense:
   saturation first occurs at the last relevant sleeve and every sleeve is
   needed.
2. Prove critical saturation forces a visible low-fold / delta-clock scaffold,
   not merely high additive energy.  In THM-400 language, the scaffold must
   expose observer-coupled nonzero augmentation, not just balanced support.
   Unit-shift AP is the warning that many folds are still insufficient if the
   clock is blocked.
3. Route multiple-of-`n` or endpoint-cover residuals to HYP-2112/HYP-2108/
   HYP-2107 terminal gates; this is where `Phi` and endpoint positivity regain
   certificate authority.
4. Treat hidden 4-term structure as a labelled middle layer.  It stays in
   Lemma-A territory until the virtual sum becomes a real runner, endpoint
   gate, residue state, or fold/shield certificate.
5. For `n=14`, use AP and `V*` as the critical-saturation calibration rows and
   unit-shift/far-shift AP as decoy controls.  The classification target is not
   "large additive structure" but "critical saturation plus terminal gate".

## Assumption Challenge

Candidate tournament vertices considered in this session included runners,
gaps, fixed circle sections, section boundaries, wall crossings, residues,
cover arcs, Fourier modes, matroid circuits, hidden sum nodes, sleeve-order
positions, proof obligations, and extremality axes.

The chosen quotient preserves the LRC predicate "which coordinate can actually
force `M` or safe measure down to the delta floor."

It destroys endpoint-owner geometry, residue languages, and component-by-
component `Phi` ramps.  Those are not discarded permanently; they re-enter only
as terminal gates through HYP-2112, HYP-2108, and HYP-2107.

The challenged assumption is that an extremal scalar statistic is automatically
a hard LRC statistic.  Shifted AP refutes this for raw 4-term energy, and
unit-shift AP refines the warning: even many visible folds are not enough if
the low clock scaffold is not critical.

## See

`04-computation/lrc_extremality_order_s586.py`,
`05-knowledge/results/lrc_extremality_order_s586.out`,
`07-reflections/lrc-extremality-order-s586.md`,
THM-400, HYP-2118, HYP-2117, HYP-2116, HYP-2115, HYP-2114, HYP-2113,
HYP-2112, HYP-2108, HYP-2107, HYP-2101, HYP-2100, HYP-2096.
