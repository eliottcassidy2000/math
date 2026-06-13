---
source: codex-2026-06-03-S586
status: exact comparative probe + proof-route synthesis
tags: [LRC, extremality, sleeve-saturation, additive-circuits, hidden-folds, n14, HYP-2119]
---

# Extremality order for the LRC additive/sleeve spine

The right lesson from this session is that extremality is plural.

Several objects can be extremal at once:

- AP is extremal for critical sleeve saturation.
- shifted AP is extremal for raw additive energy.
- hidden fixed-sum families are extremal for virtual fold multiplicity.
- dyadic AP layers are extremal for self-similar structure.
- endpoint-cover residuals are extremal for simultaneous resonance.

But these extremes do not have equal proof value.  S586 turns that into a
direct comparison.

## The calibration rows

The n=14 calibration rows are decisive.

AP has the full hard signature:

```text
M=1/14
mu_safe=0
3-term folds=36
ordered energy=1469
critical saturation=true
clock binders=[1,13]
```

`V*=(1,2,3,4,5,6,7,8,9,10,11,13,24)` has the same floor status:

```text
M=1/14
mu_safe=0
3-term folds=31
critical saturation=true
clock binders=[1,13]
```

So the AP is not the only floor row, but the floor rows share the saturation
signature: zero safe measure, saturation first at the last sleeve, and every
sleeve needed.

Now compare the decoys.

The far-shift AP at n=14 has the same ordered additive energy `1469` as AP,
but it has no visible 3-folds and is far above the floor:

```text
M=5/14
M-delta=+0.28571
mu_safe=0.11734
```

The unit-shift AP is subtler.  It keeps the same ordered energy and still has
`30` visible 3-folds, but the clock is blocked and saturation does not become
critical:

```text
M=0.125
M-delta=+0.05357
mu_safe=0.06122
```

This is the sharper warning.  It is not enough to say "visible folds cause
hardness."  The hard folds are the low folds that couple to the delta-clock
and survive the Cprime/Phi split.

The doubled-apex stress row is useful because it almost fools the criterion:

```text
(1..12,26)
M=0.07407
M-delta=+0.00265
mu_safe=0.00549
```

It is close, but still not critical.  That makes it a good stress row for a
future n=14 proof: the proof should explain the small positive safe measure,
not just classify it as "AP-like."

## What "extreme" should mean now

The old word "extreme" was overloaded.  S586 suggests the following ordered
meaning:

```text
critical sleeve saturation
> Phi / endpoint-gap terminal gates
> visible low 3-fold delta-clock folds
> circuit-free margin floors
> hidden virtual-fold pressure
> dyadic binder profile
> raw 4-term energy.
```

This order is not aesthetic.  It is a certificate order.  The top objects can
actually force or certify floor behavior.  The lower objects expose structure,
but they can be decoys unless labels are preserved until a terminal gate
consumes them.

Raw 4-term energy sits at the bottom because shifted AP has maximal AP-like
energy while becoming very safe.  Hidden virtual folds sit higher because they
remember the missing sum node: S584's hidden row
`(6,11,14,15,16,18,19,23,28)` has hidden `s=34` with multiplicity `4`, and
adding `34` lowers the augmented margin.  But the original row remains safe,
so hidden pressure is diagnostic, not terminal.

After rebasing over THM-400, the same split has a cleaner algebraic name.
Balanced relations (`epsilon=0`) are translation-invariant and observer-blind;
unbalanced relations (`epsilon != 0`) couple to the observer.  The extremality
order is therefore not just "3-term beats 4-term."  It is "observer-coupled
augmentation beats balanced support."  Raw energy is balanced; the fold and
doubling gates are hard only when their augmentation is visible to the lonely
observer.

The dyadic/fractal AP story also becomes cleaner.  The even layers reproduce
the AP pattern under doubling, but the binders at `t=1/n` live in the odd/unit
layer.  So level-4 or hyper-lacunary-looking structure is not automatically
hard.  It is hard only when it touches the delta binder layer or becomes an
endpoint/Phi terminal obstruction.

## Tournament Analysis

The Tournament Analysis uses extremality axes as vertices, not runners.

The pairwise observable is:

```text
(certificate_power, hardness_correlation, extremal_specificity,
 decoy_resistance, maturity)
```

The switch sends the arc to the larger observable, with declaration order as
the tie Hamiltonian path.  The result is transitive, with one Hamiltonian path:

```text
critical_sleeve_saturation
> Phi_endpoint_gap_terminal
> visible_low_3fold_delta_clock
> circuit_free_margin_floor
> hidden_virtual_fold_pressure
> dyadic_fractal_binder_profile
> raw_4term_energy
```

This is the useful version of transitivity for the proof search.  If
critical saturation beats visible folds and visible folds beat raw energy,
then raw energy cannot re-enter above critical saturation without carrying a
new label.  That is the hidden no-cycle fact in certificate language.

## n=14 proof angle

The n=14 target should be reframed as an extremal classification.

Do not prove: "all large additive structures are safe or AP."

Prove instead:

```text
floor row
=> critical sleeve saturation
=> low-fold / delta-clock scaffold
=> either AP/V*-type critical row or endpoint/Phi terminal positivity.
```

This separates the known calibration rows from decoys:

- AP and `V*` are allowed floor rows.
- far-shift AP is dismissed by the raw-energy decoy test.
- unit-shift AP is dismissed by the blocked-clock/non-critical saturation test.
- doubled-apex rows become stress cases for the terminal endpoint/Phi gates.
- hidden-only four-term rows stay in Lemma-A territory unless the virtual sum
  becomes visible.

The endpoint-cover circuit positivity and small-owner machinery matter exactly
where the extremality order says they should: after the low-fold/delta-clock
scaffold has reduced the problem to a terminal resonance question.

## Assumption challenge

Vertices considered: runners, gaps, fixed circle sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits,
hidden sum nodes, sleeve-order positions, proof obligations, and extremality
axes.

Chosen quotient: extremality axes plus calibrated witness rows.

Preserved predicate: which coordinate can actually force `M` or safe measure
to the delta floor.

Destroyed information: exact endpoint-owner geometry, residue languages, and
component-by-component `Phi` ramps.

That destruction is intentional only at this layer.  The terminal proof still
has to restore those labels through HYP-2112, HYP-2108, and HYP-2107.

## Next concrete lemma

The useful next lemma is:

> If a primitive n=14 row is floor-tight or measure-saturated, then its
> saturation curve is critical and its binder set contains a low unit-clock
> scaffold equivalent, after the existing n=14 quotient gates, to AP or `V*`;
> otherwise an endpoint/Phi terminal gate gives positive safe measure.

This is not proved.  But S586 makes the target much less blurry: prove
criticality, not largeness.
