# LRC14 k=8 Reflection-Fold Coordinate Resurrection

Anchors: HYP-3138, HYP-3134, HYP-3133, HYP-3132, HYP-3122, HYP-3118,
HYP-3116, HYP-3110, HYP-3085, THM-577, T1203, LTI-264, LTT-162,
OPEN-Q-108.

The useful move from the "unrelated" archive is not to import wallpaper
groups, theta functions, circuit complexity, or De Moivre as analogies.  It is
to ask the same technical question they all ask in this repository:

```text
when a symmetry quotient forgets a coordinate, what is the sidecar or adjoint
section that resurrects it before the next proof obligation needs it?
```

For the current LRC14 frontier, that question lands exactly on the `k=8`
hard row.  HYP-3132 says the resolvent quartic folds to the even biquadratic
`u^4 - 5u^2 + 4`.  HYP-3122 says the load-bearing functional is the gK8/phi4
quartic lane.  HYP-3118 says a quotient is legal only after destroyed
coordinates are repaired.  HYP-3116 says the proof still needs endpoint
`Phi/P` activation data.

The scout tests the obvious fold:

```text
q -> (q0+q6, q1+q5, q2+q4, q3).
```

The good news is strong: for primitive k=8 rows in spans `14`, `15`, and `16`,
the even fold is injective.  It has `0` collision fibers over `3431`, `6434`,
and `11432` rows.  So, at least in these bounded banks, the fold is not just a
scalar shadow.  It can serve as a finite lookup / right-adjoint candidate.

The guardrail is just as important: the best row `(0,1,2,3,4,5,6,7)` has
nonzero odd leakage

```text
(q0-q6, q1-q5, q2-q4) = (451/1470, 142/735, 131/1470).
```

So the reflection fold is not evidence that the odd coordinates vanish.  It
is evidence that they may be recoverable in the bounded-core bank from the
folded coordinate, or else must be carried as named endpoint/observer debt.

This turns the niche archive into a practical proof target:

```text
k8_reflection_fold_certificate =
  even_folded_miss_distribution
  + odd_coordinate_resurrection_table
  + de_moivre_biquadratic_resolvent_word
  + endpoint_phi_activation_status
  + observer_gluing_or_named_debt
```

The 17 wallpaper groups and 230 space groups are not proof ingredients by
themselves; they are reminders that finite quotient catalogs must carry their
stabilizer and destroyed-coordinate ledgers.  Jacobi theta is the signed
even/odd channel language.  Circuit complexity is the missing-input audit.
A000568 is the global-consistency quotient warning.  The concrete next theorem
is a k=8 fold-adjoint lemma for the bounded-core bank, tied to endpoint
`Phi/P` repair before quotienting.
