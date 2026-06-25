# LRC14, Poincare, And The Observer Frame

The Poincare cue is most useful if it is kept modest and sharp.  LRC is a
worldline problem on a cylinder: each runner is the line `x=v*t`, and loneliness
asks whether one horizontal time slice avoids every `1/14` tube around the
observer worldline.

That makes "boosts" tempting.  But the standard LRC has already fixed the
observer at speed `0`.  If we boost every worldline and keep the observer speed,
nothing changes after recentering.  If we forget the observer speed, the same
operation becomes `v -> v+c`, and the LRC predicate can change violently.  The
scout makes this exact: AP13 and GW13 have safe measure `0`, but translating
the stationary speed set by `+5` gives positive safe measure; the
observer-coupled `+5` boost stays at `0`.

So the real theorem-shape is a symmetry ledger:

```text
exact for anchored LRC:
  permutations, independent sign flips, reflection/time reversal, scale/gcd.

exact only after enlarging the state:
  common boosts carrying observer velocity.

not exact, but useful as warning lights:
  stationary speed translations, Lorentz-like velocity addition.
```

This merges cleanly with the older observer-category split.  Observer-relative
tiling is category 1; metric winding with gap widths is category 2; bare
tournament order and affine scalars are category 3.  The Poincare pass adds the
frame label that explains how one can move between categories without lying:
carry observer velocity, tube metric, and integer-lattice embedding, or name the
lost coordinate as a cocycle.

The surprising little gift is independent sign flips.  Direction is invisible
to `||v*t||`, so negative speeds are an exact kernel for the standard predicate.
But direction is very visible to winding tournaments.  That means sign flips
are not just a cute trick; they are a diagnostic.  If a tournament proof changes
when signs flip, the proof is using orientation data that the LRC predicate has
already quotiented away.  That orientation debt must be reconstructed,
annihilated, or retained.

The next concrete packet fields should be:

```text
observer_velocity_label
relative_speed_normal_form
sign_kernel_status
primitive_scale_gcd
tube_metric_label
worldline_frame_label
boost_cocycle_status
orientation_debt_for_winding
```

This does not close LRC14, but it removes a subtle source of fake motion.  A
Poincare-flavored argument is allowed only if it says whether it is acting on
the anchored LRC problem, the observer-coupled worldline groupoid, or a deformed
tube/lattice residual sector.
