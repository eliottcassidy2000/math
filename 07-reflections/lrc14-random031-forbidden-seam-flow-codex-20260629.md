# Random031 as a Forbidden Seam

HYP-3484 is the surgery refinement of the HYP-3483 recursion-flow comparator
and HYP-3482 seam atlas: it turns the random031 exception from a hard wall into
a flow-complement statement.  The row has four isolated dead islands, two
mirror pairs, and zero dead-projection edges, but unlike the HYP-3478
small-touch pockets it also has a max-delta mirror seam carrying all seven
owner labels.

The incoming HYP-3480 singleton-current audit makes the separation sharper:
the six small-touch/no-hard rows all have component-level branch-unit mirror
certificates, while random031 has `0/4` complete branch-unit components and
`0/2` mirror-unit pairs.  That pushes this row further away from a singleton
current proof and toward the seam-complement proof route.

The important experiment was surgical.  Removing the max-delta seam gates does
not change q=`14V` phase routing at all: `any_gate_hits=242` and
`no_gate_hits=40` remain fixed.  Removing the two lower-delta same-component
bypass gates changes exactly `12` witnesses: gate hits drop to `230` and
no-gate hits rise to `52`.  That is the cleanest evidence so far that the
hard pair is a forbidden seam whose complement carries phase flow.

This also clarifies the older `n+2` versus `n*2` recursion language.  The
additive face sees `C=2n-1=27`, punctures, and owner carries.  The doubling
face sees `u=2t mod 1` and routes around the seam.  Random031 is exactly where
those faces disagree: the carry view sees a charged seven-owner seam, while
the phase-doubling view only pays a lower-delta bypass toll of `6+6`.

Next useful experiments:

- Carry-lift stress: preserve or break owner classes modulo `27` and endpoint
  residues `(1,9) mod 14`, then watch whether the seam/bypass split persists.
- Bypass basis widening: add candidate secondary low-rank channels and see
  whether the `12` bypass hits remain uniquely supported.
- Puncture filling: compare HYP-3478 singleton pockets with the four
  random031 punctures and test whether isolated dead islands matter only when
  a seven-owner seam links their surrounding flow components.

The proof target should say "forbidden seam plus complement flow", not "wall".
