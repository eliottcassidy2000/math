# Polarized Cyclotomic Support And Orbit-Aware Compression

*codex-2026-06-28.  The owner asked for a new creative angle or two on the
remaining LRC14 proof targets.  This reflection records two exact bounded-bank
signals from HYP-3203.  It continues HYP-3154's Joukowski/De Moivre bridge,
HYP-3153's Lee-Yang/Worpitzky/quartic packet, and HYP-3160/HYP-3161's
ferromagnetic covariance target.*

## The Good New Angle

The seventh-cyclotomic ideal should not be read as "minimize distance to
uniform."  That is false.  Over the bounded anchored k=8 bank, AP/consec is
only rank `19` for minimum nontrivial cyclotomic energy; split-block rows are
closer to uniform but worse for coverage and covariance.

The useful object is the AP support functional:

```text
<q(E)-1/7, q_AP-1/7> <= ||q_AP-1/7||^2.
```

The exact scout verifies that this is sharp on the bounded bank, with only the
AP row and its doubled dilation maximizing it.  This looks like a real proof
target because it is linear in the miss-count vector.  It can potentially be
fed to signed SPEC, finite moment cones, or a Delsarte/MacWilliams-style
support certificate.

The conceptual correction is important: the Joukowski bridge is directional.
The AP does not minimize the cyclotomic defect; it maximizes the correct
polarized defect, the one that points toward coverage and ferromagnetic
coherence.

## The Useful Refutation

The naive compression proof is false.  If compression means "replace a larger
runner by a smaller missing runner," then the bounded k=8 bank has `19` local
no-improvement traps and greedy improvement paths get stuck on `919` states.
The doubled AP dilation is itself a trap.

This does not kill compression.  It says compression has to be theorem-facing:
first quotient or sidecar the dilation/mirror/two-block orbit, then use the AP
support functional as the Lyapunov object.  Raw coordinate-left-shift forgets
the exact information the LRC predicate needs.

## What To Try Next

1. Treat `q_AP-1/7` as a dual certificate and search for a signed-SPEC proof
   of the support inequality.
2. Classify the `19` local traps into dilation, mirror, and two-block trap
   families; each family should either be discharged by the support inequality
   or routed to a finite sidecar.
3. Test whether the same AP projection remains sharp on larger primitive banks
   after quotienting dilation.

The meta-lesson matches the controlled-forgetting ladder: cyclotomic norm and
left-compression are tempting scalar shadows.  The proof object is a function
with a named target and named lost coordinates.
