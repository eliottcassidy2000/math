# LRC14 Residual Certificate Teeth

codex-2026-06-26-S197

The user asked to keep working new LRC angles.  I pushed on the residual left
by the status gate rather than widening the full packet bank again.

S194 already isolated the important fact: the coarse ET+Henselian-unit gate has
`0` mixed boundary/open fibers, but leaves `15` route-mixed fibers containing
`38` packets.  All `38` are open.  So the next job is not to prove loneliness
again; it is to schedule which open certificate should handle each residual
packet.

The new S197 script parses the stored S194 residual ledger and tests proof
carriers on those `38` packets.  Topology alone helps but is not enough:
it leaves `3` mixed route classes.  Unit-scale alone is also not enough:
it leaves one large mixed class.  Exact `M` is surprisingly not enough as a
fallback on this stored residual ledger, leaving `2` route-mixed classes.

The useful tooth is small:

```text
topology compact signature + unit-scale tooth
```

or even the compressed version:

```text
(safe topes, quotient defect, open beta0 - closed beta0) + unit-scale tooth.
```

Both joined carriers give `21` residual fibers and `0` route-mixed fibers.

This is evidence for a proof route scheduler after status purity:

```text
status gate proves open/boundary predicate
topology+scale teeth choose q-witness versus covering/Haar/nested route
HYP-3031 zeta repair classes explain what kind of forgotten coordinate was repaired
```

Assumption challenge: the vertices are residual proof carriers, not runners,
routes, or raw scalar magnitudes.  The quotient preserves already-open status
and residual certificate scheduling.  It destroys row identity, exact ET
address, and most magnitude information.  The guardrail is that S197 parses a
stored output file; the next pass should promote `residual_topology_bucket` and
`unit_scale_tooth` into packet sidecars and rerun directly.
