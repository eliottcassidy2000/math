# LRC14 Lean Median Center-Control Frontier

Source: codex-2026-06-26-S241

I formalized the HYP-3070 median sidecar idea in Lean as
`TournamentH7.LRCMedianCenterControl`.

What is now Lean-checked:

```text
RouteLeaf has 15 leaves.
RouteTriple has 455 three-leaf packets.
RawRouteCenter is centerless for every route triple.
LegalSidecarCenter has a unique center for every route triple.
The primitive-owner split's expected center is primitive_period_router.
CenterControlPacket is now a concrete proof-bearing certificate shell.
centerControlPacket_soundness verifies that the packet numeric fields imply Mreach >= 1/14.
```

The most important theorem is not a proof of LRC14. It is the exact interface:

```text
lrc14_from_center_control :
  CenterControlCoverage ->
  CenterControlSoundness ->
  LRC14Statement

lrc14_from_center_control_coverage :
  CenterControlCoverage ->
  LRC14Statement
```

This is a good result because it forces honesty. The median/sidecar work is not
yet close to proving LRC14 by itself. After the packet shell was made concrete,
the logical frontier is:

```text
CenterControlCoverage:
  every nonzero 13-speed family has a legal proof-bearing center-control packet.

Packet numeric payload:
  each packet must carry witnessFloor >= 1/14 and witnessFloor <= Mreach v.
```

The targeted Lean check passed:

```text
lake env lean TournamentH7\LRCMedianCenterControl.lean
```

The full `lake build` remains too expensive in this checkout and timed out.
Checking `TournamentH7.lean` also stopped before the new module because an
existing aggregate dependency object, `LRCBindingPair.olean`, is missing from
`.lake/build`. So the verified claim is targeted-module verification, not a
fresh full-project build.

The closeness readout:

1. **Already solid in Lean:** the LRC predicate, denominator sieve, compact
   `Mreach -> lonely` handoff, many arithmetic kernels, and now the HYP-3070
   center-control interface shape.
2. **Still abstract/opaque:** actual HYP-2963 packet-bank objects, `shapeOf`,
   `witnessG2`, and the measure/event definitions feeding the witness route.
3. **Main missing construction:** produce non-tautological `CenterControlPacket`
   values for actual HYP-2963 rows, with the numeric fields filled by existing
   AP/GW, witness-floor, Fejer, Haar/Ramanujan, or residual discharge theorems.

The new packet record has the fields that were previously only prose:

```text
raw_route_clique_center_status
legal_sidecar_tree_center_status
median_center_expected_page
center_control_exit
finite_witness_floor
soundness_to_Mreach
```

The next useful Lean move is therefore concrete instantiation, not more
interface design: fill the record for one known AP/GW boundary row, then for one
positive residual-router family. If those two rows cannot be filled without
using `Mreach` tautologically, the sidecar vocabulary is still too vague.
