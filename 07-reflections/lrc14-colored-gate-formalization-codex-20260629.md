# LRC14 Colored Gate Formalization Reflection

The useful result of this formalization pass is negative in a productive way:
Lean forces the color story to say exactly where the proof is still missing.

HYP-3471 found a strong empirical theorem target:

```text
dead_components(row) > 0
  => rank <= 2 E/branch survivor gate.
```

HYP-3473 turns that into a formal interface rather than another paragraph of
confidence.  The producer is an abstract `DeadCoverEBranchSoundness`; the
consumer is `eBranch_gate_of_dead_cover`.  That split is important because it
prevents the current bank count `130/130` from being mistaken for a proof of
the geometric implication.

The terminal side is similarly clean.  A `ColoredGateTerminalPacket` is allowed
to choose AP84 corridor splice, AP84 color-grid placement, random031 gluing,
component conductance/Menger, owner-current, two-adic descent, signed SPEC, or
named residual debt.  But the packet is useful only if it carries the final
field:

```text
(1 : Real) / 14 <= Mreach v.
```

Then the existing skeleton gives LRC14 conditionally.  This is the right
direction of pressure: sidecars are not decorations; they are only legal when
they produce the Mreach floor.

The formal module also makes the quotient guardrail sharper.  Same-branch and
cross-branch gates are not E/branch gates, and Lean proves the separation at
the endpoint-kind level.  Raw color count is therefore a shadow.  Endpoint kind
alone, numeric mod-14 residue, AP84 four-color packets, typed residue words,
structural sidecars, and the full colored word all sit at different levels of
payload retention.  The Lean carrier ledger records that the theorem target and
the full word tie for top score, while raw color count is the sink.

The assumption challenge is now explicit.  I could have made tournament
vertices be runners, gates, components, residues, or arcs, but none of those
alone preserves the predicate that actually closes LRC14.  The chosen vertices
are proof obligations and terminal packets.  That choice destroys row geometry
and enumeration details, so the formal interface names exactly where they have
to re-enter: `DeadCoverEBranchSoundness` for the finite gate implication and
`ColoredGateGlobalCoverage` for the global route.

The next best proof-facing job is not another color census.  It is to prove the
producer.  The current/cut formulation still looks strongest: a dead component
with no low-rank E/branch leak should force impossible branch-current
divergence or a Menger cut whose terminal set is exactly the typed gate
boundary.  Once that exists, HYP-3462, HYP-3470, HYP-3461, HYP-3460, HYP-3459,
HYP-3458, HYP-3455, and HYP-3451 have a formal packet slot to fill.
