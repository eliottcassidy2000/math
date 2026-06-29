# LRC14 Random031 Spigot Stream Reflection

The spigot-algorithm prompt is useful only after translation into an exact
finite sidecar.  The imported idea is not "digits of pi" as an analogy; it is
the scheduling discipline:

```text
emit left-to-right
discard emitted output
retain only bounded deferred carry
```

For random031, the HYP-3521 terminal ledger already supplies the digit alphabet:
ordinary route certificates, free-hole single certificates, free-hole doublet
certificates, and the unique bypass owner-boundary certificate.  HYP-3522
supplies the carry law: transport owners emit first, branch-boundary owners
reduce the carry, and only `(45,173)` remains.

The new audit makes this finite.  In q-witness order:

```text
79 component events emit 77 terminal certificates.
64 ordinary components and 10 free-hole singles emit immediately.
2 doublet clusters each need one buffered component and close at the next event.
The bypass opens carry (45,147,169,173) at rank 45.
The branch-boundary lift applies at rank 46 and leaves (45,173).
```

That gives the current random031 proof a streaming invariant:

```text
emitted_prefix + predigit_buffer + owner_carry + route_sidecar_R
```

The meaningful proof reduction is not that the residual vanished.  It did not.
The reduction is that most of the finite terminal table no longer needs to be a
live proof object after emission.  The non-streamed payload is sharply named:
the residual owner pair `(45,173)` and the HYP-3513 route sidecar `R`.

## Quotient Lesson

The quotient audit is a useful canary.  `component_type`, `terminal_class`, and
even `terminal_class+cluster` all mix stream actions.  They forget whether a
free-hole doublet component is the first predigit or the closing component, and
they forget whether an ordinary component is also the branch-boundary lift that
finalizes the bypass certificate.

Adding the small spigot state makes the stream action pure:

```text
terminal_class_plus_spigot_state: mixed_fibers=0
component_index: mixed_fibers=0
```

This is exactly the kind of compression the proof needs: not raw component
identity, but a small state variable that records why a later event can still
change the emitted certificate.

## Micro Release Test

The smaller release/hold quotient test is a useful companion to the stream:
`full_spigot_packet` is the only tested carrier that emits transport, branch
lift, and free-hole caps while still holding both residual `(45,173)` and
route `R`.  The failure modes are instructive: `owner_filtration_without_route`
drops the private-firewall route, `transport_plus_bracket` drops residual and
route, and raw counts or raw seven-owner shadows either emit no certified
prefix or hide the carry split.  So the information-theoretic slogan is not
"compress harder"; it is "compress only to a sufficient statistic for terminal
emission."

## Assumption Challenge

Candidate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid
circuits, proof obligations, output certificates, carry states, quotient
guards, owner rows, free-hole packets, and stream states.

Chosen vertices are stream states with predigit/carry buffers.  They preserve
the random031 terminal predicate because each emitted token is a HYP-3521
terminal certificate, HYP-3511 doublet collapse is handled before emission, and
HYP-3522 owner carry is retained until the branch-boundary lift reduces it.

What this destroys: raw runner order, raw component index after emission, full
row-name memory, and scalar cell-count shadows.  That destruction is legal only
because the stream state keeps the predigit/carry information and because
HYP-3513 says route sidecar `R` is still required unless reconstructed.

## Next Proof Move

Turn HYP-3523 into a Lean-facing scheduler lemma:

```text
Random031StreamState =
  emitted_prefix
  + pending_free_hole_doublet?      -- size <= 1
  + owner_carry                     -- width <= 4, then residual (45,173)
  + route_sidecar_R
```

The target theorem should say every component event either emits a terminal
certificate, buffers one doublet predigit, applies the branch-boundary lift, or
leaves only the named residual pair.  That would let the random031 terminal
packet be checked by a small state machine rather than by the full component
ledger.
