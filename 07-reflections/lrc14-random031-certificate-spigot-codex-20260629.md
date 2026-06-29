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
HYP-3526 tightens that last phrase: `I/Q` are good row-free private-cut
sidecars, but they do not reconstruct the HYP-3490 route, so `R` is still part
of the terminal carry unless a separate route-reconstruction theorem is proved.
HYP-3527 then packages this as a proof contract: ordinary/free-hole emissions
are interface-ready, the bypass/route/vertical clauses retain sidecars, and the
sole open tail is the residual `(45,173)` close-tail lemma.
HYP-3528 is the executable confirmation of that package: the proof ABI is
sidecars-to-safe-emission, and the single red clause is now literally
`residual_pair_close_tail = residual (45,173) + R + no_hidden_tail_guard`.

## Owner Carry

The companion owner-carry audit checks the local emission predicate directly.
The lower bypass gates equal the pure-bypass transport word; all `12` bypass
cells have owner word `(23,93,113)`; all six mirror pairs preserve that word;
and both branch sheets have the same boundary owner word `(23,93,147,169)`
with intersection `(169,)`.  These checks justify emitting transport and
branch lift.

The result explicitly reports a typed residual:

```text
terminal_residual_closed=False
residual_after_branch=(45,173)
global_flow_minus_bypass=(45,55,147,169,173)
```

Owner `55` is the warning label.  It says the global flow union is not the
local emitter.  The proof target is no longer "make the bypass carry all seven
owners."  The target is "show no legal two-owner carry `(45,173)` survives
after ordinary-route emission, HYP-3511 doublet buffering, transport, branch
bracket, free-hole caps, HYP-3490/HYP-3513 route sidecar `R`, and the
vertical-halfturn guard are attached."

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

The companion owner audit says the same thing for owner-boundary persistence.
`flow_class`, `owner_boundary_persistence_class`, `owner_presence_word`, and
`seam_debt_word` are safe forgetful relations for the owner target.
`component_size`, `endpoint_rank_word`, `mirror_closed_shadow`,
`raw_owner_count`, and `seam_owner_count` are not, because they mix terminal
classes or erase the residual.

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

HYP-3525 is the general version of that slogan.  It should become the schema:
a visible packet may emit a proof token only when every legal hidden-tail
completion has the same target; otherwise the proof must hold route `R`, an
owner word, a branch-boundary lift, a sheet-PGF bucket, or named residual debt.
HYP-3526 is the first concrete instantiation for the route coordinate: emit
the private-firewall status through `I/Q` if desired, but keep `R` attached to
terminal certificates.

Concurrent THM-578/HYP-3529 is adjacent but not an input theorem here.  Its
doublet R-tail is the global LRC14 Obligation D, not the HYP-3511 random031
free-hole doublet buffer.  The useful transferable principle is exact
decomposition plus finite carry: close what has a local identity, bound the
remaining tail, and avoid replacing a named tail by a sharp-but-unneeded
constant.  For random031 that means the one-step free-hole doublets are closed
inside the stream, while residual `(45,173)` remains typed owner-boundary debt.

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
  + route_sidecar_R                 -- retained by HYP-3526 until reconstructed
```

The target theorem should say every component event either emits a terminal
certificate, buffers one doublet predigit, applies the branch-boundary lift, or
leaves only the named residual pair.  That would let the random031 terminal
packet be checked by a small state machine rather than by the full component
ledger.  HYP-3527 is the downstream checklist for that state machine: formalize
the three ready interfaces, keep the four carry clauses explicit, and close the
two-owner residual tail.  HYP-3528 sharpens that checklist to one red theorem:
prove `residual_pair_close_tail` and keep scalar shadows out of the proof ABI.
