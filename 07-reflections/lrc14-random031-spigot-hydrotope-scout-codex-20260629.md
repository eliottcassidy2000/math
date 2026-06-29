# Reflection: Random031 Spigot/Hydrotope Scout

The useful part of the spigot analogy is not digit computation.  It is the
discipline of irreversible emission: once a digit is emitted, the algorithm
must have a certificate that future carries cannot change it.  HYP-3524 turns
that into a proof interface for random031:

```text
seven-owner seam
  -> emit transport (23,93,113)
  -> emit bracket lift (147,169)
  -> prove residual tail (45,173)
```

This is a cleaner way to talk about HYP-3522.  The residual pair is no longer
just "whatever is left" after a set difference.  It is the terminal tail of a
monotone emitter schedule, with no duplicate emits and no tail growth.  That
is exactly the kind of finite state object that should be easier to expose to
Lean.

The hydrotope paper was useful for a different reason.  It suggests asking
which chamber/sign quotients are allowed to remember a sliced box.  On
random031, the answer is a guardrail:

```text
residue_mod14 chambers mix (45,173) with (113,147) and (147,169)
centered residues mix even more
filtration-layer weights mix heavily
owner-support-cell weights isolate transport, boundary, and residual
```

So the chamber idea does not prove the residual lemma.  It tells us which
quotients are illegal.  A chamber sign becomes proof-facing only when the
weight already retains enough owner-support geometry to distinguish the
residual pair.  Scalar sliced-box volumes are beautiful but too lossy here.

The new high-leverage target is a no-hidden-tail lemma:

```text
After transport and bracket emitters fire, every legal terminal quotient with
route sidecar R still sees residual tail exactly (45,173).
```

That is smaller than the old owner-boundary problem and sharper than raw
hydrotope or residue language.  If the next agent wants to formalize something,
the best interface is probably:

```text
EmitterState(seam_word, emitted_word, tail_word, certificate, route_sidecar)
emit_transport : State 7 -> State 4
emit_bracket   : State 4 -> State 2
close_tail     : State 2 -> terminal discharge or named debt
```

The unexpected result is that owner-support-cell weights isolate the target
subsets as singleton chambers.  That may be the computational fingerprint of a
real proof: the residual owners are not just labels; they are the only pair
with that support chamber once the random031 receiver has been fully
packetized.
