# LRC14 Resolvent-Packet Synthesis

The De Moivre-style quintic supplied by the user is useful because its
resolvent refuses to be a scalar shortcut.  The roots of the quartic are
`2,-4,8,-16`, but the visible coefficients `120` and `320` are not roots; they
are the pair and triple elementary-symmetric layers of the branch orbit.

That is exactly the safe-forgetting rule the LRC14 endgame keeps rediscovering.
Raw `p0`, raw A000568 counts, raw Lee-Yang radius, raw circuit size, and raw
De Moivre residuals all forget too much.  The live object is the packet that
retains the middle data:

```text
Q-floor constants
signed SPEC low resonance
Parseval tail
bounded-core Lee-Yang / phi4 status
far-push-out monotonicity
edge tail/tip deletion sectors
finite-address or observer-gluing exit
```

Incoming HYP-3131 is the main simplification: far elements help.  They push the
miss-PGF zeros outward, so the multi-far floor reduces to bounded-core
Lee-Yang plus a monotonicity proof.  HYP-3130 closes the Q/apex block and shows
the absolute minorant envelope cannot handle the coupling.  HYP-3129 supplies
the missing signed object: exact low resonance on the `14Z` lattice plus a
Parseval tail, already certifying `Rprime >= 0.64178` on the tested family.

So the remaining proof is smaller than the older story:

```text
covering branch S = R union 14Q
  Q-block positive       -> closed by cap constants / Lee-Yang / minorant
  far placements help    -> reduce to bounded core
  Rprime coupling        -> signed SPEC constant chase
  bounded-core bridge    -> rho>1 / phi4 implies positive coupling floor
  packet legality        -> edge tail/tip sectors or named residual debt
```

The 120/320 motif is therefore a warning label: keep the middle symmetric
payload before quotienting.  In the resolvent it is `e2/e3`; in LRC14 it is
the signed SPEC and edge-witness packet.
