# Lifted Packet Divergence

The HYP-2883 packet graph gave us a conservation law in the finite mod-7
kernel.  S102 asks what that conservation law becomes inside the real LRC
tail.

It does not lift for free.  Using the four residue-`1` core
`(1,8,15,22)`, the finite identity

```text
loop(a) + sum_b edge(a,b) = 0
```

turns into a nonzero reciprocal-lift divergence.  For the start-aligned lift,
the `H=12` L1 divergence is about `0.01934`.  Moving the pair up one 7-layer
reduces it to about `0.00610`.  That is not a proof, but it is the right kind
of object: the obstruction is now a local divergence defect, not a mysterious
global failure of signed cancellation.

This also reconciles HYP-2883 with HYP-2633.  HYP-2633 showed that finite
packet signs can flip under bounded reciprocal lifts.  S102 says the finite
packet identity is still useful if we ask it to control divergence after a
wall filtration, not raw signs packet-by-packet.

The LRC proof target should now be written as:

```text
finite packet current
  + low-height wall ledger
  + additive-frequency Abel summation
  => support-six reciprocal tail bound.
```

That is a better fit for the genuine-wide branch than a raw cluster census.
High additive energy should force the divergence into a finite wall ledger;
low additive energy should leave a spread packet graph where HYP-2636's
channelwise Cauchy/Abel estimates have room.

The incoming winding-tournament scar note gives the same object in a different
language.  Non-lonely phases have more directed 3-cycle content; they are the
cyclic/no-sink phases.  Lifted packet divergence is the reciprocal-tail shadow
of that cyclic excess.  A proof should make that sentence literal: coherent
cyclic excess is a low-height wall, and incoherent cyclic excess is an
Abel-summable packet current.
