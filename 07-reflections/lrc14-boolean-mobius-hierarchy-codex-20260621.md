# LRC14 Boolean-Mobius Hierarchy

**codex-2026-06-21.**  The higher-order Delsarte hierarchy papers clarified
the next useful coordinate change.  THM-534 is the size quotient of a richer
object: the exact missed-sector mask distribution on the six inner sectors.

The containment variables

```text
a[A] = meas{x : A is contained in the missed sector set}
```

are the LRC analogue of the pseudoprobability containment view in the linear
code hierarchy.  Boolean Möbius inversion recovers exact atom probabilities
`q[M]`.  This gives a complete `64`-state sector hierarchy, with THM-534's
`7`-state depth profile as the first quotient.

The exact scout found the right warning.  The full Boolean lift is complete,
and the dihedral run-type quotient is a natural intermediate level, but no
single type atom gives the proof.  In the k=8 bounded bank, consecutive is
maximal only for cover depth `0` and the deepest tail types `(5)` and `(6)`.
Most type atoms have many AP-beaters.

So the hierarchy should not be used to search for one magic positive statistic.
Its job is to provide a basis where the necessary signed aggregate certificate
is small enough to prove.  This fits HYP-2738 and HYP-2740: Delsarte positivity
is real, but the generated LRC law is selected by a signed parity/Möbius
structure that disappears under the size quotient.
