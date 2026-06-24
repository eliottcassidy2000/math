# LRC14 Haar-Baire Taut Boundary Finite Check Reflection

The useful surprise in S146 is that AP and Goddyn-Wong have the same boundary
owner skeleton:

```text
(1,13), (5,9), (3,11), (3,11), (5,9), (1,13)
```

with left/right orientations around the six denominator-14 boundary points.
The speed `12` is not one of the active owners on this skeleton, and neither is
the replacement `24`.  That makes the GW move a hidden transfer at the boundary
level: Haar/Baire sees the same closed endpoint set, while C27/unital sees the
`H12 -> D3` legality.

That is exactly the kind of split the proof route needs.  Haar/Baire fronts are
good for saying whether an interval has opened.  They are bad at remembering
why a hidden nonactive transfer is permitted.  So HBT*/Haar-Baire Wave should
not replace the C27/unital machinery; it should decide when to use it.

The finite evidence is also clean:

```text
one-swap add<=160: boundary-only rows are AP and 12->24 only
two-swap add<=40:  no boundary-only rows
```

The S138 two-swap rows are therefore open fronts, not new endpoint atoms.  The
first one simply inherits the petal `10->20` front; the second is a union of the
petal front and the near `12->36` front.  This makes the proof shape feel less
mysterious: after AP/GW, low-frontier perturbations split into named open
fronts with active owner pairs like `5/36`, `7/20`, or `1/26`.

The next lemma I would try to prove is boundary-owner skeleton rigidity.  A
zero-interior row with threshold support should have to preserve the AP/GW
owner skeleton; once that skeleton is fixed, the only hidden AP one-swap
replacement in the scanned range is `12->24`.  The missing unbounded step is to
show no far replacement can shadow all six skeleton gaps without either opening
an interval or becoming exactly the C27/GW transfer.
