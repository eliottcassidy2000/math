# LRC14 Borel-Baire-Haar Any-Angle Carrier Reflection

S145 gives a clean language for a split that has been showing up everywhere in
the LRC14 endgame:

```text
open witness interval
closed threshold support
boundary-only tight atom
```

AP and Goddyn-Wong are not positive-Haar rows at threshold `1/14`; they are
boundary-support rows.  The near-miss `12->36` and the C27 petals already have
positive strict Haar mass.  That matters because a proof that only tracks Haar
mass will miss the exact tight locus, while a proof that only tracks boundary
points will miss the easy open escape after perturbation.

The any-angle analogy is useful when it is made exact enough.  ANYA works with
interval nodes and taut turns; CWave works with primitive front pieces.  The LRC
analogue is not a grid path but a wavefront of unsafe arcs on the circle or on a
relation-lattice subtorus.  The proposed "Haar-Baire Wave*" carrier keeps three
labels at once:

```text
strict Haar mass
Baire interior
closed boundary support
```

That triple is the point.  It refuses the scalar collapse from "has a witness"
to one number.

The assumption challenge is also live.  The tournament vertices do not have to
be runners, arcs, or residues.  Here the useful vertices are interval fronts,
wall-crossing events, closed boundary packets, subtorus cells, and proof
obligations.  This quotient preserves the open-versus-boundary LRC predicate
and destroys runner ownership unless C27/Farey/unital labels are reattached.

The next proof move I would try is a boundary-support lemma: after the current
reductions, every threshold-safe but strict-Haar-zero row should be AP or GW.
Equivalently, covering should forbid a binding denominator `D=14q-r` with small
positive `r` unless an open safe interval appears.  If that lemma fails, the
counterexample should be a new boundary packet with exact denominator support,
which is precisely the kind of object HYP-2908/THM-572 might be able to lift.
