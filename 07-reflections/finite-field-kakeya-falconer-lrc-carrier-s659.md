# Finite-Field Kakeya/Falconer as an LRC Carrier S659

The useful import from finite-field Kakeya/Falconer is not the drama of the
big conjectures.  It is the small carrier discipline.

Kakeya says: keep a line in every direction.  Falconer says: measure the
distance fibers.  The repo translation is immediate:

```text
direction cover       -> clock/direction obligations
line concurrency      -> owner/multiplicity side channel
distance support      -> scalar safe/distance support
pinned distance       -> observer-coupled witness
unit directions       -> unit spine / direction support
finite-field character -> Paley/Fourier phase
```

The computation made the scalar-collapse warning very concrete.  Over `F_p^2`,
the odd-prime parabola Kakeya carrier has size `(p^2+2p-1)/2`; S659 verifies
this as the exact minimum for `p=3,5,7`, with exact minima `7,17,31`.  But those
same carriers already realize every quadratic distance for `p=3,5,7,11,13`.
At the distance-support level, the hard-looking object becomes boring.

The `p=5` exhaustive audit is the cleanest lesson.  Every one-line-per-direction
family, all `15625` of them, has full distance support.  Yet their union sizes
vary:

```text
17, 18, 19, 21, 25
```

and their pinned-distance minima split between `4` and `5`.  So Falconer support
is not the proof object.  The retained packet is pinned fibers plus concurrency
multiplicity.

That maps almost too neatly onto the current LRC `n=14` state.  S655 found that
the scalar `(14,21)` branch is only the shadow; the actual wall is the
off-diagonal odd pairs `(1,13),(3,11),(5,9)`.  S654 found that the carry `k` in
`v=r+27k` controls both parity and the `mod 14` obstruction.  HYP-2222 tells us
the `C=27` gcd shell is a fixed clock.

The Kakeya/Falconer model says what to do next: hold the scalar wall fixed and
jackknife the owner data.  In finite fields, hold direction coverage or full
distance support fixed and vary line intercepts; the hidden variable is
concurrency and pinned distance.  In LRC `n=14`, hold the odd wall/pair-sum
support fixed and vary `C=27` shell/carry labels; the hidden variable is
owner/carry conservativity.

So the new proof target is a no-leak lemma:

```text
odd wall fixed + C=27 shell fixed + carry/owner visible
  => AP, Vstar, nonprimitive 2AP, or strict looseness.
```

Finite-field Kakeya/Falconer does not prove that lemma.  It gives the right
mental model for why it should be narrower than a global search: direction and
distance scalars saturate early, while pinned/concurrency side channels still
carry the obstruction.

The S659 Tournament Analysis ranked pinned distance fibers and concurrency
multiplicity highest, with unit-distance direction spine and direction coverage
just behind.  Raw cardinality came last.  That is exactly the repo's recent
theme in another costume: do not trust the scalar after it starts looking
beautiful.
