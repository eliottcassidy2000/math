# LRC n=14 / Unit-Distance n=21 Difficulty Bridge S641

The clean answer is yes, but the connection is not the tempting one.

It is not:

```text
14 relates to 21, or unit-distance n=21 realizes forbidden H=21.
```

It is:

```text
27-quantum + retained section + bulk/lift side channel.
```

On the LRC side, `27` is the shell modulus `C=2n-1` at `n=14`.  This is the
first even row where the shell action does not stay transitive: the `<2,-1>`
fold leaves gcd strata `1,3,9`.  The base `Res_27` quotient is now highly
classified: AP, `V*`, and nonprimitive `2*AP` are the only floor rows in the
least-positive section.  The open seam is lift/CRT carry and owner data.

On the unit-distance side, `27` is the Moser slab edge increment.  THM-408 gives

```text
E(P_m^+) = 27m + 6,
E(P_m^-) = 27m + 3.
```

The exact `n=21` Moser row is `P_2^-`, so it has `57=2*27+3` unit edges.  Its
Hamiltonian unit spine has `20` edges and the remaining `37=C_hex(3)` edges are
bulk.  The section is the unit spine; the open seam is bulk/ear/embedding data.

That makes the real analogy pretty tight:

```text
LRC unit-shell section   <-> unit Hamiltonian spine
LRC carry/owner lift     <-> Moser bulk/ear gluing
gcd-9 high-depth packet  <-> small direction/gain packets
lift conservativity      <-> endpoint-ear conservativity
```

The useful transfer is therefore not a theorem equating the two problems.  It
is a proof architecture:

1. Name the section that makes the object tractable.
2. Name the side channel the section forgets.
3. Damage that side channel one coordinate at a time.
4. Prove every stable damaged state either routes to a known certificate or is
   one of the named residual cocycles.

For LRC `n=14`, that says: do not enumerate more raw rows until the carry
cocycle has been priced against owner and pinch labels.

For unit-distance `n=21 -> 22`, that says: do not widen the beam before the
exact `21` cores have an endpoint-compatible ear ledger recording direction
support, gain packet, spine endpoint, bulk class, and obstruction labels.

The pleasing part is that the bridge goes through an actual shared arithmetic
quantity.  The LRC modulus and the Moser slab edge quantum are both `27`.  The
danger is that the same observation can become numerology if the section and
side channel are dropped.  The whole point is to keep them.

The incoming monad S6 exhaustive H-spectrum result adds a useful warning light:
at `n=9`, the low forbidden H-values are still exactly `[7,21]`, while the old
high gaps fill in.  So `21` is a real durable tournament obstruction value.
But that should not tempt us into saying the unit-distance `21` core is an
`H=21` object.  The better reading is two-channel: tournament `21` marks the
persistent forbidden value, and LRC/UD share the `27` carrier whose side channel
explains why proving the frontier is hard.
