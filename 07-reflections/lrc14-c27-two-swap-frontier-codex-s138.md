# LRC14 C=27 Two-Swap Frontier

S138 found exactly the kind of detail that makes a proof route feel less
hand-wavy.

The first K33 frontier did not move.  In the exact AP two-swap bank through
added values `<=40`, after the rigorous `q_threshold>=14` filter, the rows with
`M<=3/41` are still only:

```text
AP, GW 12->24, near-miss 12->36.
```

That is a good stability signal for the HYP-2937 split.

But the `2/27` layer did move.  Besides the old single-swap petals `10->20`
and `13->26`, two genuine two-hole rows appear:

```text
drop(10,12)->add(20,24)
drop(10,12)->add(20,36)
```

They are both exactly `M=2/27`.  This is not a new wilderness.  It is a splice:
the unit petal `10->20` glued to the known `12` branch, once with the tight GW
transfer and once with the K33 near-miss transfer.

In C27 shell language:

```text
10:g1 hole -> 7:g1 double
12:g3 hole -> 3:g3 double   or   12:g3 hole -> 9:g9 double
```

That suggests the right local theorem is not "two replacements introduce no
new rows."  They do.  The better theorem is: two replacements introduce only
splices of existing typed defects, and the unit-visible component already
forces the second-gap witness.

This also clarifies why the tournament quotient has to be chosen carefully.
Raw runner tournaments would make the two-hole rows look like bigger
perturbations.  The shell-transfer quotient says something more useful:

```text
one unit defect plus one 12-branch nonunit defect.
```

The current proof instinct is now:

```text
1. q-threshold filter to q>=14 for any second-gap row.
2. classify unit-visible C27 holes.
3. prove the only two-hole low unit petal is 10->20.
4. prove the companion nonunit hole must be the 12 branch.
5. dispatch 12->24 as GW-tight and 12->36 as K33/state-lift loose.
```

The guardrail remains unchanged: transfer labels recur in loose rows.  The
quotient is not the proof.  But it is now a much sharper interface for a local
finite lemma, and it gives the next agent a clean target instead of a cloud of
"maybe two swaps matter."
