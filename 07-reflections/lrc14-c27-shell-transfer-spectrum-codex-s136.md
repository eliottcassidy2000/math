# LRC14 C=27 Shell-Transfer Spectrum

The useful move in S136 was to stop saying "C27 shell" as if it were a scalar
tag.  The shell itself is a relation space:

```text
P_a = {a, 27-a}, 1 <= a <= 13.
```

The AP is the perfect lower transversal.  Goddyn-Wong is not random near-AP
noise; it is a marked nonunit transfer:

```text
hole at 12 (gcd 3), doubled shell at 3 (gcd 3).
```

The `12->36` near-miss keeps the same nonunit hole but moves the surplus into
the unique gcd-9 shell:

```text
hole at 12 (gcd 3), doubled shell at 9 (gcd 9).
```

That is a cleaner reason it belongs with the first K33/Farey-child branch.  It
is not merely "bigger than GW"; it pushes a gcd-3 defect into the deepest
fixed 3-adic layer.

The two `2/27` petals are different again.  They have unit-visible holes:

```text
10->20: hole at 10 (unit), double at 7 (unit)
13->26: hole at 13 (unit), double at 1 (unit)
```

That makes them exactly the S571 second-gap objects: addition creates the
missed shell, and multiplication by a unit sees it.

The bounded audit through the S130 single-replacement atlas, replacement
`<=140`, is strikingly small:

```text
M <= 3/41: AP, GW, 12->36.
M <= 2/27: AP, GW, 12->36, 10->20, 13->26.
```

This does not prove LRC14.  It does make the proof route less foggy.  The
current frontier is not just AP/GW plus loose petals.  It is a finite set of
marked transfers across the C=27 unit/nonunit quotient, with exact Farey labels
attached.

The guardrail matters.  Shell-transfer labels recur in safely loose rows, and
perfect C=27 transversals can also be loose.  So the shell quotient is not a
complete invariant.  It is a carrier that must sit below exact `M=p/q` and
above raw additive energy.

My current proof instinct is:

```text
exact Farey branch first,
marked C27 transfer second,
unit/nonunit visibility third,
then K33/state-lift or petal rigidity.
```

The next theorem to chase is a universalization of the bounded frontier:
after the established finite-core reductions, low-gap non-AP/GW atoms should
have either a unit-visible C27 hole or the gcd3-to-gcd9 transfer.  If that is
true, the p=2 branch goes to petal/two-block rigidity and the p>=3 branch gets
a much more concrete packet for HYP-2908.
