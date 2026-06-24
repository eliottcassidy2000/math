# LRC14 C27 Unital Block Lift

The q=3 unital finally became useful once I stopped asking it to be canonical.

S107 had already warned that the Hermitian unital `2-(28,4,1)` is not a natural
`S8`-invariant block system on the `28=C(8,2)` AP8 pair slots.  That warning is
still right.  But the user asked for the more careful thing: attach AP and
Goddyn-Wong labels first, then see whether HYP-2937's marked C27 transfers lift
into 4-point unital blocks.

The clean labelling is:

```text
AP, GW, H1..H13, D1..D13.
```

That uses all 28 unital points exactly.  Calibrate one block as:

```text
{AP, GW, H12, D3}.
```

Now the Goddyn-Wong transfer `H12 -> D3` is not just near the AP/GW story; it is
the AP/GW anchor block.

The two-splice result is the best part:

```text
10+12 -> 20+24:
  petal block and AP/GW block are disjoint.

10+12 -> 20+36:
  petal block and near/K33 block share D9.
```

That turns the S138 observation into a small incidence picture:

```text
AP/GW --H12-- near/K33 --D9-- petal10.
```

This is exactly the kind of packet that might feed the HYP-2908 state-lift
endpoint.  It is not a proof, but it is no longer a metaphor.  It is an actual
`2-(28,4,1)` incidence object with AP/GW labels attached and exact transfer
pairs lifted to unique completion blocks.

The guardrail is equally important: the labelling is noncanonical.  The point
is not "the unital is the AP8 pair graph."  The point is:

```text
given a marked transfer pair, the unital supplies one 4-point completion.
```

So the next theorem target is a packet classification:

```text
low-gap non-AP/GW residual
  -> AP/GW anchor block
  or unit-petal block
  or AP/GW--near/K33--petal chain.
```

That would give the C27 branch a finite pair-completion grammar before the
terminal forbidden-H state lift.
